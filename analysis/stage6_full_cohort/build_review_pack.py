"""Package completed capture/section evidence for offline human review."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
from html import escape
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shutil
from urllib.parse import unquote, urlsplit


RNA_IMAGES = ("qc_review.png", "identity_marker_review.png", "source_identity_marker_review.png", "scrublet_scores.png")
RNA_TABLES = ("mad_thresholds.tsv", "seed_program_coverage.tsv", "marker_axis_coverage.tsv", "identity_coverage.tsv", "group_marker_evidence.tsv.gz", "scrublet.json")
SPATIAL_IMAGES = ("tissue_overview.jpg", "fixed_grid_zoom.jpg", "coordinate_overlay.png", "raw_umi.png", "detected_genes.png", "author_in_tissue_mask.png")
BOUNDARY = "Review material is complete; QC, cell identities, anatomical inclusion and endpoint eligibility require human adjudication. No final cohort lock or biological discovery is reported."


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def safe(value):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def write_table(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def table(rows, fields):
    head = "".join(f"<th>{escape(title)}</th>" for _, title in fields)
    body = "".join("<tr>" + "".join(f"<td>{escape(str(row.get(key, '')))}</td>" for key, _ in fields) + "</tr>" for row in rows)
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.targets, self.ids = [], set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        for key in ("href", "src"):
            if attrs.get(key):
                self.targets.append((tag, key, attrs[key]))


def verify_package(output):
    parser = Links()
    parser.feed((output / "index.html").read_text())
    errors, images = [], []
    from PIL import Image
    for tag, kind, value in parser.targets:
        parts = urlsplit(value)
        if parts.scheme or parts.netloc:
            errors.append(f"Unexpected network link: {value}")
            continue
        if not parts.path:
            if parts.fragment not in parser.ids:
                errors.append(f"Missing page anchor: {value}")
            continue
        target = (output / unquote(parts.path)).resolve()
        if not target.is_relative_to(output.resolve()) or not target.is_file():
            errors.append(f"Missing/nonportable target: {value}")
        elif tag == "img" and kind == "src":
            try:
                with Image.open(target) as image:
                    width, height = image.size
                    image.verify()
                if width < 1 or height < 1:
                    raise ValueError("Empty image")
                images.append({"path": str(target.relative_to(output)), "width": width, "height": height})
            except Exception as error:
                errors.append(f"Unreadable image {value}: {error}")
    result = {"status": "passed" if not errors else "failed", "local_references_checked": len(parser.targets),
              "images_checked": len(images), "image_dimensions": images, "errors": errors}
    (output / "link_image_validation.json").write_text(json.dumps(result, indent=2) + "\n")
    if errors:
        raise ValueError(errors)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    run, output = args.run.resolve(strict=True), args.output.resolve()
    intake_path = run / "config/cohort_intake.json"
    intake = read_json(intake_path)
    expected_rna = {u["unit_uid"] for name, resource in intake["resources"].items() if name != "Franzen" for u in resource["units"]}
    expected_sections = {u["library_or_section_id"] for u in intake["resources"]["Franzen"]["units"]}
    if len(expected_rna) != 109 or len(expected_sections) != 25:
        raise ValueError("This review package expects the specified 109 RNA inputs and 25 sections")
    captures = {}
    for path in sorted((run / "rna_review").rglob("summary.json")):
        value = read_json(path)
        if "record" not in value or "input_sha256" not in value:
            continue
        record = value["record"]
        uid = f"{record['resource']}::{record['library_id']}"
        if uid in captures:
            raise ValueError("Duplicate capture summary: " + uid)
        if not value.get("all_profiles_retained") or not value.get("source_counts_file_unchanged"):
            raise ValueError("Raw-data retention check failed: " + uid)
        if record.get("checkpoint_sha256") and record["checkpoint_sha256"] != value["input_sha256"]:
            raise ValueError("Prepared checkpoint hash differs from review input: " + uid)
        if record.get("cells") is not None and int(record["cells"]) != int(value["source_shape"][0]):
            raise ValueError("Prepared profile count differs from review input: " + uid)
        captures[uid] = (path, value)
    if set(captures) != expected_rna:
        raise ValueError({"missing_RNA_inputs": sorted(expected_rna - set(captures)), "unexpected_RNA_inputs": sorted(set(captures) - expected_rna)})
    spatial_path = run / "spatial_review/summary.json"
    spatial = read_json(spatial_path)
    sections = {s["section_id"]: s for s in spatial["sections"]}
    if set(sections) != expected_sections or len(spatial["sections"]) != 25:
        raise ValueError("Spatial section identities differ from locked intake")
    if spatial.get("status") != "passed_technical_checks_human_review_pending":
        raise ValueError("Spatial material has not completed technical preparation")
    if any(s.get("histology", {}).get("expert_status") != "unreviewed" for s in sections.values()):
        raise ValueError("Use a separate reviewed release; this builder preserves an unadjudicated review package")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Review output must be empty; preserve existing review decisions and choose a new directory")
    output.mkdir(parents=True, exist_ok=True)
    copied, expected_artifacts = [], {}
    artifact_index = run / "spatial_review/artifact_sha256.json"
    if artifact_index.exists():
        expected_artifacts.update({str((run / "spatial_review" / item["path"]).resolve()): item["sha256"]
                                   for item in read_json(artifact_index)["artifacts"]})
    for path, summary in captures.values():
        expected_artifacts.update({str((path.parent / name).resolve()): digest for name, digest in summary.get("artifacts", {}).items()})

    def copy(source, relative):
        target = output / relative
        if target.exists():
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        source_hash = sha256(source)
        expected = expected_artifacts.get(str(source.resolve()))
        if expected is not None and source_hash != expected:
            raise ValueError("Source evidence changed after its processing manifest: " + str(source))
        if source_hash != sha256(target):
            raise ValueError("Copied artifact differs: " + str(relative))
        copied.append({"source": str(source), "package_path": str(relative), "bytes": target.stat().st_size, "sha256": source_hash})
        return str(relative)

    def asset(source, relative):
        relative = copy(source, relative)
        return f'<a href="{escape(relative)}"><img loading="lazy" src="{escape(relative)}" alt="{escape(source.stem.replace("_", " "))}"></a>'

    def links(folder, names, relative_folder):
        result = []
        for name in names:
            if (folder / name).is_file():
                target = copy(folder / name, relative_folder / name)
                result.append(f'<a href="{escape(target)}">{escape(name)}</a>')
        return " · ".join(result)

    copy(intake_path, Path("provenance/cohort_intake.json"))
    copy(spatial_path, Path("provenance/spatial_preparation_summary.json"))
    rna_cards, capture_rows, decisions = [], [], []
    resource_totals = defaultdict(Counter)
    donor_sets = defaultdict(set)
    source_inputs = []
    for uid, (path, summary) in captures.items():
        record, folder = summary["record"], path.parent
        relative = Path("rna") / (safe(uid) + "_" + hashlib.sha256(uid.encode()).hexdigest()[:8])
        copy(path, relative / "summary.json")
        pictures = []
        for name in RNA_IMAGES:
            if (folder / name).is_file():
                pictures.append(asset(folder / name, relative / name))
            elif name != "source_identity_marker_review.png":
                raise FileNotFoundError(folder / name)
        evidence_links = links(folder, RNA_TABLES, relative)
        count = int(summary["source_shape"][0])
        flags = int(summary["mad_review_flags"]["3"])
        doublets = summary["scrublet"]
        row = {"resource": record["resource"], "donor_id": record["donor_id"], "library_id": record["library_id"],
               "source_diagnosis": record.get("diagnosis", "unknown"), "analysis_group": record.get("analysis_group", "unknown"),
               "assay": record["assay"], "fraction": record.get("fraction") or "unknown",
               "source_profiles": count, "source_features": summary["source_shape"][1], "review_flags_3mad": flags,
               "zero_count_profiles": summary["hard_fail_zero_profiles"], "doublet_status": doublets["status"],
               "doublet_unscored_profiles": doublets["unscored_profiles"], "qc_status": "not_adjudicated", "identity_status": "not_adjudicated"}
        capture_rows.append(row)
        resource_totals[record["resource"]].update({"inputs": 1, "profiles": count, "flagged_3mad": flags,
                                                  "zero": row["zero_count_profiles"], "doublets_unscored": row["doublet_unscored_profiles"]})
        donor_sets[record["resource"]].add(record["donor_id"])
        source_inputs.append({"unit_uid": uid, "source_counts_checkpoint": record["path"], "source_counts_sha256": summary["input_sha256"],
                              "preparation_script_sha256": summary["script_sha256"], "summary": str(path), "summary_sha256": sha256(path)})
        decisions.append({**row, "donor_uid": f"{record['resource']}::{record['donor_id']}", "reviewer": "", "review_date": "",
                          "qc_decision": "not_adjudicated", "identity_decision": "not_adjudicated", "endpoint_id": "", "endpoint_eligibility": "not_adjudicated",
                          "affected_profiles_or_identity": "", "decision_reason": "", "follow_up_action": ""})
        title = f"{record['resource']} · {record['library_id']} · {record['donor_id']}"
        detail = f"{row['source_diagnosis']} · {row['assay']} · fraction: {row['fraction']} · {count:,} retained profiles · {flags:,} 3-MAD review flags"
        caution = ("Check joint QC distributions and potential cell-type loss. Inspect the simulated score distribution before using a doublet flag. "
                   "Compare positive and competing identity markers with original labels. An unusual identity or mixed program requires review.")
        rna_cards.append(f'<details class="review-card" data-search="{escape(title + " " + detail)}"><summary><strong>{escape(title)}</strong><span>{escape(detail)}</span></summary>'
                         f'<p class="state">QC and identity: review pending</p><p>{escape(caution)}</p>'
                         f'<p>Doublet diagnostic: {escape(doublets["status"])}; unscored profiles: {row["doublet_unscored_profiles"]:,}. '
                         f'Zero-total profiles retained with hard-fail flags: {row["zero_count_profiles"]:,}. '
                         f'Capture provenance: {escape(record.get("capture_identity_status", "unknown"))}.</p>'
                         f'<div class="image-grid">{"".join(pictures)}</div><p class="downloads">{evidence_links}</p></details>')

    templates = output / "human_decisions"
    templates.mkdir()
    write_table(templates / "rna_capture_decisions.csv", decisions, list(decisions[0]))
    write_table(output / "rna_capture_overview.csv", capture_rows, list(capture_rows[0]))
    for path in sorted((run / "rna_review").rglob("*.tsv")):
        if path.name not in ("library_identity_coverage.tsv", "donor_identity_coverage.tsv", "library_review_summary.tsv"):
            continue
        copy(path, Path("tables") / path.relative_to(run / "rna_review"))

    spatial_cards = []
    with (run / "spatial_review/human_review/section_review.tsv").open() as handle:
        section_rows = list(csv.DictReader(handle, delimiter="\t"))
    if {r["section_id"] for r in section_rows} != expected_sections:
        raise ValueError("Anatomical review queue is not the full 25-section intake")
    for row in section_rows:
        sid, section = row["section_id"], sections[row["section_id"]]
        folder, relative = run / "spatial_review" / sid, Path("spatial") / sid
        pictures = [asset(folder / name, relative / name) for name in SPATIAL_IMAGES]
        evidence_links = links(folder, ("section_summary.json", "review_image_provenance.json", "anatomy_mask_template.json"), relative)
        row["review_image"] = str(relative / "tissue_overview.jpg")
        row["reviewer"] = row["reviewer_role"] = row["review_date"] = ""
        row["expert_status"], row["endpoint_eligibility"] = "unreviewed", "not_adjudicated"
        title = f"{sid} · {row['donor_id']} · {row['diagnosis']}"
        context = f"{row['region_label']} · block {row['block_id']} · replicate {row['replicate']} · {row['procurement']}"
        spatial_cards.append(f'<details class="review-card" data-search="{escape(title + " " + context)}"><summary><strong>{escape(title)}</strong><span>{escape(context)}</span></summary>'
                             '<p class="state">Anatomy and endpoint eligibility: unreviewed</p>'
                             '<p>Review H&amp;E and the fixed image grid first. Confirm compartment and basal-facing epithelial–stromal interface. '
                             'Mark lumens, folds, tears, tissue edges and disconnected islands. Empty mask templates mean no review has been performed.</p>'
                             f'<p>Original histology: {escape(row["histology_mapping_status"])}; image dimensions: {escape(row["histology_dimension_check"])}. '
                             f'Raw capture positions: {section["capture_positions"]:,}; author in-tissue: {section["author_in_tissue_positions"]:,}.</p>'
                             f'<div class="image-grid">{"".join(pictures)}</div><p class="downloads">{evidence_links}</p></details>')
    write_table(templates / "spatial_section_anatomy.csv", section_rows, list(section_rows[0]))
    with (run / "spatial_review/human_review/roi_inventory.tsv").open() as handle:
        roi_fields = next(csv.reader(handle, delimiter="\t"))
    write_table(templates / "roi_inventory.csv", [], roi_fields)
    contact = asset(run / "spatial_review/human_review/all25_histology_contact_sheet.jpg", Path("spatial/all25_histology_contact_sheet.jpg"))
    gene_links = []
    gene_reports = {"mapping_summary.json", "coverage_summary.json", "reference_manifest.json",
                    "batch_mapping_summary.json", "unit_feature_mapping_crosswalk.tsv",
                    "identity_marker_coverage.tsv", "identity_marker_coverage_summary.json",
                    "shared_strict_features.tsv", "shared_unique_candidates.tsv"}
    for path in sorted((run / "genes").rglob("*")):
        if not path.is_file() or (path.name not in gene_reports and path != run / "genes/RESULTS.md"):
            continue
        relative = Path("genes") / path.relative_to(run / "genes")
        target = copy(path, relative)
        gene_links.append(f'<li><a href="{escape(target)}">{escape(str(path.relative_to(run / "genes")))}</a></li>')
    resource_rows = [{"resource": name, "inputs": counts["inputs"], "source_donors": len(donor_sets[name]), "profiles": f"{counts['profiles']:,}",
                      "flags": f"{counts['flagged_3mad']:,}", "unscored": f"{counts['doublets_unscored']:,}"} for name, counts in sorted(resource_totals.items())]
    overview = table(resource_rows, [("resource", "Resource"), ("inputs", "RNA inputs"), ("source_donors", "Source donor labels"), ("profiles", "Retained profiles"), ("flags", "3-MAD review flags"), ("unscored", "Unscored doublets")])
    captured_profiles = sum(v["profiles"] for v in resource_totals.values())
    page = '''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Stage 6 · Human review</title>
<style>:root{color-scheme:light;--ink:#183041;--muted:#536b7b;--line:#d8e2e8;--paper:#f4f7f9;--accent:#125e83}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:16px/1.55 system-ui,sans-serif}main{max-width:1260px;margin:auto;padding:34px 28px 80px}header{padding:28px;background:white;border:1px solid var(--line);border-radius:16px}h1{font-size:36px;margin:8px 0}h2{margin-top:38px;font-size:25px}h3{font-size:19px}.eyebrow{color:var(--accent);font-size:13px;letter-spacing:.12em;text-transform:uppercase}.intro{max-width:900px}.state{display:inline-block;color:#754a00;background:#fff2d4;border:1px solid #ead6a3;padding:5px 10px;border-radius:8px;font-weight:600}.metrics{display:flex;gap:14px;flex-wrap:wrap;margin:22px 0}.metric{padding:14px 20px;background:#edf4f8;border-radius:10px;min-width:190px}.metric b{font-size:26px;display:block}.metric span{color:var(--muted);font-size:13px}nav{display:flex;gap:14px;flex-wrap:wrap;margin:24px 0}a{color:var(--accent);text-decoration-thickness:1px;text-underline-offset:3px}p{max-width:1000px}.note{border-left:4px solid #b87a1a;padding:12px 18px;background:#fff9ec}.table-wrap{overflow:auto;background:white;border:1px solid var(--line);border-radius:12px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:12px;text-align:left;border-bottom:1px solid var(--line);white-space:nowrap}th{background:#edf3f7}.review-card{margin:12px 0;background:white;border:1px solid var(--line);border-radius:12px;padding:0 18px}.review-card summary{cursor:pointer;padding:16px 0}.review-card summary strong,.review-card summary span{display:block}.review-card summary span{font-size:13px;color:var(--muted);margin-top:4px}.image-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.image-grid img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px}.downloads{font-size:13px;overflow-wrap:anywhere}.contact img{max-width:100%;max-height:1000px;object-fit:contain;background:white;border:1px solid var(--line)}input{width:min(640px,100%);padding:12px;border:1px solid #9eb5c4;border-radius:8px;font:inherit;margin:8px 0}button{padding:10px 16px;border:1px solid #aabdc9;background:white;border-radius:8px;cursor:pointer;color:var(--ink)}li{margin-bottom:7px}footer{margin-top:36px;color:var(--muted);font-size:13px}@media(max-width:720px){main{padding:18px 14px}header{padding:18px}h1{font-size:28px}.image-grid{grid-template-columns:1fr}}</style></head><body><main>'''
    page += f'<header><div class="eyebrow">IPF Stage 6 · Preparation evidence</div><h1>Human review: QC, cell identity and tissue anatomy</h1><p class="intro">{escape(BOUNDARY)}</p><span class="state">Human decisions pending</span>'
    page += f'<div class="metrics"><div class="metric"><b>{len(captures)} / 109</b><span>RNA inputs prepared</span></div><div class="metric"><b>{captured_profiles:,}</b><span>Original RNA profiles retained</span></div><div class="metric"><b>{len(sections)} / 25</b><span>Spatial sections prepared</span></div></div></header>'
    page += '<nav><a href="#overview">Overview</a><a href="#rna">RNA evidence</a><a href="#histology">Histology</a><a href="#genes">Gene identity</a><a href="#decisions">Decision templates</a><a href="#provenance">Provenance</a></nav>'
    page += '<section id="overview"><h2>Review the measurement, identity and usable question</h2>' + overview
    page += '<p>Source donor labels are counted within each resource. Repeated libraries, sorted fractions, sections and ROIs do not add independent people. All retained profiles remain represented; flags identify review work and are not final exclusions.</p>'
    page += '<div class="note"><strong>Interpretation limits:</strong> GSE159354 uses different FACS routes. Ordinary EPCAM+ and Triple− primary disease comparisons each start with at most three IPF and three control donors. Unsorted/special gates remain discovery and internal sensitivity inputs. GSE159585 procurement is fully confounded with disease. Franzén HC is post-mortem and IPF transplant/resection: prioritize within-IPF regional/interface questions and retain this restriction for HC comparisons.</div></section>'
    page += '<section id="decisions"><h2>Record human decisions</h2><ol><li>Review each capture’s joint QC, marker evidence and doublet diagnostics. Record retained, held or excluded profiles/groups and the reason.</li><li>Review every section’s H&amp;E, registration and tissue geometry. Inventory all appropriate and uncertain regions before examining candidate RNA.</li><li>Confirm which endpoint each donor/region can support. A missing target interface or sparse identity affects that endpoint and need not discard all donor measurements.</li></ol>'
    page += '<p>Download editable CSV templates: <a href="human_decisions/rna_capture_decisions.csv">RNA capture decisions</a> · <a href="human_decisions/spatial_section_anatomy.csv">25-section anatomy review</a> · <a href="human_decisions/roi_inventory.csv">ROI inventory</a>. Save completed copies with reviewer, date and reasons. An automated or provisional annotation cannot supply expert sign-off.</p></section>'
    page += '<label for="search"><strong>Find a resource, donor, library or section</strong></label><br><input id="search" type="search" placeholder="For example: GSE159354, N1857, IPF_1"><button id="close-all" type="button">Collapse all cards</button>'
    page += '<section id="rna"><h2>RNA QC and provisional identity</h2><p>Open a capture to view QC, fixed seed-marker groups, original-label comparisons where available, and doublet scores. Missing feature rows mean coverage is unknown; a zero observed count does not prove biological absence. Unscored doublets remain explicitly missing.</p>' + "".join(rna_cards) + '</section>'
    page += '<section id="histology"><h2>All 25 sections: candidate-independent anatomy review</h2><p>Whole-field images and fixed 3×3 zoom grids cover the entire source image. They are resized review derivatives; original full-resolution images remain on the server. Technical count maps and author in-tissue flags do not certify compartments or direct cell contact.</p><div class="contact">' + contact + '</div>' + "".join(spatial_cards) + '</section>'
    page += '<section id="genes"><h2>Gene identity and coverage</h2><p>Review source-reference uncertainty, ambiguous mappings, duplicated symbols and missing export coverage. Source counts are preserved. Symbol-only correspondences remain distinguishable from source stable-ID evidence; zero strict overlap does not imply no biological expression.</p><p class="note">RNA marker plots use original source symbols. Multiple rows with one symbol are combined only for provisional review evidence; original count matrices are unchanged. In this run, CLIC5 identity is held in GSE135893/GSE159585; GSE159354 has two RGS5 source rows, including an ID-symbol conflict. These marker signals cannot independently certify identity or an unambiguous cross-resource gene test. Review the marker coverage table before adjudication.</p><ul>' + ("".join(gene_links) if gene_links else '<li>Gene reference/mapping reports are not present in this export; this gate remains pending.</li>') + '</ul></section>'
    page += '<section id="provenance"><h2>Evidence and portability checks</h2><p><a href="rna_capture_overview.csv">Capture overview</a> · <a href="provenance/cohort_intake.json">Locked intake</a> · <a href="provenance/source_inputs.json">Source input hashes</a> · <a href="package_manifest.json">Package file hashes</a> · <a href="link_image_validation.json">Local link and image checks</a></p><p>Full counts, per-cell review tables and original histology remain in the server run. This portable package contains review images, grouped marker evidence, summaries and editable decisions; no raw H5AD or original full-resolution histology was copied.</p></section>'
    page += '<footer>Generated ' + escape(datetime.now(timezone.utc).isoformat()) + ' UTC. The package records preparation evidence, not expert acceptance, a completed Checkpoint A or validated niches.</footer>'
    page += '''</main><script>const cards=[...document.querySelectorAll('.review-card')];document.getElementById('search').addEventListener('input',event=>{const query=event.target.value.toLowerCase().trim();cards.forEach(card=>{card.hidden=!card.dataset.search.toLowerCase().includes(query);});});document.getElementById('close-all').addEventListener('click',()=>cards.forEach(card=>card.open=false));</script></body></html>'''
    (output / "index.html").write_text(page)
    (output / "provenance/source_inputs.json").write_text(json.dumps({"RNA": source_inputs, "spatial_sources": spatial["sources"]}, indent=2) + "\n")
    (output / "README.md").write_text("# Stage 6 human review package\n\nOpen `index.html` locally. The package uses local images and no network assets.\n\n" + BOUNDARY + "\n\nFull source counts, per-cell review tables and original H&E remain on the server. Use the CSV files under `human_decisions/` to record reviewer/date/endpoint-specific decisions. Keep the generated originals and save reviewed copies separately.\n")
    # Create the manifest before HTML link validation, then refresh its hashes to
    # include the validation result. Self-hashing is deliberately excluded.
    manifest_path = output / "package_manifest.json"
    manifest_path.write_text("{}\n")
    (output / "link_image_validation.json").write_text("{}\n")
    verification = verify_package(output)
    manifest = {"status": "complete_review_material_human_adjudication_pending", "boundary": BOUNDARY,
                "created_utc": datetime.now(timezone.utc).isoformat(), "server_run": str(run), "builder_sha256": sha256(Path(__file__)),
                "RNA_inputs": len(captures), "RNA_source_profiles": captured_profiles, "spatial_sections": len(sections),
                "exact_intake_ID_sets_matched": True, "local_link_and_image_checks": verification["status"],
                "expert_QC_adjudicated": False, "expert_identity_adjudicated": False, "expert_anatomy_adjudicated": False,
                "checkpoint_A_passed": False, "formal_discovery_started": False, "copied_evidence": copied,
                "files": [{"path": str(path.relative_to(output)), "bytes": path.stat().st_size, "sha256": sha256(path)}
                          for path in sorted(output.rglob("*")) if path.is_file() and path != manifest_path]}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"output": str(output), "RNA_inputs": len(captures), "RNA_profiles": captured_profiles,
                      "spatial_sections": len(sections), "portable_files": len(manifest["files"]) + 1,
                      "local_links_checked": verification["local_references_checked"], "images_checked": verification["images_checked"]}))


if __name__ == "__main__":
    main()
