"""Audit complete preparation evidence without declaring human QC or anatomy passed.

Raw count payloads are not rehashed: their parser/review hashes are linked and
HDF5 structural metadata is rechecked. Exported review tables and portable
artifacts are checked directly, including every profile's pending review state.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from html.parser import HTMLParser
import json
import math
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit


RESOURCES = {"GSE135893": 43, "GSE159354": 48, "GSE159585": 18}
SOURCE_PROFILE_COUNTS = {"GSE135893": 114396, "GSE159354": 184019, "GSE159585": 148154}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024**2), b""):
            h.update(block)
    return h.hexdigest()


def axis_hash(values):
    h = hashlib.sha256()
    for value in values:
        raw = str(value).encode()
        h.update(len(raw).to_bytes(8, "little"))
        h.update(raw)
    return h.hexdigest()


class Audit:
    def __init__(self):
        self.checks, self.limitations = [], []

    def require(self, name, condition, **context):
        self.checks.append({"check": name, "passed": bool(condition), **context})
        return bool(condition)

    def load(self, path):
        try:
            value = json.loads(path.read_text())
            self.require("readable JSON", True, path=str(path))
            return value
        except Exception as error:
            self.require("readable JSON", False, path=str(path), error=str(error))
            return {}

    def file(self, path, digest=None, size=None):
        if not self.require("artifact exists", path.is_file(), path=str(path)):
            return False
        if size is not None:
            self.require("artifact byte count", path.stat().st_size == size, path=str(path))
        if digest is not None:
            self.require("artifact SHA256", sha256(path) == digest, path=str(path))
        return True


def h5_metadata(path):
    import h5py
    with h5py.File(path, "r") as handle:
        group = handle["X"]
        shape = tuple(int(n) for n in group.attrs["shape"])
        obs_index = handle["obs"].attrs["_index"]
        var_index = handle["var"].attrs["_index"]
        return {
            "shape": shape, "nnz": int(group["data"].shape[0]),
            "integer_dtype": group["data"].dtype.kind in "iu",
            "sparse_encoding": group.attrs["encoding-type"],
            "obs_hash": axis_hash(handle["obs"][obs_index].asstr()[()]),
            "var_hash": axis_hash(handle["var"][var_index].asstr()[()]),
        }


def acquisition(audit, run):
    manifest = audit.load(run / "config/acquisition_manifest.json")
    files = manifest.get("files", [])
    audit.require("58 distinct declared original files", len(files) == 58 and len({(r["cohort"], r["filename"]) for r in files}) == 58)
    source_root = run.parents[1]
    total_bytes = 0
    for item in files:
        record = audit.load(source_root / "data/acquisition_records" / f"{item['cohort']}__{item['filename']}.json")
        path = source_root / "data/raw" / item["cohort"] / item["filename"]
        audit.require("acquisition complete and source URL/hash recorded", record.get("status") == "complete" and record.get("url") == item["url"] and re.fullmatch(r"[0-9a-f]{64}", record.get("sha256", "")) is not None, path=str(path))
        audit.file(path, size=record.get("bytes"))
        if item.get("expected_sha256"):
            audit.require("acquisition matches frozen small-source SHA", record.get("sha256") == item["expected_sha256"], path=str(path))
        total_bytes += record.get("bytes", 0)
    audit.require("58 original payloads sum to declared full acquisition", total_bytes == 3993968410, observed_bytes=total_bytes)
    return {"original_source_files": len(files), "original_source_bytes": total_bytes}


def review_table(audit, path, unit, raw, summary):
    count = total = detected_sum = zero = 0
    flags = Counter()
    states_ok = identity_ok = ids_ok = numbers_ok = doublets_ok = True
    unscored = 0
    ids = []
    try:
        with gzip.open(path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                count += 1
                ids.append(row["source_obs_name"])
                value = float(row["review_total_counts"])
                features = float(row["review_detected_source_features"])
                numbers_ok &= math.isfinite(value) and value >= 0 and value.is_integer()
                numbers_ok &= math.isfinite(features) and features >= 0 and features.is_integer()
                total += int(value)
                detected_sum += int(features)
                is_zero = row["review_zero_counts"] == "True"
                zero += is_zero
                for n in (3, 4, 5):
                    flags[str(n)] += row[f"review_any_{n}mad"] == "True"
                states_ok &= all(row[k] == "not_adjudicated" for k in (
                    "review_expert_qc_status", "review_expert_identity_status", "review_endpoint_eligibility"))
                states_ok &= row["review_automated_disposition"] == (
                    "HARD_FAIL_ZERO_COUNTS_RETAINED" if is_zero else "REVIEW_PENDING")
                identity_ok &= row["review_identity_status"] == "provisional_requires_human_review"
                ids_ok &= row["review_resource"] == raw["resource"]
                ids_ok &= row["review_library_id"] == unit["library_or_section_id"]
                ids_ok &= row["review_donor_uid"] == unit["donor_uid"]
                available = row["review_doublet_score_available"] == "True"
                unscored += not available
                threshold = summary["scrublet"].get("threshold")
                for prior in ("0.02", "0.05", "0.10"):
                    score = row[f"review_doublet_score_prior_{prior}"]
                    flag = row[f"review_doublet_flag_prior_{prior}"]
                    if available:
                        score = float(score)
                        doublets_ok &= math.isfinite(score) and 0 <= score <= 1
                        doublets_ok &= (flag == str(score > threshold)) if threshold is not None else flag == ""
                    else:
                        doublets_ok &= score == "" and flag == ""
        audit.require("every RNA profile retained in review", count == raw["cells"], path=str(path))
        audit.require("review obs IDs unique and original-order hash matched", len(set(ids)) == count and axis_hash(ids) == summary["obs_axis_sha256"], path=str(path))
        audit.require("review original counts and detected-feature mass conserved", numbers_ok and total == raw["original_count_sum"] == summary["total_counts"] and detected_sum == raw["nnz"], path=str(path))
        audit.require("review zero/3-4-5MAD counts reproduce summary", zero == summary["hard_fail_zero_profiles"] and dict(flags) == summary["mad_review_flags"], path=str(path))
        audit.require("RNA expert decisions remain pending and hard-fail profiles retained", states_ok, path=str(path))
        audit.require("RNA identity labels remain explicitly provisional", identity_ok, path=str(path))
        audit.require("RNA profile source donor/library identity matches intake", ids_ok, path=str(path))
        audit.require("doublet scores or explicitly missing values preserve diagnostic meaning", doublets_ok and unscored == summary["scrublet"].get("unscored_profiles"), path=str(path))
    except Exception as error:
        audit.require("review profile table readable", False, path=str(path), error=str(error))


def rna(audit, run, intake):
    total_profiles, raw_axes, code_hashes = 0, {}, set()
    for resource, expected_count in RESOURCES.items():
        units = {u["library_or_section_id"]: u for u in intake["resources"][resource]["units"]}
        prepared = audit.load(run / "rna_raw" / resource / "inputs.json").get("records", [])
        prepared_by_id = {r["library_id"]: r for r in prepared}
        audit.require("exact full RNA raw intake", len(prepared) == expected_count and set(prepared_by_id) == set(units), resource=resource)
        raw_summary = audit.load(run / "rna_raw" / resource / "summary.json")
        audit.require("raw parser roundtrips complete; no final eligibility", raw_summary.get("all_count_roundtrips") is True and raw_summary.get("final_eligibility") == "not_adjudicated", resource=resource)
        audit.require("source profile total matches full deposited intake", sum(r["cells"] for r in prepared) == SOURCE_PROFILE_COUNTS[resource] == raw_summary.get("profiles"), resource=resource)
        reviewed = audit.load(run / "rna_review" / resource / "summary.json")
        captures = reviewed.get("capture_summaries", [])
        by_id = {s["record"]["library_id"]: s for s in captures}
        audit.require("all RNA captures reviewed without processing quarantine", reviewed.get("status") == "review_pack_generated" and reviewed.get("completed_captures") == expected_count and reviewed.get("expected_captures") == expected_count and not reviewed.get("failures") and len(captures) == expected_count and set(by_id) == set(units), resource=resource)
        audit.require("RNA review totals match raw intake", reviewed.get("retained_profiles") == SOURCE_PROFILE_COUNTS[resource], resource=resource)
        audit.require("RNA expert/Checkpoint A/discovery states remain false", all(reviewed.get(k) is False for k in ("expert_qc_adjudicated", "expert_identity_adjudicated", "checkpoint_A_passed", "formal_discovery_started")), resource=resource)
        summaries = {}
        for path in (run / "rna_review" / resource).glob("*/summary.json"):
            value = audit.load(path)
            if value.get("record"):
                summaries[value["record"]["library_id"]] = (path, value)
        for library, raw in prepared_by_id.items():
            unit = units.get(library)
            if unit is None:
                continue
            audit.require("raw RNA identity matches locked intake", raw["donor_id"] == unit["donor_id"] and raw["diagnosis"] == unit["source_diagnosis"] and raw["analysis_group"] == unit["analysis_group"] and raw.get("fraction") == unit.get("facs_fraction"), resource=resource, library=library)
            path = Path(raw["path"])
            raw_axes[raw["feature_axis"]] = raw["feature_axis_sha256"]
            if audit.file(path):
                try:
                    h5 = h5_metadata(path)
                    audit.require("raw H5AD structural axes/nnz/integer dtype", h5["shape"] == (raw["cells"], raw["features"]) and h5["nnz"] == raw["nnz"] and h5["integer_dtype"] and h5["sparse_encoding"] == "csr_matrix", resource=resource, library=library)
                except Exception as error:
                    audit.require("raw H5AD structure readable", False, path=str(path), error=str(error))
                    h5 = None
            else:
                h5 = None
            if library not in by_id or library not in summaries:
                continue
            review = by_id[library]
            summary_path, stored = summaries[library]
            code_hashes.add(review["script_sha256"])
            audit.require("stored capture summary matches cohort review summary", stored == review, path=str(summary_path))
            audit.require("review input SHA and raw record linked", raw == review["record"] and raw["checkpoint_sha256"] == review["input_sha256"] and review.get("source_counts_file_unchanged") is True and review.get("all_profiles_retained") is True and raw.get("count_roundtrip") is True, resource=resource, library=library)
            audit.require("review shape/count sum matches prepared integer counts", review["source_shape"] == [raw["cells"], raw["features"]] and review["total_counts"] == raw["original_count_sum"], resource=resource, library=library)
            if h5:
                audit.require("raw and review obs/var identity hashes match", h5["obs_hash"] == review["obs_axis_sha256"] and h5["var_hash"] == review["var_axis_sha256"], resource=resource, library=library)
            for name, digest in review["artifacts"].items():
                audit.file(summary_path.parent / name, digest)
            review_table(audit, summary_path.parent / "obs_review.tsv.gz", unit, raw, review)
            diagnostic = review["scrublet"]
            audit.require("canonical exact-neighbor Scrublet or explicit missing diagnostic", diagnostic.get("use_approx_neighbors") is False and diagnostic.get("knn_dist_metric") == "euclidean" and diagnostic.get("status") in ("scored_threshold_requires_review", "scored_threshold_unresolved", "algorithm_failed_requires_review", "too_few_profiles_for_stable_doublet_review"), resource=resource, library=library)
            if diagnostic.get("unscored_profiles", 0) or diagnostic.get("status") != "scored_threshold_requires_review":
                audit.limitations.append({"library": f"{resource}::{library}", "doublet_diagnostic": diagnostic.get("status"), "unscored_profiles": diagnostic.get("unscored_profiles")})
            total_profiles += raw["cells"]
    for path, digest in raw_axes.items():
        audit.file(Path(path), digest)
    audit.require("109 RNA inputs total 446569 retained profiles", total_profiles == sum(SOURCE_PROFILE_COUNTS.values()), observed=total_profiles)
    return {"RNA_profiles": total_profiles, "unique_RNA_feature_axes": raw_axes, "RNA_review_script_hashes": sorted(code_hashes)}


def spatial(audit, run, intake):
    folder = run / "spatial_review"
    value = audit.load(folder / "summary.json")
    expected = {u["library_or_section_id"]: u for u in intake["resources"]["Franzen"]["units"]}
    sections = value.get("sections", [])
    audit.require("exact 25-section spatial intake", len(sections) == 25 and {s["section_id"] for s in sections} == set(expected))
    audit.require("all original histology and technical parsing complete", value.get("status") == "passed_technical_checks_human_review_pending" and value.get("whole_spaceranger_archive_crc_passed") is True and value.get("sections_completed") == 25 and value.get("original_histology_sections") == 25 and value.get("original_histology_dimension_holds") == 0)
    audit.require("global anatomy remains unreviewed", value.get("expert_review_status") == "unreviewed" and value.get("formal_anatomical_eligibility") == "pending_human_review")
    for section in sections:
        sid = section["section_id"]
        unit = expected.get(sid)
        if unit is None:
            continue
        audit.require("spatial source identity and raw count joins", section["donor_id"] == unit["donor_id"] and section["diagnosis"] == unit["source_diagnosis"] and section["capture_positions"] == 4992 and all(section.get(k) is True for k in ("all_position_and_raw_barcode_sets_equal", "filtered_count_submatrix_equal", "reload_verified", "unique_gene_ids")), section=sid)
        path = Path(section["h5ad"])
        if audit.file(path):
            try:
                h5 = h5_metadata(path)
                audit.require("spatial H5AD integer raw axes preserved", h5["shape"] == (4992, section["features"]) and h5["integer_dtype"], section=sid)
            except Exception as error:
                audit.require("spatial H5AD readable", False, section=sid, error=str(error))
        histology = section.get("histology", {})
        original = histology.get("original") or {}
        audit.require("original H&E uniquely mapped, dimension checked, expert unreviewed", histology.get("mapping", {}).get("status") == "unique_exact_token" and histology.get("dimension_status") == "dimensions_consistent_visual_review_pending" and histology.get("expert_status") == "unreviewed" and original.get("original_bytes_preserved") is True and original.get("crc_verified_on_read") is True, section=sid)
        if original.get("server_only_original_path"):
            audit.file(Path(original["server_only_original_path"]), original["sha256"], original["bytes"])
        axis = section.get("feature_axis", {})
        if axis.get("path"):
            audit.file(Path(axis["path"]), axis["sha256"])
        mask = audit.load(folder / sid / "anatomy_mask_template.json")
        audit.require("mask template unreviewed and empty means unreviewed", mask.get("review_status") == "unreviewed" and mask.get("endpoint_eligibility") == "not_adjudicated" and all(mask.get(k) == [] for k in ("artifact_polygons", "barrier_polygons", "compartment_polygons", "interface_polylines")), section=sid)
    with (folder / "human_review/section_review.tsv").open() as handle:
        queue = list(csv.DictReader(handle, delimiter="\t"))
    audit.require("all 25 anatomy queue records await human decisions", len(queue) == 25 and {r["section_id"] for r in queue} == set(expected) and all(r["expert_status"] == "unreviewed" and r["interface_status"] == "not_assessed" and r["endpoint_eligibility"] == "not_adjudicated" and not r["reviewer"] for r in queue))
    audit.require("full spatial capture and author-tissue totals", sum(s["capture_positions"] for s in sections) == 124800 and sum(s["author_in_tissue_positions"] for s in sections) == 101907)
    artifacts = audit.load(folder / "artifact_sha256.json").get("artifacts", [])
    for record in artifacts:
        path = folder / record["path"]
        # Raw H5AD content was verified at creation; do not rehash these payloads.
        audit.file(path, None if path.suffix == ".h5ad" else record["sha256"], record["bytes"])
    return {"spatial_capture_positions": sum(s["capture_positions"] for s in sections), "spatial_author_in_tissue_positions": sum(s["author_in_tissue_positions"] for s in sections)}


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links, self.ids = [], set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        for key in ("href", "src"):
            if attrs.get(key):
                self.links.append((tag, attrs[key]))


def package(audit, folder):
    manifest = audit.load(folder / "package_manifest.json")
    audit.require("portable package exact intake and pending expert gates", manifest.get("status") == "complete_review_material_human_adjudication_pending" and manifest.get("RNA_inputs") == 109 and manifest.get("RNA_source_profiles") == sum(SOURCE_PROFILE_COUNTS.values()) and manifest.get("spatial_sections") == 25 and manifest.get("exact_intake_ID_sets_matched") is True and all(manifest.get(k) is False for k in ("expert_QC_adjudicated", "expert_identity_adjudicated", "expert_anatomy_adjudicated", "checkpoint_A_passed", "formal_discovery_started")))
    for record in manifest.get("files", []):
        path = (folder / record["path"]).resolve()
        audit.require("portable artifact stays inside package", path.is_relative_to(folder.resolve()), path=str(path))
        audit.file(path, record["sha256"], record["bytes"])
    page = folder / "index.html"
    if not audit.file(page):
        return
    parser = Links()
    parser.feed(page.read_text())
    for tag, value in parser.links:
        parsed = urlsplit(value)
        if parsed.scheme or parsed.netloc:
            audit.require("offline page has no network dependency", False, target=value)
        elif parsed.path:
            target = (folder / unquote(parsed.path)).resolve()
            audit.require("offline HTML reference resolves locally", target.is_relative_to(folder.resolve()) and target.is_file(), target=value)
        else:
            audit.require("offline HTML anchor exists", parsed.fragment in parser.ids, target=value)
    text = page.read_text()
    audit.require("HTML contains 109 RNA and 25 spatial review cards", text.count('<details class="review-card"') == 134)
    link_check = audit.load(folder / "link_image_validation.json")
    audit.require("builder image verification passed", link_check.get("status") == "passed" and not link_check.get("errors") and link_check.get("images_checked", 0) >= 478)
    for name, count in (("rna_capture_decisions.csv", 109), ("spatial_section_anatomy.csv", 25), ("roi_inventory.csv", 0)):
        with (folder / "human_decisions" / name).open() as handle:
            rows = list(csv.DictReader(handle))
        audit.require("human decision template complete and unreviewed", len(rows) == count and all(not r.get("reviewer") and all(r[k] == "not_adjudicated" for k in ("endpoint_eligibility",) if k in r) for r in rows), file=name)


def genes(audit, run, intake, rna_axes):
    manifests = list((run / "genes").glob("**/reference_manifest.json"))
    audit.require("frozen gene reference manifest exists", len(manifests) == 1)
    reference_hashes = set()
    for path in manifests:
        value = audit.load(path)
        audit.require("gene identity target is fixed Ensembl116", value.get("common_identity_target") == "Homo sapiens GRCh38 Ensembl 116 (June 2026)")
        for key in ("ensembl_registry", "ensembl_genes", "hgnc_complete", "hgnc_withdrawn"):
            asset = value.get("assets", {}).get(key, {})
            audit.require("required gene reference asset recorded", bool(asset.get("path")) and bool(asset.get("sha256")), asset=key)
            if asset.get("path"):
                audit.file(path.parent / asset["path"], asset["sha256"], asset["bytes"])
        reference_hashes.add(sha256(path))
    mapping_summaries = list((run / "genes").glob("**/mapping_summary.json"))
    mapped_resources = set()
    mapped_hashes = set()
    for path in mapping_summaries:
        value = audit.load(path)
        mapped_resources.add(value.get("resource_id"))
        mapped_hashes.add(value.get("source_file_sha256"))
        audit.require("gene mapper preserves counts and links pinned references", value.get("counts_read_or_changed") is False and value.get("reference_manifest_sha256") in reference_hashes, path=str(path))
        audit.file(path.with_name("gene_mapping.tsv"), value.get("mapping_sha256"))
        if value.get("provisional_unique_rows", 0) or value.get("canonical_colliding_rows", 0):
            audit.limitations.append({"resource": value.get("resource_id"), "mapping": "Symbol-only/conflicting identities remain provisional or held; counts were not projected."})
    audit.require("all four resource gene mappings present", mapped_resources == set(RESOURCES) | {"Franzen"}, observed=sorted(str(x) for x in mapped_resources))
    audit.require("every RNA source feature axis has a mapping", set(rna_axes.values()).issubset(mapped_hashes))
    coverage = list((run / "genes").glob("**/coverage_summary.json"))
    audit.require("cross-resource coverage summary present", len(coverage) >= 1)
    for path in coverage:
        value = audit.load(path)
        audit.require("coverage uses all four resource axes", len(value.get("inputs", {})) == 4, path=str(path))
        if value.get("shared_strict_features") == 0:
            audit.limitations.append({"gene_identity": "Zero strict four-resource overlap reflects symbol-only source identities; unique candidates remain provisional, not absence of expression."})
    batch_paths = list((run / "genes").glob("**/batch_mapping_summary.json"))
    audit.require("full-unit gene mapping batch manifest exists", len(batch_paths) == 1)
    for path in batch_paths:
        batch = audit.load(path)
        audit.require("134 native units mapped without reading/changing counts", batch.get("status") == "complete_identity_metadata_human_review_pending" and batch.get("native_units_mapped") == 134 and batch.get("units_by_resource") == {**RESOURCES, "Franzen": 25} and batch.get("counts_read_or_changed") is False and batch.get("unique_axes_by_content") == 4)
        crosswalk = path.with_name("unit_feature_mapping_crosswalk.tsv")
        if not audit.file(crosswalk, batch.get("unit_crosswalk_sha256")):
            continue
        with crosswalk.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected = {u["unit_uid"]: u for resource in intake["resources"].values() for u in resource["units"]}
        audit.require("gene source-axis crosswalk covers exact full intake", len(rows) == 134 and {r["unit_uid"] for r in rows} == set(expected))
        checked_pairs = set()
        for row in rows:
            unit = expected.get(row["unit_uid"], {})
            audit.require("gene crosswalk preserves unit identity", row["donor_id"] == unit.get("donor_id") and row["library_or_section_id"] == unit.get("library_or_section_id") and row["mapping_status"] == "metadata_mapped_human_review_pending", unit=row["unit_uid"])
            pair = (row["source_axis_path"], row["mapping_path"])
            if pair in checked_pairs:
                continue
            checked_pairs.add(pair)
            audit.file(Path(row["source_axis_path"]), row["source_axis_sha256"])
            audit.file(Path(row["mapping_path"]), row["mapping_sha256"])
        for source_path, digest in batch.get("input_manifest_sha256", {}).items():
            audit.file(Path(source_path), digest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--package", type=Path)
    args = parser.parse_args()
    run = args.run.resolve(strict=True)
    output = args.output or run / "full_run_validation.json"
    audit = Audit()
    intake = audit.load(run / "config/cohort_intake.json")
    totals = {}
    if intake:
        audit.require("intake does not claim Checkpoint A/discovery", intake.get("checkpoint_A_passed") is False and intake.get("formal_discovery_started") is False)
        for task in ("acquisition", "RNA", "spatial", "genes", "package"):
            try:
                if task == "acquisition":
                    totals.update(acquisition(audit, run))
                elif task == "RNA":
                    totals.update(rna(audit, run, intake))
                elif task == "spatial":
                    totals.update(spatial(audit, run, intake))
                elif task == "genes":
                    genes(audit, run, intake, totals.get("unique_RNA_feature_axes", {}))
                else:
                    package(audit, args.package or run / "review_pack")
            except Exception as error:
                audit.require("audit stage completed", False, stage=task, error=repr(error))
    failures = [c for c in audit.checks if not c["passed"]]
    result = {
        "status": "passed_preparation_evidence_human_review_pending" if not failures else "incomplete_or_failed_preparation_evidence",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(), "run": str(run),
        "verification_depth": "Raw H5AD structural metadata, linked parser/review source hashes, profile-level review tables, rehashed small review/portable artifacts. Raw count payload hashes were verified during creation and review, not recomputed here.",
        "human_qc_passed": False, "human_identity_passed": False, "expert_anatomy_passed": False,
        "checkpoint_A_passed": False, "formal_discovery_started": False,
        "totals": totals, "checks_passed": len(audit.checks) - len(failures),
        "checks_failed": len(failures), "failed_checks": failures,
        "limitations": audit.limitations, "checks": audit.checks,
        "validator_sha256": sha256(Path(__file__)),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ("status", "checks_passed", "checks_failed", "failed_checks")}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
