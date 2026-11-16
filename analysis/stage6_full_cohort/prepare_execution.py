"""Prepare explicit full-cohort inputs from the locked intake; no expression tests."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
readiness = PROJECT / "analysis/stage6_readiness"
intake = json.loads((PROJECT / "analysis/stage6_specification/cohort_intake.json").read_text())
pilot = json.loads((PROJECT / "analysis/stage6_pilot/acquisition_manifest.json").read_text())
ledger = json.loads((readiness / "donor_library_ledger.json").read_text())["cohorts"]
files = [r for r in pilot["files"] if r["cohort"] != "GSE159354"]
for row in ledger["GSE159354"]["library_rows"]:
    url, = row["count_or_component_urls"]
    files.append({"cohort": "GSE159354", "url": url,
                  "filename": url.rsplit("/", 1)[1], "role": "native_library:" + row["gsm"]})
for filename, size in (("hs_he_fullrers_images.zip", 757472562), ("hs_image_alignment_jsons.zip", 1981556)):
    files.append({"cohort": "Franzen", "url": "https://www.ebi.ac.uk/biostudies/files/S-BSST1410/" + filename,
                  "filename": filename, "role": "original_histology_or_alignment", "expected_bytes": size})
assert len(files) == 58 and len({(r["cohort"], r["filename"]) for r in files}) == 58
(HERE / "acquisition_manifest.json").write_text(json.dumps({"scope": "Full native intake and original spatial histology; external candidate outcomes remain reserved", "files": files}, indent=2) + "\n")
(HERE / "cohort_intake.json").write_text(json.dumps(intake, indent=2) + "\n")
(HERE / "donor_library_ledger.json").write_text(json.dumps(ledger, indent=2) + "\n")
print("Prepared 58 declared source files for 134 native intake units.")
