"""Check review-image provenance and exact histology/feature identity handling."""
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import zipfile

from PIL import Image


def main():
    module_path = Path(__file__).with_name("spatial_review.py")
    spec = importlib.util.spec_from_file_location("spatial_review_fixture", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    checks = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    with tempfile.TemporaryDirectory(prefix="stage6-spatial-review-fixture-") as work:
        folder = Path(work)
        image = Image.new("RGB", (613, 409), "pink")
        buffer = io.BytesIO()
        image.save(buffer, format="PNG")
        payload = buffer.getvalue()
        archive_path = folder / "histology.zip"
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("originals/V19S23-092-A1_HE.png", payload)
            archive.writestr("originals/V19S23-092-A10_HE.png", payload)
            archive.writestr("__MACOSX/._V19S23-092-A1_HE.png", b"ignored")
            archive.writestr("other/HC_2.TH011.B0.1_HE.png", payload)
            archive.writestr("duplicate/V19S23-092-B1_HE.png", payload)
            archive.writestr("duplicate/HC_3.TH012.B0.1_HE.png", payload)
        rows = [
            {"sample_id": "V19S23-092-A1", "sample_name": "HC_1.TH010.B0.1"},
            {"sample_id": "V19S23-092-C1", "sample_name": "HC_2.TH011.B0.1"},
            {"sample_id": "V19S23-092-B1", "sample_name": "HC_3.TH012.B0.1"},
            {"sample_id": "V19S23-092-D1", "sample_name": "HC_4.TH013.B0.1"},
        ]
        with zipfile.ZipFile(archive_path) as archive:
            mapping = module.histology_members(archive, rows)
            check("exact section excludes longer identifier and Mac metadata",
                  mapping[rows[0]["sample_id"]]["candidates"] == ["originals/V19S23-092-A1_HE.png"])
            check("exact source sample-name fallback", mapping[rows[1]["sample_id"]]["status"] == "unique_exact_token")
            check("ambiguous matches retained as hold", mapping[rows[2]["sample_id"]]["status"] == "ambiguous_exact_token")
            check("missing image remains missing", mapping[rows[3]["sample_id"]]["status"] == "missing_exact_token")
            extracted = folder / "original.png"
            provenance = module.extract_image(archive, mapping[rows[0]["sample_id"]]["candidates"][0], extracted)
        check("original image bytes and SHA conserved", extracted.read_bytes() == payload and provenance["sha256"] == hashlib.sha256(payload).hexdigest())
        decoded, orientation = module.open_histology(extracted)
        check("original dimensions and orientation preserved", decoded.size == (613, 409) and orientation == 1)
        module.review_images(decoded, folder, "Fixture | unreviewed", "original fixture", side=300)
        review = json.loads((folder / "review_image_provenance.json").read_text())
        boxes = [tile["source_pixel_box_xyxy"] for tile in review["zoom_grid"]]
        area = sum((x1 - x0) * (y1 - y0) for x0, y0, x1, y1 in boxes)
        check("fixed grid covers entire non-divisible image", len(boxes) == 9 and area == 613 * 409 and boxes[0][:2] == [0, 0] and boxes[-1][2:] == [613, 409])
        check("review derives from original pixels without expert signoff", review["expert_review_status"] == "unreviewed" and review["source_width_height"] == [613, 409])
        check("overview offset is explicit", review["overview_content_offset_xy"] == [0, 100])

        class Features:
            def iterrows(self):
                # Same symbol, distinct native gene IDs must remain distinct rows.
                for gene in ("ENSG00000000001", "ENSG00000000002"):
                    yield gene, {"gene_id": gene, "gene_symbol": "DUP.SYMBOL", "genome": "GRCh38", "feature_type": "Gene Expression"}

        feature_path = folder / "features.tsv.gz"
        feature_provenance = module.export_feature_axis(Features(), feature_path)
        with gzip.open(feature_path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            features = list(reader)
        check("feature schema matches shared mapper", reader.fieldnames == module.FEATURE_FIELDS)
        check("duplicate symbols preserve distinct gene IDs and row order", [r["source_gene_id"] for r in features] == ["ENSG00000000001", "ENSG00000000002"] and [r["source_row_1based"] for r in features] == ["1", "2"])
        check("feature hash and count recorded", feature_provenance["source_features"] == 2 and feature_provenance["sha256"] == module.sha256(feature_path))

    print(json.dumps({"status": "passed", "checks_passed": len(checks), "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
