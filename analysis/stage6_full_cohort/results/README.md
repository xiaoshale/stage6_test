# Review export availability

The complete offline export (1,410 files, 539 images) is verified locally in the research workspace and at the server run, not stored in full in this technical Git repository. This checkpoint contains source code, frozen gene-reference assets, complete small execution records, full-export manifests, human decision templates and representative images.

Server run: `/home/data/t120631/project/ipf_stage6_pilot/runs/preparation_001`. Retrieve `review_export.tar.gz` and `review_export_record.json` from that directory using Termius or SCP, verify the recorded SHA-256, and open `review_pack/index.html` after extraction. The transfer record contains the exact archive hash. All source counts, per-cell tables and full-resolution H&E remain on the server.

The package manifest describes the complete export; missing image files in this Git subset are deliberate. Human QC/identity/anatomy and endpoint eligibility remain pending.
