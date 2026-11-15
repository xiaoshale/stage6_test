# Native IPF pilot continuation state

Updated: 2026-11-15. **The technical pilot is now complete across all four components.** See [final results](RESULTS.md). The interruption snapshot and once-pending instructions below are retained as historical provenance and are superseded by the completed-run report.

## Verified results

Server root: `/home/data/t120631/project/ipf_stage6_pilot`.

| Component | Retained units | Source features | Actual result |
|---|---:|---:|---|
| GSE135893 | 8,577 source-annotated cells, four donors | 33,694 | Full original MTX read; exact barcode join; source nCount/nFeature agree; H5AD roundtrip passed. Process peak RSS 6,940,749,824 bytes. |
| GSE159354 | 55,189 source-filtered cells across 16 sorted libraries, four donors | 58,302 per library | Every complete H5 parsed and per-library raw H5AD roundtrip passed. Peak RSS 631,738,368 bytes. |
| GSE159585 | 28,970 nuclei, four donors | 27,714 | Full CSV stream through EOF/gzip CRC; global 176,540-cell header/metadata match; selected raw counts and H5AD roundtrip passed. Approximately 275 seconds, peak RSS 2,071,814,144 bytes. |

These are original measurement units, not new independent donor counts or biological discoveries. No QC flag removed cells/genes. Procurement confounding in GSE159585 remains unresolved.

The fixed 17-marker, read-only technical capture audit exited successfully on the server. Its table and summary remain at `runs/pilot_001/capture_markers/`; their attempted local retrieval was blocked by permission-review timeouts. It does not assign new identities or compare disease effects.

The download, Habermann, snRNA and spatial small fixtures passed. A second spatial fixture passed after representing undefined zero-count mitochondrial percentages as missing values. The actual real spatial run has **not** been submitted yet.

## Transfer state and interruptions

- Twenty-four complete source files were declared in `config/acquisition_manifest.json`.
- Initial slow direct acquisition ran in tmux session `ipf-stage6-acquire-20261115`, then was intentionally stopped only after a complete source snRNA file had been downloaded temporarily on the Mac and uploaded through SSH.
- The relayed snRNA file is 509,536,603 bytes, SHA-256 `e6e3fe8d9b984fba7ac779d9948b84a90b659ba82b15cf5f431f234fa893590c`, equal on Mac and server. Its full source/relay record and old direct-transfer record are preserved. The old direct partial is retained and is not an analysis input.
- Resumed acquisition uses tmux `ipf-stage6-acquire-resume-002`, log `logs/acquisition_resume_002.log`. Existing complete files were rehashed; spatial continuation uses the original entity validator/If-Range.
- Last remotely verified spatial ZIP progress was **1,019,215,872 of 1,093,092,912 bytes** at 2026-11-14 23:09:44 UTC (2026-11-15 07:09:44 Asia/Shanghai). The background transfer may subsequently have finished, but this was not reverified.

Automatic permission-approval review repeatedly timed out for SSH/SCP requests. Those rejected-by-timeout requests did not execute. An unprivileged SSH check could not resolve the hostname in the restricted environment. Termius computer-use fallback first returned `cgWindowNotFound`; app inventory then hit the same approval timeout. No authentication failure or source-data failure was established. The user was asked to make the logged-in Termius terminal visible as a possible fallback; this is not a request for a password or private key.

## Exact next actions after access is restored

1. Read `logs/acquisition_resume_002.log` and `data/acquisition.json`; require all 24 complete source records and correct hashes before declaring acquisition complete. Do not restart or duplicate a still-active transfer blindly.
2. Confirm that `runs/pilot_001/Franzen/` does not already exist. Then run the already deployed, fixture-tested spatial reader:

   ```bash
   cd /home/data/t120631/project/ipf_stage6_pilot
   /home/data/t120631/project/ipf_stage6_smoke_test_20261113/.venv/bin/python -u code/run_component.py --root /home/data/t120631/project/ipf_stage6_pilot --component Franzen --run-id pilot_001
   ```

   Use a fresh log with non-clobbering redirection or a dedicated tmux session. The launcher snapshots code, verifies input hashes and applies four threads/48 GiB process virtual-memory limit. Never overwrite a partial/failed run; inspect it first and use another run ID if required.
3. Run the deployed `code/make_qc_figures.py --run-dir /home/data/t120631/project/ipf_stage6_pilot/runs/pilot_001` after all four components pass. This script has only been syntax-checked so far; verify its actual execution. It plots QC-table distributions, not expression effects.
4. Upload `analysis/stage6_pilot/verify_completed_pilot.py` from the Mac; **it is not yet deployed**. Run it with `--root /home/data/t120631/project/ipf_stage6_pilot --run-id pilot_001`. It requires all components, all 24 source files/hashes, whole-ZIP CRC, 24 H5AD sparse checkpoints, executed-code hashes and four QC figures/six overlays before writing `SUCCESS.json`.
5. Retrieve only small summaries/tables/figures, compare with the server artifact hashes, and visually inspect all six real coordinate overlays and four technical QC figures. Leave raw data and H5ADs on the server.
6. Write final `RESULTS.md`, update authoritative project phase, and push the next verified checkpoint to the existing `xiaoshale/stage6_test` repository. Do not claim formal biological analysis readiness: doublet/ambient-RNA strategy, gene harmonization, state coverage, anatomical adjudication and donor independence remain separate gates.

## Local and Git state

- Research root remains on `master`; it has extensive unrelated pre-existing changes and no research-root remote. Do not stage those changes.
- The designated test checkout is `/private/tmp/codex_stage6_test_checkout_20261113`, current `master`, origin `https://github.com/xiaoshale/stage6_test.git`.
- Verified milestone 01 was committed and pushed as `89f61d1`, adding `ipf_stage6_pilot/` with code/configuration, fixtures and the two scRNA results. Milestone 02, including the actual snRNA results and continuation state, was subsequently committed locally as `c4def83`; its push is pending tool access restoration. The exported snRNA TSV retains the server-generated CRLF line endings; the Git whitespace check passes with CR-at-EOL allowed, preserving artifact bytes.
- The first normal HTTPS push lacked a configured credential helper. It succeeded through a command-scoped GitHub CLI helper using the previously authorized xiaoshale keyring login; no persistent remote/configuration changes were made. Do not print tokens or copy authentication files into the repository.
- New scripts/records are English. Frozen Stage 5 evidence and the baseline 96-paper resource catalog remain unchanged.
