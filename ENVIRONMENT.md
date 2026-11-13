# Isolated technical-test environment

Verified on 2026-11-13. Status: environment dependency and core-import checks PASS; biological analysis is outside this check.

- Server directory: `/home/data/t120631/project/ipf_stage6_smoke_test_20261113`.
- Interpreter: Python 3.10.20, with a virtual environment based on the existing `litreview` interpreter.
- `include-system-site-packages = false`; existing environments were not modified.
- Scanpy 1.11.5, AnnData 0.11.4, SciPy 1.15.3, igraph 1.0.0 imported successfully.
- `pip check`: no broken requirements found; setup exit code 0.
- Full resolved versions are in `environment.freeze.txt`. To reproduce the resolved stack, install that file into an empty compatible Python environment.
- The declared index was PyPI; installation logs also show downloads from the server's existing USTC package mirror. No system package-source configuration was changed.
- Data, caches, logs and temporary files remain in the dedicated test directory.

This environment milestone does not certify clustering accuracy, IPF relevance, spatial analysis, GPU compatibility or large-dataset capacity. The technical pipeline result is recorded separately after completion.
