# Python STABL parity assets for installed vignettes

These files are verbatim copies of the matching files in the repository's canonical, frozen `tests/testthat/fixtures/python_parity/` directory. They are packaged under `inst/extdata` because a CRAN-installed vignette cannot assume access to its source-tree `tests/testthat` directory.

**Source:** `gregbellan/Stabl` at commit `1d07f85a13cfbecb4f08ce21075bf4fbb8e34678`, Python 3.11.9. The generator and the full frozen inputs remain in the canonical test-fixtures directory. Run `sha256sum -c manifest.sha256` from this directory (or the platform's equivalent) to verify this copied subset. The SHA-256 hashes match the canonical fixture manifest as inspected on 18 October 2026.

Do not update the bundled copies independently of the canonical fixtures. If the original algorithm contract changes, version the fixture set, compare it explicitly, and update the vignette's evidence statement. Solver-case rows in `solver_cases.csv` are **minimum gates**, not achieved scores. These fixtures never certify real-data FDP calibration.
