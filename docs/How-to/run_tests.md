# Testing suite

This page describes the current ZEPHYRUS test suite and how to run it. The suite currently consists of a single regression test that exercises the end-to-end Earth-like escape calculation, including MORS stellar evolution lookup, XUV flux conversion, and the `EL_escape` call without tidal correction. The tidal branch of `EL_escape`, the `scaling=3` radius option, and error handling on invalid inputs are not yet covered by dedicated tests.

## Running the tests

### Prerequisites

Install ZEPHYRUS with the development dependencies:

```bash
pip install -e .[develop]
```

Download the stellar evolution data (required for the regression test, which uses MORS):

```bash
export FWL_DATA=/path/to/your/data
mors download all
```

### Run the full suite

```bash
coverage run -m pytest
```

### Run a specific file

```bash
pytest tests/test_earth.py
```

### Run with verbose output

```bash
pytest -v
```

### View coverage report

```bash
coverage run -m pytest
coverage report
```

---

## Current test files

### `test_earth.py`: end-to-end EL escape regression

Runs three parametrised regression tests at different stellar ages (150, 5813, and 10020 Myr) against known-good values. For each age, the test:

1. Loads a 1 M$_\odot$ Spada track at the median rotation rate via `mors.Star(Mstar=1.0, Omega=1.0)`.
2. Reads the X-ray and EUV luminosities at the specified age from the MORS track.
3. Computes the XUV flux at Earth's semi-major axis using the inverse-square law and converts to SI units.
4. Calls `EL_escape` with Earth-like parameters ($M_\oplus$, $R_\oplus$, $e_\oplus$, $a_\oplus$), $\epsilon = 0.15$, $R_\mathrm{XUV} = R_p$, and `tidal_contribution=False`.

The four quantities checked against regression values are:

| Quantity | Description | Units |
|---|---|---|
| `Lx` | X-ray luminosity from MORS | erg s$^{-1}$ |
| `Leuv` | EUV luminosity from MORS | erg s$^{-1}$ |
| `Fxuv_star_SI` | XUV flux at Earth's orbit | W m$^{-2}$ |
| `escape` | EL escape mass-loss rate | kg s$^{-1}$ |

Tolerance: `rtol=1e-5`, `atol=0`. The test calls `mors.DownloadEvolutionTracks('Spada')` to ensure data are present before running.

---

## CI matrix

Tests run automatically on every push and pull request to `main` via GitHub Actions, as well as on manual `workflow_dispatch`:

| OS | Python versions |
|---|---|
| Ubuntu (latest) | 3.10, 3.11, 3.12 |
| macOS (latest) | 3.10, 3.11, 3.12 |

The workflow caches the Python virtual environment (keyed on `pyproject.toml` hash) and the `FWL_DATA` directory across runs to speed up CI.

Coverage is reported in the GitHub Actions summary for every run. A coverage badge is generated from the `ubuntu-latest` + Python 3.10 run on `main`.