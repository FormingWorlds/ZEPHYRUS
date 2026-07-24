# Testing suite

[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/ZEPHYRUS?label=coverage&logo=codecov&color=brightgreen)](https://app.codecov.io/gh/FormingWorlds/ZEPHYRUS)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/tests.yaml?branch=main&label=Unit%20Tests&color=brightgreen)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/nightly.yml?branch=main&label=Integration%20Tests&color=brightgreen)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/nightly.yml)
[![tests](https://img.shields.io/endpoint?url=https://proteus-framework.org/ZEPHYRUS/badges/tests-total.json)](https://proteus-framework.org/testing)

Tests verify that the code does what was written; physical correctness is judged by data, not by tests. The suite catches regressions in the energy-limited mass-loss formula, the tidal correction, the radius-scaling selection, and the MORS stellar-flux coupling, but it cannot certify that those formulae match nature. That judgement belongs to the validation runs and the published comparisons the model is built on.

This page describes the suite as a whole and how to run it.

## Test quality contract

Five layers enforce test rigor across the suite:

1. A four-marker tier scheme (`unit`, `smoke`, `integration`, `slow`) selects what runs in the PR gate versus the nightly.
2. Two validation markers (`physics_invariant`, `reference_pinned`) tag tests that carry physical meaning beyond pure code coverage.
3. A source-to-test mirroring rule pairs the physics source with a same-named test file.
4. An AST linter (`tools/check_test_quality.py`) rejects weak-test patterns on every PR.
5. A coverage gate held at the 90 % ceiling, enforced on every PR; a helper (`tools/update_coverage_threshold.py`, run manually) raises `fail_under` one-way toward that ceiling and never lets it drop.

Layers 1, 4, and 5 are blocking on PRs. Layers 2 and 3 are advisory: the linter reports gaps but does not fail the build.

## The four-marker tier scheme

Every test in the suite carries a tier marker, applied at module level (`pytestmark = [pytest.mark.X, pytest.mark.timeout(...)]`).

| Marker | What it tests | Per-test budget | CI surface |
|---|---|---|---|
| `unit` | Python logic: the closed-form mass-loss rate, the tidal factor, the radius-scaling branches, the MORS-flux hand-off with the stellar lookup mocked. No real MORS download. | < 100 ms | PR + nightly |
| `smoke` | Real dependency on a minimal input. | < 30 s | PR + nightly |
| `integration` | The real MORS-coupled escape: downloads the stellar-evolution tracks and drives the escape formula end to end. | minutes | Nightly only |
| `slow` | Long parameter sweeps and convergence studies. | up to an hour | Nightly only |
| `skip` | Placeholder, deliberately disabled. | n/a | Never |

Tests without a tier marker are invisible to CI. The PR gate runs `pytest -m "(unit or smoke) and not skip"`; the nightly runs everything in `(unit or smoke or integration or slow) and not skip`.

## Module-level marker and timeout

Every test file declares its tier and its wall-time ceiling at module top:

```python
pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| Tier | Timeout |
|---|---|
| `unit` | 30 s |
| `smoke` | 60 s |
| `integration` | 300 s |
| `slow` | 3600 s |

The timeout is a defensive ceiling, not a target. A unit test that takes 25 s of wall time has either picked the wrong tier or has a leak somewhere; the ceiling catches future regressions that introduce a hang. Per-function markers (for example `@pytest.mark.parametrize`) are additive and do not replace the module-level marker.

## Physics-invariant tiering

A unit test on a physics source (`escape.py`, `collision.py`) must assert at least one of four invariant families:

- **Conservation**: the mass-loss rate equals the deposited XUV power divided by the gravitational binding energy per unit mass, up to the geometric and efficiency factors of the formula; in a coupled step the removed mass never exceeds the available atmospheric mass.
- **Positivity or boundedness**: the escape rate is non-negative for valid inputs; the tidal factor `K_tide` lies in `(0, 1)` when the Hill radius exceeds the XUV radius (`ksi = Rhill / Rxuv > 1`), rising toward 1 as the orbit widens, and the tidal branch raises `ValueError` for `ksi <= 1`; radii, masses, and semi-major axis are strictly positive; the XUV flux and efficiency are non-negative.
- **Monotonicity or symmetry**: the rate is linear in the XUV flux at fixed geometry; it scales as `1 / Mp`; it is larger with the tidal correction than without (`K_tide < 1`) for a close-in orbit; it is zero when the flux is zero.
- **Pinned numeric value with a discrimination guard**: a closed-form value pinned via `pytest.approx`, plus explicit assertions that a wrong exponent, a wrong scaling branch, a sign flip, or a unit slip would each differ from the correct value by more than the tolerance.

For the collision law the same families read: the loss fraction is bounded in `[0, 1]` and capped at 1 for total erosion; it vanishes at grazing incidence and at zero contact speed; it never decreases with contact speed; and its closed-form pins carry guards against the wrong mass-ratio denominator, wrong exponents, and a wrong gravitational constant.

Tests that meet one or more of these are tagged `@pytest.mark.physics_invariant`. The marker is per-function, not module-level: an error-contract test (an unsupported `scaling` raising `ValueError`) does not carry it.

The physics sources in ZEPHYRUS are `escape.py` and `collision.py`. The utility sources (`__init__.py`, `constants.py`, `planets_parameters.py`) are exempt from the invariant requirement but remain subject to the anti-happy-path rules.

## Reference-pinned validation

Tests that pin behaviour against an external anchor are tagged `@pytest.mark.reference_pinned`. The anchor is one of a published benchmark, an analytical limit, or a cross-implementation check. Each physics source carries at least one such test, recorded on the [`escape`](../Validation/escape.md) and [`collision`](../Validation/collision.md) validation pages. The current anchors:

| Source | Anchor | Test |
|---|---|---|
| `escape.py` | Erkaev et al. (2007), A&A 472:329, Eq. 21: closed-form energy-limited rate for the default `scaling=2` radius term | `tests/test_escape.py::test_el_escape_scaling2_matches_erkaev2007_closed_form` |
| `escape.py` | Lehmer & Catling (2017), ApJ 845:130, Eq. 1: closed-form rate for the `scaling=3` radius term | `tests/test_escape.py::test_el_escape_scaling3_matches_lehmer_catling_closed_form` |
| `collision.py` | Kegerreis et al. (2020), ApJL 901:L31, Eq. 1: closed-form erosion fraction for identical twin bodies | `tests/test_collision.py::test_scaling_law_pins_the_kegerreis_closed_form` |
| `collision.py` | Kegerreis et al. (2020), ApJL 901:L31, Tables 1 and 2: simulated loss fractions of the SPH suite | `tests/test_collision.py::test_scaling_law_reproduces_kegerreis_table2_simulations` |

The marker is not the same thing as physical correctness: a reference-pinned test certifies that this implementation reproduces that anchor; it does not certify that the anchor is the right physics for every planetary regime.

## Source-to-test mirroring

Each source module with executable content has a same-named companion in `tests/`:

| Source | Test |
|---|---|
| `src/zephyrus/escape.py` | `tests/test_escape.py` |
| `src/zephyrus/collision.py` | `tests/test_collision.py` |
| `src/zephyrus/constants.py` | `tests/test_constants.py` |
| `src/zephyrus/planets_parameters.py` | `tests/test_planets_parameters.py` |

Cross-cutting and companion tests are the documented exception, not the rule:

- `tests/test_mors_coupling.py`: the MORS-to-escape flux hand-off with the stellar lookup mocked, so the coupling recipe runs in the fast unit tier without a download.
- `tests/test_earth.py`: an Earth-analogue regression that spans the real MORS lookup and the escape formula end to end. It carries the `integration` tier because it downloads the stellar-evolution tracks.
- `tests/test_escape_properties.py`: the Hypothesis-driven property sweeps for `escape.py`, kept in their own module so the `pytest.importorskip('hypothesis')` skip applies only to these tests and the closed-form pins in `tests/test_escape.py` still run when the develop-extra dependency is absent.
- `tests/test_collision_properties.py`: the same pattern for `collision.py`: boundedness, speed monotonicity, angle monotonicity at equal densities, and the exact reduction of the interacting mass to the interacting volume.

`__init__.py` holds only the package version string and has no dedicated test file.

## AST test-quality linter

`tools/check_test_quality.py` walks `tests/test_*.py` as an AST and enforces these rules:

| Rule | What it flags |
|---|---|
| `missing_module_pytestmark` | Test file with no module-level `pytestmark` (a tier marker is required). |
| `missing_docstring` | Test function with no docstring. |
| `single_assert` | Function with exactly one assertion (a single assert is rarely enough to discriminate the correct formula from plausible wrong ones). |
| `no_assertions` | Function with zero assertions (only valid for a test that exercises an exception path with `pytest.raises`). |
| `weak_assert_*` | Standalone `assert result is not None`, `assert result > 0`, `assert len(result) > 0`, or `assert isinstance(...)` as the sole meaningful check. A weak assertion alongside a strong primary one (the sign guard in a discrimination pattern, for example) is not flagged. |
| `float_eq_literal` | `==` adjacent to a non-zero float literal in a test body (use `pytest.approx`). |
| `missing_importorskip` | An optional dependency (`hypothesis`) imported at module top without a preceding `pytest.importorskip('<name>')`. The PR image uses `pip install --no-deps`; without `importorskip`, collection fails. |

The linter runs in two blocking modes:

- `python tools/check_test_quality.py --baseline` walks the suite and writes the per-rule violation counts to `tools/test_quality_baseline.json`. This is the floor. Regenerate it only after a deliberate reduction in violations.
- `python tools/check_test_quality.py --check` (CI mode) walks the suite, compares the current counts to the baseline, and exits non-zero if any rule's count exceeds the baseline.

The baseline ratchets one way. Override with `ZEPHYRUS_TEST_QUALITY_ALLOW_REGRESS=1` only when a new rule was added that surfaces pre-existing violations.

Two advisory modes report gaps without failing CI:

- `python tools/check_test_quality.py --reference-pinned-status`: lists physics sources whose test file has no `@pytest.mark.reference_pinned` test.
- `python tools/check_test_quality.py --physics-invariant-status`: lists physics-source tests that assert no invariant and are not tagged `@pytest.mark.physics_invariant`.

## Local commands

```console
pytest -m unit                              # fast unit tests
pytest -m smoke                             # minimal-input dependency tests
pytest -m integration                       # real MORS-coupled escape
pytest -m slow                              # long-running validation (nightly tier)
pytest -m "(unit or smoke) and not skip"    # PR-gate selection
pytest -m "not skip"                        # everything that should ever run
```

The integration tier needs the stellar-evolution data downloaded once:

```console
export FWL_DATA=/path/to/your/data
mors download all
```

Coverage:

```console
pytest --cov=zephyrus --cov-report=term -m "not skip"
pytest --cov=zephyrus --cov-report=html -m "not skip"   # htmlcov/
```

Lint and structure:

```console
bash tools/validate_test_structure.sh         # module-level marker validator
python tools/check_test_quality.py --check     # AST linter against baseline
python tools/check_test_quality.py --reference-pinned-status
python tools/check_test_quality.py --physics-invariant-status
```

## Public-facing badges versus internal taxonomy

Public-facing badges (README, project website) collapse `smoke + integration + slow` into a single `Integration Tests` category, because a four-way taxonomy is confusing to non-developer readers. The four-marker internal scheme remains for CI granularity: the PR gate runs `(unit or smoke)`, the nightly runs everything, and the test-count badges fetch the JSON files written into the documentation site during the docs deploy.

## Badge system

The documentation deploy (`.github/workflows/docs.yaml`) regenerates three JSON files, `tests-{total,unit,integration}.json`, from `pytest --collect-only` and writes them into the published site under `badges/`. Shields.io fetches them live via the endpoint URL embedded in the test-count badge. The counts refresh on every documentation build, so they track the suite without running it.

## Coverage gates

Two gates are declared in `pyproject.toml`:

| Gate | Tests included | Threshold key | Where it runs |
|---|---|---|---|
| Fast | `unit + smoke` | `[tool.zephyrus.coverage_fast].fail_under` | PR `Run unit + smoke tests` step |
| Full | `unit + smoke + integration + slow` | `[tool.coverage.report].fail_under` | Nightly `Run full suite with coverage` step |

Both gates sit at the 90 % ceiling. The `tools/update_coverage_threshold.py` helper (run manually) raises `fail_under` one-way toward the ceiling and never lets it drop; `ECOSYSTEM_CEILING = 90.0` caps it, and neither gate may be manually decreased. ZEPHYRUS has no compiled dependencies, so the unit and smoke tier already covers the whole package, which is why both gates are already at the ceiling. The PR gate has a pre-flight step that fetches the base branch's `pyproject.toml` and rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.

## PR validation pipeline

`.github/workflows/tests.yaml` runs on every push and pull request to `main`, and on manual `workflow_dispatch`. Draft PRs run only `ubuntu-latest` with Python 3.12; non-draft events run the full matrix (`ubuntu-latest`, `macos-latest` x Python 3.10, 3.11, 3.12). The step sequence:

1. **Validate test structure** (`bash tools/validate_test_structure.sh`): rejects any test file without a module-level `pytestmark`.
2. **Run test-quality lint** (`python tools/check_test_quality.py --check`): blocking; rejects regression against `tools/test_quality_baseline.json`.
3. **Pre-flight fail_under ratchet check**: rejects any PR that drops `[tool.coverage.report].fail_under` below `min(base, 90.0)`.
4. **Run unit + smoke tests**: `pytest -m "(unit or smoke) and not skip" --cov=zephyrus --cov-fail-under=${FAST_FAIL_UNDER}`, where `FAST_FAIL_UNDER` is read live from `[tool.zephyrus.coverage_fast].fail_under`.

Steps 1 to 3 are gated to `ubuntu-latest, python 3.12` to avoid redundant runs; step 4 runs across the full matrix.

Nightly (`.github/workflows/nightly.yml`) runs the full suite, uploads coverage to Codecov, and enforces `--cov-fail-under=90`.

## Canonical specification

The repository-wide rules that every PROTEUS-ecosystem submodule follows are at [proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/](https://proteus-framework.org/PROTEUS/Explanations/ecosystem_testing_standard/).
