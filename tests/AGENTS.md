# ZEPHYRUS test instructions

<!-- fwl-tests-core:begin sha256=0d210e346c572a96 -->
## Test rules shared by the PROTEUS ecosystem

Each test file starts with a module-level tier marker and a timeout. CI selects tests by marker, so a file without one runs in no CI job; the timeout stops a hang, it is not a target.

```python
pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]
```

| tier | what it tests | timeout |
|---|---|---|
| `unit` | Python logic with the heavy physics mocked; aim for < 100 ms | 30 s |
| `smoke` | the real binaries or solvers, one step, low resolution | 60 s |
| `integration` | several modules coupled | 300 s |
| `slow` | full physics validation | 3600 s |

Each test carries exactly one tier marker, so that the tier filters select it in the one CI job meant for it; a function marker for a second tier breaks that. Split a file whose tests need different tiers. `skip` excludes a test from every CI job.

Every new test covers an edge case (a boundary value, an empty input, an extreme physical parameter), exercises the error contract (a documented exception with the check that no side effect ran, a guard, a clamp, or the limit input of the formula), and asserts values that do not follow trivially from the implementation. A pinned value of 1 that every exponent reproduces checks nothing.

`python tools/check_test_quality.py --check` fails when the count of any rule rises above `tools/test_quality_baseline.json`; the offenders it prints are the first few of that rule in the tree, not necessarily yours. The rules: a file without a tier marker, a test without a docstring, a test with one assertion or none, a weak assertion as the only one (`is None`, `is not None`, `> 0`, `len(...) > 0`, `isinstance`), `==` next to a non-zero float literal, and an optional dependency imported without `pytest.importorskip`. `bash tools/validate_test_structure.sh` runs the repository's own structure check; the repository part of this file says what it checks.

Physics tests carry markers so their coverage is tracked apart from line coverage:
- `@pytest.mark.physics_invariant` on each test function that asserts a conservation law, a bound (T > 0, fractions in [0, 1]), a monotonicity or symmetry, or a pinned value with a discrimination guard. The marker goes on the function, not the module: structural tests in the same file do not carry it.
- `@pytest.mark.reference_pinned` (together with `physics_invariant`) on a test that pins a published benchmark, an analytical limit or a cross-implementation result; cite the paper, table or figure in the docstring.
- A discrimination guard asserts that the most plausible wrong formula (a missing factor, a swapped exponent, the wrong unit) gives a result outside the tolerance of the pinned value.

Floats: compare with `pytest.approx` or `np.testing.assert_allclose` and a tolerance you can justify from the method, not the one that makes the test pass.

Mocks: mock at the narrowest scope (the one external call), and return physically plausible values, so the code under test runs its real branches. Set random seeds and write files only under `tmp_path`.

A module-level constant read from an environment variable at import time does not change with `monkeypatch.setenv`; patch the constant with `monkeypatch.setattr`.
<!-- fwl-tests-core:end -->

## ZEPHYRUS specifics

Structure: `src/zephyrus/<file>.py` is tested in `tests/test_<file>.py`. The exceptions: `test_mors_coupling.py` (the MORS flux hand-off, MORS mocked, unit tier), `test_earth.py` (a real MORS lookup, integration tier), the Hypothesis sweeps in `test_escape_properties.py` and `test_collision_properties.py`, kept apart so `pytest.importorskip('hypothesis')` skips only them, and `test_nightly_data_cache.py` for `tools/nightly_data_cache.py`. `bash tools/validate_test_structure.sh` checks that every test carries exactly one of `unit`, `smoke`, `integration`, `slow` and `skip` (module, class or function); a tier marker together with `skip` fails.

CI: pull requests run `pytest -m "(unit or smoke) and not skip"` with the fast coverage gate, the structure check and `check_test_quality.py --check`, which blocks here; the nightly runs all tiers.

### Physics sources and invariants

`escape.py` and `collision.py` are the physics sources (`PHYSICS_SOURCES` in `tools/check_test_quality.py`); each has `physics_invariant` and `reference_pinned` tests and a page in `docs/Validation/`. `python tools/check_test_quality.py --reference-pinned-status` lists physics sources without a pinned test.

Invariants for `escape.py`: the rate equals the deposited XUV power over the binding energy per unit mass, up to the geometric and efficiency factors; the rate is non-negative, linear in `Fxuv`, decreasing with `Mp` and zero at `Fxuv = 0`; the rate with the tidal correction exceeds the rate without it for a close-in orbit; `K_tide` is in (0, 1) for `ksi > 1`, and `ksi <= 1` raises `ValueError`.

### Discriminating values

- Use `Rp != Rxuv`: at `Rp == Rxuv` the `scaling=2` and `scaling=3` branches give the same value, so a swapped default passes. Pin both branches and assert they differ.
- Test the tidal branch close in (`a` about 0.02 au), where `K_tide` differs from 1 by tens of percent; at 1 au it is about 0.99 and a dropped correction passes a loose tolerance. Pin `K_tide` and use the no-tidal value as the guard.
- Use two `Fxuv` values and assert the rate ratio equals the flux ratio; use two `Mp` values and assert the rate falls.
- Pin the magnitude with a scale guard (for example `1e6 < rate < 1e8`) and name the unit in the comment: a radius in cm instead of m moves the rate to about `1e12` and fails the band.

### Mocks, constants, seeds

- Unit tests mock MORS at the narrowest scope (`patch('mors.Star')`) and return a plausible `Lx`, `Leuv` pair, not a constant; assert the derived flux against a hand-computed value. `mors` is a runtime dependency: a unit test that needs it mocks it (`test_mors_coupling.py`) rather than skip; `test_nightly_data_cache.py` skips when `mors` or `fwl_io` is missing, because the tool it tests needs both.
- `hypothesis` is the one module-top optional dependency the linter knows (`OPTIONAL_DEPS`); property tests use `@settings(derandomize=True)` or a fixed `--hypothesis-seed`, because the default sequence changes between Hypothesis releases.
- `escape.py` star-imports `G`, so a test that changes `G` patches `zephyrus.escape.G` (the use site), not only `zephyrus.constants.G`.
- Test parameters are SI; parametrize ids name the physical scenario (Earth-like, close-in super-Earth, sub-Neptune).

### Docstrings and names

The file docstring names the source under test and lists the invariants it checks. Each test docstring states the physical scenario or contract clause. A comment explains why an input was chosen ("a = 0.02 au so K_tide differs from 1 by tens of percent"). Names describe behaviour (`test_el_escape_linear_in_xuv_flux`, not `test_el_escape`).

### Coverage gates

Both gates are at 90 %: the fast gate (`[tool.zephyrus.coverage_fast]`, unit and smoke, every pull request) and the full gate (`[tool.coverage.report]`, all tiers, nightly). The PR workflow rejects a change that lowers `[tool.coverage.report].fail_under` below `min(base, 90)`. `tools/update_coverage_threshold.py` raises a threshold toward 90 and never lowers it; it updates the full gate by default and the fast gate with `--target fast`.

A pull request that adds or changes more than 50 lines under `tests/` (`git diff origin/main...HEAD -- tests/`) gets an independent review of its tests before merge.
