# ZEPHYRUS Test Quality Rules

This file is the canonical deep-dive on test quality. The high-level summary lives in [`.github/copilot-instructions.md`](../../copilot-instructions.md) under "Testing Standards". The two files MUST stay in sync. If you change one, mirror the change in the other.

> **Discovery note.** ZEPHYRUS keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `zephyrus-code-review.md` explicitly so AI tooling and human readers know to load them. **When opening or editing any file under `tests/**` or `src/zephyrus/**`, read this file first.**

Sister rule files:

- [`.github/copilot-instructions.md`](../../copilot-instructions.md): high-level rules, applied repo-wide.
- [`.github/.claude/rules/zephyrus-code-review.md`](zephyrus-code-review.md): review-pass gate, domain-aware code review (tidal-correction propagation, radius-scaling exponent safety, XUV-flux coupling, PROTEUS-coupling patterns). Test-marker discipline lives there too.

ZEPHYRUS is scientific simulation code and the test suite is held to physics-grade rigor. Tests exist to catch real bugs. A test that asserts the wrong thing, or that passes for the wrong reason, is worse than no test because it generates false confidence. The rules below codify what "real test" means here.

---

## 1. Anti-happy-path rules (every new test)

Every new test function MUST include:

1. **At least one edge case**: a boundary value (`Fxuv = 0`, `e = 0`, `epsilon = 0`, `tidal_contribution` toggled, a close-in orbit where `K_tide` departs strongly from 1), an empty input, or an extreme physical parameter.
2. **At least one path that exercises the error contract**:
   - If the function under test has documented validation (`EL_escape` raises `ValueError` on an unsupported `scaling`, and on `ksi = Rhill/Rxuv <= 1` in the tidal branch), test that the error fires AND that no result is returned.
   - If the function has no validation (closed-form mathematics: the mass-loss formula), exercise the **limit-input behavior** (`Fxuv = 0` gives zero escape; `tidal_contribution=False` gives `K_tide = 1`) and assert the corresponding mathematical invariant.
   - "No validation in source therefore no error test" is not an exemption; the limit-input substitute is.
3. **Assertion values NOT trivially derivable from the implementation**: discriminating numeric pins (see Section 2 below) or property-based assertions (monotonicity, linearity, boundedness).

### Forbidden patterns

These are flagged by `tools/check_test_quality.py` and rejected at PR time.

- **Single-assert test functions**. Two or more assertions per test; the second usually pins the invariant the first hand-waves over. Exception: a single assertion of a hard-fail invariant is acceptable if the test is the only test of that invariant in the file.
- **Weak assertions when they stand alone as the sole meaningful check in the test.** The shapes are:
  - `assert result is not None`
  - `assert result > 0`
  - `assert len(result) > 0`
  - `assert isinstance(result, dict)`
  - `assert result is None` where the function returns `None` implicitly

  Required carve-out: the three-class discrimination guard (Section 2) uses `assert val > 0` as the sign-error guard and `assert lo < val < hi` as the scale-error guard alongside a primary `pytest.approx(...)` pin. Those secondary lines look like weak assertions in isolation; they are NOT flagged when paired with a stronger primary assertion in the same test. The linter applies the carve-out automatically: weak shapes are flagged only when the test has exactly one `assert` statement and that assertion is itself the weak shape.
- **Tests with no function-level docstring**. The docstring states which physical scenario or contract clause is being verified.
- **`==` adjacent to a float literal**. Use `pytest.approx(val, rel=...)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`. Comparing two floats with `==` is a known flake source even for "exact" identities like 0.0 (-0.0 vs +0.0, NaN propagation).
- **Tests asserting on a fixture's implicit default**: e.g. `assert fixture_returning_none() is None`. This is trivially true. Delete the test; do not strengthen it by adding more `is None` assertions.

---

## 2. Discriminating test values

The test contract is: a regression that introduces a plausible bug must fail the test. "Plausible bug" means off-by-one exponent, wrong sign, swapped factor of 2, missing factor of pi, dimensionally-wrong unit, **wrong radius-scaling selection**, **dropped tidal correction**. Pick input values where the wrong-formula result is far from the correct one.

### Bad / good examples

| Pattern | Bad (any-formula-passes) | Good (discriminates) |
|---|---|---|
| EL mass-loss rate | `Rp = Rxuv` with `scaling=2` (degenerate: `Rp*Rxuv**2 == Rxuv**3`, so `scaling=2` and `scaling=3` are indistinguishable) | `Rp != Rxuv` so the `scaling=2` and `scaling=3` branches give different values; pin both and assert they differ beyond tolerance |
| Tidal correction | `a` large (1 au): `K_tide ~ 0.99`, indistinguishable from 1 | Close-in orbit (`a ~ 0.02 au`) so `K_tide` departs from 1 by tens of percent and a dropped correction fails the test |
| XUV-flux linearity | Single `Fxuv` value | Two `Fxuv` values; assert the escape rate ratio equals the flux ratio (linearity), so a spurious offset term fails |
| Mass dependence | Single `Mp` | Two `Mp` values; assert escape decreases with mass (the `1/Mp` dependence), so a sign flip fails |

### Discrimination guard (REQUIRED for pinned-value tests)

When a test pins a numeric value, include explicit assertions that the wrong-formula result would differ from the correct one for **each plausible bug class**. "Each plausible bug class" means at minimum:

1. **Exponent or factor error** (off-by-one exponent, missing factor of 2 / pi). `abs(val - wrong_value)` discriminates.
2. **Sign error** (`-x` vs `+x`). `abs()` hides this; assert the sign explicitly with `val > 0`.
3. **Unit-conversion error** (m vs cm, kg vs g, erg vs W, au vs m). Pin the absolute scale with the unit named in the comment.
4. **Wrong radius-scaling / dropped-tidal selection**. When the function dispatches by `scaling`, the discrimination guard MUST include a value that distinguishes the chosen path from the sibling path. When the test exercises the tidal branch, it MUST include the no-tidal value so a dropped `K_tide` fails.

**Carve-out for conservation-style invariants.** When the primary assertion IS a conservation closure or an exact linear relation (`escape(2*F) == pytest.approx(2*escape(F))`), the relation already discriminates exponent / factor errors by construction. The exponent guard is satisfied by the relation itself; sign and scale guards remain mandatory.

Canonical pattern:

```python
def test_el_escape_scaling2_matches_closed_form():
    """Pin the EL mass-loss rate for the scaling=2 (Erkaev 2007) branch."""
    val = EL_escape(False, a, e, Mp, Ms, epsilon, Rp, Rxuv, Fxuv, scaling=2)
    expected = 4.4025157e7  # kg s-1, hand-evaluated from epsilon*pi*Rp*Rxuv**2*Fxuv/(G*Mp)
    assert val == pytest.approx(expected, rel=1e-6)
    # Wrong-scaling discrimination: scaling=3 uses Rxuv**3 instead of Rp*Rxuv**2.
    # With Rp != Rxuv the two branches differ well beyond tolerance.
    wrong_scaling3 = epsilon * np.pi * Rxuv**3 * Fxuv / (G * Mp)
    assert abs(val - wrong_scaling3) > 0.01 * expected
    # Sign guard: an energy-limited escape rate is always non-negative.
    assert val > 0
    # Scale guard: order of magnitude is 1e7 kg s-1, not 1e0 (dropped pi/epsilon)
    # or 1e14 (cm-vs-m radius slip). Pin the magnitude.
    assert 1e6 < val < 1e8
```

The guard lines are mandatory whenever the test's primary assertion is a `pytest.approx` against a hand-calculated or published value. Property-based assertions (monotonicity, linearity, boundedness) do not need a separate guard because they are already discriminating across the input space.

---

## 3. Physics-invariant assertions (tiered)

### When required

Every unit test on a **physics source** must assert at least one of the four invariants below. The physics source is:

```
src/zephyrus/escape.py
```

Per-source-file granularity: `escape.py` needs at least one `@pytest.mark.physics_invariant` test and at least one `@pytest.mark.reference_pinned` test in `tests/test_escape.py`. Granularity is per source file, not per directory.

Utility sources are exempt from the physics-invariant requirement but still subject to all anti-happy-path rules:

```
src/zephyrus/__init__.py               (version string)
src/zephyrus/constants.py              (pure physical constants, no derivation)
src/zephyrus/planets_parameters.py     (tabulated star-planet parameters)
```

### The four invariant families

1. **Conservation**
   - Energy-limited budget: the mass-loss rate equals the deposited XUV power divided by the gravitational binding energy per unit mass, up to the geometric and efficiency factors of the formula.
   - Escape-mass budget: in a coupled step the removed mass never exceeds the available atmospheric mass.
2. **Positivity / boundedness**
   - Escape rate non-negative for physically valid inputs.
   - `K_tide` in `(0, 1)` for `ksi = Rhill/Rxuv > 1`, rising toward 1 as the orbit widens; the tidal branch raises `ValueError` for `ksi <= 1`.
   - Radii, masses, semi-major axis strictly positive; XUV flux and `epsilon` non-negative.
3. **Monotonicity or symmetry**
   - Escape rate linear in `Fxuv` at fixed geometry (doubling the flux doubles the rate).
   - Escape rate decreasing with planet mass (`1/Mp`).
   - Escape rate strictly larger with the tidal correction than without (`K_tide < 1`) for a bound close-in orbit.
   - Escape rate zero when `Fxuv = 0`.
4. **Pinned numeric value with a discrimination guard**: see Section 2. This is acceptable as the sole invariant when a closed-form result or published table value is the contract.

Property-based assertions (monotonicity, linearity, boundedness) are preferred over point-value pins when both are possible. They hold for any valid input and so catch bugs across the entire input space.

### Validation certification markers

Two markers track validation quality independently of line coverage:

- **`@pytest.mark.physics_invariant`** -- this test asserts at least one of the four invariants. Tag every qualifying test in a physics-source test file. `tools/check_test_quality.py` warns when a physics-source test asserts no invariant and is not tagged.
- **`@pytest.mark.reference_pinned`** -- this test pins behavior against a **published benchmark** (paper, figure, table; cite explicitly in the test docstring), an **analytical limit** (`K_tide -> 1` as `a -> infinity`; escape linear in `Fxuv`; the closed-form EL rate at a fixed geometry), or a **cross-implementation cross-check**.
  - **Per-source-file**: `escape.py` must have at least one `reference_pinned` test in `tests/test_escape.py`. Anchor type is one of {published benchmark, analytical limit, cross-implementation cross-check}; the specific paper or limit is chosen by the test author and recorded in `docs/Validation/escape.md`.
  - **Tracking**: each physics source gets a page at `docs/Validation/<file>.md`, created when the first reference_pinned test for that source lands. The page records: the source under test, the reference cited, the test ids carrying the marker, and the date of last comparison against the source.
  - **Status report**: `python tools/check_test_quality.py --reference-pinned-status` reports the physics source files missing a `reference_pinned` test. This is the punch list for follow-up validation work.

Both markers are registered in `pyproject.toml` under `[tool.pytest.ini_options] markers`. They do not gate CI on their own; their coverage is a separate KPI surfaced in the PR summary comment.

---

## 4. Mocking discipline

- Default to `unittest.mock` for ALL external calls in unit tests: MORS lookups, file I/O, HTTP, subprocess.
- Mock at the **narrowest scope**: patch the specific function (`unittest.mock.patch('mors.Star')`), not the whole module.
- A mocked dependency MUST return **physically plausible** values. A MORS mock that returns `0.0` or `1.0` for every luminosity will mask sign / scale / dilution bugs; return a plausible `Lx`, `Leuv` pair instead.
- NEVER mock the function under test. If you're tempted to, the test is asking the wrong question.
- Smoke tests exercise the real code path with minimal input; integration and slow tests use the real MORS lookup. The rules in this file still apply to those tiers, but the mocking discipline is relaxed because the real call is the contract.

---

## 5. Optional-dependency imports

Any test that imports an optional dependency MUST call `pytest.importorskip` at module top so `pip install --no-deps` CI runs do not fail collection:

```python
import pytest

hypothesis = pytest.importorskip('hypothesis')
```

Optional deps recognized by the linter (`OPTIONAL_DEPS` constant in `tools/check_test_quality.py`):

- `hypothesis` (used in property-based / fuzz tests; lives in `[develop]` extras).

`mors` is a hard runtime dependency (`fwl-mors` in `[project.dependencies]`), so it is imported unconditionally and is NOT in `OPTIONAL_DEPS`. A unit test that would otherwise reach a real MORS lookup should mock it rather than skip on it.

The lint script enforces this. Rule key `missing_importorskip`: any module-top `import <optional_dep>` or `from <optional_dep> import ...` that is not preceded by a module-scope `pytest.importorskip('<optional_dep>')` is flagged.

---

## 6. Module-level constants and `monkeypatch`

When the source under test reads an env var or a class-level default into a **module-level constant** at import time, `monkeypatch.setenv` alone is not sufficient. The constant is frozen at the original import.

`escape.py` imports its constants via `from zephyrus.constants import *`, so `G` is bound into the `escape` namespace at import. To override `G` in a test, patch it where it is used, not only where it is defined:

```python
# Test (wrong):
monkeypatch.setattr('zephyrus.constants.G', 1.0)   # escape.py already bound the old value

# Test (right):
monkeypatch.setattr('zephyrus.escape.G', 1.0, raising=False)
```

When in doubt, patch both the definition site and the use site. The lint script does NOT currently flag this pattern (it would require source-side analysis); this is a discipline rule enforced via the >50 LOC review trigger and the recurring-trap table in Section 16.

---

## 7. Marker discipline and timeouts

### Module-level marker is mandatory

Every test file MUST begin with:

```python
import pytest

pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]
```

with budgets:

- `unit` -> `timeout(30)` (target wall-time per test is `< 100 ms`; the 30 s cap is a defensive net).
- `smoke` -> `timeout(60)` (target `< 30 s`).
- `integration` -> `timeout(300)`.
- `slow` -> `timeout(3600)`.

PR CI runs `pytest -m "(unit or smoke) and not skip"`. Tests without the tier marker are invisible to CI and shipped untested. The lint script blocks any file missing the module-level `pytestmark`.

### Per-function markers

Per-function `@pytest.mark.<tier>` markers are **additive**, not a replacement for the module-level marker. They are useful when a file's tests span multiple tiers (rare; prefer separate files).

### Timeout is a safety net, not a target

The `timeout` ceiling exists so a future regression that introduces a hang (real MORS lookup, infinite loop, network retry) surfaces as a specific-test failure rather than a generic job timeout. Current unit test wall times are 100x below the ceiling; if you find yourself needing the full 30 s for a unit test, something has gone wrong and you should reduce scope or move the test to a slower tier.

---

## 8. Float and numerical comparison

- NEVER use `==` for floats. Use `pytest.approx(val, rel=1e-5)` or `np.testing.assert_allclose(actual, expected, rtol=..., atol=...)`.
- State the tolerance rationale in a comment when the choice is non-obvious. E.g. "`rtol=1e-6` because the reference literal is hand-evaluated to seven significant figures".
- For pinned numeric values, include a **discrimination guard** (Section 2).
- For property-based assertions, use `pytest.approx` against the exact symbolic relation, with the tightest tolerance the implementation can hit (typically `rel=1e-12` for closed-form algebra; looser for MORS-coupled outputs).

---

## 9. Voice rule for test artifacts

The repo-wide voice rule (zero AI-process disclosure in any public artifact) applies to test code with the same strictness as to source. The voice rule is **scoped** to public artifacts other contributors and external readers see; it does NOT apply to the rule documents themselves (this file, `zephyrus-code-review.md`, `copilot-instructions.md`), which legitimately name the procedures they define.

In scope (the voice rule is BANNED here):

- Test-skip reasons (`@pytest.mark.skip(reason='...')`).
- Test-file docstrings.
- Test-function and test-class names.
- Test-function docstrings.
- Parametrize ids (`@pytest.mark.parametrize('name', [...], ids=[...])`).
- Log-capture assertions (regex against `caplog.records`).
- Commit messages on test-touching commits (subject AND body).
- **Pull-request titles and bodies on test-touching PRs**.
- GitHub Actions job names and step names that ship to the PR Checks tab.
- Inline source comments and docstrings on `src/zephyrus/**`.
- Log strings that ship with the repo.
- **All public-facing documentation** (anything under `docs/`, the repo README, CONTRIBUTING.md, tutorials, wiki pages). Public docs apply the rule silently; they do NOT enumerate the banned phrases, name the voice rule, advertise the existence of `.github/.claude/` rule infrastructure, or cross-reference `.github/.claude/rules/*.md` files. A user docs page that describes the testing infrastructure must do so without naming the AI-process rules that produced it.

Out of scope (these may NAME the procedures they define):

- This file (`zephyrus-tests.md`).
- `zephyrus-code-review.md`.
- `copilot-instructions.md`.

Banned phrases inside the in-scope artifacts: "audit", "review pass", "adversarial review", "Phase X" (when "X" is an AI-organized roadmap label, not a real project phase), "T1.x", "Group A/B/C/D" (when AI-organized work groups), `claude-config/...` paths, "Generated with Claude", AI-tool names, em-dashes, en-dashes (except in bibliographic page ranges within citations), process meta-commentary ("after careful analysis").

Write the OUTCOME (what the test verifies; what the PR achieves) never the PROCESS (how the rule was derived; which review caught what). First-person Tim voice. Going-forward only, no history rewrite.

---

## 10. Fixture and parameter conventions

- Use SI units in test parameters unless the function under test explicitly expects other units. `EL_escape` is SI throughout (m, kg, W m-2, s).
- Use `@pytest.mark.parametrize` when the same logic spans multiple physical regimes (Earth-like, close-in super-Earth, sub-Neptune). Each parametrize id must read like a physical scenario, not a tuple of numbers.
- Set seeds for any randomness:
  ```python
  np.random.seed(42)
  random.seed(42)
  ```
  Hypothesis tests use `@settings(derandomize=True)` or an explicit `--hypothesis-seed` to keep replays stable across versions (see Section 16 trap).
- Use `tmp_path` (pytest fixture) for temporary files. Do not produce large outputs in the test path.

---

## 11. Documentation per test

- **File-level docstring**: name the source file under test (`Tests for src/zephyrus/<file>.py`), list the invariants and contract clauses the file exercises, link to `docs/How-to/run_tests.md`. Required.
- **Function-level docstring**: state the physical scenario or contract clause in plain language. Required (lint-enforced).
- **Inline comments**: explain **why** a specific input range was chosen ("`a = 0.02 au` so `K_tide` departs from 1 by tens of percent, well above tolerance").

---

## 12. Naming

- Test names describe behavior, not the called function: `test_el_escape_linear_in_xuv_flux`, NOT `test_el_escape`.
- Test names use snake_case and read as full sentences.
- Group related tests in classes (`class TestTidalCorrection:`) when they share setup; use the class to thread a single fixture through several scenarios.
- Test file names mirror source 1:1: `src/zephyrus/<file>.py` -> `tests/test_<file>.py`. Documented exceptions to the 1:1 rule:
  - **Cross-cutting coupling tests** (`test_mors_coupling.py`, `test_earth.py`): regressions that span the MORS flux hand-off and the escape formula rather than a single source file. `test_mors_coupling.py` mocks the stellar-luminosity lookup so the coupling recipe runs in the fast `unit` tier without a download; `test_earth.py` performs a real MORS lookup end to end and carries the `integration` tier.
  - **Property-based companion** (`test_escape_properties.py`): the Hypothesis-driven sweeps for `escape.py`. They sit in their own module so the `pytest.importorskip('hypothesis')` gate covers only the property tests, leaving the closed-form pins and error-contract guards in `test_escape.py` running when the develop-extra dependency is absent.

---

## 13. Adversarial review trigger

A pull request that adds or substantially modifies **> 50 lines of test code across all its commits** triggers an independent review pass before merge. This is a discipline rule, not CI-automated: the author runs the review pass via a `code-reviewer` agent before pushing the final test-touching commit. The denominator is PR-level, not per-commit: `git diff origin/main...HEAD -- 'tests/**'` is the source of truth. Splitting one large change into 49 + 49 + 49 line commits does NOT dodge the trigger.

The reviewer's mandate:

- Cite the anti-happy-path rule (Section 1) and the discrimination-guard requirement (Section 2).
- Flag single-assert tests, weak `is not None` patterns, missing module-level marker, missing `physics_invariant` tag on a physics-source test, missing `reference_pinned` tag on a per-source benchmark test, dead tests (passes for the wrong reason), tests that mock the function under test.
- Verify discriminating values: re-derive the expected value from a plausible wrong formula (wrong scaling exponent, dropped tidal correction) and assert the test fails with that wrong formula.
- Verify physics-source coverage of the four invariants: which of the four does this test exercise? If none, why is the test in `tests/test_escape.py`?

The reviewer is a separate process from the test author. For Claude-Code workflow this means spawning a review skill or a `code-reviewer` agent with the test files in scope; the review must complete and surface findings before the test commit is pushed.

The reviewer's findings are addressed in a follow-up commit (not amended into the test commit). The follow-up subject line is in plain language describing the OUTCOME ("sharpen the tidal-branch assertion to distinguish a dropped K_tide", NOT "address review findings").

---

## 14. Tooling

The repo provides:

- `bash tools/validate_test_structure.sh` -- structural check (marker presence, file naming).
- `python tools/check_test_quality.py --check` -- CI mode: AST scan for the forbidden patterns in Section 1 and the marker requirement in Section 7. Fails the PR if violations exceed the baseline.
- `python tools/check_test_quality.py --baseline` -- after a deliberate sweep, regenerates `tools/test_quality_baseline.json`. Only run when you have intentionally reduced violations.
- `python tools/check_test_quality.py --reference-pinned-status` -- prints physics source files missing a `reference_pinned` test.
- `python tools/update_coverage_threshold.py` -- ratchet the fast PR gate upward when measured coverage exceeds the current `fail_under`. Capped at the 90% ecosystem ceiling.
- `ruff check src/ tests/ tools/` and `ruff format src/ tests/ tools/` -- run before commit.

The lint script is wired into PR CI (`tests.yaml`). The step runs in **blocking** mode: any regression above the baseline fails the PR.

---

## 15. Coverage strategy (operator's view)

ZEPHYRUS uses two coverage gates with explicit sub-targets. The fast gate is for PR cycle time; the full nightly gate is the long-running KPI.

| Gate | Tests | Target | When |
|---|---|---|---|
| Fast gate (`tool.zephyrus.coverage_fast`) | unit + smoke | **90%** (PROTEUS-ecosystem ceiling) | Every PR |
| Full gate (`tool.coverage.report`) | unit + smoke + integration + slow | **90%** | Nightly |

The ratchet is one-way (`tools/update_coverage_threshold.py`), capped at 90%. Never manually decrease the threshold. The CI guard in `.github/workflows/tests.yaml` rejects any PR that lowers `[tool.coverage.report].fail_under` below `min(base_ref, 90.0)`.

What this means for adding tests:

- A new closed-form helper in a utility source: a unit test is sufficient.
- A new function in `escape.py`: a unit test (counts toward both gates), plus a `physics_invariant` tag if it qualifies. If the function feeds a published benchmark, a `reference_pinned` test goes with it (counts toward the per-source-file inventory in `docs/Validation/escape.md`).
- A new MORS-coupled path: an integration-tier test that performs the real lookup and pins against a fixed stellar-age reference.

---

## 16. Failure modes to recognize on review

These are real patterns to watch for. The lint script catches some of them mechanically; reviewers catch the rest.

| Pattern | Example | Why it slips | Fix |
|---|---|---|---|
| **Degenerate radius geometry** | A `scaling=2` test uses `Rp == Rxuv`, so `Rp*Rxuv**2 == Rxuv**3` and the `scaling=3` branch gives the identical value. A regression that swaps the default scaling passes silently. | The test input makes the two branches indistinguishable. | Choose `Rp != Rxuv` and pin both branches; assert they differ beyond tolerance (Section 2 rule 4). |
| **Weak tidal signal** | A tidal test at `a = 1 au` gives `K_tide ~ 0.99`; a dropped `K_tide` changes the result by ~1 %, below a loose tolerance. | The geometry makes the tidal correction negligible. | Use a close-in orbit (`a ~ 0.02 au`) so `K_tide` departs from 1 by tens of percent; pin `K_tide` and include the no-tidal value as the discrimination guard. |
| **Hypothesis seed and version stability** | `@given(...)` passes on hypothesis 6.0 with the default seed; on 6.100 the strategy produces a different sequence and a previously-hidden flake surfaces. | Hypothesis seed semantics changed between minor releases; the author relied on implicit determinism. | Add `pytest.importorskip('hypothesis')` at module top. Use `@settings(derandomize=True)` or pass `--hypothesis-seed=<fixed>` in CI. Document the chosen seed in the test docstring. |
| **Unit slip masked by round numbers** | The reference literal is computed with `a` in au instead of m; the pinned value looks plausible. | The test author reused a config-unit value without converting to SI. | Pin the magnitude with an explicit scale guard (`1e6 < val < 1e8`) and name the unit in the comment. |
| **MORS mock returns constant** | A unit test mocks `mors.Star.Value` to return `1.0` for every quantity; a dilution or conversion bug produces a plausible-looking but wrong flux. | The constant mock hides the scale of the real quantity. | Return a plausible `Lx`, `Leuv` pair and assert the derived flux against a hand-computed value. |
| **Module-level constant patched at the wrong site** | `monkeypatch.setattr('zephyrus.constants.G', ...)` on a source that already star-imported `G` | The star import bound the old value into the `escape` namespace | `monkeypatch.setattr('zephyrus.escape.G', ...)` at the use site (Section 6) |
| **Optional dep imported unconditionally** | `import hypothesis` at module top | `pip install --no-deps` build skips the optional install | `pytest.importorskip('hypothesis')` at module top |
| **Stale marker after refactor** | File renamed without re-applying the module-level `pytestmark` | CI marker filter passed because of per-function markers; coverage tier became invisible | Restore module-level `pytestmark = [pytest.mark.unit, pytest.mark.timeout(30)]` |
| **Trivially-true on implicit None** | `def fixture(): pass`; `def test_x(fixture): assert fixture is None` | Fixture returned None implicitly; test passes for the wrong reason | Delete the test |

When you spot a new variant of these, add it here.

---

## 17. Sister rules (cross-link)

- `.github/copilot-instructions.md` "Testing Standards" -- the high-level summary readers without `tests/**` context see first.
- `.github/.claude/rules/zephyrus-code-review.md` "Test marker discipline" -- the review-pass gate that backs up the rules in this file. Also contains domain-aware physics checks (tidal-correction propagation, radius-scaling exponent safety, XUV-flux coupling, PROTEUS-coupling patterns) that apply when reviewing the **source** code that tests cover.

Any change to the rule set: update both files in the same commit and call out the cross-reference in the commit body.
