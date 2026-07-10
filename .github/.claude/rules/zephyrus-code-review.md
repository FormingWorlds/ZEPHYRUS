---
description: ZEPHYRUS-specific code review criteria for the atmospheric-escape module. Applies domain expertise (energy-limited mass loss, tidal correction, XUV-flux coupling, MORS coupling, PROTEUS coupling) to all code review in this repo.
---

# ZEPHYRUS Code Review Criteria

When reviewing ZEPHYRUS code (either your own or via code-reviewer agents), apply these domain-specific checks in addition to standard code quality review.

> **Discovery note.** ZEPHYRUS keeps its Claude-Code rule files under `.github/.claude/rules/` (not the conventional repo-root `.claude/`) so they can be tracked in git and shared across collaborators. Claude does NOT auto-discover them at this path; the repo-root `CLAUDE.md` (symlinked to `.github/copilot-instructions.md`) names this file and `zephyrus-tests.md` explicitly. **Before opening any review pass, read both this file and `zephyrus-tests.md`.**

## Physics plausibility

- The energy-limited mass-loss rate must be non-negative. Flag any code path where `escape_EL` could go negative for physically valid inputs (positive `epsilon`, `Fxuv`, radii, mass).
- The tidal correction factor `K_tide` must lie in `(0, 1]` for a bound orbit. `K_tide = 1 - 3/(2*ksi) + 1/(2*ksi**3)` with `ksi = Rhill/Rxuv`; it falls below zero when `ksi` is small (a close-in planet whose Hill radius approaches the XUV radius). A non-positive `K_tide` makes the escape rate diverge or flip sign. Flag any tidal path that does not guard against `ksi` approaching the `K_tide -> 0` root (`ksi ~ 1.6`).
- The Hill radius `Rhill = a*(1-e)*(Mp/(3*Ms))**(1/3)` must exceed `Rxuv` for the tidal correction to be physical (`ksi > 1`). Flag any caller that can pass a semi-major axis small enough to violate this.
- Planetary and stellar masses, radii, and semi-major axis must be strictly positive. XUV flux and escape efficiency must be non-negative. Flag any path that lets a zero or negative geometric quantity reach the division.
- Escape efficiency `epsilon` is a dimensionless factor in the literature range `0.1 < epsilon < 0.6`. Flag a hard-coded `epsilon` outside `[0, 1]`.

## Unit convention boundaries

ZEPHYRUS works in SI internally:

- **`EL_escape` inputs**: `a`, `Rp`, `Rxuv` in metres; `Mp`, `Ms` in kilograms; `Fxuv` in W m-2; `e`, `epsilon` dimensionless; `scaling` an integer.
- **Output**: mass-loss rate in kg s-1.
- **MORS coupling**: MORS returns stellar XUV luminosities `Lx`, `Leuv` in erg s-1. Converting to a flux at the planet requires `ergcm2stoWm2` (erg s-1 cm-2 -> W m-2) AND the orbital distance in centimetres (`a_au * au2cm`) so the `1 / (4*pi*a**2)` geometric dilution is dimensionally consistent. The erg-vs-W and au-vs-m/cm boundaries are the recurring traps.
- **`constants.py`**: `G` is in SI (`m3 kg-1 s-2`); `G_cgs` is the cgs sibling. Verify the SI `G` is used in `EL_escape` and the cgs one never leaks into an SI expression.

When reviewing code that crosses these boundaries (a new flux calculation, a new MORS caller, a new PROTEUS-side caller), verify the unit is correct at each conversion site.

## Radius-scaling exponent safety

`EL_escape` dispatches the radius term by the `scaling` argument:

- `scaling=2` (default): `R_cubed = Rp * Rxuv**2` (Lopez, Fortney & Miller 2012).
- `scaling=3`: `R_cubed = Rxuv**3` (Lehmer & Catling 2017).
- any other value: `ValueError`.

When the default `scaling` is changed, or a new scaling branch is added, the change has fan-out across:

1. Tests pinning a scaling-specific reference value (the `reference_pinned` escape test names the scaling in its docstring and pins the exact literal).
2. Documentation pages citing the scaling-specific number (`docs/Validation/escape.md`, the model-overview page).
3. PROTEUS-side callers that pass `scaling` explicitly or rely on the default.

Required workflow for any scaling default change:

1. `git grep` the `scaling` argument; update every reference value tied to it.
2. Update the test discrimination guards (Section 2 rule of `zephyrus-tests.md`) so the test would fail loudly under the wrong-scaling regression.
3. Update `docs/Validation/escape.md` and any page that quotes a scaling-specific number.

A PR that changes the default `scaling` but does not touch the test reference values is a red flag during review.

## Tidal-correction propagation

The tidal branch (`tidal_contribution=True`) introduces three coupled quantities: `Rhill`, `ksi`, `K_tide`. Flag any new tidal code path that:

- Computes `K_tide` without asserting `ksi > 1` (the Hill radius must exceed the XUV radius).
- Lets `K_tide` reach the escape-rate denominator without a positivity guard, so a small `ksi` produces a negative or singular escape rate.
- Reorders the `Rhill` factors in a way that changes the `(1-e)` periapsis dependence (the correction uses the periapsis distance, not the semi-major axis alone).

The unit tests of `escape.py` verify the tidal branch increases the escape rate relative to the no-tidal branch (`K_tide < 1`) and pin `K_tide` at a discriminating close-in geometry; the review's job is to make sure the positivity defense is present before the test asks it to fire.

## PROTEUS coupling patterns

ZEPHYRUS is called by PROTEUS through the escape step. Two coupling patterns need explicit care during review.

### 1. XUV-flux hand-off

PROTEUS supplies the XUV flux `Fxuv` (W m-2) at the planet, computed from the MORS stellar track and the current orbital distance, and ZEPHYRUS returns the mass-loss rate. The contract:

- `EL_escape` MUST treat `Fxuv` as an already-diluted flux at the planet; it must NOT re-apply the `1/(4*pi*a**2)` geometric factor.
- A change that folds the orbital dilution into `EL_escape` would double-count it against the PROTEUS-side caller. Flag any such change.

### 2. Escape-mass budget consistency

PROTEUS multiplies the ZEPHYRUS mass-loss rate by the timestep to remove atmospheric mass, and the coupled invariant is that the removed mass never exceeds the available atmospheric mass in a step. A regression in `EL_escape` that returns an inflated rate (wrong scaling exponent, missing `K_tide` in the denominator, dropped `epsilon`) would silently violate the PROTEUS-side budget. Flag any change to the mass-loss formula that is not accompanied by an updated discrimination guard in the escape test.

## Config mutability

Any dataclass or parameter object used to carry user input must not be mutated at runtime after IC. Flag any code that sets `config.X.Y = value` outside of config initialization. Use local variables instead.

## Cross-module constant duplication

Physical constants (`G`, `kb`, `c`) and unit conversions (`au2m`, `au2cm`, `ergcm2stoWm2`) are defined in `src/zephyrus/constants.py`. When reviewing code that uses a physical constant, check that the import is from `zephyrus.constants` and not re-derived. A new constant introduced as a literal in a body (e.g. `6.674e-11` inline for `G`) is a red flag.

## Star-import boundary

`escape.py` uses `from zephyrus.constants import *` and `from zephyrus.planets_parameters import *`. This is the existing module convention; ruff's `F403` / `F405` are ignored for this repo. When reviewing a new source file, prefer explicit imports for anything new; do not extend the star-import surface without cause, and never let a star import shadow a function-local name.

## Test marker discipline

Every test file must begin with a module-level `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<budget>)]` (unit/30 s, smoke/60 s, integration/300 s, slow/3600 s). Per-function markers are additive but do not replace the module-level marker; CI runs `pytest -m "(unit or smoke) and not skip"` and any file missing the tier marker ships untested.

## Test quality (cross-reference)

Test-content rules (anti-happy-path, discriminating-value guards, physics-invariant tiering, `physics_invariant` / `reference_pinned` certification markers, adversarial-review trigger, mocking discipline, `importorskip` + module-constant-monkeypatch traps, tidal-correction propagation, MORS-flux coupling, hypothesis seed stability) live in [`zephyrus-tests.md`](zephyrus-tests.md). When reviewing tests, apply both files: this one for marker discipline and review-pass gate, the deep-dive for the content contract.

## Sister rules (cross-link)

- [`.github/copilot-instructions.md`](../../copilot-instructions.md) "Testing Standards" -- high-level rules visible to all readers. Repo-root `CLAUDE.md` is a symlink to this file.
- [`zephyrus-tests.md`](zephyrus-tests.md) -- test quality deep-dive; the canonical source for anti-happy-path patterns and the validation certification markers.
