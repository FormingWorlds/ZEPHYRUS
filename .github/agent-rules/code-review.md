# ZEPHYRUS Code Review Criteria

The detail behind the Review section of `AGENTS.md`. Apply these domain checks in addition to general code quality review.


## Physics plausibility

- The energy-limited mass-loss rate must be non-negative. Flag any code path where `escape_EL` could go negative for physically valid inputs (positive `epsilon`, `Fxuv`, radii, mass).
- The tidal correction factor `K_tide = 1 - 3/(2*ksi) + 1/(2*ksi**3)` with `ksi = Rhill/Rxuv` factors exactly to `(ksi-1)**2 * (2*ksi + 1) / (2*ksi**3)`, so it is non-negative for every `ksi > 0` with a double (tangential) root at `ksi = 1`; it never goes negative. It lies in `(0, 1)` for `ksi > 1` (the physical regime, rising toward 1 as the orbit widens) and exceeds 1 for `ksi < 1/sqrt(3)`. Because the escape rate divides by `K_tide`, the rate diverges as `ksi -> 1`, so the energy-limited approximation is valid only for `ksi > 1`. The source raises `ValueError` for `ksi <= 1`. Flag any tidal path that computes `K_tide` without the `ksi > 1` domain guard, or that reintroduces a `K_tide -> 0` singularity into the escape-rate denominator.
- The Hill radius `Rhill = a*(1-e)*(Mp/(3*Ms))**(1/3)` must exceed `Rxuv` for the tidal correction to be physical (`ksi > 1`); the source enforces this with a `ValueError` at `ksi <= 1`. Flag any change that lets a semi-major axis small enough to violate `ksi > 1` reach the `K_tide` division without raising.
- Planetary and stellar masses, radii, and semi-major axis must be strictly positive. XUV flux and escape efficiency must be non-negative. Flag any path that lets a zero or negative geometric quantity reach the division.
- Escape efficiency `epsilon` is a dimensionless factor in the literature range `0.1 < epsilon < 0.6`. Flag a hard-coded `epsilon` outside `[0, 1]`.

## Unit convention boundaries

ZEPHYRUS works in SI internally:

- **`EL_escape` inputs**: `tidal_contribution` a bool; `a`, `Rp`, `Rxuv` in metres; `Mp`, `Ms` in kilograms; `Fxuv` in W m-2; `e`, `epsilon` dimensionless; `scaling` an integer.
- **Output**: mass-loss rate in kg s-1.
- **MORS coupling**: MORS returns stellar XUV luminosities `Lx`, `Leuv` in erg s-1. Converting to a flux at the planet requires `ergcm2stoWm2` (erg s-1 cm-2 -> W m-2) AND the orbital distance in centimetres (`a_au * au2cm`) so the `1 / (4*pi*a**2)` geometric dilution is dimensionally consistent. The erg-vs-W and au-vs-m/cm boundaries are the recurring traps.
- **`constants.py`**: `G` is in SI (`m3 kg-1 s-2`); `G_cgs` is the cgs sibling. Verify the SI `G` is used in `EL_escape` and the cgs one never leaks into an SI expression.

When reviewing code that crosses these boundaries (a new flux calculation, a new MORS caller, a new PROTEUS-side caller), verify the unit is correct at each conversion site.

## Radius-scaling exponent safety

`EL_escape` dispatches the radius term by the `scaling` argument:

- `scaling=2` (default): `R_cubed = Rp * Rxuv**2` (Watson 1981; Lammer 2003; Erkaev 2007, Eq. 21).
- `scaling=3`: `R_cubed = Rxuv**3` (Lopez, Fortney & Miller 2012; Lopez & Fortney 2013; Lehmer & Catling 2017).
- any other value: `ValueError`.

PROTEUS passes `scaling=3` explicitly (`run_zephyrus` in `src/proteus/escape/wrapper.py`), and the escape tests pass `scaling` explicitly, so a change of the default reaches only callers that omit the argument. With a default change, update the default named in the `EL_escape` docstring and in `docs/Validation/escape.md`. A new scaling branch needs its own pinned test and a wrong-scaling discrimination guard against the other branches.

## Tidal-correction propagation

The tidal branch (`tidal_contribution=True`) introduces three coupled quantities: `Rhill`, `ksi`, `K_tide`. Flag any new tidal code path that:

- Computes `K_tide` without a `ksi > 1` domain guard (the Hill radius must exceed the XUV radius).
- Lets `K_tide` reach the escape-rate denominator without that guard, so a `ksi` near 1 produces a singular (divergent) escape rate. `K_tide` is non-negative for all `ksi > 0`, so the failure mode is divergence at `ksi -> 1`, not a sign flip.
- Reorders the `Rhill` factors in a way that changes the `(1-e)` periapsis dependence (the correction uses the periapsis distance, not the semi-major axis alone).

The unit tests of `escape.py` verify the tidal branch increases the escape rate relative to the no-tidal branch (`K_tide < 1`), pin `K_tide` at a discriminating close-in geometry, and assert the `ksi <= 1` domain guard raises `ValueError`; the review's job is to make sure the domain guard is present before the test asks it to fire.

## PROTEUS coupling patterns

ZEPHYRUS is called by PROTEUS through the escape step. Two coupling patterns need explicit care during review.

### 1. XUV-flux hand-off

PROTEUS supplies the XUV flux `Fxuv` (W m-2) at the planet, computed from the MORS stellar track and the current orbital distance, and ZEPHYRUS returns the mass-loss rate. The contract:

- `EL_escape` MUST treat `Fxuv` as an already-diluted flux at the planet; it must NOT re-apply the `1/(4*pi*a**2)` geometric factor.
- A change that folds the orbital dilution into `EL_escape` would double-count it against the PROTEUS-side caller. Flag any such change.

### 2. Escape-mass budget consistency

PROTEUS multiplies the ZEPHYRUS mass-loss rate by the timestep and caps the mass removed in one step at a fraction of the escapable reservoir (`limit_escape_step`). A regression in `EL_escape` that inflates the rate (a wrong scaling exponent, a dropped `epsilon`) hits that cap and changes the coupled evolution; a dropped `K_tide` lowers the rate of close-in planets. Flag any change to the mass-loss formula that is not accompanied by an updated discrimination guard in the escape test.

## Giant-impact mass loss

`collision.mass_loss` returns the fraction of the atmosphere lost (Kegerreis et al. 2020, Eqn. 1), clipped to [0, 1]. Keep its input guards: impact parameter in [0, 1], masses, densities and radii positive and finite, collision speed non-negative and finite.

## Cross-module constant duplication

Physical constants (`G`, `kb`, `c`) and unit conversions (`au2m`, `au2cm`, `ergcm2stoWm2`) are defined in `src/zephyrus/constants.py`. When reviewing code that uses a physical constant, check that the import is from `zephyrus.constants` and not re-derived. A new constant introduced as a literal in a body (e.g. `6.674e-11` inline for `G`) is a red flag.

## Star-import boundary

`escape.py` uses `from zephyrus.constants import *` and `from zephyrus.planets_parameters import *`. This is the existing module convention; ruff's `F403` / `F405` are ignored for this repo. When reviewing a new source file, prefer explicit imports for anything new; do not extend the star-import surface without cause, and never let a star import shadow a function-local name.

