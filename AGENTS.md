# ZEPHYRUS agent instructions

ZEPHYRUS computes atmospheric escape for PROTEUS: the energy-limited mass-loss rate driven by stellar XUV irradiation (`escape.py`, `EL_escape`) and the atmospheric mass lost in a giant impact (`collision.py`, `mass_loss`). Before a first edit:

- Tests for `src/zephyrus/<file>.py` go in `tests/test_<file>.py`; the test rules are in `tests/AGENTS.md`.
- `escape.py` and `collision.py` are the physics sources; `constants.py` and `planets_parameters.py` are utilities.
- The escape rate is non-negative, linear in `Fxuv`, decreasing with planet mass, and zero at `Fxuv = 0`; PROTEUS removes `rate * dt` from the atmosphere, so an inflated rate breaks its mass budget.
- These commands decide whether a change is ready (CI runs the tests, the structure check, the test-quality lint and the agent-file check):

```bash
pytest -m "(unit or smoke) and not skip"
ruff check <changed files> && ruff format <changed files>
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
python tools/agents/check_agents_md.py
```

<!-- fwl-core:begin sha256=c232bac058cc749d -->
## PROTEUS ecosystem rules

PROTEUS couples separate module repositories into one model of the evolution of rocky planets and their atmospheres. Its config selects one module per process; the config values are the lowercase names:

- `interior_energetics.module`: `aragog`, `spider` (or the built-in `boundary`); `interior_struct.module`: `zalmoxis`, `spider`
- `atmos_clim.module`: `agni`, `janus`, both with SOCRATES spectral radiative transfer
- `outgas.module`: `calliope`, `atmodeller`; `atmos_chem.module`: `vulcan`
- `star.module`: `mors`; `escape.module`: `zephyrus`, `boreas`; `orbit.module` (tides): `obliqua`, `lovepy`
- fwl-io downloads the reference data the modules use

A module change to anything PROTEUS calls or reads (a function signature, a config key, an output column, a unit) can break the coupled model while the module's own tests pass. Search `src/proteus/` in PROTEUS for the name before you change it, and name the affected PROTEUS call sites in the pull request.

### Physics and numerics

- Units differ between modules: MORS and VULCAN use cgs; CALLIOPE and atmodeller work with pressures in bar; AGNI takes surface pressure in bar at setup and holds Pa inside; PROTEUS keeps pressure in bar and time in years in `hf_row` and SI otherwise, and its config reference pages state each key's unit (for example stellar mass in M_sun, initial partial surface pressures in bar, stellar age in Gyr). State the unit of every physical quantity in its docstring and convert explicitly at the boundary: a unit mismatch between two modules passes the tests of both and shows only in the coupled run.
- A physics test must fail for the most plausible wrong formula: check a conservation law, a bound, a monotonicity, or a published or analytical value, and assert that the wrong result falls outside the tolerance.
- Take each physical constant from one source per repository (in Python `scipy.constants` or the module's constants file). Two retyped values of one constant differ at round-off and hide real differences between code paths.
- Do not loosen a solver tolerance, a conservation check or a clamp to make a run or a test pass. Find the cause; a loosened check also hides the next defect.

### Branches and pull requests

- Work on a feature branch `<initials>/<short-description>`; `main` changes only through a reviewed pull request.
- Fill every section of the repository's pull-request template, where it has one.
- Before you push, run the checks listed at the top of this `AGENTS.md` and in `tests/AGENTS.md`.

### Where knowledge goes

Rules for every contributor go in the `AGENTS.md` files, the reason for a change in its commit message and pull-request description, and the scientific validation of a module in its `docs/Validation/` pages, where the repository has them. Do not add memory, notes or decision-log files to the repository: nobody maintains them, and they go stale.
<!-- fwl-core:end -->

<!-- fwl-voice:begin sha256=a943ab6c93ddad24 -->
### Commit messages and public text

Commit messages, pull-request text, code comments, docstrings and test names describe the change and the current state of the code. They name no tool used to write the change and carry no tool-attribution trailer.
<!-- fwl-voice:end -->

## Environment

`pip install -e ".[develop]"` and `pre-commit install -f`; there is nothing to compile. `fwl-mors` is a runtime dependency; only the integration tier (`test_earth.py`) needs its stellar data (`FWL_DATA` set, `mors download all`). The pre-commit hook runs `ruff check --fix`; run `ruff format` yourself on the files you change. `constants.py` and `planets_parameters.py` keep hand-aligned tables that `ruff format` would rewrite; leave their layout. Docs build with `zensical serve` after `pip install -e '.[docs]'`.

## Physics and coupling contract

- Units are SI throughout: `EL_escape` takes `a`, `Rp`, `Rxuv` in m, `Mp`, `Ms` in kg, `Fxuv` in W m-2, `e` and `epsilon` dimensionless, and returns kg s-1. MORS returns `Lx`, `Leuv` in erg s-1: convert with `ergcm2stoWm2` and the orbital distance in cm (`a_au * au2cm`). The erg against W and au against m or cm conversions are where errors enter.
- `Fxuv` arrives from PROTEUS already diluted to the planet (`src/proteus/escape/wrapper.py`, `run_zephyrus`); `EL_escape` must not apply `1 / (4 pi a**2)` again.
- `scaling=2` (default) uses `Rp * Rxuv**2`, `scaling=3` uses `Rxuv**3`, any other value raises `ValueError`. Changing the default changes the pinned values in the escape tests, `docs/Validation/escape.md` and every PROTEUS caller that relies on the default: update all three together.
- Tidal branch: `ksi = Rhill / Rxuv` with `Rhill = a (1 - e) (Mp / (3 Ms))**(1/3)`, and `K_tide = (ksi - 1)**2 (2 ksi + 1) / (2 ksi**3)`. `K_tide` is in (0, 1) for `ksi > 1`, and the rate divides by it, so it diverges as `ksi` approaches 1. The source raises `ValueError` for `ksi <= 1`; every tidal path keeps that guard, and the periapsis factor `(1 - e)` stays in `Rhill`.
- `collision.py` (Kegerreis et al. 2020, Eqn. 1) raises `ValueError` for an impact parameter outside [0, 1], a non-positive or non-finite mass, density or radius, and a negative collision speed.
- Constants and conversions (`G`, `kb`, `au2m`, `au2cm`, `ergcm2stoWm2`) come from `zephyrus.constants`; `G` is SI and `G_cgs` must not enter an SI expression. `escape.py` star-imports `constants` and `planets_parameters` (ruff `F403`, `F405` ignored); new code imports names explicitly.

## Review

Check each change against these points; `.github/agent-rules/code-review.md` has the detail.

- The escape rate stays non-negative for valid inputs; geometric quantities stay strictly positive before any division; `epsilon` stays in [0, 1].
- Units at the MORS and PROTEUS boundaries; no second orbital dilution of `Fxuv`.
- A formula change comes with an updated discrimination guard in the escape tests (wrong scaling, dropped `K_tide`, dropped `epsilon`).
- The `ksi > 1` guard is in place before `K_tide` is used.
- A change of the default `scaling` updates the pinned values and `docs/Validation/escape.md`.
- No retyped constant literals; no run-time mutation of parameter objects.
- Tests follow `tests/AGENTS.md`.
