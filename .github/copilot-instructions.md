# ZEPHYRUS short instructions

The full instructions are in `AGENTS.md` (repository rules and review points) and `tests/AGENTS.md` (test rules). This file repeats what a review or a chat needs when it cannot read them.

## Commands

```bash
pip install -e ".[develop]"
pytest -m "(unit or smoke) and not skip"
ruff check <changed files> && ruff format <changed files>
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
```

## Review checklist

- Units are SI: `EL_escape` takes m, kg and W m-2 and returns kg s-1; MORS luminosities in erg s-1 are converted with `ergcm2stoWm2` and the distance in cm.
- `Fxuv` from PROTEUS is already diluted to the planet; `EL_escape` does not divide by `4 pi a**2` again.
- The escape rate stays non-negative, linear in `Fxuv` and decreasing with `Mp`; geometric quantities are positive before a division.
- The tidal branch keeps the `ksi > 1` guard before `K_tide` is used, and the `(1 - e)` periapsis factor in `Rhill`.
- A change of the default `scaling` updates the pinned escape values and `docs/Validation/escape.md`.
- Constants come from `zephyrus.constants`; `G_cgs` never enters an SI expression.
- Tests: `src/zephyrus/<file>.py` is tested in `tests/test_<file>.py`; each file has `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<s>)]`; each test has a docstring, at least 2 assertions, an edge case and the error path; floats are compared with a tolerance; physics tests carry `physics_invariant`, and pinned values have sign, scale and wrong-formula guards.
- Commit messages and pull-request text describe the change, with no tool attribution.
