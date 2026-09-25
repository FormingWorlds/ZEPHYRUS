# ZEPHYRUS short instructions

The full instructions are in `AGENTS.md` (repository rules and review points) and `tests/AGENTS.md` (test rules). This file repeats what a review or a chat needs when it cannot read them.

## Commands

```bash
pip install -e ".[develop]"
pytest -m "(unit or smoke) and not skip" --cov=zephyrus --cov-fail-under=90
ruff check <changed files> && ruff format <changed files>  # not constants.py, planets_parameters.py
bash tools/validate_test_structure.sh
python tools/check_test_quality.py --check
python tools/agents/check_agents_md.py
```

## Review checklist

- Units are SI: `EL_escape` takes m, kg and W m-2 and returns kg s-1; MORS luminosities in erg s-1 become a flux by dividing by `4 pi a**2` with `a` in cm and multiplying by `ergcm2stoWm2`.
- `Fxuv` from PROTEUS is already diluted to the planet; `EL_escape` does not divide by `4 pi a**2` again.
- The escape rate stays non-negative, linear in `Fxuv` and decreasing with `Mp`; geometric quantities are positive before a division.
- The tidal branch keeps the `ksi > 1` guard before `K_tide` is used, and the `(1 - e)` periapsis factor in `Rhill`.
- PROTEUS passes `scaling=3` and the tests pass `scaling` explicitly; a change of the default updates the `EL_escape` docstring and `docs/Validation/escape.md`, which name it.
- Constants come from `zephyrus.constants`; `G_cgs` never enters an SI expression.
- Tests: `src/zephyrus/<file>.py` is tested in `tests/test_<file>.py`; each file has `pytestmark = [pytest.mark.<tier>, pytest.mark.timeout(<s>)]`; each test has a docstring, at least 2 assertions, an edge case and the error path; floats are compared with a tolerance; a test that asserts a physical invariant carries `physics_invariant`, a test against a published or analytical value also carries `reference_pinned`, and pinned values have sign, scale and wrong-formula guards.
- Commit messages and pull-request text describe the change, with no tool attribution.
