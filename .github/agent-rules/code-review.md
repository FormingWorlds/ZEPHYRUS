# ZEPHYRUS review detail

The detail behind the Review section of `AGENTS.md`. Apply these domain checks in addition to general code quality review.

## Escape efficiency

`epsilon` is dimensionless; published values lie between about 0.1 and 0.6. `EL_escape` does not check it, so flag a hard-coded `epsilon` outside [0, 1].

## Radius scaling

- `scaling=2` (default): `R_cubed = Rp * Rxuv**2` (Watson 1981; Lammer 2003; Erkaev 2007, Eq. 21).
- `scaling=3`: `R_cubed = Rxuv**3` (Lopez, Fortney & Miller 2012; Lopez & Fortney 2013; Lehmer & Catling 2017).
- A new scaling branch needs its own pinned test and a wrong-scaling discrimination guard against the other branches.

## Escape in the coupled model

A formula change that moves the rate (a wrong scaling exponent, a dropped `epsilon` or `K_tide`) changes the coupled evolution in PROTEUS: an inflated rate runs into the step cap named in `AGENTS.md`, and a dropped `K_tide` lowers the rate of close-in planets. Flag a change to the mass-loss formula that comes without an updated discrimination guard in the escape tests.

## Giant-impact mass loss

`collision.py` provides `impact_loss`, `mass_loss` (Kegerreis et al. 2020), and `mass_loss_roche2026` (Roche et al. 2026), returning fractional loss in [0, 1]. Keep input guards: scalar inputs, impact parameter in [0, 1], strictly positive finite masses and radii, sub-luminal collision speed 0 <= v_c < c, and envelope mass fraction f_atm in [0, 1) with f_atm = 0 returning zero. For roche2026, M_t is the refractory mass while mutual escape speed v_esc uses total target mass M_t / (1 - f_atm).

## Star imports

`escape.py` uses `from zephyrus.constants import *` and `from zephyrus.planets_parameters import *` (ruff `F403` and `F405` are ignored for this). New code imports names explicitly; do not extend the star-import surface, and do not define a module-level name in `escape.py` that a star import would rebind.
