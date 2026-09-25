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

`collision.mass_loss` returns the fraction of the atmosphere lost (Kegerreis et al. 2020, Eqn. 1), clipped to [0, 1]. Keep its input guards: impact parameter in [0, 1]; masses, densities and radii positive and finite; collision speed non-negative and finite.

## Star imports

`escape.py` uses `from zephyrus.constants import *` and `from zephyrus.planets_parameters import *` (ruff `F403` and `F405` are ignored for this). New code imports names explicitly; do not extend the star-import surface, and do not define a module-level name in `escape.py` that a star import would rebind.
