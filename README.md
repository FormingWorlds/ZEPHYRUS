![ZEPHYRUS banner](https://raw.githubusercontent.com/FormingWorlds/ZEPHYRUS/main/docs/assets/ZEPHYRUS_logo_white.png#gh-light-mode-only)
![ZEPHYRUS banner](https://raw.githubusercontent.com/FormingWorlds/ZEPHYRUS/main/docs/assets/ZEPHYRUS_logo_black.png#gh-dark-mode-only)

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/docs.yaml?branch=main&label=Docs)](https://proteus-framework.org/ZEPHYRUS/)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/ZEPHYRUS?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/ZEPHYRUS)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/nightly.yml)

**ZEPHYRUS** is the atmospheric escape module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It computes the rate at which an exoplanet loses its atmosphere to space under high-energy stellar irradiation, and the fraction of its atmosphere removed by a giant impact.

Given a planet's mass and radius, the radius at which XUV radiation is absorbed, the stellar XUV flux at the planet, and an escape efficiency, ZEPHYRUS returns the energy-limited mass-loss rate. An optional tidal correction after Erkaev et al. (2007) enhances the rate for close-in or eccentric orbits, and two forms of the XUV-heated cross section are available. The stellar XUV flux is typically supplied by [MORS](https://proteus-framework.org/MORS/) evolution tracks. For giant impacts during accretion, ZEPHYRUS returns the eroded fraction of the target's atmosphere from the collision speed, masses, densities, and impact angle, following the scaling law of Kegerreis et al. (2020).

Named after Zephyrus, the Greek god of the west wind and a messenger of spring. Pronounced *ZEF-ir-us*.

## Escape model

Energy-limited (EL) escape: the deposited XUV power divided by the gravitational binding energy of the escaping gas, with an optional tidal enhancement `K_tide`. The XUV cross section scales as `Rp * Rxuv**2` (default) or `Rxuv**3`, and the tidal factor applies only outside the Roche lobe (`R_Hill / R_XUV > 1`).

## Documentation

Full documentation is at **[proteus-framework.org/ZEPHYRUS](https://proteus-framework.org/ZEPHYRUS/)**, including:

- [Getting started](https://proteus-framework.org/ZEPHYRUS/getting_started.html): installation and the quickest path to a first run.
- [Tutorial](https://proteus-framework.org/ZEPHYRUS/Tutorials/first_run.html): a first escape calculation for the Sun-Earth system.
- [How-to guides](https://proteus-framework.org/ZEPHYRUS/How-to/installation.html): install, run the tests, build the documentation.
- [Explanations](https://proteus-framework.org/ZEPHYRUS/Explanations/model.html): model overview, coupling to PROTEUS, limitations, and the testing suite.
- [API reference](https://proteus-framework.org/ZEPHYRUS/Reference/api/index.html) and [validation anchors](https://proteus-framework.org/ZEPHYRUS/Validation/escape.html): every public function with NumPy-style docstrings, plus the reference-pinned test inventory.

## Installation

```console
pip install fwl-zephyrus
```

Or, for development:

```console
git clone https://github.com/FormingWorlds/ZEPHYRUS.git
cd ZEPHYRUS
pip install -e .[develop,docs]
```

The `docs` extra pulls in [Zensical](https://zensical.org/) so you can build this documentation locally with `zensical serve`.

### Input data

ZEPHYRUS uses [MORS](https://proteus-framework.org/MORS/) for the stellar XUV flux, which needs a set of stellar evolution tracks from the [OSF repository](https://osf.io/9u3fb/). Set the `FWL_DATA` environment variable to the directory where the data should live, then download the tracks:

```console
export FWL_DATA=/your/local/path/FWL_DATA   # add to ~/.bashrc to persist
mors download all
```

## Quick start

```python
import numpy as np
import mors
from zephyrus.escape import EL_escape
from zephyrus.constants import au2cm, au2m, ergcm2stoWm2
from zephyrus.planets_parameters import Me, Re, a_earth, e_earth

# XUV flux at Earth from a solar-mass star, over its evolution (MORS tracks)
star = mors.Star(Mstar=1.0, Omega=1.0)
Lxuv = star.Tracks['Lx'] + star.Tracks['Leuv']          # erg s-1
Fxuv = Lxuv / (4 * np.pi * (a_earth * au2cm) ** 2)      # erg s-1 cm-2
Fxuv_SI = Fxuv * ergcm2stoWm2                           # W m-2

# Energy-limited mass-loss rate [kg s-1] across the star's history.
# a, e and Ms only enter the tidal branch (tidal_contribution=True).
mdot = EL_escape(
    tidal_contribution=False,
    a=a_earth * au2m, e=e_earth,
    Mp=Me, Ms=1.0, epsilon=0.15,
    Rp=Re, Rxuv=Re, Fxuv=Fxuv_SI,
)   # aligned with star.Tracks['Age']
```

See the [first-run tutorial](https://proteus-framework.org/ZEPHYRUS/Tutorials/first_run.html) for the full walkthrough.

## License

[Apache License 2.0](LICENSE.md). ZEPHYRUS is part of the [PROTEUS framework](https://proteus-framework.org/).
