![ZEPHYRUS banner](https://raw.githubusercontent.com/FormingWorlds/ZEPHYRUS/main/docs/assets/ZEPHYRUS_logo_white.png#gh-light-mode-only)
![ZEPHYRUS banner](https://raw.githubusercontent.com/FormingWorlds/ZEPHYRUS/main/docs/assets/ZEPHYRUS_logo_black.png#gh-dark-mode-only)

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Docs](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/docs.yaml?branch=main&label=Docs)](https://proteus-framework.org/ZEPHYRUS/)
[![codecov](https://img.shields.io/codecov/c/github/FormingWorlds/ZEPHYRUS?label=coverage&logo=codecov)](https://app.codecov.io/gh/FormingWorlds/ZEPHYRUS)
[![Unit Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/tests.yaml?branch=main&label=Unit%20Tests)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/tests.yaml)
[![Integration Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/nightly.yml?branch=main&label=Integration%20Tests)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/nightly.yml)

**ZEPHYRUS** is the atmospheric escape module of the [PROTEUS](https://proteus-framework.org/PROTEUS) coupled atmosphere-interior evolution framework. It computes the rate at which a rocky or sub-Neptune planet loses its atmosphere to space, how that loss is shared between chemical species, and the fraction of the atmosphere a giant impact removes.

Named after Zephyrus, the Greek god of the west wind and a messenger of spring. Pronounced *ZEF-ir-us*.

## Escape model

The escape model is the regime dispatcher, `zephyrus.dispatch`. It takes one planetary state (mass, radius, orbit, stellar fluxes, and an atmosphere profile), decides which escape physics applies, and returns a regime label, the mass-loss rate, the rate of each element, flags, and diagnostics. The six regimes are boil-off, three limits on an XUV-driven hydrodynamic wind (energy-limited, recombination-limited, and photon-limited), hydrostatic Jeans escape with diffusion-limited supply, and Roche-lobe overflow. A hydrodynamic wind is split over species by a fractionation closure, so heavy species can lag the flow or stay behind.

`zephyrus.escape.EL_escape` applies the energy-limited rate alone, with an optional tidal correction after Erkaev et al. (2007), to every state without checking which regime applies. For giant impacts, `zephyrus.collision.impact_loss` returns the eroded fraction of the atmosphere from the scaling laws of Kegerreis et al. (2020) and Roche et al. (2026).

## Documentation

Full documentation is at **[proteus-framework.org/ZEPHYRUS](https://proteus-framework.org/ZEPHYRUS/)**, including:

- [Getting started](https://proteus-framework.org/ZEPHYRUS/getting_started.html): installation and the quickest path to a first run.
- [Tutorials](https://proteus-framework.org/ZEPHYRUS/Tutorials/dispatch.html): dispatching escape regimes on synthetic atmospheres, and a first energy-limited calculation for the Sun-Earth system.
- [How-to guides](https://proteus-framework.org/ZEPHYRUS/How-to/installation.html): install, run the tests, build the documentation.
- [Explanations](https://proteus-framework.org/ZEPHYRUS/Explanations/model.html): model overview, the escape regimes, fractionation, giant impacts, coupling to PROTEUS, limitations, and the testing suite.
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
import math

from scipy.constants import Stefan_Boltzmann as sigma_sb

from zephyrus import EscapeInputs, dispatch
from zephyrus.constants import au2m
from zephyrus.planets_parameters import Ls, Me, Ms, Re
from zephyrus.profiles import isothermal_profile

# An Earth-mass, Earth-radius planet with a pure CO2 atmosphere at 0.0775 au
a = 0.0775 * au2m
T_eq = (Ls / (16 * math.pi * sigma_sb * a**2)) ** 0.25        # 1000 K
state = EscapeInputs(
    M_p=Me, R_p=Re, M_star=Ms, a=a, e=0.0, T_eq=T_eq,
    F_xuv=10.0,                                  # XUV flux at the planet [W m-2]
    F_bol=Ls / (4 * math.pi * a**2),             # bolometric flux [W m-2]
    F_int=1.0,                                   # interior heat flux [W m-2]
    kappa_photo=0.01,                            # grey photospheric opacity [m2 kg-1]
    profile=isothermal_profile(Me, Re, T_eq, {'CO2': 1.0}, 1e7, 1e-5),
)
result = dispatch(state)
print(result.regime, f'{result.mdot:.3e} kg/s')   # hydrodynamic:EL 2.348e+06 kg/s
print({el: f'{r:.3e}' for el, r in result.per_species.items()})
```

The [dispatcher tutorial](https://proteus-framework.org/ZEPHYRUS/Tutorials/dispatch.html) builds on this state step by step, and the [first-run tutorial](https://proteus-framework.org/ZEPHYRUS/Tutorials/first_run.html) applies `EL_escape` along a MORS stellar history.

## License

[Apache License 2.0](LICENSE.md). ZEPHYRUS is part of the [PROTEUS framework](https://proteus-framework.org/).
