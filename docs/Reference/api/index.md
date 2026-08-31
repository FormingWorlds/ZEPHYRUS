# API overview

This is an overview of ZEPHYRUS' API for the user's reference. If you want to understand the underlying model, please visit [the ZEPHYRUS model](../../Explanations/model.md) and its topic pages. <br>

| Module | Description |
|---|---|
| **Escape prescriptions** | |
| [`zephyrus.escape`](escape.md) | `EL_escape`: energy-limited atmospheric mass-loss rate |
| [`zephyrus.dispatcher`](dispatcher.md) | `dispatch`: the escape-regime dispatcher (one call, one regime, one rate) |
| [`zephyrus.boiloff`](boiloff.md) | Bolometrically driven boil-off with Bondi and luminosity caps |
| [`zephyrus.nozzle`](nozzle.md) | Roche-lobe overflow: the tidally driven L1 nozzle transfer rate |
| [`zephyrus.hydrodynamic`](hydrodynamic.md) | Energy-limited and radiation-recombination-limited rates |
| [`zephyrus.hydrostatic`](hydrostatic.md) | Per-species Jeans escape with the diffusion-limited supply cap |
| [`zephyrus.collision`](collision.md) | `mass_loss`: fractional atmospheric loss in a giant impact |
| **Supporting physics** | |
| [`zephyrus.profiles`](profiles.md) | Atmosphere-profile container and the escape working levels |
| [`zephyrus.thermostat`](thermostat.md) | Wind temperature from local heating against radiative cooling |
| [`zephyrus.knudsen`](knudsen.md) | Collision cross sections and the sonic-point Knudsen switch |
| [`zephyrus.fractionation`](fractionation.md) | Simultaneous N-species fractionation closure |
| [`zephyrus.diagnostics`](diagnostics.md) | Regime diagnostics reported beside every dispatch verdict |
| **Reference data** | |
| [`zephyrus.diffusion`](diffusion.md) | Binary diffusion coefficient library with provenance classes |
| [`zephyrus.atomic_data`](atomic_data.md) | Cooling data tables and closed-form rate coefficients |
| [`zephyrus.composition`](composition.md) | Element masses, formula parsing, and composition handling |
| [`zephyrus.constants`](constants.md) | Physical constants and unit conversions (SI and CGS) |
| [`zephyrus.planets_parameters`](planet_parameters.md) | Sun, Earth, Jupiter, and TOI-561 system reference values |


The **source tree** is given by:

```
src/zephyrus
    ├── atomic_data.py          # Cooling data tables and rate coefficients
    ├── boiloff.py              # Bolometrically driven boil-off escape
    ├── collision.py            # Giant-impact atmospheric mass loss (mass_loss)
    ├── composition.py          # Element masses, formula parsing, atomization
    ├── constants.py            # Physical constants and unit conversions
    ├── diagnostics.py          # Regime diagnostics beside every verdict
    ├── diffusion.py            # Binary diffusion coefficients with provenance
    ├── dispatcher.py           # The escape-regime dispatcher (dispatch)
    ├── escape.py               # Energy-limited atmospheric escape (EL_escape)
    ├── fractionation.py        # N-species fractionation closure
    ├── hydrodynamic.py         # EL and radiation-recombination-limited rates
    ├── hydrostatic.py          # Jeans escape with diffusion-limited supply
    ├── __init__.py             # Package entry point and top-level exports
    ├── knudsen.py              # Cross sections and the Knudsen switch
    ├── nozzle.py               # Roche-lobe overflow through the L1 nozzle
    ├── planets_parameters.py   # Sun, Earth, Jupiter, TOI-561 reference values
    ├── profiles.py             # Profile container and escape working levels
    └── thermostat.py           # Wind-temperature thermostat
```
