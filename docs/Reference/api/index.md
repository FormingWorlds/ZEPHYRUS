# API overview

This is an overview of ZEPHYRUS' API for the user's reference. If you want to understand the underlying model, please visit the [model overview](../../Explanations/model.md). <br>

| Module | Description |
|---|---|
| **Physical model** | |
| [`zephyrus.escape`](escape.md) | `EL_escape`: energy-limited atmospheric mass-loss rate |
| **Reference values** | |
| [`zephyrus.constants`](constants.md) | Physical constants and unit conversions (SI and CGS) |
| [`zephyrus.planets_parameters`](planet_parameters.md) | Sun, Earth, and TOI-561 system reference values |


The **source tree** is given by:

```
src/zephyrus
    ├── constants.py            # Physical constants and unit conversions
    ├── escape.py               # Energy-limited atmospheric escape (EL_escape)
    ├── __init__.py             # Package entry point
    └── planets_parameters.py   # Sun, Earth, and TOI-561 reference values
```