# API overview

This is an overview of ZEPHYRUS' API for the user's reference. If you want to understand the underlying model, please visit the [model overview](../../Explanations/model.md).

| Module | Description |
|---|---|
| **User-facing classes** | |
| [`mors.star`](star.md) | `Star` class: full rotation and XUV evolutionary tracks for a single star |
| [`mors.cluster`](cluster.md) | `Cluster` class: vectorised `Star` evolution for an array of stars |
| **Physical model** | |
| [`mors.physicalmodel`](physicalmodel.md) | Torque equations, X-ray/EUV/Ly-α emission, HZ boundaries, ODE right-hand side |
| [`mors.rotevo`](rotevo.md) | ODE solvers (RK4, RKF, Rosenbrock) and `FitRotation` bisection |
| **Stellar structure** | |
| [`mors.stellarevo`](stellarevo.md) | Spada et al. (2013) grid loading, interpolation, and `StarEvo` class |
| [`mors.baraffe`](baraffe.md) | Baraffe et al. (2015) track loading and interpolation |
| **Spectral synthesis** | |
| [`mors.spectrum`](spectrum.md) | `Spectrum` class: continuous stellar spectra and band integration |
| [`mors.synthesis`](synthesis.md) | Historical spectrum scaling from a modern reference spectrum |
| **Utilities** | |
| [`mors.miscellaneous`](miscellaneous.md) | `Load`, `ModelCluster`, `ActivityLifetime`, `IntegrateEmission` |
| [`mors.parameters`](parameters.md) | Default parameter dictionary, `NewParams`, `PrintParams` |
| [`mors.constants`](constants.md) | Physical constants and solar reference values (CGS) |
| [`mors.data`](data.md) | Stellar evolution track downloads (Zenodo + OSF) |
| [`mors.cli`](cli.md) | Command-line interface (`mors download`, `mors env`) |

---

## Source tree

```
src/mors
├── __init__.py           # Public API, re-exports all user-facing functions and classes
├── constants.py          # Physical constants and solar reference values (CGS units)
├── parameters.py         # Default model parameter dictionary and NewParams/PrintParams helpers
├── data.py               # Stellar evolution track downloads (Zenodo + OSF)
├── data/
│   └── ModelDistribution.dat  # Built-in 1 Myr rotation distribution (masses + Omega)
├── stellarevo.py         # Spada et al. (2013) grid loading, interpolation, and StarEvo class
├── baraffe.py            # Baraffe et al. (2015) track loading and interpolation
├── physicalmodel.py      # Torques, X-ray/EUV/Ly-α emission, HZ boundaries, ODE right-hand side
├── rotevo.py             # ODE solvers (RK4, RKF, Rosenbrock) and FitRotation bisection
├── star.py               # Star class: full rotation + XUV evolutionary tracks for one star
├── cluster.py            # Cluster class: vectorised Star evolution for an array of stars
├── spectrum.py           # Spectrum class: continuous stellar spectra and band integration
├── synthesis.py          # Historical spectrum scaling from a modern reference spectrum
├── miscellaneous.py      # Load, ModelCluster, ActivityLifetime, IntegrateEmission utilities
├── logs.py               # Logging setup (setup_logger)
└── cli.py                # Command-line interface (mors download, mors env)
```