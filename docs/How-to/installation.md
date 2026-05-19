# Installation

!!! note "PROTEUS users"
    The standard way of installing ZEPHYRUS is within the PROTEUS Framework, as described in the [PROTEUS installation guide](https://proteus-framework.org/PROTEUS/installation.html#9-install-submodules-as-editable). When installed as part of PROTEUS, ZEPHYRUS is set up automatically alongside all other modules. The standalone instructions below are only needed if you want to use ZEPHYRUS independently of PROTEUS.

!!! info "Prerequisites"
    - **Python** ≥ 3.11
    - **pip** (`python -m pip --version`)
    - **Git**: only needed for the developer install (`git --version`)
    - **Internet access**: required once to download the stellar evolution tracks used by MORS
    - *(Optional)* **Conda** or **venv** : recommended to isolate the installation

---

## Standard install

ZEPHYRUS is available on [PyPI](https://pypi.org/project/fwl-zephyrus/). Install it with:

```sh
pip install fwl-zephyrus
```

ZEPHYRUS uses [`fwl-mors`](https://pypi.org/project/fwl-mors/) to compute the XUV stellar flux incident on the planet, which is the main driver of energy-limited escape. Once installed, download the required stellar evolution data:

```sh
mors download all
```

That's it. Jump to [first run](../Tutorials/first_run.md) to run your first escape computation.

---

## Developer install

Use this route if you want to modify the source code or contribute to ZEPHYRUS.

### 1. Create an isolated environment (recommended)

=== "Conda"

    ```sh
    conda create -n zephyrus python=3.11 -y
    conda activate zephyrus
    ```

=== "venv"

    ```sh
    python -m venv .venv
    source .venv/bin/activate   # Windows: .venv\Scripts\activate
    ```

### 2. Clone the repository

```sh
git clone git@github.com:FormingWorlds/ZEPHYRUS.git
cd ZEPHYRUS
```

### 3. Install in editable mode

```sh
pip install -e .
```

This installs ZEPHYRUS and all its dependencies (see [`pyproject.toml`](https://github.com/FormingWorlds/ZEPHYRUS/blob/main/pyproject.toml)) while keeping the source directory live: any edits you make are immediately reflected without reinstalling.

---

## Stellar evolution data

ZEPHYRUS relies on `fwl-mors` for the stellar XUV history, which in turn requires a set of pre-computed stellar evolution tracks. After installation, download them with:

```sh
mors download all
```

This fetches both the [Spada](https://zenodo.org/records/15729101) and [Baraffe](https://zenodo.org/records/15729114) track sets from two Zenodo records, with an automatic fallback to [OSF](https://osf.io/9u3fb/) if Zenodo is unavailable. If you only need one set, you can download them individually:

```sh
mors download spada
mors download baraffe
```

### Data location

By default, MORS stores data according to the [XDG Base Directory specification](https://specifications.freedesktop.org/basedir-spec/latest/). To check where data will be stored on your system:

```sh
mors env
```

### Changing the data directory

Set the `FWL_DATA` environment variable to redirect MORS to a different location:

```sh
export FWL_DATA=/path/to/your/data
```

To make this permanent, add that line to your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.) and reload it:

```sh
echo 'export FWL_DATA=/path/to/your/data' >> ~/.bashrc
source ~/.bashrc
```

---

## Verifying the installation

After installing, run the following to confirm everything works:

```python
from zephyrus.escape import EL_escape
from zephyrus.planets_parameters import Me, Re, Fxuv_earth_today

mdot = EL_escape(
    tidal_contribution=False,
    a=0.0, e=0.0,                  # unused when tidal_contribution=False
    Mp=Me, Ms=0.0,                 # Ms unused when tidal_contribution=False
    epsilon=0.15,
    Rp=Re, Rxuv=Re,
    Fxuv=Fxuv_earth_today,
)
print(f"Earth-like EL mass-loss rate: {mdot:.3e} kg/s")
```

You should see a mass-loss rate printed without errors. If you run into any issues, feel free to [contact the developers](../Community/contact.md) or open an issue on [GitHub](https://github.com/FormingWorlds/ZEPHYRUS/issues).