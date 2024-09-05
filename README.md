[![Documentation Status](https://readthedocs.org/projects/fwl-zephyrus/badge/?version=latest)](https://fwl-zephyrus.readthedocs.io/en/latest/?badge=latest)
![Coverage](https://gist.githubusercontent.com/lsoucasse/152250f71914339d24537977e64aba55/raw/covbadge_zephyrus.svg)

# ZEPHYRUS

Zephyrus provides models to compute atmospheric escape of exoplanets.

### Documentation
https://fwl-zephyrus.readthedocs.io

### Contributors
* Emma Postolec (e.n.postolec[at]rug.nl)
* Tim Lichtenberg (tim.lichtenberg@rug.nl)
* Laurent Soucasse (l.soucasse@esciencecenter.nl)

### Repository structure
* `README.md`       - This file
* `docs/`           - Documentation sources
* `src/zephyrus/`   - Zephyrus sources
* `examples/`       - Typical use scripts
* `tests/`          - Zephyrus tests
* `tools/`          - Useful tools

### Installation instructions
1. Basic install
```console
pip install fwl-zephyrus
```
2. Developper install with code sources
```console
git clone git@github.com:FormingWorlds/ZEPHYRUS.git
cd ZEPHYRUS
pip install -e .
```
3. Download input data
The `fwl-mors` python package is used to compute the XUV stellar incoming flux on the planet, which most escape models rely on. This package requires a set of stellar evolution data, stored in the [OSF repository](https://osf.io/9u3fb/). To download the data follow the following steps.
    * Set the environment variable FWL_DATA to define where the data files will be stored
        * `export FWL_DATA=...`
        * This can be done permanently by entering this line into you `~/.profile` file.
    * Run the following command to download all evolution track data
        * `mors download all`

### Run instructions
In the examples folder you can find python scripts showing typical usecases/workflows of escape computation with Zephyrus.
