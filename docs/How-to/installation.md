# Installation

1. Basic install
```console
pip install fwl-zephyrus
```
2. Developer install with code sources
```console
git clone git@github.com:FormingWorlds/ZEPHYRUS.git
cd ZEPHYRUS
pip install -e .
```
3. Download input data
The `fwl-mors` python package is used to compute the XUV stellar incoming flux on the planet, which most escape models rely on. This package requires a set of stellar evolution data, stored in the [OSF repository](https://osf.io/9u3fb/). To download the data follow the following steps.
    * Set the environment variable FWL_DATA to define where the data files will be stored
        * `export FWL_DATA=...`
        * This can be done permanently by entering this line into your `~/.bashrc` file.
    * Run the following command to download all evolution track data
        * `mors download all`

### Run instructions
In the example folder, you can find python scripts showing typical usecases/workflows of escape computation with Zephyrus.
