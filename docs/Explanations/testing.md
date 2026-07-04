# Testing suite

![coverage](https://gist.githubusercontent.com/lsoucasse/152250f71914339d24537977e64aba55/raw/covbadge_zephyrus.svg)
[![Tests](https://img.shields.io/github/actions/workflow/status/FormingWorlds/ZEPHYRUS/tests.yaml?branch=main&label=Tests)](https://github.com/FormingWorlds/ZEPHYRUS/actions/workflows/tests.yaml)
[![tests](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/FormingWorlds/ZEPHYRUS/badges/tests-total.json)](https://proteus-framework.org/testing)

The badge shows the number of tests in the suite. It reads a small JSON file on the repository's `badges` branch that the `Refresh test count badges` workflow regenerates whenever the test suite or source changes on `main`, so the count stays current without hand editing. The same file backs the module's row on the central [PROTEUS testing dashboard](https://proteus-framework.org/testing).

The suite is a regression net: a passing run confirms that recent changes have not perturbed locked behaviour, not that the behaviour is itself physically correct. Physical correctness is judged separately, against analytic limits, benchmark codes, and published references.

## Badge system

`tests-total.json` sits at the root of the `badges` branch in the shields.io endpoint schema. Publishing to a dedicated branch keeps the generated file off `main`, where direct pushes need review. The `Refresh test count badges` workflow regenerates the count whenever the test suite or source changes on `main` and publishes it there; shields.io fetches it live.

For how to run the suite, see [Run tests](../How-to/run_tests.md).
