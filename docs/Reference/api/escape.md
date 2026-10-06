# zephyrus.escape

`zephyrus.escape` is the released entry point for energy-limited escape: it re-exports `EL_escape`, whose implementation lives in [`zephyrus.hydrodynamic`](hydrodynamic.md) beside the other hydrodynamic rate prescriptions. Existing code importing `from zephyrus.escape import EL_escape` keeps working unchanged.

::: zephyrus.hydrodynamic.EL_escape
    options:
      show_source: true
