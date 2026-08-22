# Validation: `src/zephyrus/dispatcher.py`

This page tracks the `@pytest.mark.reference_pinned` test that anchors the assembled dispatcher of `zephyrus.dispatcher`.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_dispatcher.py::test_routing_hydrodynamic_and_el_candidate_matches_el_escape` | Cross-implementation check against `zephyrus.escape.EL_escape` (Erkaev et al. 2007 form, scaling = 2) | The dispatcher's energy-limited candidate equals the released public entry point evaluated with the same efficiency, radii, flux, and tidal factor to $10^{-9}$ relative, and the dispatched rate is the min(EL, RR) winner named by the sub-label. |

## Notes

The totality contract (200 random physically posed inputs returning one label, finite non-negative rates, and per-species sums equal to the bulk rate) and the boxedness of the diagnostics container are asserted as physics invariants in the same file; every branch of the routing is exercised by dedicated scenario tests.

## Anchor type

Cross-implementation cross-check against the package's released energy-limited contract.

Date of last comparison against the sources: 2026-08-20.
