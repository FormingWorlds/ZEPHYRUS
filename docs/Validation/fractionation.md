# Validation: `src/zephyrus/fractionation.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the N-species closure of `zephyrus.fractionation` against every special case of the escape-fractionation lineage it generalizes.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_fractionation.py::test_ternary_deuterium_reductions` | Gu & Chen (2023), Eqs. 4, 8, 9, and 12; Cherubim & Wordsworth (2024), ApJ 967, 139, Eq. 11 | The H-He-D system reproduces both branches and both critical rates exactly, including the helium admixture factor on the deuterium threshold. |
| `tests/test_fractionation.py::test_two_majors_trace_minor_relations` | Odert et al. (2018), Icarus 307, 327, Eq. 5; Zahnle et al. (1990), Icarus 84, 502, Eqs. 35 and 36; Zahnle & Kasting (1986), Icarus 68, 462, Eq. 36 | Entrained trace minors follow the Odert relation (equal to Zahnle et al. Eq. 35) to 1e-11 relative, the closure clamps exactly to zero where the printed formula goes negative, the limiting flux follows Zahnle et al. Eq. 36, and the earlier 1986 drag-deficit weighting is demonstrably NOT reproduced (the adjudication between the two printed variants). |
| `tests/test_fractionation.py::test_first_entrainment_with_two_retained_heavies` | Zahnle et al. (1990), Icarus 84, 502, Eq. 42 | The first-entrainment threshold with two retained heavy backgrounds matches the printed expression and is sharp. |
| `tests/test_fractionation.py::test_zk23_nontrace_ternary_relations` | Zahnle & Kasting (2023), GeCoA 361, 228, Eqs. 19 and 20 | The non-trace H-O-CO2 relations hold: Eq. 19 in closed form to machine precision and through the solver's own bisected activation threshold, Eq. 20 at every flux inside the two-species band. |
| `tests/test_fractionation.py::test_chassefiere_prescribed_flux_partition` | Chassefiere (1996), Icarus 124, 537, Eqs. 1, 6, and 7 | The prescribed-total-flux binary partition holds to machine precision at every draw, and his crossover-mass dropout test coincides exactly with the closure's own threshold. |
| `tests/test_fractionation.py::test_hunten_anchors_earth_mars_venus` | Hunten, Pepin & Walker (1987), Icarus 69, 532 (worked Earth, Mars, and Venus anchors) | The printed numerical anchors are reproduced from the threshold relation and through the solver's own activation threshold within 1%. |
| `tests/test_fractionation_ensembles.py::test_binary_limit_matches_closed_form_over_random_draws` | Cherubim & Wordsworth (2024), ApJ 967, 139, Eqs. 7 to 9 | 200 random binaries match the closed-form partition on both branches to 1e-12 relative, with flux continuity at the crossover and exact mass conservation. |

## Notes

The randomized ensemble suite additionally verifies global properties across activation thresholds (conservation, non-negativity, monotonicity, active-set growth, two-sided continuity, piecewise linearity), active-set uniqueness by Karush-Kuhn-Tucker enumeration against a brute-force oracle, and the stability of composition fixed points under well-mixed-layer relaxation.

## Anchor type

Published benchmarks (exact reductions to six printed formulations spanning 1987 to 2024) plus a brute-force cross-implementation oracle.

Date of last comparison against the sources: 2026-08-20.
