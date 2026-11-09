# Validation: `src/zephyrus/fractionation.py`

This page tracks the `@pytest.mark.reference_pinned` tests that anchor the N-species closure of `zephyrus.fractionation` against every special case of the escape-fractionation lineage it generalises.

| Test id | Reference | Scope |
|---|---|---|
| `tests/test_fractionation.py::test_ternary_deuterium_reductions` | Gu & Chen (2023)[^cite-gu2023], Eqs. (4), (8), (9), and (12); Cherubim et al. (2024)[^cite-cherubim2024], ApJ 967, 139, Eq. (11) | The H-He-D system reproduces both branches and both critical rates, including the helium admixture factor on the deuterium threshold. |
| `tests/test_fractionation.py::test_two_majors_trace_minor_relations` | Odert et al. (2018)[^cite-odert2018], Icarus 307, 327, Eq. (5); Zahnle et al. (1990)[^cite-zahnle1990], Icarus 84, 502, Eqs. (35) and (36); Zahnle & Kasting (1986)[^cite-zahnle1986], Icarus 68, 462, Eq. (36) | Entrained trace minors follow the Odert relation (equal to Zahnle et al. Eq. (35)) to $10^{-11}$ relative, the closure clamps to zero where the printed formula goes negative, the limiting flux follows Zahnle et al. Eq. (36), and the earlier 1986 drag-deficit weighting is not reproduced (the adjudication between the two printed variants). |
| `tests/test_fractionation.py::test_first_entrainment_with_two_retained_heavies` | Zahnle et al. (1990), Icarus 84, 502, Eq. (42) | The first-entrainment threshold with two retained heavy backgrounds matches the printed expression and is sharp. |
| `tests/test_fractionation.py::test_zk23_nontrace_ternary_relations` | Zahnle & Kasting (2023)[^cite-zahnle2023], GeCoA 361, 228, Eqs. (19) and (20) | The non-trace H-O-CO2 relations hold: Eq. (19) in closed form to machine precision and through the solver's own bisected activation threshold, Eq. (20) at every flux inside the two-species band. |
| `tests/test_fractionation.py::test_chassefiere_prescribed_flux_partition` | Chassefière (1996)[^cite-chassefiere1996], Icarus 124, 537, Eqs. (1), (6), and (7) | The prescribed-total-flux binary partition holds to machine precision at every draw, and his crossover-mass dropout test coincides with the closure's own threshold. |
| `tests/test_fractionation.py::test_hunten_anchors_earth_mars_venus` | Hunten et al. (1987)[^cite-hunten1987], Icarus 69, 532 (worked Earth, Mars, and Venus anchors) | The printed numerical anchors are reproduced from the threshold relation and through the solver's own activation threshold within 1%. |
| `tests/test_fractionation_ensembles.py::test_binary_limit_matches_closed_form_over_random_draws` | Cherubim et al. (2024), ApJ 967, 139, Eqs. (7) to (9) | 200 random binaries match the closed-form partition on both branches to $10^{-12}$ relative, with flux continuity at the crossover and exact mass conservation. |

## Notes

The randomised ensemble suite additionally verifies global properties across activation thresholds (conservation, non-negativity, monotonicity, active-set growth, two-sided continuity, piecewise linearity), active-set uniqueness by Karush-Kuhn-Tucker enumeration against a brute-force oracle, and the stability of composition fixed points under well-mixed-layer relaxation.

## Anchor type

Published benchmarks (exact reductions to six printed formulations spanning 1987 to 2024) plus a brute-force cross-implementation oracle.

Date of last comparison against the sources: 2026-08-20.

[^cite-gu2023]: P.-G. Gu, H. Chen, *[Deuterium Escape on Photoevaporating Sub-Neptunes](https://doi.org/10.3847/2041-8213/acee01)*, The Astrophysical Journal Letters, 953, L27, 2023. [SciX](https://scixplorer.org/abs/2023ApJ...953L..27G/abstract).
[^cite-cherubim2024]: C. Cherubim, R. Wordsworth, R. Hu, E. Shkolnik, *[Strong Fractionation of Deuterium and Helium in Sub-Neptune Atmospheres along the Radius Valley](https://doi.org/10.3847/1538-4357/ad3e77)*, The Astrophysical Journal, 967, 139, 2024. [SciX](https://scixplorer.org/abs/2024ApJ...967..139C/abstract).
[^cite-odert2018]: P. Odert, H. Lammer, N. V. Erkaev, A. Nikolaou, H. I. M. Lichtenegger, C. P. Johnstone, K. G. Kislyakova, M. Leitzinger, N. Tosi, *[Escape and fractionation of volatiles and noble gases from Mars-sized planetary embryos and growing protoplanets](https://doi.org/10.1016/j.icarus.2017.10.031)*, Icarus, 307, 327-346, 2018. [SciX](https://scixplorer.org/abs/2018Icar..307..327O/abstract).
[^cite-zahnle1990]: K. Zahnle, J. F. Kasting, J. B. Pollack, *[Mass fractionation of noble gases in diffusion-limited hydrodynamic hydrogen escape](https://doi.org/10.1016/0019-1035(90)90050-J)*, Icarus, 84, 502-527, 1990. [SciX](https://scixplorer.org/abs/1990Icar...84..502Z/abstract).
[^cite-zahnle1986]: K. J. Zahnle, J. F. Kasting, *[Mass fractionation during transonic escape and implications for loss of water from Mars and Venus](https://doi.org/10.1016/0019-1035(86)90051-5)*, Icarus, 68, 462-480, 1986. [SciX](https://scixplorer.org/abs/1986Icar...68..462Z/abstract).
[^cite-zahnle2023]: K. J. Zahnle, J. F. Kasting, *[Elemental and isotopic fractionation as fossils of water escape from Venus](https://doi.org/10.1016/j.gca.2023.09.023)*, Geochimica et Cosmochimica Acta, 361, 228-244, 2023. [SciX](https://scixplorer.org/abs/2023GeCoA.361..228Z/abstract).
[^cite-chassefiere1996]: E. Chassefière, *[Hydrodynamic Escape of Oxygen from Primitive Atmospheres: Applications to the Cases of Venus and Mars](https://doi.org/10.1006/icar.1996.0229)*, Icarus, 124, 537-552, 1996. [SciX](https://scixplorer.org/abs/1996Icar..124..537C/abstract).
[^cite-hunten1987]: D. M. Hunten, R. O. Pepin, J. C. G. Walker, *[Mass fractionation in hydrodynamic escape](https://doi.org/10.1016/0019-1035(87)90022-4)*, Icarus, 69, 532-549, 1987. [SciX](https://scixplorer.org/abs/1987Icar...69..532H/abstract).
