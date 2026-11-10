# The ZEPHYRUS model

ZEPHYRUS computes the atmospheric mass loss of rocky and sub-Neptune exoplanets. It runs standalone or as the escape module of the [PROTEUS](https://proteus-framework.org) coupled atmosphere and interior framework, where it is called at every time step with the current planetary state and returns the mass-loss rate that depletes the volatile inventory.

!!! info "Two entry points"
    The escape model of ZEPHYRUS is the regime dispatcher, `zephyrus.dispatch`, which decides which escape physics applies to a state before computing a rate. The older entry point `zephyrus.escape.EL_escape` applies one prescription, the energy-limited rate, to every state without that decision. Energy-limited escape is also one of the dispatcher's regimes, `hydrodynamic:EL`, so the name appears in both places.

Mass loss happens through two physically distinct channels, and ZEPHYRUS models both:

1. **Continuous thermal escape**: the steady outflow or evaporation of the upper atmosphere, driven by stellar irradiation and by the planet's own heat. This is one model with several regimes, described below.
2. **Impact-driven erosion**: the impulsive removal of atmosphere by a single giant collision during accretion, a separate channel with its own prescription (see [giant impacts](impacts.md)).

## Continuous escape: one model, six regimes

Which physics carries the continuous loss depends on how tightly the atmosphere is bound, how strongly it is irradiated, and how collisional its outer layers are. Applying a prescription outside its regime gives rates that are wrong by orders of magnitude, so the dispatcher decides the regime of each atmospheric state before choosing a rate. Every state receives one of six regime labels:

- `boiloff`: the atmosphere is so loosely bound that it flows out on the planet's own thermal energy, faster than stellar XUV heating could drive it. This is typical of young, hot, hydrogen-rich planets fresh out of the nebula. The rate is the isothermal Parker wind of Owen & Wu (2016)[^cite-owenwu2016] at the wind temperature of Misener et al. (2025)[^cite-misener2025], and the state is recognised by the restricted Jeans parameter of Fossati et al. (2017)[^cite-fossati2017].
- `hydrodynamic:EL`: a fluid wind driven by stellar XUV heating, with the rate set by the energy budget (Watson et al. 1981[^cite-watson1981], in the tidally corrected form of Erkaev et al. 2007[^cite-erkaev2007]). The energy-limited rate is the smallest of the three hydrodynamic limits here.
- `hydrodynamic:RR`: the same fluid wind, with the rate held below the energy limit because the absorbed energy goes into ionisation that recombination re-radiates (the recombination-limited rate of Murray-Clay et al. 2009[^cite-murrayclay2009] is the smallest).
- `hydrodynamic:PL`: the same fluid wind in a well shallow enough that each ionising photon carries more energy than lifting one particle out costs, so the number of photons rather than their energy sets the rate (the photon-limited rate of Owen & Alvarez 2016[^cite-owenalvarez2016] is the smallest).
- `hydrostatic`: the gas is too rarefied to sustain a fluid wind, and escape proceeds particle by particle from the exosphere (Jeans escape, with the kinetic enhancement of Volkov et al. 2011b[^cite-volkov2011b]), species by species, capped by how fast diffusion resupplies each species (Yelle 2024[^cite-yelle2024]). Whether the wind is collisional enough to exist is decided by the Knudsen number at its sonic point, after Chatterjee & Pierrehumbert (2026)[^cite-chatterjee2026].
- `roche_overflow`: either the tidally driven transfer through the inner Lagrange point (Jackson et al. 2017[^cite-jackson2017]) is the largest rate and is dispatched, or the flow reaches the planet's Hill sphere and the label replaces that of the regime that produced the rate, which is then a lower limit (the overflow geometries of Owen & Jackson 2012[^cite-owen2012]). `diagnostics['roche']['rate_branch']` says which reading applies.

The logic reduces to four questions:

```mermaid
flowchart TD
    IN(["Planet state + atmosphere profile"]) --> Q1{"Is the atmosphere<br/>loosely bound?<br/>(restricted Jeans parameter<br/>below threshold)"}
    Q1 -- yes --> QB{"Does the XUV wind<br/>outrun the bolometric one?"}
    QB -- yes --> HD
    QB -- no --> BO["BOIL-OFF<br/>bolometric wind"]
    Q1 -- no --> Q2{"Is an XUV wind collisional<br/>at its sonic point?<br/>(Knudsen number below<br/>threshold)"}
    Q2 -- yes --> HD["HYDRODYNAMIC WIND<br/>label EL, RR, or PL,<br/>whichever rate is smallest"]
    Q2 -- no --> HS["HYDROSTATIC<br/>per-species Jeans escape"]
    BO --> Q3{"Does the flow reach<br/>the Hill sphere?"}
    HD --> Q3
    HS --> Q3
    Q3 -- yes --> RO["ROCHE OVERFLOW"]
    Q3 -- no --> OUT(["Regime label + bulk rate<br/>+ per-species rates"])
    RO --> OUT
    class BO,HD,HS,RO regime
    class Q1,QB,Q2,Q3 decision
```

The figure shows the logic, not the full machinery. Each branch carries its own rate physics, limits, and checks, and three refinements are left out:

- An unstable exosphere lets the wind rate compete with the Jeans rate.
- A setting can let the bolometric wind compete after boil-off ends.
- The transfer through L1 competes wherever it applies, and takes the overflow label with its own rate when it wins.

The [escape regimes](regimes.md) section walks through every step with its equations and thresholds, and its admission table is the reference for the routing.

The regime boundaries are not sharp lines in nature. Each threshold carries a physical band: the collisionality threshold spans a factor of 30 across heating geometries, and the boil-off threshold a factor of about two across the literature. Beside every verdict, ZEPHYRUS reports how close the state sat to each boundary and what the label would have been at the edges of the band.

Whichever regime sets the bulk rate, the loss is also split over chemical species. Hydrodynamic winds fractionate: heavy species lag the outflow through diffusive drag and can drop out of it entirely, which the N-species closure of the [fractionation](fractionation.md) page resolves (Attia & Lichtenberg 2026[^cite-attia2026]). The other regimes split the rate by reservoir mass fractions, and the hydrostatic regime is natively per-species.

## The energy-limited entry point

`zephyrus.escape.EL_escape` is the energy-limited rate on its own, with its two radius scalings and its tidal correction, defined on the [energy-limited escape](energy_limited.md) page. It does not check whether the state sustains a collisional XUV wind, which is the regime where the rate applies and which dominates the loss during roughly the first $10^6$ to $10^8$ years of a close-in planet's life. A study that wants one prescription applied everywhere can use it.

## The impact channel

A giant collision removes part of the target's atmosphere in a single event, on a timescale unrelated to the continuous escape. ZEPHYRUS computes the eroded fraction with the scaling law of Kegerreis et al. (2020)[^cite-kegerreis2020] for thin atmospheres or that of Roche et al. (2026)[^cite-roche2026] for envelopes of 1% to 20% of the planet mass, selected through `zephyrus.collision.impact_loss`. The [giant impacts](impacts.md) page defines both laws and their fitted domains. The impact channel is not part of the continuous regimes: a label is reserved for it, and the caller applies impact erosion as a discrete event.

## Reading guide

- [Escape regimes](regimes.md): how the dispatcher chooses a rate, the notation, and one page per regime with its equations, thresholds, and bands.
- [Dispatching a regime](../Tutorials/dispatch.md): the dispatcher driven end to end on synthetic atmospheres, boundary crossings included.
- [Fractionation](fractionation.md): how a wind splits over species, and when heavy species drop out.
- [Energy-limited escape](energy_limited.md): the single-prescription entry point, its radius scaling, and its tidal correction.
- [Giant impacts](impacts.md): the erosion scaling laws and their fitted domains.
- [Coupling to PROTEUS](proteus.md): configuration keys, the per-time-step sequence, and reservoir bookkeeping.
- [Limitations](limitations.md): what each entry point does not model, and what that implies for results.
- [Troubleshooting the dispatcher](../How-to/troubleshooting.md): what to look at when a flag fires or a label surprises you.
- [Parameter reference](../Reference/parameters.md), [dispatch results](../Reference/results.md), and [API reference](../Reference/api/index.md).

[^cite-owenwu2016]: J. E. Owen, Y. Wu, *[Atmospheres of Low-mass Planets: The “Boil-off”](https://doi.org/10.3847/0004-637X/817/2/107)*, The Astrophysical Journal, 817, 107, 2016. [SciX](https://scixplorer.org/abs/2016ApJ...817..107O/abstract).

[^cite-misener2025]: W. Misener, M. Schulik, H. E. Schlichting, J. E. Owen, *[Blowin’ in the Nonisothermal Wind: Core-powered Mass Loss with Hydrodynamic Radiative Transfer](https://doi.org/10.3847/1538-4357/ada777)*, The Astrophysical Journal, 980, 152, 2025. [SciX](https://scixplorer.org/abs/2025ApJ...980..152M/abstract).

[^cite-fossati2017]: L. Fossati, N. V. Erkaev, H. Lammer, P. E. Cubillos, P. Odert, I. Juvan, K. G. Kislyakova, M. Lendl, D. Kubyshkina, S. J. Bauer, *[Aeronomical constraints to the minimum mass and maximum radius of hot low-mass planets](https://doi.org/10.1051/0004-6361/201629716)*, Astronomy & Astrophysics, 598, A90, 2017. [SciX](https://scixplorer.org/abs/2017A%26A...598A..90F/abstract).

[^cite-watson1981]: A. J. Watson, T. M. Donahue, J. C. G. Walker, *[The dynamics of a rapidly escaping atmosphere: Applications to the evolution of Earth and Venus](https://doi.org/10.1016/0019-1035(81)90101-9)*, Icarus, 48, 150-166, 1981. [SciX](https://scixplorer.org/abs/1981Icar...48..150W/abstract).

[^cite-erkaev2007]: N. V. Erkaev, Yu. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).

[^cite-murrayclay2009]: R. A. Murray-Clay, E. I. Chiang, N. Murray, *[Atmospheric Escape From Hot Jupiters](https://doi.org/10.1088/0004-637X/693/1/23)*, The Astrophysical Journal, 693, 23-42, 2009. [SciX](https://scixplorer.org/abs/2009ApJ...693...23M/abstract).

[^cite-owenalvarez2016]: J. E. Owen, M. A. Alvarez, *[UV Driven Evaporation of Close-in Planets: Energy-limited, Recombination-limited, and Photon-limited Flows](https://doi.org/10.3847/0004-637X/816/1/34)*, The Astrophysical Journal, 816, 34, 2016. [SciX](https://scixplorer.org/abs/2016ApJ...816...34O/abstract).

[^cite-volkov2011b]: A. N. Volkov, O. J. Tucker, J. T. Erwin, R. E. Johnson, *[Kinetic simulations of thermal escape from a single component atmosphere](https://doi.org/10.1063/1.3592253)*, Physics of Fluids, 23, 066601, 2011. [SciX](https://scixplorer.org/abs/2011PhFl...23f6601V/abstract).

[^cite-yelle2024]: R. V. Yelle, *[Diffusion limited escape of hydrogen from Mars](https://doi.org/10.1016/j.icarus.2024.116099)*, Icarus, 416, 116099, 2024. [SciX](https://scixplorer.org/abs/2024Icar..41616099Y/abstract).

[^cite-chatterjee2026]: R. D. Chatterjee, R. T. Pierrehumbert, *[Novel Physics of Escaping Secondary Atmospheres May Shape the Cosmic Shoreline](https://doi.org/10.3847/1538-4357/ae2ffa)*, The Astrophysical Journal, 998, 236, 2026. [SciX](https://scixplorer.org/abs/2026ApJ...998..236C/abstract).

[^cite-jackson2017]: B. Jackson, P. Arras, K. Penev, S. Peacock, P. Marchant, *[A New Model of Roche Lobe Overflow for Short-period Gaseous Planets and Binary Stars](https://doi.org/10.3847/1538-4357/835/2/145)*, The Astrophysical Journal, 835, 145, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...835..145J/abstract).

[^cite-owen2012]: J. E. Owen, A. P. Jackson, *[Planetary evaporation by UV & X-ray radiation: basic hydrodynamics](https://doi.org/10.1111/j.1365-2966.2012.21481.x)*, Monthly Notices of the Royal Astronomical Society, 425, 2931-2947, 2012. [SciX](https://scixplorer.org/abs/2012MNRAS.425.2931O/abstract).

[^cite-attia2026]: M. Attia, T. Lichtenberg, *[Atmospheric escape fractionates secondary but not primary atmospheres](https://doi.org/10.48550/arXiv.2608.30106)*, arXiv e-prints, arXiv:2608.30106, 2026. [SciX](https://scixplorer.org/abs/2026arXiv260830106A/abstract).

[^cite-kegerreis2020]: J. A. Kegerreis, V. R. Eke, D. C. Catling, R. J. Massey, L. F. A. Teodoro, K. J. Zahnle, *[Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density](https://doi.org/10.3847/2041-8213/abb5fb)*, The Astrophysical Journal Letters, 901, L31, 2020. [SciX](https://scixplorer.org/abs/2020ApJ...901L..31K/abstract).

[^cite-roche2026]: M. J. Roche, S. J. Lock, P. J. Carter, Z. M. Leinhardt, *[Giant Impacts Preferentially Remove Low-mass Atmospheres: A Generalised Scaling Law for Impact-driven Atmospheric Loss](https://doi.org/10.3847/2041-8213/aeb26a)*, The Astrophysical Journal Letters, 1011, L52, 2026. [SciX](https://scixplorer.org/abs/2026ApJ..1011L..52R/abstract). Dataset: [Zenodo](https://doi.org/10.5281/zenodo.23192423).
