# Energy-limited escape

Energy-limited (EL) escape is one of the three limits on the hydrodynamic wind of the [escape-regime dispatcher](regimes.md), valid when the atmosphere sustains a collisional, XUV-driven fluid wind. It is also available on its own through the entry point `zephyrus.escape.EL_escape`, which applies it to every state without checking the regime. This page defines the prescription itself, its radius scaling, and its tidal correction.

The physical idea is an energy budget. The stellar X-ray and extreme-ultraviolet (XUV) flux absorbed high in the atmosphere heats the gas. If a fixed fraction of that power goes into lifting gas out of the planet's gravitational well, the mass-loss rate follows from dividing the absorbed power by the escape energy per unit mass[^cite-watson1981][^cite-lammer2003]. The rate is computed as

$$\dot{M}_\mathrm{EL} = \frac{\epsilon\,\pi\,R^3\,F_\mathrm{XUV}}{G\,M_\mathrm{p}\,K_\mathrm{tide}} \tag{1}$$

where $\epsilon$ is the escape efficiency (the fraction of intercepted XUV power converted into work against gravity), $R^3$ is the radius term selected below, $F_\mathrm{XUV}$ is the XUV flux received at the planet's orbit (supplied by [MORS](https://proteus-framework.org/MORS) in coupled runs), $G$ is the gravitational constant, $M_\mathrm{p}$ is the planetary mass, and $K_\mathrm{tide}$ is the tidal correction of Eq. (2).

Canonical efficiencies for rocky planets lie between 0.1 and 0.3. The dispatcher's settings accept any $\epsilon \in (0, 1]$ and reject values outside it, while `EL_escape` evaluates the formula for whatever efficiency it is given. Hydrodynamic simulations find that for strongly bound planets the effective efficiency collapses far below the canonical band, reaching of order $10^{-2}$ above a threshold gravitational potential near $\log_{10}(G M_\mathrm{p} / R_\mathrm{p}) \approx 12.9$ to $13.2$ in cgs units[^cite-caldiroli2022]. The threshold is on the untidal specific binding energy, which is the convention the screen in the diagnostics uses and the one Salz et al. (2016)[^cite-salz2016] share, and the tidal factor enters the rate, not the potential the threshold is quoted against. The dispatcher offers that fitted efficiency as an option, while the entry point itself treats $\epsilon$ as an input.

Two radii enter the problem. $R_\mathrm{p}$ is the planetary (interior) radius. $R_\mathrm{XUV}$ is the radius at which the atmosphere becomes optically thick to XUV photons. In PROTEUS it is recomputed at each time step from the atmospheric structure at the pressure set by `escape.zephyrus.Pxuv` (default $5 \times 10^{-5}$ bar). The dispatcher takes it instead at its photospheric working level, 20 mbar by default, the photosphere-type level of Baumeister et al. (2023)[^cite-baumeister2023]. Either level is a bookkeeping radius for the intercepting area: the XUV heating is actually deposited, and the wind launched, at the far lower pressure (of order a nanobar) where the gas first reaches unit optical depth to ionising photons[^cite-lopez2017], and the [escape-regime dispatcher](regimes.md) locates that launching level on the profile when it needs it.

## Radius scaling

The cubic radius term in Eq. (1) is selected at runtime by the `scaling` argument of `EL_escape`:

| `scaling` | $R^3$ | Description |
|---|---|---|
| `2` | $R_\mathrm{p}\,R^2_\mathrm{XUV}$ | Default: the XUV-absorbing cross section paired with the potential measured at the surface radius, the form of Erkaev et al. (2007)[^cite-erkaev2007] |
| `3` | $R^3_\mathrm{XUV}$ | All three powers taken at the XUV radius, the single-radius form of Lopez et al. (2012)[^cite-lopez2012] and Lehmer & Catling (2017)[^cite-lehmer2017] |

Both forms reduce to $R_\mathrm{p}^3$ when $R_\mathrm{XUV} = R_\mathrm{p}$, the conservative lower bound adopted by Luger & Barnes (2015)[^cite-luger2015] and Moore et al. (2023)[^cite-moore2023]. Allowing $R_\mathrm{XUV} > R_\mathrm{p}$ increases the effective XUV-absorbing area and therefore the escape rate.

## Tidal correction

When the `tidal_contribution` flag is `True`, the effective potential barrier is reduced by the host star's tidal field, following Erkaev et al. (2007), their Eq. (17)[^cite-erkaev2007]:

$$K_\mathrm{tide} = 1 - \frac{3}{2\xi} + \frac{1}{2\xi^3}, \qquad \xi = \frac{R_\mathrm{Hill}}{R} \tag{2}$$

with the periapsis Hill radius

$$R_\mathrm{Hill} = a\,(1-e)\,\left(\frac{M_\mathrm{p}}{3\,M_\star}\right)^{1/3} \tag{3}$$

where $a$ is the semi-major axis, $e$ the orbital eccentricity, and $M_\star$ the stellar mass. The radius $R$ in $\xi$ is the one that appears linearly in the $R^3$ term of Eq. (1): $R_\mathrm{p}$ for `scaling=2` (the convention of Erkaev et al. 2007, whose own $\xi$ is the Roche-lobe distance over the planetary radius) and $R_\mathrm{XUV}$ for `scaling=3`, where the XUV radius is the only radius in the problem.

Factoring the numerator gives $K_\mathrm{tide} = (\xi - 1)^2\,(2\xi + 1) / (2\xi^3)$, non-negative for every $\xi > 0$ with a double root at $\xi = 1$. In the physical regime $\xi > 1$ it lies in $(0, 1)$, rising toward 1 for $\xi \gg 1$ (the escape level deep inside the Hill sphere) and falling toward 0 as the atmosphere expands toward the Roche lobe at $\xi = 1$. Because the rate divides by $K_\mathrm{tide}$, the tidal correction enhances escape and diverges at the Roche lobe, so the tidally corrected rate is defined only for $\xi > 1$: `EL_escape` raises a `ValueError` for $\xi \le 1$, where the atmosphere reaches the Roche lobe and the energy-limited approximation no longer applies. When `tidal_contribution` is `False`, $K_\mathrm{tide} = 1$. The [dispatcher](regimes.md) handles the $\xi \le 1$ geometry itself, as the `roche_overflow` label.

## When the prescription applies

The EL form is appropriate in the high-irradiation, collisional-wind regime that dominates the loss during roughly the first $10^6$ to $10^8$ years of a close-in planet's evolution[^cite-watson1981][^cite-lammer2003]. Outside it, at lower XUV flux or for a less extended atmosphere, particle-by-particle (nonthermal and Jeans) escape becomes comparable or dominant and the bulk EL prescription no longer applies. The [dispatcher](regimes.md) classifies each state before choosing a rate. When using `EL_escape` alone, verify that the XUV-driven loss genuinely dominates (for scale, present-day nonthermal rates for an Earth-mass planet are of order $10^4$ to $10^5$ kg s$^{-1}$, Kislyakova et al. 2014[^cite-kislyakova2014]).

Bulk removal is the second assumption. When the escaping particle flux drops below the critical value needed to drag a heavy species along, the outflow fractionates: hydrogen escapes preferentially and the residual atmosphere is enriched in heavy species[^cite-wordsworth2018][^cite-cherubim2024]. For scale, the critical flux for water in a hydrogen background is about $1.9 \times 10^{5}$ kg s$^{-1}$ (Yoshida et al. 2022[^cite-yoshida2022]). `EL_escape` removes everything in bulk, and the [fractionation](fractionation.md) page describes the closure that resolves the partition when the regime framework confirms a wind.

[^cite-watson1981]: A. J. Watson, T. M. Donahue, J. C. G. Walker, *[The dynamics of a rapidly escaping atmosphere: Applications to the evolution of Earth and Venus](https://doi.org/10.1016/0019-1035(81)90101-9)*, Icarus, 48, 150-166, 1981. [SciX](https://scixplorer.org/abs/1981Icar...48..150W/abstract).

[^cite-lammer2003]: H. Lammer, F. Selsis, I. Ribas, E. F. Guinan, S. J. Bauer, W. W. Weiss, *[Atmospheric Loss of Exoplanets Resulting from Stellar X-Ray and Extreme-Ultraviolet Heating](https://doi.org/10.1086/380815)*, The Astrophysical Journal Letters, 598, L121-L124, 2003. [SciX](https://scixplorer.org/abs/2003ApJ...598L.121L/abstract).

[^cite-erkaev2007]: N. V. Erkaev, Yu. N. Kulikov, H. Lammer, F. Selsis, D. Langmayr, G. F. Jaritz, H. K. Biernat, *[Roche lobe effects on the atmospheric loss from “Hot Jupiters”](https://doi.org/10.1051/0004-6361:20066929)*, Astronomy & Astrophysics, 472, 329-334, 2007. [SciX](https://scixplorer.org/abs/2007A%26A...472..329E/abstract).

[^cite-lopez2012]: E. D. Lopez, J. J. Fortney, N. Miller, *[How Thermal Evolution and Mass-loss Sculpt Populations of Super-Earths and Sub-Neptunes: Application to the Kepler-11 System and Beyond](https://doi.org/10.1088/0004-637X/761/1/59)*, The Astrophysical Journal, 761, 59, 2012. [SciX](https://scixplorer.org/abs/2012ApJ...761...59L/abstract).

[^cite-lehmer2017]: O. R. Lehmer, D. C. Catling, *[Rocky Worlds Limited to ∼1.8 Earth Radii by Atmospheric Escape during a Star’s Extreme UV Saturation](https://doi.org/10.3847/1538-4357/aa8137)*, The Astrophysical Journal, 845, 130, 2017. [SciX](https://scixplorer.org/abs/2017ApJ...845..130L/abstract).

[^cite-lopez2017]: E. D. Lopez, *[Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor](https://doi.org/10.1093/mnras/stx1558)*, Monthly Notices of the Royal Astronomical Society, 472, 245-253, 2017. [SciX](https://scixplorer.org/abs/2017MNRAS.472..245L/abstract).

[^cite-baumeister2023]: P. Baumeister, N. Tosi, C. Brachmann, J. L. Grenfell, L. Noack, *[Redox state and interior structure control on the long-term habitability of stagnant-lid planets](https://doi.org/10.1051/0004-6361/202245791)*, Astronomy & Astrophysics, 675, A122, 2023. [SciX](https://scixplorer.org/abs/2023A%26A...675A.122B/abstract).

[^cite-caldiroli2022]: A. Caldiroli, F. Haardt, E. Gallo, R. Spinelli, I. Malsky, E. Rauscher, *[Irradiation-driven escape of primordial planetary atmospheres. II. Evaporation efficiency of sub-Neptunes through hot Jupiters](https://doi.org/10.1051/0004-6361/202142763)*, Astronomy & Astrophysics, 663, A122, 2022. [SciX](https://scixplorer.org/abs/2022A%26A...663A.122C/abstract).

[^cite-salz2016]: M. Salz, P. C. Schneider, S. Czesla, J. H. M. M. Schmitt, *[Energy-limited escape revised. The transition from strong planetary winds to stable thermospheres](https://doi.org/10.1051/0004-6361/201527042)*, Astronomy & Astrophysics, 585, L2, 2016. [SciX](https://scixplorer.org/abs/2016A%26A...585L...2S/abstract).

[^cite-luger2015]: R. Luger, R. Barnes, *[Extreme Water Loss and Abiotic O₂ Buildup on Planets Throughout the Habitable Zones of M Dwarfs](https://doi.org/10.1089/ast.2014.1231)*, Astrobiology, 15, 119-143, 2015. [SciX](https://scixplorer.org/abs/2015AsBio..15..119L/abstract).

[^cite-moore2023]: K. Moore, N. B. Cowan, C.-É. Boukaré, *[The role of magma oceans in maintaining surface water on rocky planets orbiting M-dwarfs](https://doi.org/10.1093/mnras/stad3138)*, Monthly Notices of the Royal Astronomical Society, 526, 6235-6247, 2023. [SciX](https://scixplorer.org/abs/2023MNRAS.526.6235M/abstract).

[^cite-kislyakova2014]: K. G. Kislyakova, C. P. Johnstone, P. Odert, N. V. Erkaev, H. Lammer, T. Lüftinger, M. Holmström, M. L. Khodachenko, M. Güdel, *[Stellar wind interaction and pick-up ion escape of the Kepler-11 “super-Earths”](https://doi.org/10.1051/0004-6361/201322933)*, Astronomy & Astrophysics, 562, A116, 2014. [SciX](https://scixplorer.org/abs/2014A%26A...562A.116K/abstract).

[^cite-wordsworth2018]: R. D. Wordsworth, L. K. Schaefer, R. A. Fischer, *[Redox Evolution via Gravitational Differentiation on Low-mass Planets: Implications for Abiotic Oxygen, Water Loss, and Habitability](https://doi.org/10.3847/1538-3881/aab608)*, The Astronomical Journal, 155, 195, 2018. [SciX](https://scixplorer.org/abs/2018AJ....155..195W/abstract).

[^cite-cherubim2024]: C. Cherubim, R. Wordsworth, R. Hu, E. Shkolnik, *[Strong Fractionation of Deuterium and Helium in Sub-Neptune Atmospheres along the Radius Valley](https://doi.org/10.3847/1538-4357/ad3e77)*, The Astrophysical Journal, 967, 139, 2024. [SciX](https://scixplorer.org/abs/2024ApJ...967..139C/abstract).

[^cite-yoshida2022]: T. Yoshida, N. Terada, M. Ikoma, K. Kuramoto, *[Less Effective Hydrodynamic Escape of H₂-H₂O Atmospheres on Terrestrial Planets Orbiting Pre-main-sequence M Dwarfs](https://doi.org/10.3847/1538-4357/ac7be7)*, The Astrophysical Journal, 934, 137, 2022. [SciX](https://scixplorer.org/abs/2022ApJ...934..137Y/abstract).
