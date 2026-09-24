# Energy-limited escape

Energy-limited (EL) escape is the default rate prescription of ZEPHYRUS: the one PROTEUS consumes at each coupled time step, through the released entry point `zephyrus.escape.EL_escape`. Within the [escape-regime framework](regimes.md) it is one of the two hydrodynamic limits, valid when the atmosphere sustains a collisional, XUV-driven fluid wind; this page defines the prescription itself, its radius scaling, and its tidal correction.

The physical idea is an energy budget. The stellar X-ray and extreme-ultraviolet (XUV) flux absorbed high in the atmosphere heats the gas; if a fixed fraction of that power goes into lifting gas out of the planet's gravitational well, the mass-loss rate follows from dividing the absorbed power by the escape energy per unit mass [^watson][^lammer2003]. The rate is computed as

$$\dot{M}_\mathrm{EL} = \frac{\epsilon\,\pi\,R^3\,F_\mathrm{XUV}}{G\,M_\mathrm{p}\,K_\mathrm{tide}} \tag{1}$$

where $\epsilon$ is the escape efficiency (the fraction of intercepted XUV power converted into work against gravity), $R^3$ is the radius term selected below, $F_\mathrm{XUV}$ is the XUV flux received at the planet's orbit (supplied by [MORS](https://proteus-framework.org/MORS) in coupled runs), $G$ is the gravitational constant, $M_\mathrm{p}$ is the planetary mass, and $K_\mathrm{tide}$ is the tidal correction of Eq. (2). Canonical efficiencies for rocky planets lie between 0.1 and 0.3, and ZEPHYRUS accepts any $\epsilon \in (0, 1]$; hydrodynamic simulations find that for strongly bound planets the effective efficiency collapses far below the canonical band, reaching of order $10^{-2}$ above a threshold gravitational potential near $\log_{10}(G M_\mathrm{p} / R_\mathrm{p}) \approx 12.9$ to $13.2$ in cgs units [^caldiroli]. The threshold is on the untidal specific binding energy, which is the convention the screen in the diagnostics uses and the one Salz et al. (2016) share; the tidal factor enters the rate, not the potential the threshold is quoted against. The regime framework offers that fitted efficiency as an option; the entry point itself treats $\epsilon$ as an input.

Two radii enter the problem, and keeping them apart matters. $R_\mathrm{p}$ is the planetary (interior) radius. $R_\mathrm{XUV}$ is the radius at which the atmosphere becomes optically thick to XUV photons; in PROTEUS it is recomputed at each time step from the atmospheric structure at a configured reference pressure, by default 20 mbar following the photosphere-type level of Baumeister et al. (2023) [^baumeister]. That level is a bookkeeping radius for the intercepting area: the XUV heating is actually deposited, and the wind launched, at the far lower pressure (of order a nanobar) where the gas first reaches unit optical depth to ionizing photons [^lopez2017]; the [escape-regime framework](regimes.md) locates that launching level on the profile when it needs it.

## Radius scaling

The cubic radius term in Eq. (1) is selected at runtime by the `scaling` argument of `EL_escape`:

| `scaling` | $R^3$ | Description |
|---|---|---|
| `2` | $R_\mathrm{p}\,R^2_\mathrm{XUV}$ | Default; the XUV-absorbing cross section paired with the potential measured at the surface radius, the form of Erkaev et al. (2007) [^erkaev] |
| `3` | $R^3_\mathrm{XUV}$ | All three powers taken at the XUV radius, the single-radius form of Lopez, Fortney & Miller (2012) [^lfm2012] and Lehmer & Catling (2017) [^lehmer] |

Both forms reduce to $R_\mathrm{p}^3$ when $R_\mathrm{XUV} = R_\mathrm{p}$, the conservative lower bound adopted by Luger & Barnes (2015) [^luger] and Moore et al. (2023) [^moore]. Allowing $R_\mathrm{XUV} > R_\mathrm{p}$ increases the effective XUV-absorbing area and therefore the escape rate.

## Tidal correction

When the `tidal_contribution` flag is `True`, the effective potential barrier is reduced by the host star's tidal field, following Erkaev et al. (2007), their Eq. 17 [^erkaev]:

$$K_\mathrm{tide} = 1 - \frac{3}{2\xi} + \frac{1}{2\xi^3}, \qquad \xi = \frac{R_\mathrm{Hill}}{R} \tag{2}$$

with the periapsis Hill radius

$$R_\mathrm{Hill} = a\,(1-e)\,\left(\frac{M_\mathrm{p}}{3\,M_\star}\right)^{1/3} \tag{3}$$

where $a$ is the semi-major axis, $e$ the orbital eccentricity, and $M_\star$ the stellar mass. The radius $R$ in $\xi$ is the one that appears linearly in the $R^3$ term of Eq. (1): $R_\mathrm{p}$ for `scaling=2` (the convention of Erkaev et al. 2007, whose own $\xi$ is the Roche-lobe distance over the planetary radius) and $R_\mathrm{XUV}$ for `scaling=3`, where the XUV radius is the only radius in the problem.

Factoring the numerator gives $K_\mathrm{tide} = (\xi - 1)^2\,(2\xi + 1) / (2\xi^3)$, non-negative for every $\xi > 0$ with a double root at $\xi = 1$. In the physical regime $\xi > 1$ it lies in $(0, 1)$, rising toward 1 for $\xi \gg 1$ (the escape level deep inside the Hill sphere) and falling toward 0 as the atmosphere expands toward the Roche lobe at $\xi = 1$. Because the rate divides by $K_\mathrm{tide}$, the tidal correction enhances escape and diverges at the Roche lobe, so the tidally corrected rate is defined only for $\xi > 1$: `EL_escape` raises a `ValueError` for $\xi \le 1$, where the atmosphere reaches the Roche lobe and the energy-limited approximation no longer applies. When `tidal_contribution` is `False`, $K_\mathrm{tide} = 1$. The [regime framework](regimes.md) handles the $\xi \le 1$ geometry itself, as the `roche_overflow` label.

## When the prescription applies

The EL form is appropriate in the high-irradiation, collisional-wind regime that dominates the loss during roughly the first $10^6$ to $10^8$ years of a close-in planet's evolution [^watson][^lammer2003]. Outside it, at lower XUV flux or for a less extended atmosphere, particle-by-particle (nonthermal and Jeans) escape becomes comparable or dominant and the bulk EL prescription no longer applies. The [regime framework](regimes.md) classifies each state before choosing a rate; when using `EL_escape` alone, verify that the XUV-driven loss genuinely dominates (for scale, present-day nonthermal rates for an Earth-mass planet are of order $10^7$ to $10^8$ g s$^{-1}$; Kislyakova et al. 2014 [^kislyakova]).

Bulk removal is the second assumption. When the escaping particle flux drops below the critical value needed to drag a heavy species along, the outflow fractionates: hydrogen escapes preferentially and the residual atmosphere is enriched in heavy species [^wordsworth2018][^cherubim2024]. For scale, the critical flux for water in a hydrogen background is about $1.9 \times 10^{8}$ g s$^{-1}$ (Yoshida et al. 2022 [^yoshida]). `EL_escape` removes everything in bulk; the [fractionation](fractionation.md) page describes the closure that resolves the partition when the regime framework confirms a wind.

---

[^watson]: Watson, A. J., Donahue, T. M., & Walker, J. C. G. (1981). The dynamics of a rapidly escaping atmosphere: applications to the evolution of Earth and Venus. *Icarus, 48*(2), 150–166. https://doi.org/10.1016/0019-1035(81)90101-9

[^lammer2003]: Lammer, H., Selsis, F., Ribas, I., et al. (2003). Atmospheric loss of exoplanets resulting from stellar X-ray and extreme-ultraviolet heating. *The Astrophysical Journal, 598*(2), L121–L124. https://doi.org/10.1086/380815

[^erkaev]: Erkaev, N. V., Kulikov, Y. N., Lammer, H., et al. (2007). Roche lobe effects on the atmospheric loss from "Hot Jupiters". *Astronomy & Astrophysics, 472*(1), 329–334. https://doi.org/10.1051/0004-6361:20066929

[^lfm2012]: Lopez, E. D., Fortney, J. J., & Miller, N. (2012). How thermal evolution and mass-loss sculpt populations of super-Earths and sub-Neptunes: application to the Kepler-11 system and beyond. *The Astrophysical Journal, 761*(1), 59.

[^lehmer]: Lehmer, O. R., & Catling, D. C. (2017). Rocky worlds limited to 1.8 Earth radii by atmospheric escape during a star's extreme UV saturation. *The Astrophysical Journal, 845*(2), 130.

[^lopez2017]: Lopez, E. D. (2017). Born dry in the photoevaporation desert: Kepler's ultra-short-period planets formed water-poor. *Monthly Notices of the Royal Astronomical Society, 472*(1), 245–253.

[^baumeister]: Baumeister, P., Tosi, N., Brachmann, C., Grenfell, J. L., & Noack, L. (2023). Redox state and interior structure control on the long-term habitability of stagnant-lid planets. *Astronomy & Astrophysics, 675*, A122. https://doi.org/10.1051/0004-6361/202245791

[^caldiroli]: Caldiroli, A., Haardt, F., Gallo, E., Spinelli, R., Malsky, I., & Rauscher, E. (2022). Irradiation-driven escape of primordial planetary atmospheres II. Evaporation efficiency of sub-Neptunes through hot Jupiters. *Astronomy & Astrophysics, 663*, A122. https://doi.org/10.1051/0004-6361/202142763

[^luger]: Luger, R., & Barnes, R. (2015). Extreme water loss and abiotic O$_2$ buildup on planets throughout the habitable zones of M dwarfs. *Astrobiology, 15*(2), 119–143. https://doi.org/10.1089/ast.2014.1231

[^moore]: Moore, K., Cowan, N. B., & Boukaré, C.-É. (2023). The role of magma oceans in maintaining surface water on rocky planets orbiting M-dwarfs. *Monthly Notices of the Royal Astronomical Society, 526*(4), 6235–6249. https://doi.org/10.1093/mnras/stad3138

[^kislyakova]: Kislyakova, K. G., Johnstone, C. P., Odert, P., et al. (2014). Stellar wind interaction and pick-up ion escape of the Kepler-11 "super-Earths". *Astronomy & Astrophysics, 562*, A116. https://doi.org/10.1051/0004-6361/201322933

[^wordsworth2018]: Wordsworth, R. D., Schaefer, L. K., & Fischer, R. A. (2018). Redox evolution via gravitational differentiation on low-mass planets: implications for abiotic oxygen, water loss, and habitability. *The Astronomical Journal, 155*(5), 195. https://doi.org/10.3847/1538-3881/aab608

[^cherubim2024]: Cherubim, C., Wordsworth, R., Hu, R., & Shkolnik, E. (2024). Strong Fractionation of Deuterium and Helium in Sub-Neptune Atmospheres along the Radius Valley. *The Astrophysical Journal, 967*(2), 139. https://doi.org/10.3847/1538-4357/ad3e77

[^yoshida]: Yoshida, T., Terada, N., Ikoma, M., & Kuramoto, K. (2022). Less Effective Hydrodynamic Escape of H$_2$–H$_2$O Atmospheres on Terrestrial Planets Orbiting Pre-main-sequence M Dwarfs. *The Astrophysical Journal, 934*(2), 137. https://doi.org/10.3847/1538-4357/ac7be7
