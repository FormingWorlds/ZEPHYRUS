# ZEPHYRUS model overview

ZEPHYRUS computes the bulk hydrodynamic atmospheric escape rate driven by stellar XUV irradiation for rocky exoplanets coupled to the [PROTEUS](https://proteus-framework.org) interior–atmosphere framework. It implements an energy-limited (EL) formalism following Watson et al. (1981) [^watson] and Lopez & Fortney (2013) [^lopez]. ZEPHYRUS is called at each PROTEUS time step with the current planetary radius and mass, the stellar XUV flux supplied by [MORS](https://proteus-framework.org/MORS) , and the escape radius computed from the atmospheric structure produced by AGNI or JANUS. The mass-loss rate it returns is distributed across atmospheric species according to their elemental mass mixing ratios, so the atmosphere is depleted in bulk without elemental fractionation.

A model parameter reference can be found [here](../Reference/parameters.md).

## Energy-limited escape

The mass-loss rate is computed by `escape.EL_escape` as

$$\dot{M}_\mathrm{EL} = \frac{\epsilon\,\pi\,R^3_\mathrm{XUV}\,F_\mathrm{XUV}}{G\,M_p\,K_\mathrm{tide}} \tag{1}$$

where $\epsilon$ is the escape efficiency factor (`epsilon`), $R_\mathrm{XUV}$ is the planetary radius at which the atmosphere becomes optically thick to stellar XUV photons (`Rxuv`), $F_\mathrm{XUV}$ is the XUV flux received at the planet (`Fxuv`) supplied by MORS, $M_p$ is the planetary mass (`Mp`), $G$ is the gravitational constant, and $K_\mathrm{tide}$ is the tidal correction factor described below. The efficiency $\epsilon$ quantifies the fraction of incident XUV energy that is converted into work against gravity to drive the outflow; canonical values for rocky planets lie in the range $0.1 \leq \epsilon \leq 0.3$, although ZEPHYRUS accepts any $\epsilon \in (0, 1]$.

### Radius scaling

The cubic radius factor in the numerator of Eq. (1) is selected at runtime by the `scaling` argument of `escape.EL_escape`:

| `scaling` | Expression | Description |
|---|---|---|
| `2` | $R_p\,R^2_\mathrm{XUV}$ | Default; XUV-absorbing cross-section weighted by surface radius |
| `3` | $R^3_\mathrm{XUV}$ | All three powers taken at the XUV radius |

Both forms reduce to $R_p^3$ when $R_\mathrm{XUV} = R_p$, which is the conservative lower bound on the mass-loss rate adopted by Luger & Barnes (2015) [^luger] and Moore et al. (2023) [^moore]. Allowing $R_\mathrm{XUV} > R_p$ increases the effective XUV-absorbing area and therefore the escape rate. In PROTEUS, $R_\mathrm{XUV}$ is recomputed at each time step from the atmospheric pressure–temperature profile at a user-specified reference pressure $P_\mathrm{XUV}$.

### Tidal correction $K_\mathrm{tide}$

When the `tidal_contribution` flag is `True`, the effective gravitational potential is reduced by the host star's tidal field through the factor

$$K_\mathrm{tide} = 1 - \frac{3}{2\xi} + \frac{1}{2\xi^3}, \qquad \xi = \frac{R_\mathrm{Hill}}{R_\mathrm{XUV}} \tag{2}$$

with the Hill radius

$$R_\mathrm{Hill} = a\,(1-e)\,\left(\frac{M_p}{3\,M_\star}\right)^{1/3} \tag{3}$$

where $a$ is the planetary semi-major axis, $e$ is the orbital eccentricity, and $M_\star$ is the stellar mass. $K_\mathrm{tide} \to 1$ for $\xi \gg 1$ (i.e., when the XUV radius lies well inside the Hill sphere) and decreases as the atmosphere expands toward the Roche lobe, enhancing escape. When `tidal_contribution` is `False`, $K_\mathrm{tide} = 1$ is enforced.

---

## Coupling to PROTEUS

ZEPHYRUS treats atmospheric escape as a bulk process: at each PROTEUS time step the total mass-loss rate from Eq. (1) is partitioned across atmospheric species in proportion to their elemental mass mixing ratios as computed by CALLIOPE. No elemental fractionation between light and heavy species is imposed in the outflow itself; however, because only outgassed volatiles are subject to escape while dissolved species remain in the magma ocean reservoir, escape fractionates the planet's *total* (interior + atmosphere) volatile budget over time, preferentially retaining species that are highly soluble in silicate melts (e.g. H$_2$O, S$_2$).


## Regime of validity

The EL formalism is appropriate in the high-irradiation, hydrodynamic regime that dominates atmospheric loss during the first $\sim 10^6$–$10^8$ yr of evolution for close-in rocky planets [^watson][^lammer2003]. Outside this regime—at lower XUV fluxes or for less extended atmospheres—non-thermal escape (Jeans escape, ion pickup, charge exchange) becomes comparable to or exceeds the hydrodynamic rate, and the bulk EL prescription no longer applies. ZEPHYRUS does not currently include these processes; users should verify that the integrated XUV-driven loss exceeds non-thermal estimates (e.g. $\sim 10^7$–$10^8$ g s$^{-1}$ for an Earth-mass planet; Kislyakova et al. 2014 [^kislyakova]) before interpreting model outputs.

Similarly, the bulk-removal assumption breaks down when the hydrodynamic particle flux drops below the critical flux required to drag heavy species against gravity, at which point compositional fractionation in the outflow becomes significant [^wordsworth2018][^cherubim2024]. Following Yoshida et al. (2022) [^yoshida], the critical flux for H$_2$O in an H$_2$ background is $\approx 1.9 \times 10^{8}$ g s$^{-1}$.

---

[^watson]: Watson, A. J., Donahue, T. M., & Walker, J. C. G. (1981). The dynamics of a rapidly escaping atmosphere: applications to the evolution of Earth and Venus. *Icarus, 48*(2), 150–166. https://doi.org/10.1016/0019-1035(81)90101-9


[^lopez]: Lopez, E. D., & Fortney, J. J. (2013). The role of core mass in controlling evaporation: the Kepler radius distribution and the Kepler-36 density dichotomy. *The Astrophysical Journal, 776*(1), 2. https://doi.org/10.1088/0004-637X/776/1/2

[^luger]: Luger, R., & Barnes, R. (2015). Extreme water loss and abiotic O$_2$ buildup on planets throughout the habitable zones of M dwarfs. *Astrobiology, 15*(2), 119–143. https://doi.org/10.1089/ast.2014.1231

[^moore]: Moore, K., Cowan, N. B., & Boukaré, C.-É. (2023). Atmospheric retention on Earth-sized planets around M dwarfs through magma-ocean outgassing. *Monthly Notices of the Royal Astronomical Society, 526*(4), 6235–6249. https://doi.org/10.1093/mnras/stad3098

[^lammer2003]: Lammer, H., Selsis, F., Ribas, I., et al. (2003). Atmospheric loss of exoplanets resulting from stellar X-ray and extreme-ultraviolet heating. *The Astrophysical Journal, 598*(2), L121–L124. https://doi.org/10.1086/380815

[^kislyakova]: Kislyakova, K. G., Johnstone, C. P., Odert, P., et al. (2014). Stellar wind interaction and pick-up ion escape of the Kepler-11 "super-Earths". *Astronomy & Astrophysics, 562*, A116. https://doi.org/10.1051/0004-6361/201322933

[^wordsworth2018]: Wordsworth, R. D., Schaefer, L. K., & Fischer, R. A. (2018). Redox evolution via gravitational differentiation on low-mass planets: implications for abiotic oxygen, water loss, and habitability. *The Astronomical Journal, 155*(5), 195. https://doi.org/10.3847/1538-3881/aab608

[^cherubim2024]: Cherubim, C., Wordsworth, R., Hu, R., & Shkolnik, E. (2024). Strong water loss and atmospheric escape from oxygen-rich rocky exoplanets. *The Astrophysical Journal, 967*(2), 139. https://doi.org/10.3847/1538-4357/ad3e6a

[^yoshida]: Yoshida, T., Terada, N., Ikoma, M., & Kuramoto, K. (2022). Effect of radiative cooling on the heating efficiency and mass-loss rate of hot Jupiters in hydrodynamic escape. *The Astrophysical Journal, 934*(2), 137. https://doi.org/10.3847/1538-4357/ac7be7