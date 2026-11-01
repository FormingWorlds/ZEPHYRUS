# Roche overflow

A planet close to its star sits inside a Roche lobe, the region where its own gravity holds gas against the star's tidal field. The lobe touches the lobe of the star at the inner Lagrange point L1, the saddle of the combined potential on the line between the two bodies. Gas that reaches L1 can flow to the star through the narrow cross section of the saddle, the L1 nozzle. The label `roche_overflow` has two readings:

- **A transfer rate.** The flow through L1 is the fifth candidate of the dispatcher. When it is admitted and is the largest rate, it is dispatched with the label `roche_overflow`.
- **A rename.** After the rate is chosen, the Roche screen tests whether the flow of the winning branch reaches beyond the Hill radius. If it does, the state is renamed `roche_overflow` and keeps its rate.

`diagnostics['roche']['rate_branch']` says which reading applies. Symbols shared across the section are in the [notation tables](../regimes.md#notation).

## The transfer through L1

The rate is that of Jackson et al. (2017) [^jackson17], in the mass-transfer tradition of Ritter (1988) [^ritter], rebuilt to hold at planetary mass ratios. In their reading, Roche-lobe overflow and evaporative escape are one unbound hydrodynamic outflow seen at two separations, distinguished by whether the photosphere nearly coincides with the lobe.

**The launch level.** The flow is isothermal at the launch level's temperature $T_\mathrm{launch}$ and mean particle mass $\mu_\mathrm{launch}$. The setting `nozzle_temperature` chooses the level:

- `'photospheric'` (default): the photospheric working level, so $(R_\mathrm{launch}, \rho_\mathrm{launch}, T_\mathrm{launch}, \mu_\mathrm{launch}) = (R_\mathrm{ph}, \rho_\mathrm{ph}, T_\mathrm{ph}, \mu_\mathrm{ph})$, as in the source.
- `'wind'`: the wind base at the wind's own temperature and mean mass, $(R_\mathrm{base},\ P_\mathrm{base}\, \mu_\mathrm{w} m_\mathrm{p} / (k_\mathrm{B} T_\mathrm{w}),\ T_\mathrm{w},\ \mu_\mathrm{w} m_\mathrm{p})$.

**The geometry.** At an orbital separation $d$, with mass ratio $q = M_\mathrm{p}/M_\star$,

$$v_\mathrm{th} = \sqrt{\frac{k_\mathrm{B} T_\mathrm{launch}}{\mu_\mathrm{launch}}}, \qquad \Omega^2 = \frac{G (M_\mathrm{p} + M_\star)}{d^3}, \qquad A(q) = 4 + \frac{b_1}{b_2 + q^{1/3} + q^{-1/3}} \tag{1}$$

where $v_\mathrm{th}$ is the isothermal sound speed, $\Omega$ the orbital frequency, and $A$ the dimensionless curvature of the potential at L1 in the fit of their Eq. (10), with $b_1 = 2 \cdot 3^{2/3}$ and $b_2 = b_1/4 - 2$, accurate to 0.3% for all mass ratios. The lobe is the Eggleton (1983) volume-equivalent radius as printed in their Section 2.1 [^eggleton],

$$R_\mathrm{lobe} \;=\; d\, \frac{0.49\, q^{2/3}}{0.6\, q^{2/3} + \ln(1 + q^{1/3})} \tag{2}$$

**The barrier.** The barrier is the difference between two values of their volume-averaged Roche potential, their Eq. (14), written $\Psi$ here (they write $\Phi$):

$$\begin{aligned} \Psi(r) \;=\;& -\left(\frac{G M_\star}{d} + \frac{G M_\star^2}{2 d (M_\mathrm{p} + M_\star)}\right) \\ &- \frac{G M_\mathrm{p}}{r}\Bigg[1 + \frac{M_\mathrm{p} + M_\star}{3 M_\mathrm{p}} \left(\frac{r}{d}\right)^3 \\ &\qquad + \frac{4}{45}\, \frac{(M_\mathrm{p} + M_\star)^2 + 9 M_\star^2 + 3 M_\star (M_\mathrm{p} + M_\star)}{M_\mathrm{p}^2} \left(\frac{r}{d}\right)^6\Bigg] \end{aligned} \tag{3}$$

This is the potential of the equipotential enclosing the volume of a sphere of radius $r$, with $\Psi_\mathrm{L1} = \Psi(R_\mathrm{lobe})$ and $\Psi_\mathrm{launch} = \Psi(R_\mathrm{launch})$.

**The rate.** A Bernoulli integral from the launch level, where the flow is slow, to L1, where it is sonic, gives the density at L1 (their Eqs. (11) to (13), the source of the factor $e^{-1/2}$). The gas crosses the elliptical nozzle that the curvature at L1 admits (their Eqs. (8) and (9)) at the sound speed. The transfer rate is their Eq. (3):

$$\begin{gathered} \dot{M}_\mathrm{L1}(d) \;=\; e^{-1/2}\, \rho_\mathrm{launch}\, \exp\!\left(-\frac{\Psi_\mathrm{L1} - \Psi_\mathrm{launch}}{v_\mathrm{th}^2}\right) v_\mathrm{th}\; A_\mathrm{noz}, \\ A_\mathrm{noz} = \frac{2\pi\, v_\mathrm{th}^2}{\Omega^2 \sqrt{A (A - 1)}} \end{gathered} \tag{4}$$

with $A_\mathrm{noz}$ the nozzle area. The expansion behind Eq. (3) converges only inside the lobe, so a launch level at or beyond the lobe, $R_\mathrm{launch} \ge R_\mathrm{lobe}$, sets the exponential to 1. That is their lobe-filling case (their Figure 5), flagged `nozzle_saturated`, and a lower bound on the transfer.

### Where the transfer applies

The rate is admitted on an orbital phase only where two conditions hold.

**The sonic radius reaches L1.** The nozzle is the flow's narrowest point only where the isothermal sonic radius of the launch gas reaches L1, so that no spherical transonic wind fits inside the lobe:

$$R_\mathrm{s,L1} = \frac{G M_\mathrm{p}}{2\, v_\mathrm{th}^2} \;\ge\; R_\mathrm{L1}, \qquad R_\mathrm{L1} = d\left(\varepsilon - \frac{\varepsilon^2}{3} - \frac{\varepsilon^3}{9}\right), \qquad \varepsilon = \left(\frac{q}{3}\right)^{1/3} \tag{5}$$

with $R_\mathrm{L1}$ the small-mass-ratio expansion of the distance to L1, which is the Hill radius at leading order. Where the sonic radius sits inside it, the gas becomes supersonic in a spherical wind first, and the wind branches describe it. Jackson et al. draw the same comparison qualitatively in their Figure 9. Turning it into a criterion on the candidate is this module's choice.

**The gas at L1 is a fluid.** Eq. (5) alone holds for almost every cool, bound planet, whose sonic radius is enormous. The rate assumes a collisional isothermal column all the way to L1, so the density it implies there, $n_\mathrm{L1} = \dot{M}_\mathrm{L1} / (v_\mathrm{th} A_\mathrm{noz} \mu_\mathrm{launch})$, must make the gas a fluid on the nozzle's own scale:

$$\mathrm{Kn}_\mathrm{L1} \;=\; \frac{1}{\sqrt{2}\, \sigma_\mathrm{C}\, n_\mathrm{L1}\, H_\mathrm{L1}} \;\le\; \mathrm{Kn}_\mathrm{c}, \qquad H_\mathrm{L1} = \sqrt{\frac{A_\mathrm{noz}}{2\pi}} \tag{6}$$

where $\sigma_\mathrm{C}$ is the mixture cross section of the [collisionality switch](collisionality.md) for the launch composition at $T_\mathrm{launch}$, $\mathrm{Kn}_\mathrm{c}$ is its threshold (`kn_crit`, without the hysteresis window), and $H_\mathrm{L1}$ is the geometric mean of the nozzle's two widths, since their Eqs. (8) and (9) give the area as $2\pi$ times their product. This is the sonic-point test of the collisionality switch moved to L1, which Jackson et al. do not state. A planet whose gas at L1 is collisionless keeps the label of its own branch.

### The orbit average

The source treats a circular, synchronously rotating planet and has no eccentric formulation. Each orbital phase is therefore evaluated with the circular formula at its own separation $d(E) = a (1 - e \cos E)$, $E$ being the eccentric anomaly, and the result is averaged in time, with Kepler's $\mathrm{d}t \propto (1 - e \cos E)\, \mathrm{d}E$, over the arc where the transfer applies:

$$\langle \dot{M}_\mathrm{L1} \rangle \;=\; \frac{\sum_k (1 - e \cos E_k)\; \chi_k\; \dot{M}_\mathrm{L1}\!\left(d(E_k)\right)}{\sum_k (1 - e \cos E_k)}, \qquad E_k = \frac{2\pi (k + 1/2)}{64} \tag{7}$$

The sums run over 64 orbital phases, numbered $k = 0$ to 63 and spaced evenly in eccentric anomaly at the midpoints $E_k$. The weight $1 - e \cos E_k$ turns equal steps in eccentric anomaly into equal steps in time, so a phase near apoapsis, where the planet moves slowly, counts for more. $\chi_k$ is 1 on the phases that satisfy Eqs. (5) and (6) and 0 elsewhere. On a circular orbit the average is the instantaneous rate. The applicable arc surrounds periapsis, because the distance to L1 grows with separation while the sonic radius does not. The flags `nozzle_orbit_averaged` and `nozzle_partial_orbit` mark an eccentric result, and the average without $\chi_k$ is kept in the diagnostics so the rate stays comparable with the source's published values.

### Energy

The transfer rate has no energy limit. Like the source, it assumes the gas keeps one temperature all the way to L1, and it does not model the heating that would keep it there. To show where that assumption is strained, the diagnostics report the power the flow would need, to lift the gas over the barrier and to accelerate it to the sound speed, beside the two sources that could supply it: the planet's interior luminosity and the starlight the planet intercepts. Where the required power exceeds both, the gas cannot stay at that temperature, and the transfer rate is likely an overestimate. The required power is computed for the transfer rate whether or not it won, so a state where it lost still shows what it would have cost. The [limitations page](../limitations.md) lists the other conventions and omissions of this rate.

## The Roche screen

Every candidate except the transfer through L1 assumes the flow is bound to the planet. After the rate is chosen, the flow radius of the branch that produced it is tested against the periapsis Hill radius of [Eq. (1) of the overview](../regimes.md#tidal-geometry):

$$\xi_\mathrm{flow} \;=\; \frac{R_\mathrm{Hill}}{R_\mathrm{flow}}, \qquad R_\mathrm{flow} = \begin{cases} R_\mathrm{B}, & \text{boil-off rate or residual} \\ \max\left(R_\mathrm{XUV},\ R_\mathrm{s}\right), & \text{hydrodynamic rate} \\ R_\mathrm{exo}, & \text{hydrostatic rate} \end{cases} \tag{8}$$

with $R_\mathrm{B}$ from [boil-off](boiloff.md), Eq. (2), $R_\mathrm{s}$ from [the hydrodynamic wind](wind.md), Eq. (4), and $R_\mathrm{exo}$ from [hydrostatic escape](hydrostatic.md), Eq. (2). When the transfer through L1 wins, $R_\mathrm{flow}$ keeps the value of the branch it beat.

A state with $\xi_\mathrm{flow} \le 1$, or with $\xi \le 1$ for the planet itself, is spilling over the gravitational boundary rather than escaping through a bound outflow, and it is renamed `roche_overflow`. `zephyrus.dispatcher.roche_label` implements the screen. The second test decides alone only when the flow radius lies inside the planet. No profile that `isothermal_profile` builds gets there, but the test keeps a planet larger than its own Hill sphere from being reported as bound.

**The rate stands.** The screen renames a state and never changes its rate. The radius it tests belongs to the branch that won the rate comparison, so both sides of the screen's boundary hold the same branch. Reporting that branch's own rate keeps the dispatched rate continuous across the boundary, while substituting another branch's formula would make it jump by orders of magnitude at a line the physics puts nowhere in particular. When the screen renames a bound branch, the label means that the flow reaches the lobe, and the rate beside it is the bound-flow estimate, a lower limit on what tides would do.

**Which geometry fired.** Owen & Jackson (2012) separate two overflow geometries [^oj12], and the subflag `roche_subflag` says which one applies:

- `dynamical`: the atmosphere itself reaches its Roche lobe, or the planet fills its Hill sphere. The outer extent of the atmosphere, the larger of the profile-top radius and $R_\mathrm{exo}$, is reported as `r_atmosphere` in `diagnostics['roche']` and compared with the periapsis lobe radius of Eq. (2), `r_lobe` in `diagnostics['nozzle']`. This takes precedence whichever candidate carries the rate.
- `no_transonic`: the atmosphere sits inside its lobe, but its sonic surface would have to sit outside it, so no transonic solution exists. Owen & Jackson describe this as a narrow band and suggest a subsonic wind out to the lobe.
- `neither`: the label was won by the transfer rate while neither the atmosphere nor the flow radius reaches out.

The comparison for `dynamical` uses the Roche lobe, the critical surface, which sits about 0.70 of the way out to the Hill radius. The subflag is worth reading before a label is trusted. On the boil-off branch the tested radius is a sonic radius that grows with the Jeans parameter, so a tightly bound heavy atmosphere can push it several Hill radii out while the atmosphere itself sits deep inside, at a rate with no numerical content. `r_atmosphere` and `diagnostics['rate_floor']` separate that case from a real overflow ([troubleshooting](../../How-to/troubleshooting.md)).

**Near misses.** A state with $\xi_\mathrm{flow} < 1.5$ that the screen did not rename raises the flag `near_roche`, because the tidal factor inflates the energy-limited rate steeply there. The flag reports, and never modifies, the rate.

## Measured behaviour

Numbers measured with this code, on the states named:

- An Earth-mass steam atmosphere at 1000 K with no XUV puts $4 \times 10^{-28}$ particles per cubic metre at L1, a Knudsen number of $6 \times 10^{37}$, and returns `hydrostatic`.
- The fluid condition of Eq. (6) is a criterion boundary. On an 8 Earth-mass, 3 Earth-radius hydrogen envelope at $F_\mathrm{XUV} = 10^{-4}$ W m⁻², the dispatched rate drops by a factor of 113 where it fails, at 0.0526 au.
- Under `nozzle_temperature = 'wind'`, the rate holds to 0.5% across a factor of four in launch radius along the wind's isothermal column, and moves by a factor of eleven along the profile's colder one. On a mildly inverted profile, the rate varies by a factor of 2.4 across three decades of launch level, against 1.14 on an isothermal one.

[^jackson17]: Jackson, B., Arras, P., Penev, K., Peacock, S., & Marchant, P. (2017). A new model of Roche lobe overflow for short-period gaseous planets and binary stars. *The Astrophysical Journal, 835*(2), 145. https://doi.org/10.3847/1538-4357/835/2/145

[^ritter]: Ritter, H. (1988). Turning on and off mass transfer in cataclysmic binaries. *Astronomy & Astrophysics, 202*, 93.

[^eggleton]: Eggleton, P. P. (1983). Approximations to the radii of Roche lobes. *The Astrophysical Journal, 268*, 368. https://doi.org/10.1086/160960

[^oj12]: Owen, J. E., & Jackson, A. P. (2012). Planetary evaporation by UV and X-ray radiation: basic hydrodynamics. *Monthly Notices of the Royal Astronomical Society, 425*(4), 2931. https://doi.org/10.1111/j.1365-2966.2012.21481.x
