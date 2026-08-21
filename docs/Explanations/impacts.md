# Giant impacts

A giant collision removes part of the target planet's atmosphere in a single event: a shock launched by the impact accelerates atmosphere past the escape velocity, locally near the impact site and globally through ground motion. This is the second mass-loss channel of ZEPHYRUS, physically and numerically separate from the continuous escape of the [regime framework](regimes.md): it is applied per collision rather than per time step, and its rate question ("what fraction is lost in this event") replaces the continuous channel's ("how fast is mass leaving"). A regime label is reserved for routing impacts through the same interface in the future; today the caller invokes the channel directly.

ZEPHYRUS computes the eroded fraction with `zephyrus.collision.mass_loss`, which implements the scaling law of Kegerreis et al. (2020), their Eq. 1 [^kegerreis]:

$$X \;=\; \min\!\left\{0.64 \left[ \left(\frac{v_\mathrm{c}}{v_\mathrm{esc}}\right)^2 \left(\frac{M_\mathrm{i}}{M_\mathrm{tot}}\right)^{1/2} \left(\frac{\rho_\mathrm{i}}{\rho_\mathrm{t}}\right)^{1/2} f_M(b) \right]^{0.65},\; 1\right\} \tag{1}$$

where $X$ is the fraction of the target's atmosphere removed (capped at 1 for total erosion), subscript $\mathrm{i}$ denotes the impactor and $\mathrm{t}$ the target, $v_\mathrm{c}$ is the speed at first contact, $M_\mathrm{tot} = M_\mathrm{i} + M_\mathrm{t}$ is the total mass, $\rho_\mathrm{i}$ and $\rho_\mathrm{t}$ are the bulk densities of the atmosphere-free bodies, and $b \equiv \sin\beta$ is the dimensionless impact parameter for impact angle $\beta$ (0 head-on, 1 fully grazing). The prefactor and exponent are least-squares fits to the paper's suite of 259 smoothed-particle-hydrodynamics simulations, each fitted with an uncertainty of 0.01. The mutual escape speed of the pair at contact is

$$v_\mathrm{esc} \;=\; \sqrt{\frac{2\,G\,(M_\mathrm{t} + M_\mathrm{i})}{R_\mathrm{t} + R_\mathrm{i}}} \tag{2}$$

with $R_\mathrm{t}$ and $R_\mathrm{i}$ the body radii at the base of any atmosphere, and $f_M(b)$ is the fractional interacting mass of the pair (their Eq. B1), built from density-weighted spherical caps of common height $d = (R_\mathrm{t} + R_\mathrm{i})(1 - b)$:

$$f_M \;=\; \frac{\rho_\mathrm{t}\, V^\mathrm{cap}_\mathrm{t} + \rho_\mathrm{i}\, V^\mathrm{cap}_\mathrm{i}}{\rho_\mathrm{t}\, V_\mathrm{t} + \rho_\mathrm{i}\, V_\mathrm{i}}, \qquad V^\mathrm{cap}_\mathrm{t,i} = \frac{\pi}{3}\, d^2 \left(3 R_\mathrm{t,i} - d\right) \tag{3}$$

where $V_\mathrm{t,i}$ are the full body volumes. At equal bulk densities $f_M$ reduces to the fractional interacting volume of their Eq. B2. The common-height caps are a linearized bookkeeping: outside the fitted geometry, for a much denser and much smaller impactor near head-on, the raw $f_M$ can leave $[0, 1]$ and vary non-monotonically with $b$, so ZEPHYRUS clamps $f_M$ to $[0, 1]$. Within the fitted domain the clamp never engages.

Three input conventions follow the paper and must be honored by the caller: $v_\mathrm{c}$ is the speed at first contact, not the relative speed at infinity; the masses and radii exclude any atmosphere, with radii taken at the base of the atmosphere; and the densities are bulk values of the atmosphere-free bodies.

## Fitted domain and accuracy

The law is constrained by simulations spanning target masses of roughly 0.3 to 3 Earth masses, impactor masses down to about 0.05 Earth masses, bulk densities from about half to double Earth's, contact speeds of 1 to 3 $v_\mathrm{esc}$, all impact angles, and thin atmospheres of order 1 percent of the planet mass. The median deviation of the simulations from the law is 9 percent, rising to about 20 percent for slow, head-on impacts, whose outcomes are chaotic. The loss depends only mildly on the atmosphere mass in this thin-atmosphere regime, with a factor of 10 less atmosphere increasing the eroded fraction by roughly 10 percent; substantially thicker atmospheres, which can cushion the impactor, fall outside the law's domain. The function evaluates the law for any physically valid inputs and does not warn when masses, densities, or speeds leave the fitted ranges; staying inside them is the caller's responsibility, and the [limitations page](limitations.md) lists what the channel leaves out (impactor-side atmosphere, volatile delivery, and mantle stripping among them).

The returned fraction applies to the target's atmosphere as a whole. Consistent with the bulk-removal treatment of the continuous channel's unfractionated splits, the caller partitions the lost mass across atmospheric species without elemental fractionation.

---

[^kegerreis]: Kegerreis, J. A., Eke, V. R., Catling, D. C., Massey, R. J., Teodoro, L. F. A., & Zahnle, K. J. (2020). Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density. *The Astrophysical Journal Letters, 901*(2), L31. https://doi.org/10.3847/2041-8213/abb5fb
