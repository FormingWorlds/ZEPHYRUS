# Fractionation

A hydrodynamic wind does not carry every species equally. Light species stream out; heavier ones are dragged along through collisions, lag the flow, and below a species-specific threshold flux they stop escaping altogether while continuing to exert drag on everything that still escapes. Over time this fractionates the atmosphere, enriching it in heavy species, which is one of the main observable signatures escape leaves behind. When the [regime framework](regimes.md) confirms a hydrodynamic wind, ZEPHYRUS partitions the bulk rate over species with a simultaneous N-species closure (Attia & Lichtenberg 2026, in prep. [^attia]); this page describes what the closure solves, where its coefficients come from, and which regimes it applies to.

## The problem and the closure

The classic treatment is the two-species problem of Hunten, Pepin & Walker (1987) [^hunten]: a light major species escaping through one heavy species, with a crossover mass separating dragged-along from left-behind. Real atmospheres carry many species at once, and each pair interacts through its own binary diffusion coefficient, so the two-species answer cannot just be applied pairwise. The closure generalizes the constant-composition treatment of the multispecies wind equations of Zahnle, Kasting & Pollack (1990) [^z90] to N species escaping simultaneously.

The solved system couples the species drift velocities. For each escaping species $j$, the drag exerted by every other species balances its weight surplus,

$$\sum_{i\,\mathrm{escaping}} \frac{X_i\,(w_i - w_j)}{b_{ij}} \;-\; w_j \sum_{k\,\mathrm{retained}} \frac{X_k}{b_{jk}} \;=\; \frac{m_j\, g}{k_\mathrm{B} T} - \frac{1}{\bar{H}} \tag{1}$$

where $X_i$ is the mole fraction of species $i$, $w_i$ its escape velocity scale (the number flux is $\Phi_i = X_i w_i$), $b_{ij}$ the binary diffusion parameter of the pair, $m_j$ the particle mass, $g$ the gravity at the wind base, $T$ the wind temperature, and $\bar{H}$ the one density scale height that every escaping gas shares, itself an unknown of the solve rather than an input. The system closes with the mass constraint that the per-species fluxes carry the bulk rate the regime framework dispatched, $\sum_j m_j X_j \Phi_j = \phi$.

Which species escape is part of the solution, not an input. A heavy species whose settling under gravity beats the drag the outflow can exert on it drops out of the escaping set and moves to the retained set, where it still appears in the drag sums of Eq. (1). The solver finds the unique partition into escaping and retained species for which every escaping species has a positive flux and every retained species genuinely cannot be lifted; each heavy species therefore has a threshold bulk flux at which it starts to escape, and below the lowest threshold only the lightest species leaves. The returned per-species rates are non-negative and sum to the bulk rate at machine precision.

The closure reproduces, as exact special cases, the published treatments it generalizes: the two-species crossover of Hunten et al. (1987) in the form of Cherubim et al. (2024) [^cherubim], the three-species deuterium system of Gu & Chen (2023) [^guchen], the trace-minor relations of Odert et al. (2018) [^odert] and Zahnle et al. (1990), the non-trace three-species relations of Zahnle & Kasting (2023) [^zk23], and the prescribed-flux partition of Chassefière (1996) [^chassefiere], along with the worked Earth, Mars, and Venus numbers of Hunten et al. (1987). The test suite asserts every one of these reductions, and the general formulation is the subject of Attia & Lichtenberg (2026, in prep.) [^attia].

## Coefficients and their provenance

Everything species-dependent enters through the binary diffusion parameters $b_{ij}$, and no compilation measures every pair, so each pair carries a provenance class that travels with the result. Measured rows come from the compilations of Zahnle & Kasting (1986) [^zk86] and Zahnle & Kasting (2023) [^zk23], which trace to the reference measurements of Marrero & Mason (1972) [^marrero], with the noble-gas rows of Sasaki & Nakazawa (1988) [^sasaki] verified against them. Pairs no compilation prints are built by the reduced-mass and kinetic-diameter scaling rule of Zahnle & Kasting (2023), validated in and out of sample against the printed entries, and land in a wider uncertainty class. Rock-forming species (Na, Mg, Si, Fe) sit in the widest class of all, and results involving them carry a dedicated flag: no measured coefficient exists for any of their pairs, and sodium and magnesium ionize at the temperatures where rock vapor exists while these are neutral-gas coefficients.

## Where it applies

The closure evaluates at the XUV wind base on the atomized composition (molecules are photodissociated well below the launching level, so the escaping gas is atomic). It applies only where a hydrodynamic branch produced the rate; the other branches split their rates differently:

| Branch that produced the rate | Per-species split |
|---|---|
| `hydrodynamic:EL`, `hydrodynamic:RR` | The N-species closure at the wind base (this page); with fractionation disabled, reservoir mass fractions |
| `boiloff` | Reservoir mass fractions (no fractionation: the flow is fast and bulk) |
| `hydrostatic` | Natively per-species: each species carries its own Jeans flux and supply cap (see [escape regimes](regimes.md)) |

The split follows the branch and not the label, which matters under `roche_overflow`: that label renames a state without changing its rate, so the split is whatever the branch named in `diagnostics['roche']['rate_branch']` would have produced.

---

[^attia]: Attia, M., & Lichtenberg, T. (2026). In preparation.

[^hunten]: Hunten, D. M., Pepin, R. O., & Walker, J. C. G. (1987). Mass Fractionation in Hydrodynamic Escape. *Icarus, 69*, 532–549.

[^z90]: Zahnle, K., Kasting, J. F., & Pollack, J. B. (1990). Mass Fractionation of Noble Gases in Diffusion-Limited Hydrodynamic Hydrogen Escape. *Icarus, 84*(2), 502–527.

[^zk86]: Zahnle, K. J., & Kasting, J. F. (1986). Mass Fractionation during Transonic Escape and Implications for Loss of Water from Mars and Venus. *Icarus, 68*(3), 462–480.

[^zk23]: Zahnle, K. J., & Kasting, J. F. (2023). Elemental and isotopic fractionation as fossils of water escape from Venus. *Geochimica et Cosmochimica Acta, 361*, 228–244.

[^marrero]: Marrero, T. R., & Mason, E. A. (1972). Gaseous diffusion coefficients. *Journal of Physical and Chemical Reference Data, 1*(1), 3–118.

[^sasaki]: Sasaki, S., & Nakazawa, K. (1988). Origin of isotopic fractionation of terrestrial Xe: hydrodynamic fractionation during escape of the primordial H$_2$-He atmosphere. *Earth and Planetary Science Letters, 89*(3-4), 323–334.

[^cherubim]: Cherubim, C., Wordsworth, R., Hu, R., & Shkolnik, E. (2024). Strong Fractionation of Deuterium and Helium in Sub-Neptune Atmospheres along the Radius Valley. *The Astrophysical Journal, 967*(2), 139. https://doi.org/10.3847/1538-4357/ad3e77

[^guchen]: Gu, P.-G., & Chen, H. (2023). Deuterium Escape on Photoevaporating Sub-Neptunes. *The Astrophysical Journal Letters, 953*(2), L27. https://doi.org/10.3847/2041-8213/acee01

[^odert]: Odert, P., et al. (2018). Escape and fractionation of volatiles and noble gases from Mars-sized planetary embryos and growing protoplanets. *Icarus, 307*, 327–346.

[^chassefiere]: Chassefière, E. (1996). Hydrodynamic Escape of Oxygen from Primitive Atmospheres: Applications to the Cases of Venus and Mars. *Icarus, 124*, 537–552.
