# Limitations

ZEPHYRUS implements the **energy-limited (EL) approximation** to hydrodynamic atmospheric escape, given by Eq. (1) of the [model overview](model.md), the **giant-impact erosion scaling law** of Eq. (4), and the **[regime dispatcher](dispatcher.md)**, which classifies an atmosphere into the boil-off, hydrodynamic, or hydrostatic regime before choosing a rate and partitions it over species. All are deliberate simplifications of much richer physical problems. The most important regimes and processes the package does not cover, and which entry point covers what, are summarised below.

---

## What ZEPHYRUS *does* model

Three entry points. The first is bulk hydrodynamic escape driven by stellar XUV irradiation, in the energy-limited approximation, with an optional tidal correction (Eq. 2 of the [model overview](model.md)). The tidal correction is defined only outside the Roche lobe, where $\xi > 1$ with $\xi$ measured from the radius the `scaling` argument selects; ZEPHYRUS raises an error for $\xi \le 1$, at which point the atmosphere reaches the Roche lobe and the energy-limited approximation no longer holds. The mass-loss rate is partitioned across atmospheric species in proportion to their elemental mass mixing ratios. The second is the fraction of the target's atmosphere eroded by a single giant impact, from a fitted power law in the collision speed, mass ratio, density ratio, and impact angle (Eq. 4 of the [model overview](model.md)). The third is the [regime dispatcher](dispatcher.md), which adds bolometrically driven boil-off, a radiation-recombination-limited cap with a radiatively cooled wind temperature, per-species Jeans escape with a diffusion-limited supply, a collisionality switch between the fluid and kinetic regimes, Roche-overflow handling, and N-species fractionation of the hydrodynamic outflow.

Everything below is **not modelled**, or modelled only by the dispatcher entry point where stated.

---

## Giant-impact erosion

The collision channel is a single fitted power law, not an impact simulation, and inherits the scope of the simulation suite behind it:

- **Thin atmospheres only.** The fit covers atmospheres of order 1 percent of the planet mass. A substantially thicker envelope cushions the impactor and alters its trajectory, and the eroded fraction is no longer described by the law.
- **Target-side loss only.** The law returns what the target's atmosphere loses. Any atmosphere the impactor itself carries, and any volatile delivery from the impactor into the merged body, is outside the function; the underlying paper shows that in slow, grazing collisions with an atmosphere-hosting impactor the target can retain about 85 percent of the two bodies' combined initial atmospheres, so treating the impactor as ballastless is a caller-side assumption, not a property of the collision.
- **No mantle or core erosion.** Violent impacts also strip silicate and metal mass; the law tracks only the atmospheric fraction.
- **Chaotic regime scatter.** Slow, head-on collisions produce chaotic fall-back and sloshing; the fit carries about 20 percent scatter there, against 9 percent overall.
- **Linearised interacting-mass geometry.** The common-height cap construction behind $f_M(b)$ misbehaves for a much denser, much smaller impactor near head-on, outside the fitted density ratios; ZEPHYRUS clamps $f_M$ to $[0, 1]$ in that corner rather than extrapolating the artifact.
- **Fit-domain extrapolation is unflagged.** The function evaluates the power law for any physically valid inputs; it does not warn when masses, densities, or speeds leave the fitted ranges listed in the [model overview](model.md). Staying inside them is the caller's responsibility.

---

## Other hydrodynamic regimes

The EL approximation assumes a fixed fraction $\epsilon$ of absorbed XUV energy goes into driving the outflow. This breaks down in several ways:

- **Radiative cooling is ignored by `EL_escape`.** Atomic line cooling, molecular emission, and ionisation losses can divert XUV energy away from heating the bulk gas, reducing the effective $\epsilon$. `EL_escape` treats $\epsilon$ as a constant input rather than computing it self-consistently, and setting $\epsilon = 1$ in particular is a non-physical upper limit on the mass-loss rate. The dispatcher partially addresses this: its wind-temperature thermostat balances photoionization heating against four radiative cooling channels, and its radiation-recombination cap bounds the rate where recombination radiates the energy away, but the energy-limited efficiency itself remains an input there too (fixed, or the fitted efficiency of Caldiroli et al. 2022).
- **Fractionation in the outflow is not captured by `EL_escape`.** When the particle flux drops below the critical value required to drag heavy species along, the outflow becomes compositionally fractionated: hydrogen escapes preferentially and the residual atmosphere is enriched in heavy species. `EL_escape` removes everything in bulk. The dispatcher implements this through the N-species fractionation closure, with per-species dropout thresholds and provenance-classed diffusion coefficients.
- **$\epsilon$ is held constant in time.** In reality the efficiency evolves with planet mass, radius, and incident flux. Fixed-$\epsilon$ models can overestimate mass loss at late times; the dispatcher's fitted-efficiency mode captures the potential dependence but not a time dependence beyond it.

---

## Non-hydrodynamic escape

These processes operate on a different physical basis (kinetic rather than fluid). Jeans escape is now implemented: the dispatcher's hydrostatic branch evaluates per-species Jeans effusion with a kinetic enhancement factor and a diffusion-limited supply cap, on an extended upper structure. The remaining kinetic channels are not modelled by any entry point:

- Ion pickup
- Charge exchange
- Photochemical escape
- Sputtering
- Polar wind / unmagnetised ion outflow

For present-day Earth and Venus these mechanisms dominate over hydrodynamic escape, with total non-thermal rates around $\sim 10^3$ g s$^{-1}$; many orders of magnitude below the EL rates ZEPHYRUS produces during the early evolution phase. Because they are absent, the dispatcher's hydrostatic heavy-element rates are lower limits, and it flags every hydrostatic result accordingly; grid points where the neutral and plasma escape-temperature conventions disagree are flagged as contested, with both branch rates reported, because the unmodelled ion physics decides them.

---

## Other escape drivers

**Core-powered mass loss** is not implemented in `EL_escape`. This mechanism is driven by the planet's own internal heat and dominates for low-gravity planets at high equilibrium temperatures (~500–2000 K) over $\sim 10^9$ yr timescales. The dispatcher's bolometric branch covers it: boil-off below the activation threshold, and a luminosity-capped residual past it, which represents the long tail without adjudicating the open dispute over how long it survives.

---

## Stellar XUV uncertainties

The XUV flux $F_\mathrm{XUV}$ that enters Eq. (1) of the [model overview](model.md) carries large intrinsic uncertainties from the underlying stellar evolution model:

- Saturation timescales for the stellar XUV phase can vary from ~10 to ~300 Myr for G stars and up to ~1 Gyr for fully convective M dwarfs, depending on initial rotation.
- The integrated XUV flux, and therefore the integrated mass loss, can vary by factors of $\sim 2–10$ between standard stellar evolution prescriptions.
- The ISM absorbs stellar XUV emission, so observational anchors on young-star XUV luminosities are themselves uncertain.

Because of these uncertainties, the mass-loss rates computed by ZEPHYRUS should generally be treated as an upper bound.

---

## Atmospheric chemistry

- **No photochemistry.** Hazes, aerosols, and photochemically-produced species are not tracked in the coupled framework.
- **$R_\mathrm{XUV}$ is set by a single reference pressure** $P_\mathrm{XUV}$ specified in the config.

---

## Practical implications

For users:

- Avoid $\epsilon > 0.3$ for rocky planets unless you have a specific reason. $\epsilon \approx 0.15$ is the conservative baseline.
- For close-in M-dwarf planets where elemental fractionation is expected to matter, ZEPHYRUS bulk rates are a lower bound on the change in atmospheric mean molecular weight. The actual atmosphere should become heavier faster than the model predicts.