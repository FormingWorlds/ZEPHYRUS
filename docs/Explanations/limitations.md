# Limitations

ZEPHYRUS implements the **energy-limited (EL) approximation** to hydrodynamic atmospheric escape, given by Eq. (1) of the [model overview](model.md), and the **giant-impact erosion scaling law** of Eq. (4). Both are deliberate simplifications of much richer physical problems. The most important regimes and processes the model does not cover are summarised below.

---

## What ZEPHYRUS *does* model

Two channels. The first is bulk hydrodynamic escape driven by stellar XUV irradiation, in the energy-limited approximation, with an optional tidal correction (Eq. 2 of the [model overview](model.md)). The tidal correction is defined only outside the Roche lobe, where the Hill-to-XUV radius ratio $\xi > 1$; ZEPHYRUS raises an error for $\xi \le 1$, at which point the atmosphere reaches the Roche lobe and the energy-limited approximation no longer holds. The mass-loss rate is partitioned across atmospheric species in proportion to their elemental mass mixing ratios. The second channel is the fraction of the target's atmosphere eroded by a single giant impact, from a fitted power law in the collision speed, mass ratio, density ratio, and impact angle (Eq. 4 of the [model overview](model.md)).

Everything below is **not modelled.**

---

## Giant-impact erosion

The collision channel is a single fitted power law, not an impact simulation, and inherits the scope of the simulation suite behind it:

- **Thin atmospheres only.** The fit covers atmospheres of order 1 percent of the planet mass. A substantially thicker envelope cushions the impactor and alters its trajectory, and the eroded fraction is no longer described by the law.
- **Target-side loss only.** The law returns what the target's atmosphere loses. Any atmosphere the impactor itself carries, and any volatile delivery from the impactor into the merged body, is outside the function; the underlying paper shows that slow, grazing, atmosphere-hosting impactors can deliver most of their atmosphere, so treating the impactor as ballastless is a caller-side assumption, not a property of the collision.
- **No mantle or core erosion.** Violent impacts also strip silicate and metal mass; the law tracks only the atmospheric fraction.
- **Chaotic regime scatter.** Slow, head-on collisions produce chaotic fall-back and sloshing; the fit carries about 20 percent scatter there, against 9 percent overall.
- **Linearised interacting-mass geometry.** The common-height cap construction behind $f_M(b)$ misbehaves for a much denser, much smaller impactor near head-on, outside the fitted density ratios; ZEPHYRUS clamps $f_M$ to $[0, 1]$ in that corner rather than extrapolating the artifact.
- **Fit-domain extrapolation is unflagged.** The function evaluates the power law for any physically valid inputs; it does not warn when masses, densities, or speeds leave the fitted ranges listed in the [model overview](model.md). Staying inside them is the caller's responsibility.

---

## Other hydrodynamic regimes

The EL approximation assumes a fixed fraction $\epsilon$ of absorbed XUV energy goes into driving the outflow. This breaks down in several ways:

- **Radiative cooling is ignored.** Atomic line cooling, molecular emission, and ionisation losses can divert XUV energy away from heating the bulk gas, reducing the effective $\epsilon$. ZEPHYRUS treats $\epsilon$ as a constant input rather than computing it self-consistently. Setting $\epsilon = 1$ in particular is a non-physical upper limit on the mass-loss rate.
- **Fractionation in the outflow is not captured.** When the particle flux drops below the critical value required to drag heavy species along, the outflow becomes compositionally fractionated: hydrogen escapes preferentially and the residual atmosphere is enriched in heavy species. ZEPHYRUS removes everything in bulk. Fractionation will be implemented in the future.
- **$\epsilon$ is held constant in time.** In reality the efficiency evolves with planet mass, radius, and incident flux. Fixed-$\epsilon$ models can overestimate mass loss at late times.

---

## Non-hydrodynamic escape

These processes operate on a different physical basis (kinetic rather than fluid) and are neglected because they are subdominant in the high-XUV regime that ZEPHYRUS targets:

- Jeans escape
- Ion pickup
- Charge exchange
- Photochemical escape
- Sputtering
- Polar wind / unmagnetised ion outflow

For present-day Earth and Venus these mechanisms dominate over hydrodynamic escape, with total non-thermal rates around $\sim 10^3$ g s$^{-1}$; many orders of magnitude below the EL rates ZEPHYRUS produces during the early evolution phase.

---

## Other escape drivers

**Core-powered mass loss** is not implemented. This mechanism is driven by the planet's own internal heat and dominates for low-gravity planets at high equilibrium temperatures (~500–2000 K) over $\sim 10^9$ yr timescales. It is complementary to XUV-driven escape rather than competing with it.

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