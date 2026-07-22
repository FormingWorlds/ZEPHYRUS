# Limitations

ZEPHYRUS implements the **energy-limited (EL) approximation** to hydrodynamic atmospheric escape, given by Eq. (1) of the [model overview](model.md). This is a deliberate simplification of a much richer physical problem. The most important regimes and processes the model does not cover are summarised below.

---

## What ZEPHYRUS *does* model

A single regime: bulk hydrodynamic escape driven by stellar XUV irradiation, in the energy-limited approximation, with an optional tidal correction (Eq. 2 of the [model overview](model.md)). The tidal correction is defined only outside the Roche lobe, where the Hill-to-XUV radius ratio $\xi > 1$; ZEPHYRUS raises an error for $\xi \le 1$, at which point the atmosphere reaches the Roche lobe and the energy-limited approximation no longer holds. The mass-loss rate is partitioned across atmospheric species in proportion to their elemental mass mixing ratios.

Everything below is **not modelled.**

---

## Other hydrodynamic regimes

The EL approximation assumes a fixed fraction $\epsilon$ of absorbed XUV energy goes into driving the outflow. This breaks down in several ways:

- **Radiative cooling is ignored.** Atomic line cooling, molecular emission, and ionisation losses can divert XUV energy away from heating the bulk gas, reducing the effective $\epsilon$. ZEPHYRUS treats $\epsilon$ as a constant input rather than computing it self-consistently. Setting $\epsilon = 1$ in particular is a non-physical upper limit on the mass-loss rate.
- **Fractionation in the outflow is not captured.** When the particle flux drops below the critical value required to drag heavy species along, the outflow becomes compositionally fractionated: hydrogen escapes preferentially and the residual atmosphere is enriched in heavy species. ZEPHYRUS removes everything in bulk. Fractionation will be implemented in the future.
- **$\epsilon$ is held constant in time.** In reality the efficiency evolves with planet mass, radius, and incident flux. Fixed-$\epsilon$ models can overestimate mass loss at late times.
- **Thermally driven outflows are diagnosed, not modelled.** An envelope reaching towards its Bondi radius loses mass through a wind powered by its own thermal energy rather than by XUV heating, and there the energy-limited rate is not an upper bound: it can fall short of direct hydrodynamic calculations by a factor of a few tens. `boiloff.py` reports whether a planet is in that regime and how much of its envelope would survive, but ZEPHYRUS does not compute the thermally driven mass-loss rate itself. The isothermal Parker-wind rate is the natural next step.
- **The tidal correction is clamped near the Roche lobe.** $K_\mathrm{tide}$ has a double root at $\xi = 1$ and the rate divides by it, so the rate grows without bound as the atmosphere approaches its Roche lobe; $\xi = 1.1$ already inflates it 83-fold. A floor keeps the result finite, but a rate computed against that floor marks the edge of validity rather than a prediction, and such a planet should be treated with the boil-off diagnostics instead.

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