# Limitations

ZEPHYRUS implements the **energy-limited (EL) approximation** to hydrodynamic atmospheric escape, given by Eq. (1) of the [model overview](model.md), and the **giant-impact erosion scaling laws** of Kegerreis et al. (2020) and Roche et al. (2026). Both are deliberate simplifications of much richer physical problems. The most important regimes and processes the model does not cover are summarised below.

---

## What ZEPHYRUS *does* model

Two channels. The first is bulk hydrodynamic escape driven by stellar XUV irradiation, in the energy-limited approximation, with an optional tidal correction (Eq. 2 of the [model overview](model.md)). The tidal correction is defined only outside the Roche lobe, where the Hill-to-XUV radius ratio $\xi > 1$; ZEPHYRUS raises an error for $\xi \le 1$, at which point the atmosphere reaches the Roche lobe and the energy-limited approximation no longer holds. The mass-loss rate is partitioned across atmospheric species in proportion to their elemental mass mixing ratios. The second channel is the fraction of the target's atmosphere eroded by a single giant impact, evaluated by empirical scaling laws (Eq. 4 of the [model overview](model.md) for Kegerreis et al. 2020, and the Roche et al. 2026 formulation).

Everything below is **not modelled.**

---

## The impact channel (`collision.impact_loss`, `collision.mass_loss`, `collision.mass_loss_roche2026`)

The channel evaluates empirical scaling laws (Kegerreis et al. 2020 for thin atmospheres, and Roche et al. 2026 for envelope mass fractions from 1% to 20%), rather than simulating collisions directly, and inherits the scope of the simulation suites behind them (see the [model overview](model.md) for fitted domains):

- **Atmospheric mass limits.** The `kegerreis2020` fit applies to thin atmospheres of order 1% of the target mass. The `roche2026` fit extends calibration over $f_\mathrm{atm} \in [0.01, 0.20]$, but envelopes outside that range require extrapolation. An airless target ($f_\mathrm{atm} = 0$) returns $X = 0$ immediately, while any non-zero $f_\mathrm{atm} < 10^{-6}$ is evaluated at the $10^{-6}$ stability clamp (creating a jump at zero and a plateau below $10^{-6}$).
- **Target-side loss only.** Any atmosphere the impactor carries, and any volatile delivery into the merged body, is outside the function. The underlying simulations of Kegerreis et al. (2020) show a slow grazing collision with an atmosphere-hosting impactor can leave the target with about 85% of the two bodies' combined atmospheres, so treating the impactor as bare is a caller-side assumption.
- **Simulation-to-law scatter.** For Kegerreis et al. (2020), the scatter between simulations and the scaling law is a relative deviation: the median relative deviation is 9%, rising to about 20% relative deviation for slow, head-on collisions where outcomes are chaotic. For Roche et al. (2026), misfits are reported as absolute differences on the loss fraction ($|X_\mathrm{sim} - X_\mathrm{calc}|$): over their 790-row fitting set (296 of the 300 new simulations plus 494 runs at $f_\mathrm{atm} = 0.05$ from Roche et al. 2025; Fig. 3), the scaling law has mean and maximum absolute misfits of 0.0396 and 0.2224 respectively.
- **Physical simplifications.** The Roche et al. (2026, Sect. 4.1) scaling law does not account for pre-impact planetary rotation, surface liquid water oceans, thermal evolution, or a core mass fraction other than about 0.3 (Sect. 2.1). Magma-envelope miscibility at high pressures without a sharp boundary can cause the scaling law to overestimate the loss of massive envelopes on young planets (Sect. 4.1). The law covers the immediate shock- and vapour-plume-driven loss only; it neglects later thermally driven loss (an outflow driven by heat from the post-impact interior, Biersteker & Schlichting 2021), so for primordial H2-He envelopes the total loss can be higher (Roche et al. 2026, Sect. 4.1); that later loss becomes negligible for envelopes of higher mean molecular weight. All simulations assume H2-He envelopes; heavier atmospheres ($\mathrm{CO}, \mathrm{CO}_2$) are less susceptible to shock-driven removal, so shock-driven loss is an upper limit for a given atmosphere mass (Sect. 4.3); for PROTEUS atmospheres of higher mean molecular weight ($\mathrm{H}_2\mathrm{O}, \mathrm{CO}_2, \mathrm{O}_2$), the H2-He loss fractions likewise serve as upper limits.
- **Fit-domain flags.** `collision.impact_loss` records out-of-range parameters in `flags` for the `roche2026` law (`'f_atm'`, `'M_t_earth'`, `'gamma'`, `'b'`, `'R_ratio'`, `'v_ratio'`). It also records diagnostic flags `'v_sub_escape'` for speeds below 0.99 mutual escape speed and `'X_FF_zero_energy'` when the zero-energy far-field loss exceeds zero (which can occur inside the fitted range). The direct functions `collision.mass_loss` and `collision.mass_loss_roche2026` evaluate without warning flags; staying inside fitted bounds is the caller's responsibility.

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