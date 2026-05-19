# First run

This tutorial walks through your first ZEPHYRUS computation: the energy-limited atmospheric mass-loss history of an Earth-mass planet around a Sun-like star. By the end you will have produced a figure showing how the escape rate evolves over the lifetime of the system.

If you haven't installed ZEPHYRUS yet, follow the [installation guide](../How-to/installation.md) first and make sure you've run `mors download all` to download the stellar evolution tracks.

!!! info "What you'll do"
    - Load a stellar evolution track with MORS
    - Derive the XUV flux incident on a planet from the stellar luminosity
    - Call `EL_escape` for an evolutionary track
    - Plot the result with secondary units

---

## The full script

The complete example lives at `examples/demo_earth/demo_earth.py` in the repository. Start by creating the output directory it writes to:

```sh
mkdir -p output
```

!!! warning "Output directory"
    The script writes to `output/` relative to the **current working directory.** If the directory doesn't exist, `plt.savefig` will fail. Make sure you run the example from the same level as the output directory, or change the path. 

Then run it:

```sh
python examples/demo_earth/demo_earth.py
```

After a few seconds you should see `output/demo_earth_escape_vs_time_MORS.pdf`.

The rest of this page steps through what the script does.

---

## Step 1: Imports and system parameters

```python
import numpy as np
import matplotlib.pyplot as plt
import mors

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.escape import EL_escape
```

The two `import *` statements pull in the unit conversions (`au2m`, `au2cm`, `ergcm2stoWm2`, `s2yr`) and the Earth/Sun reference values (`Me`, `Re`, `a_earth`, `e_earth`) used throughout the script.

We then define the system:

```python
M_star                   = 1.0              # Star mass in solar mass                [M_sun, kg]
Omega_sun                = 1.0              # Solar rotation rate                    [rad s-1]
semi_major_axis_cm       = a_earth*au2cm    # Planetary semi-major axis              [cm]

tidal_contribution       = False            # Tidal correction factor                [dimensionless]
semi_major_axis          = a_earth*au2m     # Planetary semi-major axis              [m]
eccentricity             = e_earth          # Planetary eccentricity                 [dimensionless]
M_planet                 = Me               # Planetary mass                         [kg]
epsilon                  = 0.15             # Escape efficiency factor               [dimensionless]  
R_earth                  = Re               # Planetary radius                       [m]
Rxuv                     = Re               # XUV planetary radius                   [m]
```

A few choices worth noting:

- `epsilon = 0.15` is the conservative value adopted by Kasting & Pollack (1983) [^kp83] and used in many subsequent rocky-planet studies. The full plausible range is roughly $0.1$–$0.6$.
- `Rxuv = Re` sets the XUV-absorbing radius equal to the planetary radius. This is a lower bound on the mass-loss rate; allowing $R_\mathrm{XUV} > R_p$ would increase escape.
- `tidal_contribution = False` ignores the Roche-lobe enhancement of escape. At 1 au this is a negligible correction, but for close-in planets it should be enabled.

---

## Step 2: Compute the XUV flux at the planet

```python
star            = mors.Star(Mstar=M_star, Omega=Omega_sun)      # Load the stellar evolution tracks from MORS
Age_star        = star.Tracks['Age']                            # Stellar age                                      [Myr]
Lxuv_star       = star.Tracks['Lx']+star.Tracks['Leuv']         # XUV luminosity                                   [erg s-1]
Fxuv_star       = Lxuv_star/(4*np.pi*(semi_major_axis_cm)**2)   # XUV flux                                         [erg s-1 cm-2]
Fxuv_star_SI    = Fxuv_star*ergcm2stoWm2                        # XUV flux    
```

`mors.Star` loads the full rotational and high-energy evolutionary track for a $1\,M_\odot$ star rotating at the solar rate. The X-ray and EUV luminosities are summed to give the total XUV luminosity, which is then converted to an irradiation flux at the planet's orbital distance using the inverse-square law. The final SI conversion uses the `ergcm2stoWm2` factor from `zephyrus.constants`.

`Age_star`, `Lxuv_star`, and `Fxuv_star_SI` are all arrays of the same length — one entry per timestep of the MORS track.

---

## Step 3: Compute the escape rate

```python
escape = EL_escape(
    tidal_contribution, semi_major_axis, eccentricity,
    M_planet, M_star, epsilon,
    R_earth, Rxuv, Fxuv_star_SI,
)
```
!!! warning "SI units"
    `EL_escape` expects all inputs in **SI** units. The XUV flux must be in W m$^{-2}$, which is why the script applies `ergcm2stoWm2` after MORS returns cgs luminosities.

Because `Fxuv_star_SI` is a NumPy array, `escape` comes out as an array of the same length: the mass-loss rate in kg s$^{-1}$ at each age. All other arguments are scalars, so the planet itself is static. Only the stellar flux varies.

---

## Step 4: Plot

```python
fig, ax1 = plt.subplots(figsize=(10, 8))

ax1.loglog(Age_star, escape, '-', color='orange', label='MORS stellar evolution tracks')
ax1.set_xlabel('Time [Myr]', fontsize=15)
ax1.set_ylabel(r'Mass loss rate [kg $s^{-1}$]', fontsize=15)
ax1.set_title('Zephyrus : EL escape for Sun-Earth system', fontsize=15)
ax1.grid(alpha=0.4)
ax1.legend()
ax1.set_yscale('log')

ax2 = ax1.twinx()
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0]/ s2yr) / Me,(ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'Mass loss rate [$M_{\oplus}$ $yr^{-1}$]', fontsize=15)

plt.savefig('output/demo_earth_escape_vs_time_MORS.pdf', dpi=180)
```

The right-hand axis converts the SI mass-loss rate (kg s$^{-1}$) into Earth masses per year ($M_\oplus$ yr$^{-1}$) so you can read off how much of the planet is lost per unit time in more intuitive units. The conversion uses `s2yr` and `Me`, both from the ZEPHYRUS imports.

You should see a curve that peaks near $5 \times 10^5$ kg s$^{-1}$ at 1 Myr, drops by roughly a factor of 30 over the first ~30 Myr, plateaus through ~30–200 Myr, and then declines slowly to a few × 10$^3$ kg s$^{-1}$ by the end of the main sequence. On the right axis these correspond to $\sim 10^{-12}$ down to $\sim 10^{-13}\,M_\oplus$ yr$^{-1}$.

---

## A second run: a close-in planet with tidal correction

The same machinery handles much more extreme cases. Copy `demo_earth.py` to a new file and change three lines:

```python
semi_major_axis_cm = 0.05 * au2cm   # was: a_earth * au2cm
semi_major_axis    = 0.05 * au2m    # was: a_earth * au2m
tidal_contribution = True           # was: False
```

This places the Earth-mass planet at 0.05 au; well inside the Hill regime where tidal enhancement matters. 

Also make sure you use a stellar mass in kg when calculating the escape rate:

```python
M_star_kg           = 1.0  * Ms     # Star mass in kg              

# Adjust the escape call
escape = EL_escape(
    tidal_contribution, semi_major_axis, eccentricity,
    M_planet, M_star_kg, epsilon,
    R_earth, Rxuv, Fxuv_star_SI,
)
```

!!! warning "Stellar mass units"
    When tides are enabled, the stellar mass needs to be in kg, just like `M_planet`. MORS, however, expects the stellar mass in solar units, which is why an extra variable called `M_star_kg` is required.

Re-running produces a curve that is roughly $(1/0.05)^2 = 400$ times higher than the 1 AU case at every age, plus an additional boost from $K_\mathrm{tide} < 1$.

Things to try from here:

- **Sweep $\epsilon$**: loop `epsilon` over `[0.1, 0.3, 0.5, 1.0]` and plot all four curves on the same axes. 
- **Sweep semi-major axis**: keep $\epsilon = 0.15$ fixed and try `[0.05, 0.1, 0.5, 1.0]` au. The integral of each curve over the stellar lifetime gives the total atmospheric mass lost; compare it to one Earth atmosphere ($M_\mathrm{atm,\oplus} \approx 5.15 \times 10^{18}$ kg, available as `Me_atm`).
- **Vary $R_\mathrm{XUV}$**: try `Rxuv = 1.5 * Re` to see how an extended XUV-absorbing region amplifies the mass-loss rate. 

[^kp83]: Kasting, J. F., & Pollack, J. B. (1983). Loss of water from Venus. I. Hydrodynamic escape of hydrogen. *Icarus, 53*(3), 479–508. https://doi.org/10.1016/0019-1035(83)90212-9