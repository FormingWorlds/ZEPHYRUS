import numpy as np
from zephyrus.constants import G, s2yr, M_H, M_O
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O
from zephyrus.fractionation import Number_flux
from zephyrus.fractionation import Diffusion_flux, Phi_critical
from zephyrus.fractionation import Scale_height_single_species, Scale_height_binary_mixture, Molar_concentration_binary_mixture
from zephyrus.fractionation import Acceleration_of_gravity
import matplotlib.pyplot as plt

# --- fixed planetary & species parameters ---
Te   = 2880.0                       # Temperature [K]
epsilon = 1.0                       # heating efficiency
n_H  = 1.0e20                       # mol H
n_O  = 1.0e10                       # mol O
Mp, Rp = Me, Re                     # Earth analog
ge     = Acceleration_of_gravity(Mp=Mp, Rp=Rp)  # g [m/s²]
b_HO   = b_H_O(Te)                  # binary diffusion coefficient

# scale heights & concentrations (constant w/ T, g, M)
H_H = Scale_height_single_species(T=Te, g=ge, M_i=M_H)
H_O = Scale_height_single_species(T=Te, g=ge, M_i=M_O)
X_H, X_O = Molar_concentration_binary_mixture(n_light=n_H, n_heavy=n_O)
H_mix = Scale_height_binary_mixture(T=Te, g=ge,
                                    X_light=X_H, M_light=M_H,
                                    X_heavy=X_O, M_heavy=M_O)

# --- XUV flux evolution function ---
def Fxuv(t, F0, t_sat=5e8, beta=-1.23):
    """t in years; F0 in W/m²; returns W/m²"""
    if t < t_sat:
        return F0
    return F0 * (t / t_sat)**beta

# --- time grid ---
t_max   = 5e9      # years
n_steps = 500
times   = np.linspace(0.0, t_max, n_steps)  # in years
F0      = 1e6    # your chosen initial XUV flux [W/m²]

# prepare storage
Phi_vals         = np.zeros_like(times)
Phi_diff_H_vals  = np.zeros_like(times)
Phi_diff_O_vals  = np.zeros_like(times)
Phi_crit_vals    = np.zeros_like(times)
Phi_light_vals   = np.zeros_like(times)
Phi_heavy_vals   = np.zeros_like(times)
Phi_dl_vals      = np.zeros_like(times)
Phi_df_vals      = np.zeros_like(times)
Phi_c_vals       = np.zeros_like(times)
Phi_l_vals       = np.zeros_like(times)
Phi_h_vals       = np.zeros_like(times)
Phi_tot_vals     = np.zeros_like(times)

for i, t in enumerate(times):
    # 1) XUV at this time
    Fx = Fxuv(t, F0)

    # 2) Translate XUV power flux to particle flux Φ [particles/m²/s]
    #    using ε⋅F_XUV / (4 ⋅ g)
    Phi_t = epsilon * Fx / (4.0 * ge)
    Phi_vals[i] = Phi_t

    # 3) diffusion fluxes
    Phi_diff_H_vals[i] = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_H)
    Phi_diff_O_vals[i] = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_O)

    # 4) critical flux
    Phi_crit_vals[i] = Phi_critical(T=Te, g=ge,
                                    M_light=M_H, M_heavy=M_O,
                                    b=b_HO, x_light=X_H,
                                    n_light=n_H, n_heavy=n_O)

    # 5) upwelling number fluxes
    Phi_Lt, Phi_Hv = Number_flux(Phi_t, Phi_crit_vals[i],
                                 n_light=n_H, n_heavy=n_O,
                                 M_light=M_H, M_heavy=M_O,
                                 X_light=X_H, X_heavy=X_O,
                                 Phi_diffusion_light=Phi_diff_H_vals[i],
                                 Phi_diffusion_heavy=Phi_diff_O_vals[i])
    Phi_light_vals[i] = Phi_Lt
    Phi_heavy_vals[i] = Phi_Hv

    # 6) full fractionation solver
    fd = Fractionation_binary_mixture(n_light=n_H, n_heavy=n_O,
                                      Mp=Mp, Rp=Rp, b=b_HO,
                                      T=Te, M_light=M_H, M_heavy=M_O,
                                      Phi=Phi_t)
    Phi_dl_vals[i], Phi_df_vals[i], Phi_c_vals[i], \
    Phi_l_vals[i],  Phi_h_vals[i],  Phi_tot_vals[i] = fd

# --- now you have time‐series arrays:
# times, Phi_vals, Phi_light_vals, Phi_heavy_vals, Phi_tot_vals, etc.
# You can plot or save to CSV as you wish.
import pandas as pd
df = pd.DataFrame({
    'time_yr':         times,
    'Phi_input':       Phi_vals,
    'Phi_diff_H':      Phi_diff_H_vals,
    'Phi_diff_O':      Phi_diff_O_vals,
    'Phi_crit':        Phi_crit_vals,
    'Phi_light':       Phi_light_vals,
    'Phi_heavy':       Phi_heavy_vals,
    'Phi_dl':          Phi_dl_vals,
    'Phi_df':          Phi_df_vals,
    'Phi_c':           Phi_c_vals,
    'Phi_l':           Phi_l_vals,
    'Phi_h':           Phi_h_vals,
    'Phi_total':       Phi_tot_vals,
})
df.to_csv('fractionation_time_evolution.csv', index=False)
print("Done – results written to fractionation_time_evolution.csv")

# 2) Plot input particle flux Φ vs time
plt.figure()
plt.plot(df['time_yr'], df['Phi_input'])
plt.xlabel('Time [years]')
plt.ylabel('Φ input [particles/m²/s]')
plt.title('Input Particle Flux vs Time')
plt.tight_layout()
plt.show()

# 3) Plot light and heavy species number fluxes vs time
plt.figure()
plt.plot(df['time_yr'], -df['Phi_light'], label='Light Species Flux')
plt.plot(df['time_yr'], df['Phi_heavy'], label='Heavy Species Flux')
plt.yscale('log')
plt.xlabel('Time [years]')
plt.ylabel('Number Flux [particles/m²/s]')
plt.title('Number Flux of Light and Heavy Species over Time')
plt.legend()
plt.tight_layout()
plt.show()

# 4) Plot total fractionation flux vs time
plt.figure()
plt.plot(df['time_yr'], df['Phi_total'])
plt.xlabel('Time [years]')
plt.ylabel('Total Fractionation Flux [particles/m²/s]')
plt.title('Total Fractionation Flux vs Time')
plt.tight_layout()
plt.show()