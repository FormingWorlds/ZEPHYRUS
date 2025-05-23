import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from zephyrus.constants import kb, G, M_H, M_O, s2yr
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture, Planetary_surface_area, Molar_concentration_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O

# Efficiency and planet parameters
epsilon = 0.1
avogadro = 6.022e23  # Avogadro's number [particles/mole]

# Load PROTEUS data
df = pd.read_csv('/Users/emmapostolec/Downloads/runtime_helpfile.csv', sep='\t')
time = df['Time'].values            # [yr]
Teq = df['T_eqm'].values * 1000     # [K]
M_H_prot = df['H_kg_atm'].values    # [kg]
M_O_prot = df['O_kg_atm'].values    # [kg]
Fxuv = df['F_xuv'].values           # [W/m2]
M_atm_prot = df['M_atm'].values     # [kg]
Mp = df['M_planet'].values          # [kg]
Rp = df['R_obs'].values             # [m]
escape_EL = df['esc_rate_total'].values  # [kg/s]

# Pre-allocate arrays for Zephyrus outputs
n_steps = len(time)

y_H = np.zeros(n_steps)
y_O = np.zeros(n_steps)
y_H_prot = np.zeros(n_steps)
y_O_prot = np.zeros(n_steps)
x_H = np.zeros(n_steps)
x_O = np.zeros(n_steps)
x_H_prot = np.zeros(n_steps)
x_O_prot = np.zeros(n_steps)
Phi_c = np.zeros(n_steps)
Phi_l = np.zeros(n_steps)
Phi_h = np.zeros(n_steps)
Phi_t = np.zeros(n_steps)
Phi_H = np.zeros(n_steps)
Phi_O = np.zeros(n_steps)
Phi_tot = np.zeros(n_steps)
Phi_Collin_EUV = np.zeros(n_steps)
Phi_EL_model = np.zeros(n_steps)

mass_loss_H = np.zeros(n_steps)
mass_loss_O = np.zeros(n_steps)
mass_loss_tot = np.zeros(n_steps)
mass_loss_Collin_EUV = np.zeros(n_steps)
mass_loss_EL = np.zeros(n_steps)
M_atm = np.zeros(n_steps)

# Initial conditions

y_H[0] = M_H_prot[0] / M_H # Convert PROTEUS masses to moles
y_O[0] = M_O_prot[0] / M_O

x_H[0], x_O[0] = Molar_concentration_binary_mixture(n_light=y_H[0], n_heavy=y_O[0]) # Molar fractions

Phi_Collin_EUV[0] = epsilon * Fxuv[0] / (4 * (G * Mp[0] / Rp[0]**2)) # Collin_Collin_EUV collisional flux
# Fractionation fluxes
Phi_c[0], Phi_l[0], Phi_h[0], Phi_t[0] = Fractionation_binary_mixture(
    X_light=x_H[0], X_heavy=x_O[0], Mp=Mp[0], Rp=Rp[0],
    b=b_H_O(Teq[0]), T=Teq[0], mu_light=M_H/6.022e23, mu_heavy=M_O/6.022e23,
    M_light=M_H, M_heavy=M_O, Phi=Phi_Collin_EUV[0])
# Store model fluxes
Phi_H[0] = Phi_l[0] * (M_H/6.022e23)
Phi_O[0] = Phi_h[0] * (M_O/6.022e23)
Phi_tot[0] = Phi_t[0] * (x_H[0]* (M_H/6.022e23) + x_O[0]* (M_O/6.022e23))
Phi_EL_model[0] = epsilon * np.pi * Rp[0]**3 * Fxuv[0] / (G * Mp[0])
# Mass losses
area0 = Planetary_surface_area(Rp[0])
mass_loss_Collin_EUV[0] = area0 * Phi_Collin_EUV[0] * time[0] * s2yr
mass_loss_H[0] = area0 * Phi_H[0] * time[0] * s2yr
mass_loss_O[0] = area0 * Phi_O[0] * time[0] * s2yr
mass_loss_tot[0] = area0 * Phi_tot[0] * time[0] * s2yr
mass_loss_EL[0] = area0 * Phi_EL_model[0] * time[0] * s2yr
M_atm[0] = M_atm_prot[0] - mass_loss_Collin_EUV[0]

# Time integration
for i in range(1, n_steps):
    # timestep in seconds
    delta_t = (time[i] - time[i-1]) * s2yr

    # Collin_EUV collisional flux
    Phi_Collin_EUV[i] = epsilon * Fxuv[i] / (4 * (G * Mp[i] / Rp[i]**2))
    # Energy-limited loss
    Phi_EL_model[i] = epsilon * np.pi * Rp[i]**3 * Fxuv[i] / (G * Mp[i])

    # Update moles by subtracting loss
    area = Planetary_surface_area(Rp[i])
    y_H[i] = y_H[i-1] - Phi_l[i-1] * area * delta_t
    y_O[i] = y_O[i-1] - Phi_h[i-1] * area * delta_t

    # Update molar fractions
    x_H[i] = y_H[i] / (y_H[i] + y_O[i])
    x_O[i] = y_O[i] / (y_H[i] + y_O[i])

    # Fractionation
    Phi_c[i], Phi_l[i], Phi_h[i], Phi_t[i] = Fractionation_binary_mixture(
        X_light=x_H[i], X_heavy=x_O[i], Mp=Mp[i], Rp=Rp[i],
        b=b_H_O(Teq[i]), T=Teq[i], mu_light=M_H/6.022e23, mu_heavy=M_O/6.022e23,
        M_light=M_H, M_heavy=M_O, Phi=Phi_Collin_EUV[i]
    )
    # Convert to mass flux
    Phi_H[i] = Phi_l[i] * (M_H/6.022e23)
    Phi_O[i] = Phi_h[i] * (M_O/6.022e23)
    Phi_tot[i] = Phi_t[i] * (x_H[i]* (M_H/6.022e23) + x_O[i]* (M_O/6.022e23))

    # Mass losses per interval
    mass_loss_Collin_EUV[i] = area * Phi_Collin_EUV[i] * delta_t
    mass_loss_H[i] = area * Phi_H[i] * delta_t
    mass_loss_O[i] = area * Phi_O[i] * delta_t
    mass_loss_tot[i] = area * Phi_tot[i] * delta_t
    mass_loss_EL[i] = area * Phi_EL_model[i] * delta_t



    # Atmospheric mass
    M_atm[i] = M_atm_prot[i] - mass_loss_Collin_EUV[i]

print(x_H)
print(y_O)

# Plot comparison
fig, axs = plt.subplots(3, 2, figsize=(12, 8))
axs = axs.flatten()

# 1. Fluxes
axs[0].plot(time, Phi_Collin_EUV, label='Phi_Collin_EUV (Proteus)', color='green', linestyle='-')
axs[0].plot(time, Phi_c, label='Critical flux (Zephyrus)', color='orange', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_H, label='H flux (Zephyrus)', color='dodgerblue', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_O, label='O flux (Zephyrus)', color='red', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_tot, label='Total flux (Zephyrus)', color='purple', linestyle=':', linewidth=2)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlabel('Time [yr]', fontsize = 14)
axs[0].set_ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=14)
axs[0].legend()

# 2. Molar fractions
axs[1].plot(time, (M_H_prot/M_H)/(M_H_prot/M_H + M_O_prot/M_O), label='x_H (Proteus)', color='dodgerblue', linestyle='-')
axs[1].plot(time, (M_O_prot/M_O)/(M_H_prot/M_H + M_O_prot/M_O), label='x_O (Proteus)', color='red', linestyle='-')
axs[1].plot(time, x_H, label='x_H (Zephyrus)', color='dodgerblue', linestyle=':', linewidth=2)
axs[1].plot(time, x_O, label='x_O (Zephyrus)', color='red', linestyle=':', linewidth=2)
axs[1].set_xscale('log')
axs[1].set_xlabel('Time [yr]', fontsize = 14); axs[1].set_ylabel(r'$\chi_i^{atm}$ [%]', fontsize = 14)
axs[1].legend()

# 3. Mass losses
axs[2].plot(time, mass_loss_Collin_EUV/Me, color='green', label='Collin_EUV loss (Proteus)', linestyle='-')
axs[2].plot(time, mass_loss_EL/Me, color='pink', label='EUV loss (Proteus)', linestyle='-')
axs[2].plot(time, mass_loss_H/Me,color='dodgerblue',label='H loss (Zephyrus)', linestyle=':', linewidth=2)
axs[2].plot(time, mass_loss_O/Me, color='red', label='O loss (Zephyrus)', linestyle=':', linewidth=2)
axs[2].plot(time, mass_loss_tot/Me,color='purple', label='Total loss (Zephyrus)', linestyle=':', linewidth=2)
axs[2].set_xscale('log'); axs[2].set_yscale('log')
axs[2].set_xlabel('Time [yr]', fontsize = 14); axs[2].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize = 14)
axs[2].legend()

# 4. Number of moles
axs[3].plot(time, M_H_prot/M_H, label=r'$N_H$ (Proteus)', color='dodgerblue', linestyle='-')
axs[3].plot(time, M_O_prot/M_O, label=r'$N_O$ (Proteus)', color='red', linestyle='-')
axs[3].plot(time, y_H/avogadro, label=r'$N_H$ (Zephyrus)', color='dodgerblue', linestyle=':', linewidth=2)
axs[3].plot(time, y_O/avogadro, label=r'$N_O$ (Zephyrus)', color='red', linestyle=':', linewidth=2)
axs[3].plot(time, M_H_prot/M_H + M_O_prot/M_O, label=r'$N_{total}$ (Proteus)', color='purple', linestyle='-')
axs[3].plot(time, M_H_prot/M_H + M_O_prot/M_O, label=r'$N_{total}$ (Zephyrus)', color='purple', linestyle=':', linewidth=2)
axs[3].set_xscale('log')
axs[3].set_yscale('log')        
axs[3].set_xlabel('Time [yr]', fontsize = 14); axs[3].set_ylabel(r'$N_i$ [moles]', fontsize = 14)
axs[3].legend()

# 5. Atmospheric mass
axs[4].plot(time, M_atm_prot/Me, label='M_atm (Proteus)', color='darkblue', linestyle='-')
axs[4].plot(time, M_atm/Me, label='M_atm (Zephyrus)', color='darkblue', linestyle=':', linewidth=2)
axs[4].set_xscale('log')
axs[4].set_yscale('log')            
axs[4].set_xlabel('Time [yr]', fontsize = 14); axs[4].set_ylabel(r'$M_{atm}$ [M$_\oplus$]', fontsize = 14)
axs[4].legend()
            

plt.tight_layout()
plt.savefig('output/evolution_comparison_chatgpt.png', dpi=300)

