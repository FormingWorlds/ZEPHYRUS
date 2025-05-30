import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from zephyrus.constants import kb, G, M_H, M_O, s2yr
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture, Planetary_surface_area, Molar_concentration_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O
from zephyrus.escape import EL_escape

# Initial condition
epsilon = 0.1  

# Binary mixture of H/O
M_H  = 0.001008 # molar mass of H [kg/mol]
M_O = 0.015999 # molar mass of O [kg/mol]
avogadro = 6.022e23 # Avogadro's number [particles/mole]
mu_H = M_H/avogadro # atomic mass of H [kg/atom]
mu_O = M_O/avogadro # atomic mass of O [kg/atom]

# Load and extract CSV data from the case 834 from the 7 parameter grid on Habrok
df = pd.read_csv('/Users/emmapostolec/Downloads/runtime_helpfile.csv', sep='\t')
time_proteus    = df['Time'].values          # [yr]
a_proteus       = df['semimajorax'].values   # [m]
Teq_proteus     = df['T_eqm'].values         # [K]
mass_H_proteus  = df['H_kg_atm'].values      # [kg]
mass_O_proteus  = df['O_kg_atm'].values      # [kg]
Fxuv_proteus    = df['F_xuv'].values         # [W/m2]
Mstar_proteus = df['M_star'].values          # [kg]
M_atm_proteus   = df['M_atm'].values         # [kg]
Mp_proteus      = df['M_planet'].values      # [kg]
Rp_proteus      = df['R_obs'].values         # [m]
Rxuv_proteus    = df['R_xuv'].values         # [m]
escape_EL_proteus = df['esc_rate_total'].values      # [kg/s]

# Initializaton of arrays
N_light = np.zeros(len(time_proteus))   
N_heavy = np.zeros(len(time_proteus))
x_light = np.zeros(len(time_proteus))
x_heavy = np.zeros(len(time_proteus))
N_total = np.zeros(len(time_proteus))
Phi_light = np.zeros(len(time_proteus))
Phi_heavy = np.zeros(len(time_proteus)) 
Phi_total = np.zeros(len(time_proteus))
Phi_critical = np.zeros(len(time_proteus))
Phi_EUV_Collin = np.zeros(len(time_proteus))
Phi_EL = np.zeros(len(time_proteus))
Phi_EL_kg = np.zeros(len(time_proteus))
M_atm = np.zeros(len(time_proteus))
mass_loss_H = np.zeros(len(time_proteus))
mass_loss_O = np.zeros(len(time_proteus))
mass_loss_total = np.zeros(len(time_proteus))
mass_loss_EUV = np.zeros(len(time_proteus))


# Loop over time steps
# before the loop, seed your “working” mole‐counts from the arrays
n_light = None
n_heavy = None

for i, t in enumerate(time_proteus):
    if i == 0:
        # initial moles & concentrations
        N_light[0] = mass_H_proteus[0] / M_H
        N_heavy[0] = mass_O_proteus[0] / M_O
        x_light[0], x_heavy[0] = Molar_concentration_binary_mixture(n_light=N_light[0], n_heavy=N_heavy[0])

        # store them in working vars
        n_light = N_light[0]
        n_heavy = N_heavy[0]

        # initial EUV‐driven mass flux
        Phi_EUV_Collin[0] = (
            epsilon * Fxuv_proteus[0]
        ) / (4 * (G * Mp_proteus[0] / Rp_proteus[0]**2))

        # fractionation at t=0
        Phi_c, Phi_l, Phi_h, Phi_tot = Fractionation_binary_mixture(
            X_light=x_light[0], X_heavy=x_heavy[0],
            Mp=Mp_proteus[0], Rp=Rp_proteus[0],
            b=b_H_O(Teq_proteus[0]), T=Teq_proteus[0],
            mu_light=mu_H, mu_heavy=mu_O,
            M_light=M_H, M_heavy=M_O,
            Phi=Phi_EUV_Collin[0]
        )

        # seed loss arrays at index 0 (note the [0] on every array you reference)
        area0 = Planetary_surface_area(Rp_proteus[0])
        mass_loss_H[0]     = area0 * Phi_l      * mu_H * time_proteus[0] / s2yr
        mass_loss_O[0]     = area0 * Phi_h      * mu_O * time_proteus[0] / s2yr
        mass_loss_total[0] = area0 * Phi_tot    * (x_light[0]*mu_H + x_heavy[0]*mu_O) * time_proteus[0] / s2yr
        mass_loss_EUV[0]   = area0 * Phi_EUV_Collin[0] * time_proteus[0] / s2yr
        M_atm[0]           = M_atm_proteus[0] - mass_loss_EUV[0]

        # seed the rest of your flux arrays
        Phi_light[0]    = Phi_l * mu_H
        Phi_heavy[0]    = Phi_h * mu_O
        Phi_total[0]    = Phi_tot * (x_light[0]*mu_H + x_heavy[0]*mu_O)
        Phi_critical[0] = Phi_c
        Phi_EL[0]       = EL_escape(
            tidal_contribution=False,
            a=a_proteus[0], e=0.0,
            Mp=Mp_proteus[0], Ms=Mstar_proteus[0], epsilon=epsilon,
            Rp=Rp_proteus[0], Rxuv=Rxuv_proteus[0],
            Fxuv=Fxuv_proteus[0], scaling=2
        )

    else:
        # Time step 
        delta_t = t/len(time_proteus) # timestep [s]

        # subtract working mole‐counts
        n_light -= Phi_l * Planetary_surface_area(Rp_proteus[i]) * delta_t
        n_heavy -= Phi_h * Planetary_surface_area(Rp_proteus[i]) * delta_t
        # update the arrays
        N_light[i] = n_light
        N_heavy[i] = n_heavy
        x_light[i] = n_light / (n_light + n_heavy)
        x_heavy[i] = n_heavy / (n_light + n_heavy)

        # Phi
        Phi_EUV_Collin_i = (epsilon * Fxuv_proteus[i]) / (4 * (G * Mp_proteus[i] / Rp_proteus[i]**2))  # mass flux [kg/m2/s]
        Phi_EL_i = EL_escape(
            tidal_contribution=False, a=a_proteus[i], e=0.0,
            Mp=Mp_proteus[i], Ms=Mstar_proteus[i], epsilon=epsilon,
            Rp=Rp_proteus[i], Rxuv=Rxuv_proteus[i], Fxuv=Fxuv_proteus[i], scaling=2
        )
        
        # Binary diffusion coefficient for H and O
        b_HO_i = b_H_O(Teq_proteus[i])  # Binary diffusion coefficient for H and O

        # Mass loss in the atmosphere
        mass_loss_EUV[i] = Phi_EUV_Collin[i] * Planetary_surface_area(Rp_proteus[i]) * delta_t / s2yr # Mass loss due to EUV [kg]
        M_atm[i] = M_atm_proteus[i] - mass_loss_EUV[i] # Mass of the atmosphere [kg]

        # Fractionation fluxes
        Phi_c, Phi_l, Phi_h, Phi_tot = Fractionation_binary_mixture(X_light=x_light[i], X_heavy=x_heavy[i], Mp=Mp_proteus[i], Rp=Rp_proteus[i], b=b_HO_i, T=Teq_proteus[i], mu_light=mu_H, mu_heavy=mu_O, M_light=M_H, M_heavy=M_O, Phi=Phi_EUV_Collin_i)
        
        Phi_critical[i] = Phi_c             # Critical mass flux [kg/m2/s]
        Phi_light[i] = Phi_l * mu_H         # H flux [kg/m2/s]
        Phi_heavy[i] = Phi_h * mu_O         # O flux [kg/m2/s]
        Phi_total[i] = Phi_tot * (x_light[i]*mu_H + x_heavy[i]*mu_O)    # Total flux [kg/m2/s]
        Phi_EUV_Collin[i] = Phi_EUV_Collin_i    # Mass flux  [kg/m2/s]
        Phi_EL[i] = Phi_EL_i # mass loss [kg/s]

        mass_loss_H[i] = Planetary_surface_area(Rp_proteus[i]) * Phi_light[i] * mu_H * delta_t / s2yr  # Mass loss of H [kg]
        mass_loss_O[i] = Planetary_surface_area(Rp_proteus[i]) * Phi_heavy[i] * mu_O * delta_t  / s2yr # Mass loss of O [kg]
        mass_loss_total[i] = Planetary_surface_area(Rp_proteus[i]) * Phi_total[i] * (x_light[i]*mu_H + x_heavy[i]*mu_O) * delta_t / s2yr # Total mass loss [kg]
        Phi_EL_kg[i] =  Phi_EL_i / Rp_proteus[i]**2 # Mass loss due to EL [kg]


print('Phi_Critical:', Phi_critical)

# Plot the results
fig, axs = plt.subplots(2, 2, figsize=(10, 6))
axs = axs.flatten()

# 1. Evolution of the mass flux
axs[0].plot(time_proteus, Phi_EUV_Collin, label=r'$\Phi_{EUV}$ (IsoFATE)', color='purple', linestyle='-')
axs[0].plot(time_proteus, Phi_light, label='H flux Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
axs[0].plot(time_proteus, Phi_heavy, label='O flux Zephyrus', color='red', linestyle=':', linewidth=2)
axs[0].plot(time_proteus, Phi_total, label='Total flux Zephyrus', color='purple', linestyle=':', linewidth=2)
axs[0].plot(time_proteus, Phi_critical, label='Critical flux', color='orange', linestyle='--')
axs[0].set_xlabel('Time [yr]', fontsize=14)
axs[0].set_ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=14)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend()

# 2. Evolution of the number of moles
axs[1].plot(time_proteus, mass_H_proteus/M_H, label=r'$N_H$ PROTEUS', color='dodgerblue', linestyle='-')
axs[1].plot(time_proteus, mass_O_proteus/M_O, label=r'$N_O$ PROTEUS', color='red', linestyle='-')
axs[1].plot(time_proteus, N_light, label=r'$N_H$ Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
axs[1].plot(time_proteus, N_heavy, label=r'$N_O$ Zephyrus', color='red', linestyle=':', linewidth=2)   
axs[1].set_xlabel('Time [yr]', fontsize=14)
axs[1].set_ylabel(r'$N_i$ [moles]', fontsize=14)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend()

# 3. Evolution of the mass loss
axs[2].plot(time_proteus, mass_loss_H/Me, label='H loss Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
axs[2].plot(time_proteus, mass_loss_O/Me, label='O loss Zephyrus', color='red', linestyle=':', linewidth=2)
axs[2].plot(time_proteus, mass_loss_total/Me, label='Total loss Zephyrus', color='purple', linestyle=':', linewidth=2)
axs[2].plot(time_proteus, mass_loss_EUV/Me, label='EUV loss IsoFATE', color='green', linestyle='-')
axs[2].plot(time_proteus, Phi_EL_kg/Me, label='EL loss Zephyrus', color='pink', linestyle='-')
axs[2].set_xlabel('Time [yr]', fontsize=14)
axs[2].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=14)
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].legend()


# 4. Evolution of the mole fraction
axs[3].plot(time_proteus, mass_H_proteus/M_H*100, label=r'$x_H$ PROTEUS', color='dodgerblue')
axs[3].plot(time_proteus, mass_O_proteus/M_O*100, label=r'$x_O$ PROTEUS', color='red')
axs[3].plot(time_proteus, x_light*100, label=r'$x_H$ Zephyrus', color='dodgerblue')
axs[3].plot(time_proteus, x_heavy*100, label=r'$x_O$', color='red')
axs[3].set_xlabel('Time [yr]', fontsize=14)
axs[3].set_ylabel(r'$\chi_i^{atm}$ [%]', fontsize=14)
axs[3].set_xscale('log')
axs[3].legend()


# Adjust layout and save
plt.tight_layout()
plt.savefig('plot/comparison_proteus_4panels.pdf', dpi=300)
plt.close()
