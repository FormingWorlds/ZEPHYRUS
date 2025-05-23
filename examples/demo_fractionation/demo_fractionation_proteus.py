import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from zephyrus.constants import kb, G, M_H, M_O, s2yr
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture, Planetary_surface_area, Molar_concentration_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O

epsilon = 0.1  
Mp = Me
Rp = Re

M_H  = 0.001008 # molar mass of H [kg/mol]
M_O = 0.015999 # molar mass of O [kg/mol]
avogadro = 6.022e23 # Avogadro's number [particles/mole]
mu_H = M_H/avogadro # atomic mass of H [kg/atom]
mu_O = M_O/avogadro # atomic mass of O [kg/atom]

# Load and extract CSV data from the case 834 from the 7 parameter grid on Habrok
df = pd.read_csv('/Users/emmapostolec/Downloads/runtime_helpfile.csv', sep='\t')
time_proteus   = df['Time'].values          # [yr]
Teq_proteus   = df['T_eqm'].values *1000          # [K]
mass_H_proteus = df['H_kg_atm'].values     # [kg]
mass_O_proteus = df['O_kg_atm'].values      # [kg]
Fxuv_proteus   = df['F_xuv'].values          # [W/m2]
M_atm_proteus   = df['M_atm'].values          # [kg]
Mp_proteus   = df['M_planet'].values          # [kg]
Rp_proteus   = df['R_obs'].values          # [m]
escape_EL_proteus = df['esc_rate_total'].values      # [kg/s]

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
for i,t in enumerate(time_proteus):

    # Initial conditions
    if i == 0:
        # Initial conditions at time step 0
        y_light_0 = mass_H_proteus[0]/M_H # Initial numebr of moles of H [moles]
        y_heavy_0 = mass_O_proteus[0]/M_O  # Initial number of moles of O [kg]
        x_light_0, x_heavy_0 = Molar_concentration_binary_mixture(n_light=y_light_0, n_heavy=y_heavy_0)  # Initial molar concentration of H [1]
        Phi_EUV_Collin_0 = (epsilon * Fxuv_proteus[0]) / (4 * (G * Mp_proteus[0] / Rp_proteus[0]**2))  # mass flux [kg/m2/s]
        Phi_c_0, Phi_l_0, Phi_h_0, Phi_tot_0 = Fractionation_binary_mixture(X_light=x_light_0, X_heavy=x_heavy_0, Mp=Mp_proteus[0], Rp=Rp_proteus[0], b=b_H_O(Teq_proteus[0]), T=Teq_proteus[0], mu_light=mu_H, mu_heavy=mu_O, M_light=M_H, M_heavy=M_O, Phi=Phi_EUV_Collin_0)
        mass_loss_H_0 = Planetary_surface_area(Rp[0]) * Phi_l_0 * mu_H * time_proteus[0] / s2yr  # Mass loss of H [kg]
        mass_loss_O_0 = Planetary_surface_area(Rp[0]) * Phi_h_0 * mu_O * time_proteus[0]  / s2yr # Mass loss of O [kg]
        mass_loss_total_0 = Planetary_surface_area(Rp[0]) * Phi_tot_0 * time_proteus[0] / s2yr # Total mass loss [kg]
        mass_loss_EUV_0 = Planetary_surface_area(Rp[0]) * Phi_EUV_Collin_0 * time_proteus[0] / s2yr # Mass loss due to EUV [kg]
        M_atm_0 = M_atm_proteus[0]  - mass_loss_EUV_0 # Mass of the atmosphere [kg]

    else:
        # Time step 
        delta_t = t/len(time_proteus) # timestep [s]

        # Phi
        Phi_EUV_Collin_i = (epsilon * Fxuv_proteus[i]) / (4 * (G * Mp_proteus[i] / Rp_proteus[i]**2))  # mass flux [kg/m2/s]
        Phi_EL_i = (epsilon * np.pi * Rp_proteus[i]**3 * Fxuv_proteus[i]) / (G * Mp_proteus[i])  # mass loss rate [kg/s]
        
        # Binary diffusion coefficient for H and O
        b_HO_i = b_H_O(Teq_proteus[i])  # Binary diffusion coefficient for H and O

        # Number of moles and molar concentration of H and O
        y1 -= Phi_l*Planetary_surface_area(Rp[i])*delta_t
        y2 -= Phi_h*Planetary_surface_area(Rp[i])*delta_t
        x1 = y1/(y1 + y2)
        x2 = y2/(y1 + y2)

        # Mass loss in the atmosphere
        mass_loss_EUV[i] = Phi_EUV_Collin[i] * Planetary_surface_area(Rp[i])  * delta_t / s2yr # Mass loss due to EUV [kg]
        M_atm[i] = M_atm_proteus[i] - mass_loss_EUV[i] # Mass of the atmosphere [kg]

        # Fractionation fluxes
        Phi_c, Phi_l, Phi_h, Phi_tot = Fractionation_binary_mixture(X_light=x1, X_heavy=x2, Mp=Mp[i], Rp=Rp[i], b=b_HO_i, T=Teq_proteus[i], mu_light=mu_H, mu_heavy=mu_O, M_light=M_H, M_heavy=M_O, Phi=Phi_EUV_Collin_i)
        

        N_light[i] = x1  
        N_heavy[i] = x2
        x_light[i] = x1
        x_heavy[i] = x2
        Phi_critical[i] = Phi_c             # Critical mass flux [kg/m2/s]
        Phi_light[i] = Phi_l * mu_H         # H flux [kg/m2/s]
        Phi_heavy[i] = Phi_h * mu_O         # O flux [kg/m2/s]
        Phi_total[i] = Phi_tot * (x_light[i]*mu_H + x_heavy[i]*mu_O)    # Total flux [kg/m2/s]
        Phi_EUV_Collin[i] = Phi_EUV_Collin_i    # Mass flux  [kg/m2/s]
        Phi_EL[i] = Phi_EL_i # mass loss [kg/s]
        mass_loss_H[i] = Planetary_surface_area(Rp[i]) * Phi_light[i] * delta_t / s2yr  # Mass loss of H [kg]
        mass_loss_O[i] = Planetary_surface_area(Rp[i]) * Phi_heavy[i] * delta_t  / s2yr # Mass loss of O [kg]
        mass_loss_total[i] = Planetary_surface_area(Rp[i]) * Phi_total[i] * delta_t / s2yr # Total mass loss [kg]
        
        Phi_EL_kg[i] = Planetary_surface_area(Rp[i]) * Phi_EL[i] / Rp[i] * delta_t / s2yr # Mass loss due to EL [kg]

# Plot the results

fig, axs = plt.subplots(2, 2, figsize=(10, 6))
axs = axs.flatten()
#plt.title('Output PROTEUS simulation', fontsize=16)

# 1. Evolution of the number of moles
axs[1].plot(time, N_light, label=r'$N_H$', color='dodgerblue')
axs[1].plot(time, N_heavy, label=r'$N_O$', color='red')   
axs[1].plot(time, N_total, label=r'$N_{Total}$', color='purple')
axs[1].set_xlabel('Time [yr]', fontsize=14)
axs[1].set_ylabel(r'$N_i$ [moles]', fontsize=14)
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].legend()

# 2. Evolution of the mass flux
axs[0].plot(time, Phi_EUV_Collin, label=r'$\Phi_{EUV}$ (IsoFATE)', color='purple', linestyle='-')
axs[0].plot(time, Phi_light, label='H flux Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_heavy, label='O flux Zephyrus', color='red', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_total, label='Total flux Zephyrus', color='purple', linestyle=':', linewidth=2)
axs[0].plot(time, Phi_critical, label='Critical', color='orange', linestyle='--')
axs[0].set_xlabel('Time [yr]', fontsize=14)
axs[0].set_ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=14)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].legend()

# 3. Evolution of the mole fraction
axs[3].plot(time, x_light*100, label=r'$x_H$', color='dodgerblue')
axs[3].plot(time, x_heavy*100, label=r'$x_O$', color='red')
axs[3].set_xlabel('Time [yr]', fontsize=14)
axs[3].set_ylabel(r'$\chi_i^{atm}$ [%]', fontsize=14)
axs[3].set_xscale('log')
axs[3].legend()

# 4. Evolution of the mass loss
axs[2].plot(time, mass_loss_H/Me, label='H loss Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
axs[2].plot(time, mass_loss_O/Me, label='O loss Zephyrus', color='red', linestyle=':', linewidth=2)
axs[2].plot(time, mass_loss_total/Me, label='Total loss Zephyrus', color='purple', linestyle=':', linewidth=2)
axs[2].plot(time, mass_loss_EUV/Me, label='EUV loss IsoFATE', color='green', linestyle='-')
axs[2].plot(time, Phi_EL_kg/Me, label='EL loss Zephyrus', color='pink', linestyle='-')
axs[2].set_xlabel('Time [yr]', fontsize=14)
axs[2].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=14)
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].legend()

# Adjust layout and save
plt.tight_layout()
plt.savefig('output/evolution_summary.png', dpi=300)
plt.close()
