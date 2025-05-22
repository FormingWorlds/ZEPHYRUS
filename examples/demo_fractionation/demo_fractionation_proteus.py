import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from zephyrus.constants import kb, G, M_H, M_O, s2yr
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture, Planetary_surface_area
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
time   = df['Time'].values          # [yr]
Teq   = df['T_eqm'].values *1000          # [K]
mass_H = df['H_kg_atm'].values     # [kg]
mass_O = df['O_kg_atm'].values      # [kg]
Fxuv   = df['F_xuv'].values          # [W/m2]
Matm   = df['M_atm'].values          # [kg]
Mp   = df['M_planet'].values          # [kg]
Rp   = df['R_obs'].values          # [m]
escape_EL = df['esc_rate_total'].values      # [kg/s]

N_light = np.zeros(len(time))   
N_heavy = np.zeros(len(time))
x_light = np.zeros(len(time))
x_heavy = np.zeros(len(time))
N_total = np.zeros(len(time))
Phi_light = np.zeros(len(time))
Phi_heavy = np.zeros(len(time)) 
Phi_total = np.zeros(len(time))
Phi_critical = np.zeros(len(time))
Phi_EUV_Collin = np.zeros(len(time))
Phi_EL = np.zeros(len(time))
Phi_EL_kg = np.zeros(len(time))
mass_loss_H = np.zeros(len(time))
mass_loss_O = np.zeros(len(time))
mass_loss_total = np.zeros(len(time))
mass_loss_EUV = np.zeros(len(time))


for i,t in enumerate(time):
    Phi_EUV_Collin_i = (epsilon * Fxuv[i]) / (4 * (G * Mp[i] / Rp[i]**2))  # [particles/m2/s]
    Phi_EL_i = (epsilon * np.pi * Rp[i]**3 * Fxuv[i]) / (G * Mp[i])  # [kg/s]
    b_HO_i = b_H_O(Teq[i])  # Binary diffusion coefficient for H and O
    n_H_i = mass_H[i] / (M_H)  # Number of H particles in the atmosphere [mol]
    n_O_i = mass_O[i] / (M_O)  # Number of O particles in the atmosphere [mol]
    Phi_dl, Phi_df, Phi_c, Phi_l, Phi_h, Phi_tot = Fractionation_binary_mixture(n_light=n_H_i, n_heavy=n_O_i, Mp=Mp[i], Rp=Rp[i], b=b_HO_i, T=Teq[i], M_light=M_H, M_heavy=M_O, Phi=Phi_EUV_Collin_i, mu_light=mu_H, mu_heavy=mu_O)
    delta_t = t/len(time) # timestep [s]

    N_light[i] = n_H_i  
    N_heavy[i] = n_O_i
    x_light[i] = n_H_i / (n_H_i + n_O_i)  # Mole fraction of H
    x_heavy[i] = n_O_i / (n_H_i + n_O_i)  # Mole fraction of O
    N_total[i] = n_H_i + n_O_i
    Phi_critical[i] = Phi_c             # Critical flux [kg/m2/s]
    Phi_light[i] = Phi_l * mu_H         # H flux [kg/m2/s]
    Phi_heavy[i] = Phi_h * mu_O         # O flux [kg/m2/s]
    Phi_total[i] = Phi_tot * (x_light[i]*mu_H + x_heavy[i]*mu_O)    # Total flux [kg/m2/s]
    Phi_EUV_Collin[i] = Phi_EUV_Collin_i    # Mass flux  [kg/m2/s]
    Phi_EL[i] = Phi_EL_i # mass loss [kg/s]
    mass_loss_H[i] = Planetary_surface_area(Rp[i]) * Phi_light[i] * delta_t / s2yr  # Mass loss of H [kg]
    mass_loss_O[i] = Planetary_surface_area(Rp[i]) * Phi_heavy[i] * delta_t  / s2yr # Mass loss of O [kg]
    mass_loss_total[i] = Planetary_surface_area(Rp[i]) * Phi_total[i] * delta_t / s2yr # Total mass loss [kg]
    mass_loss_EUV[i] = Planetary_surface_area(Rp[i]) * Phi_EUV_Collin[i] * delta_t / s2yr # Mass loss due to EUV [kg]
    Phi_EL_kg[i] = Planetary_surface_area(Rp[i]) * Phi_EL[i] * delta_t / s2yr # Mass loss due to EL [kg]

# Plot the results

# Evolution of the number of moles
plt.figure(figsize=(10, 6))
plt.plot(time, N_light, label=r'$N_H$', color='dodgerblue')
plt.plot(time, N_heavy, label=r'$N_O$', color='red')   
plt.plot(time, N_total, label=r'$N_{Total}$', color='purple')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'$N_i$ [moles]', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()    
plt.savefig('output/evolution_mol.png', dpi=300)
plt.close()

# Evolution of the mass flux
plt.figure(figsize=(10, 6))
plt.plot(time, Phi_EUV_Collin, label=r'$\Phi_{EUV}$ (IsoFATE)', color='purple', linestyle='-')
plt.plot(time, Phi_light, label='H flux Zephyrus', color='dodgerblue', linestyle=':', linewidth=2)
plt.plot(time, Phi_heavy, label='O flux Zephyrus', color='red', linestyle=':', linewidth=2)
plt.plot(time, Phi_total, label='Total flux Zephyrus', color='purple', linestyle=':',linewidth=2)
plt.plot(time, Phi_critical, label='Critical', color='orange', linestyle='--')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()    
plt.savefig('output/evolution_mass_flux_all.png', dpi=300)
plt.close()

# Evolution of the mole fraction
plt.figure(figsize=(10, 6))
plt.plot(time, x_light*100, label=r'$x_H$', color='dodgerblue')
plt.plot(time, x_heavy*100, label=r'$x_O$', color='red')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'$\chi_i^{atm}$', fontsize=16)
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()
plt.savefig('output/evolution_mole_fraction.png', dpi=300)
plt.close()

# Evolution of the mass loss
plt.figure(figsize=(10, 6))
plt.plot(time, mass_loss_H/Me, label='H loss Zephyrus', color='dodgerblue')
plt.plot(time, mass_loss_O/Me, label='O loss Zephyrus', color='red')
plt.plot(time, mass_loss_total/Me, label=r'Total loss Zephyrus', color='purple')
#plt.plot(time, escape_EL, label='EL escape Zephyrus ', color='orange', linestyle='--')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=16)            
plt.xscale('log')
plt.yscale('log')
plt.grid(alpha=0.2)
plt.legend()
plt.savefig('output/evolution_mass_loss.png', dpi=300)
plt.close()


fig, axs = plt.subplots(2, 2, figsize=(10, 6))
axs = axs.flatten()
plt.title('Output PROTEUS simulation', fontsize=16)

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
