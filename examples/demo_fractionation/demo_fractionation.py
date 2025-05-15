import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from zephyrus.constants import kb, G, M_H, M_O
from zephyrus.planets_parameters import Me, Re
from zephyrus.fractionation import Fractionation_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O

epsilon = 0.1  
Mp = Me
Rp = Re

# Load and extract CSV data from the case 834 from the 7 parameter grid on Habrok
df = pd.read_csv('/Users/emmapostolec/Downloads/runtime_helpfile.csv', sep='\t')
time   = df['Time'].values          # [yr]
Teq   = df['T_eqm'].values           # [K]
mass_H = df['H_kg_atm'].values      # [kg]
mass_O = df['O_kg_atm'].values      # [kg]
Fxuv   = df['F_xuv'].values          # [W/m2]
Matm   = df['M_atm'].values          # [kg]

N_light = np.zeros(len(time))   
N_heavy = np.zeros(len(time))
N_total = np.zeros(len(time))
Phi_light = np.zeros(len(time))
Phi_heavy = np.zeros(len(time)) 
Phi_total = np.zeros(len(time))
Phi_diffusion_light = np.zeros(len(time))
Phi_diffusion_heavy = np.zeros(len(time))
Phi_diffusion_critical = np.zeros(len(time))
Phi_EUV_Collin = np.zeros(len(time))

for i,t in enumerate(time):
    Phi_EUV_Collin_i = (epsilon * Fxuv[i]) / (4 * (G * Mp / Rp**2))  # [particles/m2/s]
    b_HO_i = b_H_O(Teq[i])  # Binary diffusion coefficient for H and O
    n_H_i = mass_H[i] / (M_H)  # Number of H particles in the atmosphere [mol]
    n_O_i = mass_O[i] / (M_O)  # Number of O particles in the atmosphere [mol]
    Phi_dl, Phi_df, Phi_c, Phi_l, Phi_h, Phi_tot = Fractionation_binary_mixture(n_light=n_H_i, n_heavy=n_O_i, Mp=Mp, Rp=Rp, b=b_HO_i, T=Teq[i], M_light=M_H, M_heavy=M_O, Phi=Phi_EUV_Collin_i)
    
    N_light[i] = n_H_i  
    N_heavy[i] = n_O_i
    N_total[i] = n_H_i + n_O_i
    Phi_diffusion_light[i] = Phi_dl
    Phi_diffusion_heavy[i] = Phi_df
    Phi_diffusion_critical[i] = Phi_c
    Phi_light[i] = Phi_l
    Phi_heavy[i] = Phi_h
    Phi_total[i] = Phi_tot
    Phi_EUV_Collin[i] = Phi_EUV_Collin_i

# Plot the results

# Evolution of the number of moles
plt.figure(figsize=(10, 6))
plt.plot(time, N_light, label='H', color='dodgerblue')
plt.plot(time, N_heavy, label='O', color='tomato')   
plt.plot(time, N_total, label='Total', color='purple')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel('Number of particles [moles]', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()    
plt.savefig('output/evolution_mol.png', dpi=300)
plt.close()

# Evolution of the mass flux
plt.figure(figsize=(10, 6))
#plt.plot(time, Phi_EUV_Collin, label='EUV', color='black', linestyle='-')
plt.plot(time, Phi_light, label='H', color='lightskyblue', linestyle='-')
plt.plot(time, Phi_heavy, label='O', color='salmon', linestyle='-')
plt.plot(time, Phi_total, label='Total', color='purple', linestyle='-')
plt.plot(time, Phi_diffusion_light, label=r'$H_{diff}$', color='lightskyblue', linestyle=':')
plt.plot(time, Phi_diffusion_heavy, label=r'$O_{diff}$', color='salmon', linestyle=':')
plt.plot(time, Phi_diffusion_critical, label='Critical', color='red', linestyle='--')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'Number flux $\Phi$ [particles m$^{-2}$ s$^{-1}$]', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()    
plt.savefig('output/evolution_number_flux_all.png', dpi=300)
plt.close()

plt.figure(figsize=(10, 6))
plt.plot(time, Phi_EUV_Collin, label='EUV', color='black', linestyle='-')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
plt.yscale('log')
plt.xscale('log')
plt.grid(alpha=0.2)
plt.legend()    
plt.savefig('output/evolution_mass_flux.png', dpi=300)
plt.close()