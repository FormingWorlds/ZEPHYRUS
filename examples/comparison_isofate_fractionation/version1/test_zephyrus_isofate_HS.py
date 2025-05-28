from zephyrus_binary import isocalc_zephyrus
from isofate_binary import isocalc
from isofunks import *
from constants import *
import numpy as np
import matplotlib.pyplot as plt

# Initialazition parameters
f_atm = 0.0001
Mp = 1.0 * Me
Mstar = 1.0 * Ms
Rstar = 1.0 * Rs
Tstar = 5600 # [K]
T = 1000 # K
d = 1.0 * au2m
L = Luminosity(Rstar, Tstar)
Fp = Insolation(L, d)
F0 = Fp*10**(-3.5)*(Mstar/Ms) 


# Run IsoFate for a binary mixture
sol_isofate = isocalc(f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', species = 'H/S', rad_evol = True,
mu = mu_solar, eps = 0.15, activity = 'medium', flux_model = 'power law', stellar_type = 'G',
Rp_override = False, t_sat = 50e6, f_atm_final = 'null', n_TO_final = 'null', 
n_steps = int(1e5), t0 = 1e6, rho_rcb = 1.0, Johnson = False, RR = True, f_pred = False, thermal = True, 
beta = -1.23)

sol_zephyrus = isocalc_zephyrus(f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', species = 'H/S', rad_evol = True,
mu = mu_solar, eps = 0.15, activity = 'medium', flux_model = 'power law', stellar_type = 'G',
Rp_override = False, t_sat = 50e6, f_atm_final = 'null', n_TO_final = 'null', 
n_steps = int(1e5), t0 = 1e6, rho_rcb = 1.0, Johnson = False, RR = True, f_pred = False, thermal = True, 
beta = -1.23)


# Define solutions for plotting
t_a_isofate = sol_isofate['time']
rp_a_isofate = sol_isofate['rp']            # total planet radius [m]
menv_a_isofate = sol_isofate['menv']        # atmospheric mass [kg]
fenv_a_isofate = sol_isofate['fenv']        # atmospheric mass fraction [ndim]
mloss_a_isofate = sol_isofate['mloss']      # mass loss per timestep [kg]
phi_a_isofate = sol_isofate['phi']          # mass flux [kg/2/m2]
phic_a_isofate = sol_isofate['phic']        # critical mass flux [kg/s/m2]
N1_a_isofate = sol_isofate['N1']            # light species number [particles]
N2_a_isofate = sol_isofate['N2']            # heavy species number [particles]
x1_a_isofate = sol_isofate['x1']            # light species molar concentration [ndim]
x2_a_isofate = sol_isofate['x2']            # heavy species molar concentration [ndim]
dh_a_isofate = sol_isofate['X2']            # N2/N1 mole ratio [ndim]
dhfinal_a_isofate = sol_isofate['X2_final'] # final N2/N1 mole ratio [ndim]
Phi1_a_isofate= sol_isofate['Phi1']         # light species number flux [particles/s/m2]
Phi2_a_isofate = sol_isofate['Phi2']        # heavy species number flux [particles/s/m2]  

t_a_zephyrus = sol_zephyrus['time']
rp_a_zephyrus = sol_zephyrus['rp']            # total planet radius [m]
menv_a_zephyrus = sol_zephyrus['menv']        # atmospheric mass [kg]
fenv_a_zephyrus = sol_zephyrus['fenv']        # atmospheric mass fraction [ndim]
mloss_a_zephyrus = sol_zephyrus['mloss']      # mass loss per timestep [kg]
phi_a_zephyrus = sol_zephyrus['phi']          # mass flux [kg/2/m2]
phic_a_zephyrus = sol_zephyrus['phic']        # critical mass flux [kg/s/m2]
N1_a_zephyrus = sol_zephyrus['N1']            # light species number [particles]
N2_a_zephyrus = sol_zephyrus['N2']            # heavy species number [particles]
x1_a_zephyrus = sol_zephyrus['x1']            # light species molar concentration [ndim]
x2_a_zephyrus = sol_zephyrus['x2']            # heavy species molar concentration [ndim]
dh_a_zephyrus = sol_zephyrus['X2']            # N2/N1 mole ratio [ndim]
dhfinal_a_zephyrus = sol_zephyrus['X2_final'] # final N2/N1 mole ratio [ndim]
Phi1_a_zephyrus = sol_zephyrus['Phi1']        # light species number flux [particles/s/m2]
Phi2_a_zephyrus = sol_zephyrus['Phi2']        # heavy species number flux [particles/s/m2]

mlossl_a_isofate = sol_isofate['mlossl']      # mass loss per timestep [kg]
mlossh_a_isofate = sol_isofate['mlossh']      # mass loss per timestep [kg]
mlossl_a_zephyrus = sol_zephyrus['mlossl']      # mass loss per timestep [kg]
mlossh_a_zephyrus = sol_zephyrus['mlossh']      # mass loss per timestep [kg]

# # Plot
# # Fig 1: Atmospheric mass loss - IsoFATE only
# plt.figure(figsize=(14, 8))
# plt.plot(t_a_isofate*s2yr, phi_a_isofate[0],  linestyle='-',  label='Total flux',       color='darkviolet')
# plt.plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux',    color='orange')
# plt.plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-',  label='H flux',           color='dodgerblue')
# plt.plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_O, linestyle='-',  label='O flux',           color='red')
# cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
# if len(cross_isofate[0]) != 0:
#     plt.axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral', label='End fractionation ?')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/mass_flux_isofate.png', dpi=300)
# plt.close()

# # Fig 1bis: Atmospheric mass loss - IsoFATE vs Zephyrus
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, phi_a_isofate[0],  linestyle='-',  label='IsoFATE',       color='black')
# plt.plot(t_a_zephyrus[1:]*s2yr, phic_a_zephyrus[0][1:], linestyle=':', label='Zephyrus', color='black', linewidth=2)

# plt.plot(t_a_isofate*s2yr, phi_a_isofate[0],  linestyle='-',  label='Total flux',       color='darkviolet')
# plt.plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux',    color='orange')
# plt.plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-',  label='H flux',           color='dodgerblue')
# plt.plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_O, linestyle='-',  label='O flux',           color='red')
# cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
# if len(cross_isofate[0]) != 0:
#     plt.axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
# plt.plot(t_a_zephyrus*s2yr, phi_a_zephyrus[0], linestyle=':', color='darkviolet', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, phic_a_zephyrus[0], linestyle=':', color='orange', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, Phi1_a_zephyrus[0]*mu_H, linestyle=':', color='dodgerblue', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, Phi2_a_zephyrus[0]*mu_O, linestyle=':',color='red', linewidth=2)
# cross_zephyrus = np.where(np.abs(phi_a_zephyrus[0] - phic_a_zephyrus[0]) < 1e-15)
# if len(cross_zephyrus[0]) != 0:
#     plt.axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/mass_flux_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 2: Moles of H and O
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, N1_a_isofate[0]/avogadro, linestyle='-', label=r'$N_H$ IsoFATE', color='dodgerblue')
# plt.plot(t_a_isofate*s2yr, N2_a_isofate[0]/avogadro, linestyle='-', label=r'$N_O$ IsoFATE', color='red')
# plt.plot(t_a_zephyrus*s2yr, N1_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_H$ Zephyrus', color='dodgerblue', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, N2_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_O$ Zephyrus', color='red', linewidth=2)
# cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
# if len(cross_isofate[0]) != 0:
#     plt.axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
# cross_zephyrus = np.where(np.abs(phi_a_zephyrus[0] - phic_a_zephyrus[0]) < 1e-15)
# if len(cross_zephyrus[0]) != 0:
#     plt.axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral', label='End fractionation ?')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel('Number of particles [moles]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/moles_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 3: Radius
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, rp_a_isofate[0]/Re, linestyle='-', label=r'IsoFATE', color='green')
# plt.plot(t_a_zephyrus*s2yr, rp_a_zephyrus[0]/Re, linestyle=':', label=r'Zephyrus', color='green', linewidth=2)
# plt.xscale('log')   
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'Radius [R$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/radius_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 4: M_atmosphere
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, menv_a_isofate[0]/Me, linestyle='-', label=r'IsoFATE', color='darkblue')
# plt.plot(t_a_zephyrus*s2yr, menv_a_zephyrus[0]/Me, linestyle=':', label=r'Zephyrus', color='darkblue', linewidth=2)
# plt.xscale('log')     
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'M$_{ATM}$ [M$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/mass_atmosphere_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 5: Mass Loss
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, mloss_a_isofate[0]/Me, linestyle='-', label='IsoFATE', color='crimson')
# plt.plot(t_a_zephyrus*s2yr, mloss_a_zephyrus[0]/Me, linestyle=':', label='Zephyrus', color='crimson', linewidth=2)
# plt.xscale('log')   
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/mass_loss_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 6: Atmospheric mass fraction 
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, fenv_a_isofate[0]*100, linestyle='-', label=r'IsoFATE', color='orange')
# plt.plot(t_a_zephyrus*s2yr, fenv_a_zephyrus[0]*100, linestyle=':', label=r'Zephyrus', color='gold', linewidth=2)
# plt.xscale('log')   
# plt.xlabel('Time [yr]', fontsize=16)    
# plt.ylabel(r'f$_{ATM}$ [%]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/fenv_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 8 :  species/H
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, x1_a_isofate[0]*100, linestyle='-', label=r'$x_H$ IsoFATE', color='dodgerblue')
# plt.plot(t_a_isofate*s2yr, x2_a_isofate[0]*100, linestyle='-', label=r'$x_O$ IsoFATE', color='red')
# plt.plot(t_a_zephyrus*s2yr, x1_a_zephyrus[0]*100, linestyle=':', label=r'$x_H$ Zephyrus', color='dodgerblue', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, x2_a_zephyrus[0]*100, linestyle=':', label=r'$x_O$ Zephyrus', color='red', linewidth=2)
# plt.xscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'$\chi_{i}^{ATM}$ [%]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/species_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Fig 9 : mole ratio
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, (N2_a_isofate[0]/N1_a_isofate[0])/OtoH_protosolar, linestyle='-', label=r'IsoFATE', color='maroon')
# plt.plot(t_a_zephyrus*s2yr, (N2_a_zephyrus[0]/N1_a_zephyrus[0])/OtoH_protosolar, linestyle=':', label=r'Zephyrus', color='maroon', linewidth=2)
# plt.xscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'O/H [Solar]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.savefig('plots_HO/mole_ratio_comparison_zephyrus.png', dpi=300)
# plt.close()

# # Figure 10: Mass loss per timestep
# plt.figure(figsize=(10, 6))
# plt.plot(t_a_isofate*s2yr, mloss_a_isofate[0]/Me, linestyle='-', label='IsoFATE', color='green')
# plt.plot(t_a_zephyrus*s2yr, mloss_a_zephyrus[0]/Me, linestyle=':', label='Zephyrus', color='green', linewidth=2)
# plt.plot(t_a_isofate*s2yr, mlossl_a_isofate[0]/Me, linestyle='-', label='H loss IsoFATE', color='dodgerblue')
# plt.plot(t_a_isofate*s2yr, mlossh_a_isofate[0]/Me, linestyle='-', label='O loss IsoFATE', color='red')
# plt.plot(t_a_zephyrus*s2yr, mlossl_a_zephyrus[0]/Me, linestyle=':', label='H loss Zephyrus', color='dodgerblue', linewidth=2)
# plt.plot(t_a_zephyrus*s2yr, mlossh_a_zephyrus[0]/Me, linestyle=':', label='O loss Zephyrus', color='red', linewidth=2)       
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.legend()    
# plt.savefig('plots_HO/mass_loss_per_timestep_comparison_zephyrus.png', dpi=300)
# plt.close()


# Fig 7: All panels in one figure
fig, axs = plt.subplots(4, 2, figsize=(10, 10))
axs = axs.flatten()  # Flatten for easier indexing

# Panel 1: Flux
axs[0].plot(t_a_isofate*s2yr, phi_a_isofate[0], linestyle='-', label='IsoFATE', color='black')
axs[0].plot(t_a_zephyrus*s2yr, phic_a_zephyrus[0], linestyle=':', label='Zephyrus', color='black', linewidth=2)
axs[0].plot(t_a_isofate*s2yr, phi_a_isofate[0], linestyle='-', label='Total flux', color='darkviolet')
axs[0].plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux', color='orange')
axs[0].plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-', label='H flux', color='dodgerblue')
axs[0].plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_S, linestyle='-', label='S flux', color='peru')
cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
if len(cross_isofate[0]) != 0:
    axs[0].axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
axs[0].plot(t_a_zephyrus*s2yr, phi_a_zephyrus[0], linestyle=':', color='darkviolet', linewidth=2)
axs[0].plot(t_a_zephyrus*s2yr, phic_a_zephyrus[0], linestyle=':', color='orange', linewidth=2)
axs[0].plot(t_a_zephyrus*s2yr, Phi1_a_zephyrus[0]*mu_H, linestyle=':', color='dodgerblue', linewidth=2)      
axs[0].plot(t_a_zephyrus*s2yr, Phi2_a_zephyrus[0]*mu_S, linestyle=':',color='peru', linewidth=2)
cross_zephyrus = np.where(np.abs(phi_a_zephyrus[0] - phic_a_zephyrus[0]) < 1e-15)
if len(cross_zephyrus[0]) != 0:
    axs[0].axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_ylim(1e-13, 1e-9)
axs[0].set_xlabel('Time [yr]', fontsize=14)
axs[0].set_ylabel(r'$\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=14)
axs[0].legend()

# Panel 2: Moles of H and O
axs[1].plot(t_a_isofate*s2yr, N1_a_isofate[0]/avogadro, linestyle='-', label=r'$N_H$ IsoFATE', color='dodgerblue')
axs[1].plot(t_a_isofate*s2yr, N2_a_isofate[0]/avogadro, linestyle='-', label=r'$N_S$ IsoFATE', color='peru')
axs[1].plot(t_a_zephyrus*s2yr, N1_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_H$ Zephyrus', color='dodgerblue', linewidth=2)
axs[1].plot(t_a_zephyrus*s2yr, N2_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_S$ Zephyrus', color='peru', linewidth=2)
if len(cross_isofate[0]) != 0:
    axs[1].axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')#, label='End fractionation ?')
if len(cross_zephyrus[0]) != 0:
    axs[1].axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xlabel('Time [yr]', fontsize=14)
axs[1].set_ylabel(r'$N_i$ [Moles]', fontsize=14)
axs[1].legend()

# Panel 3: Radius
axs[4].plot(t_a_isofate*s2yr, rp_a_isofate[0]/Re, label=r'IsoFATE', linestyle='-', color='green')
axs[4].plot(t_a_zephyrus*s2yr, rp_a_zephyrus[0]/Re, label=r'Zephyrus', linestyle=':', color='green', linewidth=2)
axs[4].set_xscale('log')
axs[4].set_yscale('log')
axs[4].set_xlabel('Time [yr]', fontsize=14)
axs[4].set_ylabel(r'Radius [R$_\oplus$]', fontsize=14)
axs[4].legend()

# Panel 3: Species mole fraction
axs[3].plot(t_a_isofate*s2yr, x1_a_isofate[0]*100, linestyle='-', label=r'$x_H$ IsoFATE', color='dodgerblue')
axs[3].plot(t_a_isofate*s2yr, x2_a_isofate[0]*100, linestyle='-', label=r'$x_S$ IsoFATE', color='peru')          
axs[3].plot(t_a_zephyrus*s2yr, x1_a_zephyrus[0]*100, linestyle=':', label=r'$x_H$ Zephyrus', color='dodgerblue', linewidth=2)
axs[3].plot(t_a_zephyrus*s2yr, x2_a_zephyrus[0]*100, linestyle=':', label=r'$x_S$ Zephyrus', color='peru', linewidth=2)
axs[3].set_xscale('log')
# axs[3].set_yscale('log')
axs[3].set_xlabel('Time [yr]', fontsize=14)
axs[3].set_ylabel(r'$\chi_{i}^{ATM}$ [%]', fontsize=14)
axs[3].legend()

# Panel 4: Mole ratio vs solar
#axs[2].plot(t_a_isofate*s2yr, dh_a_isofate[0]/dh_a_isofate[0,0], linestyle='-', label=r'IsoFATE', color='maroon')       
#axs[2].plot(t_a_zephyrus*s2yr, dh_a_zephyrus[0]/dh_a_zephyrus[0,0], linestyle=':', label=r'Zephyrus', color='maroon', linewidth=2)
axs[2].plot(t_a_isofate*s2yr, (N2_a_isofate[0]/N1_a_isofate[0])/StoH_protosolar*100, linestyle='-', label=r'IsoFATE', color='maroon')
axs[2].plot(t_a_zephyrus*s2yr, (N2_a_zephyrus[0]/N1_a_zephyrus[0])/StoH_protosolar*100, linestyle=':', label=r'Zephyrus', color='maroon', linewidth=2)
axs[2].set_xscale('log')
# axs[2].set_yscale('log')
axs[2].set_xlabel('Time [yr]', fontsize=14)
axs[2].set_ylabel(r'S/H [Solar]', fontsize=14)
axs[2].legend()

# Panel 6: Atmospheric Mass Fraction
axs[5].plot(t_a_isofate*s2yr, fenv_a_isofate[0]*100, label=r'IsoFATE', linestyle='-', color='gold')
axs[5].plot(t_a_zephyrus*s2yr, fenv_a_zephyrus[0]*100, label=r'Zephyrus', linestyle=':', color='gold', linewidth=2)
axs[5].set_xscale('log')
axs[5].set_xlabel('Time [yr]', fontsize=14)
axs[5].set_ylabel(r'f$_{ATM}$ [%]', fontsize=14)
axs[5].legend()

# Panel 7: M_atmosphere
axs[6].plot(t_a_isofate*s2yr, menv_a_isofate[0]/Me, label=r'IsoFATE', linestyle='-', color='darkblue')
axs[6].plot(t_a_zephyrus*s2yr, menv_a_zephyrus[0]/Me, label=r'Zephyrus', linestyle=':', color='darkblue', linewidth=2)
axs[6].set_xscale('log')
axs[6].set_yscale('log')
axs[6].set_xlabel('Time [yr]', fontsize=14)
axs[6].set_ylabel(r'M$_{ATM}$ [M$_\oplus$]', fontsize=14)
axs[6].legend()

# Panel 8: Mass Loss
axs[7].plot(t_a_isofate*s2yr, mloss_a_isofate[0]/Me, linestyle='-', label='IsoFate', color='green')
axs[7].plot(t_a_zephyrus*s2yr, mloss_a_zephyrus[0]/Me, linestyle=':', label='Zephyrus', color='green', linewidth=2)
axs[7].plot(t_a_isofate*s2yr, mlossl_a_isofate[0]/Me, linestyle='-', label='H loss IsoFATE', color='dodgerblue')
axs[7].plot(t_a_isofate*s2yr, mlossh_a_isofate[0]/Me, linestyle='-', label='S loss IsoFATE', color='peru')
axs[7].plot(t_a_zephyrus*s2yr, mlossl_a_zephyrus[0]/Me, linestyle=':', label='H loss Zephyrus', color='dodgerblue', linewidth=2)
axs[7].plot(t_a_zephyrus*s2yr, mlossh_a_zephyrus[0]/Me, linestyle=':', label='S loss Zephyrus', color='peru', linewidth=2)
axs[7].set_xscale('log')
axs[7].set_yscale('log')
axs[7].set_ylim(1e-9, 1e-7)
axs[7].set_xlabel('Time [yr]', fontsize=14)
axs[7].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=14)
axs[7].legend()


# Box input parameters
# Add a box with input parameters
param_text = (
    f"$f_{{atm}}$ = {f_atm*100} % \n"
    f"$M_p$ = {Mp/Me:.2f} $M_\\oplus$\n"
    f"$T_{{p}}$ = {T} K\n"
    f"$d$ = {d/au2m:.2f} AU\n"
    f"$M_{{star}}$ = {Mstar/Ms:.2f} $M_\\odot$\n"
    f"$R_{{star}}$ = {Rstar/Rs:.2f} $R_\\odot$\n"
    f"$T_{{star}}$ = {Tstar} K")
fig.text(0.61, 0.302, param_text, fontsize=10, va='bottom', ha='left',
         bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="black", alpha=0.8))

fig.suptitle('Comparison fractionation IsoFATE vs Zephyrus : Binary mixture H-S', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig('plots_HS/8panels_comparison_zephyrus_HS.png', dpi=300)
plt.close()

