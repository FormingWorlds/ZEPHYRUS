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
T = 280 # K
d = 1.0 * au2m
L = Luminosity(Rstar, Tstar)
Fp = Insolation(L, d)
F0 = Fp*10**(-3.5)*(Mstar/Ms) 


# Run IsoFate for a binary mixture
sol_isofate = isocalc(f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', species = 'H/O', rad_evol = True,
mu = mu_solar, eps = 0.15, activity = 'medium', flux_model = 'power law', stellar_type = 'G',
Rp_override = False, t_sat = 50e6, f_atm_final = 'null', n_TO_final = 'null', 
n_steps = int(1e5), t0 = 1e6, rho_rcb = 1.0, Johnson = False, RR = True, f_pred = False, thermal = True, 
beta = -1.23)

sol_zephyrus = isocalc_zephyrus(f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', species = 'H/O', rad_evol = True,
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
Phi1_a_isofate= sol_isofate['Phi1']        # light species number flux [particles/s/m2]
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
Phi1_a_zephyrus = sol_zephyrus['Phi1']        # light species number flux [particles/s/m2]
Phi2_a_zephyrus = sol_zephyrus['Phi2']        # heavy species number flux [particles/s/m2]

print('IsoFate:')
print(phic_a_isofate[0][0], phic_a_isofate[0][1])
print('ZEPHYRUS:')
print(phic_a_zephyrus[0][0], phic_a_zephyrus[0][1])

# # Plot

# Fig 1: Atmospheric mass loss
plt.figure(figsize=(10, 6))
plt.plot(t_a_isofate*s2yr, phi_a_isofate[0],  linestyle='-',  label='Total flux',       color='darkviolet')
plt.plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux',    color='red')
plt.plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-',  label='H flux',           color='lightblue')
plt.plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_O, linestyle='-',  label='O flux',           color='red')
cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
if len(cross_isofate[0]) != 0:
    plt.axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
plt.tight_layout()
plt.legend()
plt.show()

# Fig 1: Atmospheric mass loss
plt.figure(figsize=(10, 6))
plt.plot(t_a_isofate*s2yr, phi_a_isofate[0],  linestyle='-',  label='Total flux',       color='darkviolet')
plt.plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux',    color='red')
plt.plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-',  label='H flux',           color='lightblue')
plt.plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_O, linestyle='-',  label='O flux',           color='red')
cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
if len(cross_isofate[0]) != 0:
    plt.axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
plt.plot(t_a_zephyrus[1:]*s2yr, phi_a_zephyrus[0][1:], linestyle=':', label='Total flux ZEPHYRUS', color='darkviolet')
plt.plot(t_a_zephyrus[1:]*s2yr, phic_a_zephyrus[0][1:], linestyle=':', label='Critical flux ZEPHYRUS', color='red')
plt.plot(t_a_zephyrus[1:]*s2yr, Phi1_a_zephyrus[0][1:]*mu_H, linestyle=':', label='H flux ZEPHYRUS', color='lightblue') 
plt.plot(t_a_zephyrus[1:]*s2yr, Phi2_a_zephyrus[0][1:]*mu_O, linestyle=':', label='O flux ZEPHYRUS', color='red')
cross_zephyrus = np.where(np.abs(phi_a_zephyrus[0] - phic_a_zephyrus[0]) < 1e-15)
if len(cross_zephyrus[0]) != 0:
    plt.axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Time [yr]', fontsize=16)
plt.ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=16)
plt.tight_layout()
plt.legend()
plt.show()

# # Fig 2: Moles of H and O
# plt.figure(figsize=(10, 6))
# plt.plot(t_a*s2yr, N1_a[0]/avogadro, linestyle='-', label=r'$N_H$', color='lightblue')
# plt.plot(t_a*s2yr, N2_a[0]/avogadro, linestyle='-', label=r'$N_O$', color='red')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel('Number of particles [moles]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.show()

# # Fig 3: Radius
# plt.figure(figsize=(10, 6))
# plt.plot(t_a*s2yr, rp_a[0]/Re, linestyle='-', color='green')
# plt.xscale('log')   
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'Radius [R$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.show()

# # Fig 4: M_atmosphere
# plt.figure(figsize=(10, 6))
# plt.plot(t_a*s2yr, menv_a[0]/Me, linestyle='-', color='darkblue')
# plt.xscale('log')     
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'M$_{ATM}$ [M$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.show()

# # Fig 3: Mass Loss
# plt.figure(figsize=(10, 6))
# plt.plot(t_a*s2yr, mloss_a[0]/Me, linestyle='-', label='Mass loss', color='crimson')
# plt.xscale('log')   
# plt.yscale('log')
# plt.xlabel('Time [yr]', fontsize=16)
# plt.ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=16)
# plt.tight_layout()
# plt.legend()
# plt.show()

# # Fig 4: Atmospheric mass fraction 
# plt.figure(figsize=(10, 6))
# plt.plot(t_a*s2yr, fenv_a[0]*100, linestyle='-', color='gold')
# plt.xscale('log')   
# plt.xlabel('Time [yr]', fontsize=16)    
# plt.ylabel(r'f$_{ATM}$ [%]', fontsize=16)
# plt.tight_layout()
# plt.show()

fig, axs = plt.subplots(3, 2, figsize=(10, 8))
axs = axs.flatten()  # Flatten for easier indexing

# Panel 1: Flux
axs[0].plot(t_a_isofate*s2yr, phi_a_isofate[0], linestyle='-', label='Total flux', color='darkviolet')
axs[0].plot(t_a_isofate*s2yr, phic_a_isofate[0], linestyle='--', label='Critical flux', color='red')
axs[0].plot(t_a_isofate*s2yr, Phi1_a_isofate[0]*mu_H, linestyle='-', label='H flux', color='lightblue')
axs[0].plot(t_a_isofate*s2yr, Phi2_a_isofate[0]*mu_O, linestyle='-', label='O flux', color='red')
cross_isofate = np.where(np.abs(phi_a_isofate[0] - phic_a_isofate[0]) < 1e-15)
if len(cross_isofate[0]) != 0:
    axs[0].axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')

axs[0].plot(t_a_zephyrus[1:]*s2yr, phi_a_zephyrus[0][1:], linestyle=':', label='Total flux ZEPHYRUS', color='darkviolet')
axs[0].plot(t_a_zephyrus[1:]*s2yr, phic_a_zephyrus[0][1:], linestyle=':', label='Critical flux ZEPHYRUS', color='red')
axs[0].plot(t_a_zephyrus[1:]*s2yr, Phi1_a_zephyrus[0][1:]*mu_H, linestyle=':', label='H flux ZEPHYRUS', color='lightblue')      
axs[0].plot(t_a_zephyrus[1:]*s2yr, Phi2_a_zephyrus[0][1:]*mu_O, linestyle=':', label='O flux ZEPHYRUS', color='red')
cross_zephyrus = np.where(np.abs(phi_a_zephyrus[0] - phic_a_zephyrus[0]) < 1e-15)
if len(cross_zephyrus[0]) != 0:
    axs[0].axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')

axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlabel('Time [yr]', fontsize=12)
axs[0].set_ylabel(r'Mass flux $\Phi$ [kg m$^{-2}$ s$^{-1}$]', fontsize=12)
axs[0].legend()

# Panel 2: Moles of H and O
axs[1].plot(t_a_isofate*s2yr, N1_a_isofate[0]/avogadro, linestyle='-', label=r'$N_H$', color='lightblue')
axs[1].plot(t_a_isofate*s2yr, N2_a_isofate[0]/avogadro, linestyle='-', label=r'$N_O$', color='red')
if len(cross_isofate[0]) != 0:
    axs[1].axvline(t_a_isofate[cross_isofate][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')
if len(cross_zephyrus[0]) != 0:
    axs[1].axvline(t_a_zephyrus[cross_zephyrus][0]*s2yr, ls = '--', lw = 0.5, color = 'coral')

axs[1].plot(t_a_zephyrus*s2yr, N1_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_H$ ZEPHYRUS', color='lightblue')
axs[1].plot(t_a_zephyrus*s2yr, N2_a_zephyrus[0]/avogadro, linestyle=':', label=r'$N_O$ ZEPHYRUS', color='red')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xlabel('Time [yr]', fontsize=12)
axs[1].set_ylabel('Number of particles [moles]', fontsize=12)
axs[1].legend()

# Panel 3: Radius
axs[2].plot(t_a_isofate*s2yr, rp_a_isofate[0]/Re, linestyle='-', color='green')
axs[2].plot(t_a_zephyrus*s2yr, rp_a_zephyrus[0]/Re, linestyle=':', color='green')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].set_xlabel('Time [yr]', fontsize=12)
axs[2].set_ylabel(r'Radius [R$_\oplus$]', fontsize=12)

# Panel 4: M_atmosphere
axs[3].plot(t_a_isofate*s2yr, menv_a_isofate[0]/Me, linestyle='-', color='darkblue')
axs[3].plot(t_a_zephyrus*s2yr, menv_a_zephyrus[0]/Me, linestyle=':', color='darkblue')
axs[3].set_xscale('log')
axs[3].set_yscale('log')
axs[3].set_xlabel('Time [yr]', fontsize=12)
axs[3].set_ylabel(r'M$_{ATM}$ [M$_\oplus$]', fontsize=12)

# Panel 5: Mass Loss
axs[4].plot(t_a_isofate*s2yr, mloss_a_isofate[0]/Me, linestyle='-', label='Mass loss', color='crimson')
axs[4].plot(t_a_zephyrus*s2yr, mloss_a_zephyrus[0]/Me, linestyle=':', label='Mass loss ZEPHYRUS', color='crimson')
axs[4].set_xscale('log')
axs[4].set_yscale('log')
axs[4].set_xlabel('Time [yr]', fontsize=12)
axs[4].set_ylabel(r'$\Delta$ mass [M$_\oplus$]', fontsize=12)
axs[4].legend()

# Panel 6: Atmospheric Mass Fraction
axs[5].plot(t_a_isofate*s2yr, fenv_a_isofate[0]*100, linestyle='-', color='gold')
axs[5].plot(t_a_zephyrus*s2yr, fenv_a_zephyrus[0]*100, linestyle=':', color='gold')
axs[5].set_xscale('log')
axs[5].set_xlabel('Time [yr]', fontsize=12)
axs[5].set_ylabel(r'f$_{ATM}$ [%]', fontsize=12)

plt.tight_layout()
plt.show()