import numpy as np
import matplotlib.pyplot as plt
import mors

from zephyrus.constants import *

########################### Initialization #####################################

Omega_sun        = 1.0                                             # Solar rotation rate        [Omega_sun, rad s-1]
Star_masses      = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]                  # Star mass in solar mass    [M_sun,     kg]

colors           = ['brown', 'red', 'orange', 'gold', 'yellow', 'palegoldenrod']
labels           = [f'm={mass} M$_{{\odot}}$' for mass in Star_masses]

########################### Extract XUV luminosity from MORS data #####################################

star_data = []
lxuv_data = []

for mass in Star_masses:   
    print(r'Mstar = ', mass, 'Msun')
    star = mors.Star(Mstar=mass, Omega=Omega_sun)  # extract luminosities using the mors.Star() function
    star_data.append({
        'mass': mass,
        'Age': star.Tracks['Age'],
        'Lbol': star.Tracks['Lbol'],
        'Lx': star.Tracks['Lx'],
        'Leuv': star.Tracks['Leuv'],
        'Lxuv': star.Tracks['Lx'] + star.Tracks['Leuv']
    })

    ages = star.Tracks['Age']
    lxuv_for_mass = [mors.Lxuv(Mstar=mass, Age=age, Omega=Omega_sun)['Lxuv'] for age in ages] # extract luminosities using the mors.Lxuv() function
    lxuv_data.append(lxuv_for_mass)

######### Plot #####################################

# Bolometric luminosity vs time
plt.figure(figsize=(10, 8))
for data, color, label in zip(star_data, colors, labels):
    plt.loglog(data['Age'], data['Lbol'], color=color, label=f'{label}')
plt.ylabel(r'Bolometric luminosity L$_{\mathrm{bol}}$ [erg $s^{-1}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title(f'Stellar Evolution Tracks from MORS (Omega = {Omega_sun} $\Omega_{{\odot}}$)', fontsize=15)
plt.axvline(x=4603, color='dimgrey', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e30, 'Today', color='dimgrey', rotation=90, verticalalignment='bottom')
plt.savefig('MORS_Stellar_evolution_Lbol_vs_time_masses_using_Star.pdf')

# XUV luminosity vs time
plt.figure(figsize=(10, 8))
plt.loglog(star_data[0]['Age'], star_data[0]['Lxuv'], color='black', label='Extracted using Star()')
plt.loglog(star_data[0]['Age'], lxuv_data[0], ':', color='black', label='Extracted using Lxuv()')
for data, lxuv_for_mass, color, label in zip(star_data, lxuv_data, colors, labels):
    plt.loglog(data['Age'], data['Lxuv'], color=color, label=f'{label}')
    plt.loglog(data['Age'], lxuv_for_mass, ':', color=color)
plt.ylabel(r'XUV luminosity L$_{\mathrm{XUV}}$ [erg $s^{-1}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title(f'Stellar Evolution Tracks from MORS (Omega = {Omega_sun} $\Omega_{{\odot}}$)', fontsize=15)
plt.axvline(x=4603, color='dimgrey', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e30, 'Today', color='dimgrey', rotation=90, verticalalignment='bottom')
plt.savefig('MORS_Stellar_evolution_Lxuv_vs_time_masses_using_Star_and_Lxuv.pdf')

