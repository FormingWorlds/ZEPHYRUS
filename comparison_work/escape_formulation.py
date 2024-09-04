import numpy as np
import matplotlib.pyplot as plt

import mors

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.escape import *

########################### Path to directories ###############################

path_plot           = '../plots/comparison_work/'


########################### Initialization #####################################

Omega_sun          = 1.0                                        # Solar rotation rate        [Omega_sun, rad s-1]

planet_masses      = [0.5,1.0,2.0,3.0,4.0,5.0]                # Planet masses               [Me, kg]
escape_today       = []                                       # Mass Loss Rate              [kg s-1]
escape_evolution   = []                                       # Mass Loss Rate              [kg s-1]
epsilon            = 0.15                                     # Escape efficiency factor    [dimensionless]

models             = ['Zephyrus', 'Attia+2021', 'IsoFATE']
colors             = ['steelblue', 'darkorange', 'green']


########################### Extract XUV luminosity from MORS data #####################################

Sun         = mors.Star(Mstar=1.0, Omega=Omega_sun)                         # Extract luminosities using the mors.Star() function
Sun_age     = Sun.Tracks['Age']
Sun_Lxuv    = Sun.Tracks['Lx'] + Sun.Tracks['Leuv']                         # XUV luminosity            [erg s-1]
Sun_Fxuv    = (Sun_Lxuv/(4 * np.pi * a_earth*au2cm **2)) * ergcm2stoWm2     # XUV flux                  [W m-2]


########################### Escape computations ####################################

for mass in planet_masses : 
    zephyrus    = EL_escape(False,a_earth*au2m,e_earth,mass*Me,Ms,epsilon,Re,Re,Fxuv_earth_today)
    attia_2021  = (dMdt_EL_Attia2021(False,a_earth*au2cm,e_earth,mass*Me*1e3,Ms*1e3,Re*1e2,(Fxuv_earth_today/ergcm2stoWm2)*4*np.pi*a_earth*au2cm**2,Fxuv_earth_today/ergcm2stoWm2,Rxuv=Re*1e2,epsilon=epsilon)/1e3)
    isofate     = (dMdt_Cherubim2024(mass*Me,Re,epsilon,Fxuv_earth_today)) * (4 * np.pi * Re ** 2) 

    escape_today.append({
        'Zephyrus': zephyrus, 
        'Attia+2021': attia_2021, 
        'IsoFATE': isofate 
    })

for fxuv in Sun_Fxuv : 
    zephyrus    = EL_escape(False,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv)
    attia_2021  = (dMdt_EL_Attia2021(False,a_earth*au2cm,e_earth,Me*1e3,Ms*1e3,Re*1e2,(fxuv/ergcm2stoWm2)*4*np.pi*a_earth*au2cm**2,Fxuv_earth_today/ergcm2stoWm2,Rxuv=Re*1e2,epsilon=epsilon)/1e3)
    isofate     = (dMdt_Cherubim2024(Me,Re,epsilon,fxuv)) * (4 * np.pi * Re ** 2) 

    escape_evolution.append({
        'Zephyrus': zephyrus, 
        'Attia+2021': attia_2021, 
        'IsoFATE': isofate 
    })


########################### Plots ####################################

# Escape vs masses
fig, ax1 = plt.subplots(figsize=(10, 7))
for model, color in zip(models, colors):
    escape_rates = [data[model] for data in escape_today]
    ax1.plot(planet_masses, escape_rates, 'o', color=color, label=model)
    ax1.plot(planet_masses, escape_rates, color=color)
ax1.set_xlabel('Planet mass [M$_{\oplus}$]', fontsize=15)
ax1.set_ylabel(r'Mass loss rate [kg $s^{-1}$]', fontsize=15)
ax1.set_title('Comparison of EL escape for Sun-Earth system', fontsize=15)
ax1.grid(alpha=0.4)
ax1.legend()
ax1.set_yscale('log')
ax2 = ax1.twinx()
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0] / s2yr) / Me, (ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'Mass loss rate [M$_{\oplus}$ $yr^{-1}$]', fontsize=15)
plt.savefig(path_plot+'Escape_vs_masses_3_models.pdf', dpi=180)

# Escape vs time
fig, ax1 = plt.subplots(figsize=(10, 7))
for model, color, label in zip(models, colors, models):
    escape_rates = [data[model] for data in escape_evolution]
    ax1.loglog(Sun_age, escape_rates, color=color, label=label)
ax1.axvline(x=age_earth/1e6, color='dimgrey', linestyle='--', linewidth=1)
ax1.text(3600, 1.1e5, 'Earth today', color='dimgrey', rotation=90, verticalalignment='bottom')
ax1.set_xlabel('Time [Myr]', fontsize=15)
ax1.set_ylabel(r'Mass loss rate [kg $s^{-1}$]', fontsize=15)
ax1.set_title('Comparison of EL escape for Sun-Earth system', fontsize=15)
ax1.grid(alpha=0.4)
ax1.legend(loc= 'upper right')
ax1.set_yscale('log')
ax2 = ax1.twinx()
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0] / s2yr) / Me, (ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'Mass loss rate [M$_{\oplus}$ $yr^{-1}$]', fontsize=15)
textstr = r'$\epsilon$ = 0.15' '\n' r'$R_p = R_{\mathrm{XUV}} = R_{\oplus}$' '\n' r'$M_p = M_{\oplus}$' '\n' r'F$_{\mathrm{XUV}}$ = MORS' '\n' r'a = a$_{\mathrm{Earth}}$' '\n' r'e = e$_{\mathrm{Earth}}$'
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax1.text(1.2, 5e4, textstr, fontsize=14,verticalalignment='top', bbox=props)
plt.savefig(path_plot+'Escape_vs_time_3_models.pdf', dpi=180)
