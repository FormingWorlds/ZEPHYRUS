import numpy as np
import matplotlib.pyplot as plt
import mors

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.escape import *
from zephyrus.XUV_flux import *

########################### Initialization #####################################

simulation_time         = np.logspace(6,10, 100)/s2yr               # Simulation time [s]
simulation_time_Myr     = np.logspace(6,10, 100)/1e6                # Simulation time [Myr]

Omega_sun               = 1.0                                       # Solar rotation rate        [Omega_sun, rad s-1]

escape_evolution        = []                                        # Mass Loss Rate              [kg s-1]
epsilon                 = 0.15                                      # Escape efficiency factor    [dimensionless]


########################### Extract Fxuv from different models #####################################

# Baraffe+2015
mors.DownloadEvolutionTracks('Baraffe')
baraffe = mors.BaraffeTrack(1.0)
Fxuv_Baraffe = [baraffe.BaraffeSolarConstant(t*s2yr, 1.0)/1.e3 for t in simulation_time]

# Johnstone+2021
Fxuv_johnstone2021          = Fxuv_Johnstone_Sun(simulation_time, a_earth*au2cm)                      # [W m-2]

# MORS
Sun                         = mors.Star(Mstar=1.0, Omega=Omega_sun)                         # Extract luminosities using the mors.Star() function
Sun_age                     = Sun.Tracks['Age']                                             # Age of the Sun            [Myr]
Sun_Lxuv                    = Sun.Tracks['Lx'] + Sun.Tracks['Leuv']                         # XUV luminosity            [erg s-1]
Sun_Fxuv                    = (Sun_Lxuv/(4 * np.pi * a_earth*au2cm **2)) * ergcm2stoWm2     # XUV flux                  [W m-2]

# Ribas+2005
vectorized_Fxuv             = np.vectorize(Fxuv)
Fxuv_ribas                  = vectorized_Fxuv(simulation_time, Fxuv_earth_10Myr, t_sat=50e6, beta=-1.23)


########################### Escape computations ####################################

baraffe_escape          = [EL_escape('no',a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_Baraffe]
johnstone_escape        = [EL_escape('no',a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_johnstone2021]
mors_escape             = [EL_escape('no',a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Sun_Fxuv]
ribas_escape            = [EL_escape('no',a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_ribas]

########################### Plots ####################################

# Escape vs time
fig, ax1 = plt.subplots(figsize=(10, 7))
ax1.loglog(simulation_time_Myr, baraffe_escape, color='steelblue', label='Baraffe+2015 ($F_{bol}$/$10^3$)')
ax1.loglog(simulation_time_Myr, johnstone_escape, color='orange', label='Johnstone+2021')
ax1.loglog(Sun_age, mors_escape, color='gold', label='MORS')
ax1.loglog(simulation_time_Myr, ribas_escape, color='green', label='Ribas+2005')
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
textstr = r'$\epsilon$ = 0.15' '\n' r'$R_p = R_{\mathrm{XUV}} = R_{\oplus}$' '\n' r'$M_p = M_{\oplus}$' '\n' r'a = a$_{\mathrm{Earth}}$' '\n' r'e = e$_{\mathrm{Earth}}$'
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax1.text(1.2, 4e4, textstr, fontsize=14,verticalalignment='top', bbox=props)
plt.savefig('output/Escape_vs_time_3_Fxuv_models.pdf', dpi=180)
