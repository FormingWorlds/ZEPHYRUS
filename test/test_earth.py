import numpy as np
import matplotlib.pyplot as plt

import mors

import os
import sys
zephyrus_dir = os.path.dirname('../src/zephyrus/')
sys.path.extend([zephyrus_dir])
from constants import *
from planets_parameters import *
from escape import EL_escape

########################### Path to directories ###############################

path_plot                       = '../plots/test_planets/'

########################### Initialization #####################################

M_sun                    = 1.0              # Star mass in solar mass                [M_sun,     kg]
Omega_sun                = 1.0              # Solar rotation rate                    [Omega_sun, rad s-1]
semi_major_axis_cm       = a_earth*au2cm    # Planetary semi-major axis              [cm]

tidal_contribution       = 'no'             # Tidal correction factor                [dimensionless]
semi_major_axis          = a_earth*au2m     # Planetary semi-major axis              [m]
eccentricity             = e_earth          # Planetary eccentricity                 [dimensionless]
M_earth                  = Me               # Planetary mass                         [kg]
epsilon                  = 0.15             # Escape efficiency factor               [dimensionless]  
Rxuv                     = Re               # XUV planetary radius                   [m]


########################### Fxuv computations #####################################

star            = mors.Star(Mstar=M_sun, Omega=Omega_sun)       # Load the stellar evolution tracks from MORS
Age_star        = star.Tracks['Age']                            # Stellar age                                      [Myr]
Lxuv_star       = star.Tracks['Lx']+star.Tracks['Leuv']         # XUV luminosity                                   [erg s-1]
Fxuv_star       = Lxuv_star/(4*np.pi*(semi_major_axis_cm)**2)   # XUV flux                                         [erg s-1 cm-2]
Fxuv_star_SI    = Fxuv_star*ergcm2stoWm2                        # XUV flux                                         [W m-2]


########################### EL escape computations #####################################

escape = EL_escape(tidal_contribution, semi_major_axis, eccentricity, M_earth, M_sun, epsilon, Rxuv, Fxuv_star_SI)   # Compute EL escape     [kg s-1]


########################### Plot #####################################

fig, ax1 = plt.subplots(figsize=(10, 8))
ax1.loglog(Age_star, escape, '-', color='orange', label='MORS stellar evolution tracks')
ax1.set_xlabel('Time [Myr]', fontsize=15)
ax1.set_ylabel(r'Mass loss rate [kg $s^{-1}$]', fontsize=15)
ax1.set_title('Zephyrus : EL escape for Sun-Earth system', fontsize=15)
ax1.grid(alpha=0.4)
ax1.legend()
ax1.set_yscale('log')
ax2 = ax1.twinx()
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0]/ s2yr) / Me,(ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'Mass loss rate [$M_{\oplus}$ $yr^{-1}$]', fontsize=15)
plt.savefig(path_plot+'test_earth_Escape_vs_time_MORS.pdf', dpi=180)


