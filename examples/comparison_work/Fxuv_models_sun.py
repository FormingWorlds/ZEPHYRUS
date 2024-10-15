import numpy as np
import matplotlib.pyplot as plt
import mors

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.XUV_flux import *

########################### Initialization #####################################

time_simulation             = np.logspace(6,10, 100)/s2yr           # Simulation time                     [s]
time_simulation_Myr         = time_simulation*s2yr/1e6              # Simulation time                     [Myr]
initial_incident_flux       = 14.67                                 # Fxuv received on Earth at 10 Myr    [W.m-2]
present_day_earth_XUV_flux  = 4.64e-3                               # Fxuv received on Earth today        [W.m-2]
age_earth_year              = 4.543e9                               # Age of the Earth                    [yr]
eps                         = 0.15                                  # Escape efficiency factor            [dimensionless]

########################### Flux computation ####################################

## Fxuv for different models (functions from IsoFATE code)
vectorized_Fxuv                                     = np.vectorize(Fxuv)
incident_xuv_flux                                   = vectorized_Fxuv(time_simulation, initial_incident_flux, t_sat=50e6, beta=-1.23)
vectorized_Fxuv_Ribas                               = np.vectorize(Fxuv_Ribas)
incident_xuv_flux_Ribas                             = vectorized_Fxuv_Ribas(time_simulation)
#vectorized_Fxuv_Johnstone                           = np.vectorize(Fxuv_Johnstone)
#incident_xuv_flux_Johnstone, luminosity_Johnstone   = vectorized_Fxuv_Johnstone(time_simulation, au2m, 'G5')
vectorized_SF                                       = np.vectorize(Fxuv_SF)
incident_xuv_flux_SF                                = vectorized_SF(time_simulation,au2m)
vectorized_hazmat                                   = np.vectorize(Fxuv_hazmat)
incident_xuv_flux_hazmat                            = vectorized_hazmat(time_simulation,au2m,'medium')

## Bolometric flux from Baraffe+2015
mors.DownloadEvolutionTracks('Baraffe')
baraffe = mors.BaraffeTrack(1.0)
incident_xuv_flux_Baraffe = [baraffe.BaraffeSolarConstant(t*s2yr, 1.0)/ergcm2stoWm2/1e3 for t in time_simulation]

## Compute the Fxuv tracks from MORS using mors.Star() or mors.Lxuv() to extract Fxuv
star_data               = mors.Star(Mstar=1.0, Omega=1.0)
Star_age                = star_data.Tracks['Age']
Star_lxuv               = star_data.Tracks['Lx']+star_data.Tracks['Leuv']
Star_Fxuv               = (star_data.Tracks['Lx']+star_data.Tracks['Leuv'])/(4*np.pi*(a_earth*au2cm)**2)
ages                    = star_data.Tracks['Age']
Lxuv_lxuv               = [mors.Lxuv(Mstar=1.0, Age=age, Omega=1.0)['Lxuv'] for age in ages]
Lxuv_Fxuv               = [lxuv/(4*np.pi*(a_earth*au2cm)**2) for lxuv in Lxuv_lxuv]

## Load the Fxuv from Johnstone+2021 file : 1.0 Msun, OmegaSun = 1.0
Fxuv_johnstone2021      = Fxuv_Johnstone_Sun(time_simulation, a_earth*au2cm)                                            # [W m-2]

########################### Plot ####################################

plt.figure(figsize=(14, 12))

plt.loglog(time_simulation_Myr,incident_xuv_flux/ergcm2stoWm2, color='green', linestyle='-', label='XUV model 1 : IsoFATE : adapted from Ribas+2005 (consistent for early M dwarfs)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_Ribas/ergcm2stoWm2, color='yellowgreen', linestyle='-', label='XUV model 1 : IsoFATE : adapted from Ribas+2005')
#plt.loglog(time_simulation_Myr,incident_xuv_flux_Johnstone/ergcm2stoWm2,color = 'yellow', label = 'XUV model 2 : IsoFATE : Johnstone+2021 (G5 star, 1.0 Msun, 50 Percentile)')
plt.loglog(time_simulation_Myr,Fxuv_johnstone2021/ergcm2stoWm2, color='gold', linestyle='-', label='XUV model 2 : Johnstone+2021 (1.0 Msun, 1.0 OmegaSun, downloaded from paper)')
plt.loglog(Star_age,Star_Fxuv, color='orange', linestyle='-', label='XUV model 2 : Extracted from FWL-Mors using Star() class')
plt.loglog(Star_age,Lxuv_Fxuv, color='orange', linestyle='--', label='XUV model 2 : Extracted from FWL-Mors using Lxuv() function')
plt.loglog(time_simulation_Myr, incident_xuv_flux_Baraffe, color = 'steelblue', label = 'XUV model 3 : Baraffe+2015 (F$_{bol}$/10$^3$, 1.0 Msun)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_SF/ergcm2stoWm2, color='deeppink', linestyle='-', label='XUV model 4 : IsoFATE : Sanz-Forcada+2011 (for M to F stars)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_hazmat/ergcm2stoWm2, color='purple', linestyle='-', label='XUV model 5 : IsoFATE : HAZMAT program')
plt.axvline(x=4543, color='grey', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e2, 'Today', color='grey', rotation=90, verticalalignment='bottom')
plt.loglog(age_earth_year/1e6,present_day_earth_XUV_flux/ergcm2stoWm2, 'o', color = 'black', label = 'Earth today : t = 4 543 Myr')
plt.ylabel(r'XUV Flux at 1 AU [erg $s^{-1}$ $cm^{-2}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title('Comparison of XUV stellar evolution tracks from IsoFate, FWL-Mors and Baraffe+2015', fontsize=15)

plt.text(1.5, 1.6e4, 'XUV model 1', color='green', verticalalignment='bottom')
plt.text(1.1e1, 2e2, 'XUV model 2', color='gold', verticalalignment='bottom')
plt.text(1.1e1, 1.1e3, 'XUV model 3', color='steelblue', verticalalignment='bottom')
plt.text(1.5, 2e5, 'XUV model 4', color='deeppink', verticalalignment='bottom')
plt.text(1.5, 2.8e2, 'XUV model 5', color='purple', verticalalignment='bottom')

plt.savefig('output/Comparison_Fxuv_models_MORS_IsoFate_Baraffe_SUN.pdf',dpi=180)


