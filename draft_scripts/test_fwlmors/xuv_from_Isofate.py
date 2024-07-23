import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
# To import functions from mors
mors_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/MORS-master/src/mors')
sys.path.append(mors_dir)
constants_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/constants.py')
sys.path.append(constants_dir)
planets_parameters_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/planets_parameters.py')
sys.path.append(planets_parameters_dir)
from planets_parameters import *
import mors
from constants import *
xuv_flux_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/XUV_flux.py')
sys.path.append(xuv_flux_dir)
from XUV_flux import *

## Path to directories
data_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/'
plot_save_path = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/plots/'

# ##### Initial parameters ##### 
time_simulation_array = np.logspace(6,10, 100)/s2yr                # Simulation time [s]
time_simulation_Myr = time_simulation_array*s2yr/1e6                  # in years (not Gyr)
initial_incident_flux = 14.67                                    # Fxuv received on Earth at t = 0.01 Gyr = 10e6 years = 10 Myr -> see Fig 9. Wordsworth+18 [W.m-2]
present_day_earth_XUV_flux = 4.64e-3                                # Fxuv received on Earth today [W.m-2]
age_earth_year = 4.543e9                                         # Age of the Earth [years]
eps = 0.15  

## Compute flux from Isofate functions 
vectorized_Fxuv = np.vectorize(Fxuv)
vectorized_Fxuv_Ribas = np.vectorize(Fxuv_Ribas)
vectorized_Fxuv_Johnstone = np.vectorize(Fxuv_Johnstone)
vectorized_SF = np.vectorize(Fxuv_SF)
vectorized_hazmat = np.vectorize(Fxuv_hazmat)
vectorized_Fxuv_Baraffe_Sun = np.vectorize(Fxuv_Baraffe_Sun)
# Calculate all values at once
incident_xuv_flux = vectorized_Fxuv(time_simulation_array, initial_incident_flux, t_sat=50e6, beta=-1.23)
incident_xuv_flux_Ribas = vectorized_Fxuv_Ribas(time_simulation_array)
incident_xuv_flux_Johnstone, luminosity_Johnstone = vectorized_Fxuv_Johnstone(time_simulation_array, au2m, 'G5')
incident_xuv_flux_SF = vectorized_SF(time_simulation_array,au2m)
incident_xuv_flux_hazmat = vectorized_hazmat(time_simulation_array,au2m,'medium')

incident_xuv_flux_Baraffe, luminosity_Baraffe = vectorized_Fxuv_Baraffe_Sun(time_simulation_array, au2m)


## Plot
plt.figure(figsize=(10, 8))

plt.loglog(time_simulation_Myr,incident_xuv_flux/ergcm2stoWm2, color='green', linestyle='-', label='XUV model 1 : adapted from Ribas+2005 (consistent for early M dwarfs)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_Ribas/ergcm2stoWm2, color='yellowgreen', linestyle='-', label='XUV model 1 : adapted from Ribas+2005')
plt.loglog(time_simulation_Myr,incident_xuv_flux_Johnstone/ergcm2stoWm2,color = 'yellow', label = 'XUV model 2 : Johnstone+2021 (G5 star, 1.0 Msun, 50 Percentile)')
#plt.loglog(time_simulation_Myr,(incident_xuv_flux_Baraffe/ergcm2stoWm2)/1e3, color = 'steelblue', label = 'XUV model 3 : Baraffe+2015 (F$_{bol}$/10$^3$, 1.0 Msun, NOT IN ISOFATE)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_SF/ergcm2stoWm2, color='deeppink', linestyle='-', label='XUV model 4 : Sanz-Forcada+2011 (for M to F stars)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_hazmat/ergcm2stoWm2, color='purple', linestyle='-', label='XUV model 5 : HAZMAT program')

plt.loglog(age_earth_year/1e6,present_day_earth_XUV_flux/ergcm2stoWm2, 'o', color = 'black', label = 'Earth today : 4 543 Myr')

plt.ylabel(r'XUV Flux at 1 AU [erg $s^{-1}$ $cm^{-2}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title('XUV stellar evolution tracks from IsoFate', fontsize=15)

# Annotate specific points and lines 
plt.axvline(x=4543, color='black', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e4, 'Today', rotation=90, verticalalignment='bottom')

#plt.text(1.5, 1.6e4, 'XUV model 1', color='green', verticalalignment='bottom')
#plt.text(1.5, 4e3, 'XUV model 2', color='gold', verticalalignment='bottom')

plt.savefig(plot_save_path+'Comparison_Fxuv_models_MORS_IsoFate.png',dpi=180)
plt.show()

