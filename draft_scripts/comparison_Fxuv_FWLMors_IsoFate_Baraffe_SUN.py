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
vectorized_Fxuv_Baraffe_Sun = np.vectorize(Fbol_Baraffe_Sun)
# Calculate all values at once
incident_xuv_flux = vectorized_Fxuv(time_simulation_array, initial_incident_flux, t_sat=50e6, beta=-1.23)
incident_xuv_flux_Ribas = vectorized_Fxuv_Ribas(time_simulation_array)
incident_xuv_flux_Johnstone, luminosity_Johnstone = vectorized_Fxuv_Johnstone(time_simulation_array, au2m, 'G5')
incident_xuv_flux_SF = vectorized_SF(time_simulation_array,au2m)
incident_xuv_flux_hazmat = vectorized_hazmat(time_simulation_array,au2m,'medium')
incident_xuv_flux_Baraffe, luminosity_Baraffe = vectorized_Fxuv_Baraffe_Sun(time_simulation_array, au2m)

## Functions
# Function to read data from a text file
def read_data(filepath):
    data = np.loadtxt(filepath, skiprows=1)
    ages = data[:, 0]
    xuv_luminosity = data[:, 1]
    return ages, xuv_luminosity

## Open the files and get age + Lxuv
# with Lxuv function
filename_Lxuv = 'Lxuv_vs_age_Star_1.0Msun_Omega_1.0.txt'
filepath_Lxuv = os.path.join(data_dir, filename_Lxuv)
Lxuv_age, Lxuv_lxuv = read_data(filepath_Lxuv)
Lxuv_Fxuv = Lxuv_lxuv/(4*np.pi*(a_earth*au2cm)**2)
# with Star class
filename_Star = 'Star_1.0Msun_Myr_Omega_1.0.pickle'
filepath_Star = os.path.join(data_dir, filename_Star)
Star = mors.Load(filepath_Star)
Star_age = Star.Tracks['Age']
Star_lxuv = Star.Tracks['Lx']+Star.Tracks['Leuv']
Star_Fxuv = (Star.Tracks['Lx']+Star.Tracks['Leuv'])/(4*np.pi*(a_earth*au2cm)**2)
# Johnstone+2021 : 1.0 Msun, OmegaSun = 1.0
path_to_file_johnstone_2021 = '/Users/emmapostolec/Downloads/RotationXUVTracks/TrackGrid_MstarOmega0/1p0Msun_1p0OmegaSun_basic.dat'
data_johnstone_2021 = np.loadtxt(path_to_file_johnstone_2021, unpack = True)
age_johnstone_2021 = data_johnstone_2021[0] # [Myr]
L_XUV_johnstone_2021 = (data_johnstone_2021[3] + data_johnstone_2021[4] + data_johnstone_2021[5] + data_johnstone_2021[6]) # [erg.s-1]
F_XUV_johnstone_2021 = L_XUV_johnstone_2021/(4*np.pi*(a_earth*au2cm)**2)
path_to_file_johnstone_2021_2 = '/Users/emmapostolec/Downloads/RotationXUVTracks/TrackGrid_MstarPercentile/1p0Msun_2percentile_basic.dat'
data_johnstone_2021_2 = np.loadtxt(path_to_file_johnstone_2021_2, unpack = True)
age_johnstone_2021_2 = data_johnstone_2021_2[0] # [Myr]
L_XUV_johnstone_2021_2 = (data_johnstone_2021_2[3] + data_johnstone_2021_2[4] + data_johnstone_2021_2[5] + data_johnstone_2021_2[6]) # [erg.s-1]
F_XUV_johnstone_2021_2 = L_XUV_johnstone_2021_2/(4*np.pi*(a_earth*au2cm)**2)

# PROTEUS simulations
save_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations_from_PROTEUS_simulation_for_Earth/'
df = pd.read_csv(save_dir+'EL_escape_only_xuv_flux_from_PROTEUS_simulation.txt', sep='\t')
df_all = pd.read_csv(save_dir+'EL_escape_all_wavelength_flux_from_PROTEUS_simulation.txt', sep='\t')
time_step = df['Time step [Myr]']
xuv_flux = df['XUV Flux [W.m-2]']
xuv_luminosity_proteus = xuv_flux*4*np.pi*(a_earth*au2cm)**2
mass_loss_rate = df['Mass loss rate (EL) [kg.s-1]']
time_step_all = df_all['Time step [Myr]']
bolometric_xuv_flux = df_all['Bolometric XUV Flux [W.m-2]']
mass_loss_rate_all = df_all['Mass loss rate (EL) [kg.s-1]']


## Plot
plt.figure(figsize=(14, 12))

plt.loglog(time_simulation_Myr,incident_xuv_flux/ergcm2stoWm2, color='green', linestyle='-', label='IsoFATE : XUV model 1 : adapted from Ribas+2005 (consistent for early M dwarfs)')
plt.loglog(time_simulation_Myr,incident_xuv_flux_Ribas/ergcm2stoWm2, color='yellowgreen', linestyle='-', label='IsoFATE : XUV model 1 : adapted from Ribas+2005')

plt.loglog(time_simulation_Myr,incident_xuv_flux_Johnstone/ergcm2stoWm2,color = 'yellow', label = 'IsoFATE : XUV model 2 : Johnstone+2021 (G5 star, 1.0 Msun, 50 Percentile)')
plt.loglog(age_johnstone_2021,F_XUV_johnstone_2021, color='gold', linestyle='-', label='XUV model 2 : Johnstone+2021 (1.0 Msun, 1.0 OmegaSun, downloaded from paper)')
plt.loglog(Star_age,Star_Fxuv, color='orange', linestyle='-', label='XUV model 2 : Extracted from FWL-Mors using Star() class')
plt.loglog(Lxuv_age,Lxuv_Fxuv, color='orange', linestyle='--', label='XUV model 2 : Extracted from FWL-Mors using Lxuv() function')
#plt.loglog(age_johnstone_2021_2,F_XUV_johnstone_2021_2, color='peru', linestyle='-', label='XUV model 2 : Johnstone+2021 (1.0 Msun, 1.29 OmegaSun)')

plt.loglog(time_simulation_Myr,(incident_xuv_flux_Baraffe/ergcm2stoWm2)/1e3, color = 'steelblue', label = 'XUV model 3 : Baraffe+2015 (F$_{bol}$/10$^3$, 1.0 Msun)')

plt.loglog(time_simulation_Myr,incident_xuv_flux_SF/ergcm2stoWm2, color='deeppink', linestyle='-', label='IsoFATE : XUV model 4 : Sanz-Forcada+2011 (for M to F stars)')

plt.loglog(time_simulation_Myr,incident_xuv_flux_hazmat/ergcm2stoWm2, color='purple', linestyle='-', label='IsoFATE : XUV model 5 : HAZMAT program')

plt.loglog(time_step,xuv_flux/ergcm2stoWm2, color='red', linestyle='-', label='PROTEUS simulation using old Mors version')

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
plt.text(7e1, 4e0,  'PROTEUS simulation', color='red', verticalalignment='bottom')


plt.savefig(plot_save_path+'Comparison_Fxuv_models_MORS_IsoFate_Baraffe_SUN.png',dpi=180)
plt.show()

