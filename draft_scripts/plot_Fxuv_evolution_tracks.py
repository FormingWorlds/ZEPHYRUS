import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
mors_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/MORS-master/src/mors')
sys.path.append(mors_dir)
constants_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/constants.py')
sys.path.append(constants_dir)
planets_parameters_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/planets_parameters.py')
sys.path.append(planets_parameters_dir)
xuv_flux_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/XUV_flux.py')
sys.path.append(xuv_flux_dir)
fxuv_functions_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/tools/Fxuv_functions.py')
sys.path.append(fxuv_functions_dir)
import mors
from constants import *
from planets_parameters import *
from XUV_flux import *
from Fxuv_functions import *


##### Path to the directory #####
data_fwlmors_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/'
data_isofate_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/IsoFATE_Fxuv_tracks/'
data_baraffe_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Baraffe_Fxuv_tracks/'
plot_dir         = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/plots/test_fluxes_evolution'

##### Extract the time and flux #####
time_Fxuv, flux_Fxuv = read_flux_data(data_isofate_dir+'Fxuv_Sun-Earth_system.dat')
time_Fxuv_Ribas, flux_Fxuv_Ribas = read_flux_data(data_isofate_dir+'Fxuv_Ribas_Sun-Earth_system.dat')
time_Fxuv_Johnstone, flux_Fxuv_Johnstone = read_flux_data(data_isofate_dir+'Fxuv_Johnstone_Sun-Earth_system.dat')
time_Fxuv_SF, flux_Fxuv_SF = read_flux_data(data_isofate_dir+'Fxuv_SF_Sun-Earth_system.dat')
time_Fxuv_hazmat, flux_Fxuv_hazmat = read_flux_data(data_isofate_dir+'Fxuv_hazmat_Sun-Earth_system.dat')
time_Fbol_Baraffe, flux_Fbol_Baraffe = read_flux_data(data_baraffe_dir+'Fbol_Baraffe_Sun-Earth_system.dat')

Lxuv_filename = 'Lxuv_vs_age_Star_1.0Msun_Omega_1.0.txt'
Lxuv_age, Lxuv_Lxuv, Lxuv_Fxuv = LoadLxuvfluxes(data_fwlmors_dir+Lxuv_filename, a_earth*au2cm)
Star_filename = 'Star_1.0Msun_Myr_Omega_1.0.pickle'
Star_age, Star_Lxuv, Star_Fxuv = LoadStarfluxes(data_fwlmors_dir+Star_filename, a_earth*au2cm)


##### Plot #####
plt.figure(figsize=(12, 10))

plt.loglog(time_Fxuv*s2yr/1e6, flux_Fxuv, color='green', linestyle='-', label='IsoFATE : XUV model 1 : adapted from Ribas+2005 (consistent for early M dwarfs)')
plt.loglog(time_Fxuv_Ribas*s2yr/1e6, flux_Fxuv_Ribas, color='yellowgreen', linestyle='-', label='IsoFATE : XUV model 1 : adapted from Ribas+2005')
plt.loglog(time_Fxuv_Johnstone*s2yr/1e6, flux_Fxuv_Johnstone,color = 'yellow', label = 'IsoFATE : XUV model 2 : Johnstone+2021 (G5 star, 1.0 Msun, 50 Percentile)')
#plt.loglog(age_johnstone_2021*s2yr/1e6,F_XUV_johnstone_2021, color='gold', linestyle='-', label='XUV model 2 : Johnstone+2021 (1.0 Msun, 1.0 OmegaSun, downloaded from paper)')
#plt.loglog(Star_age,Star_Fxuv*s2yr/1e6, color='orange', linestyle='-', label='XUV model 2 : Extracted from FWL-Mors using Star() class')
#plt.loglog(Lxuv_age,Lxuv_Fxuv*s2yr/1e6, color='orange', linestyle='--', label='XUV model 2 : Extracted from FWL-Mors using Lxuv() function')
#plt.loglog(age_johnstone_2021_2*s2yr/1e6,F_XUV_johnstone_2021_2, color='peru', linestyle='-', label='XUV model 2 : Johnstone+2021 (1.0 Msun, 1.29 OmegaSun)')
plt.loglog(time_Fbol_Baraffe*s2yr/1e6, flux_Fbol_Baraffe/1e3, color = 'steelblue', label = 'XUV model 3 : Baraffe+2015 (F$_{bol}$/10$^3$, 1.0 Msun)')
plt.loglog(time_Fxuv_SF*s2yr/1e6, flux_Fxuv_SF, color='deeppink', linestyle='-', label='IsoFATE : XUV model 4 : Sanz-Forcada+2011 (for M to F stars)')
plt.loglog(time_Fxuv_hazmat*s2yr/1e6, flux_Fxuv_hazmat, color='purple', linestyle='-', label='IsoFATE : XUV model 5 : HAZMAT program')

plt.text(1.5, 1.6e1, 'XUV model 1', color='green', verticalalignment='bottom')
plt.text(1.1e1, 2e-1, 'XUV model 2', color='yellow', verticalalignment='bottom')
plt.text(1e1, 1.1e0, 'XUV model 3', color='steelblue', verticalalignment='bottom')
plt.text(1.5, 2e2, 'XUV model 4', color='deeppink', verticalalignment='bottom')
plt.text(1.5, 2.8e-1, 'XUV model 5', color='purple', verticalalignment='bottom')

plt.axvline(x=4543, color='grey', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e-1, 'Today', color='grey', rotation=90, verticalalignment='bottom')
plt.loglog(age_earth/1e6,Fxuv_earth_today, 'o', color = 'black', label = 'Earth today : t = 4 543 Myr')

plt.ylabel(r'XUV Flux at 1 AU [W $m^{-2}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title('Comparison of XUV stellar evolution tracks from IsoFate and Baraffe+2015', fontsize=15)

plt.savefig(plot_dir+'Comparison_Fxuv_models_IsoFate_Baraffe_SUN.pdf',dpi=180)
plt.show()

