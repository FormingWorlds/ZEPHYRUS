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

## Initial parameters for Sun-Earth system    
time_simulation_array = np.logspace(6,10,1000)/s2yr    # Simulation time : 1 Myr - 1000 Myr     [s]

## Directory to save the files
data_isofate_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/IsoFATE_Fxuv_tracks/'
data_baraffe_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Baraffe_Fxuv_tracks/'

## Compute Fxuv for different model and store it in a .dat file
compute_and_store_fluxes(data_isofate_dir+'Fxuv_Sun-Earth_system.dat', time_simulation_array, Fxuv, Fxuv_earth_10Myr, t_sat=50e6, beta=-1.23)
compute_and_store_fluxes(data_isofate_dir+'Fxuv_Ribas_Sun-Earth_system.dat', time_simulation_array, Fxuv_Ribas)
compute_and_store_fluxes(data_isofate_dir+'Fxuv_Johnstone_Sun-Earth_system.dat', time_simulation_array, Fxuv_Johnstone, au2m, 'G5')
compute_and_store_fluxes(data_isofate_dir+'Fxuv_SF_Sun-Earth_system.dat', time_simulation_array, Fxuv_SF, au2m)
compute_and_store_fluxes(data_isofate_dir+'Fxuv_hazmat_Sun-Earth_system.dat', time_simulation_array, Fxuv_hazmat, au2m, 'medium')
compute_and_store_fluxes(data_baraffe_dir+'Fbol_Baraffe_Sun-Earth_system.dat', time_simulation_array, Fbol_Baraffe_Sun, au2m)