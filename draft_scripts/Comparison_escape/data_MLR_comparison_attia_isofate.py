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
escape_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/escape.py')
sys.path.append(escape_dir)
xuv_flux_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/XUV_flux.py')
sys.path.append(xuv_flux_dir)
fxuv_functions_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/tools/Fxuv_functions.py')
sys.path.append(fxuv_functions_dir)
isofate_function_escape = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/IsoFATE-main/isofunks.py')
sys.path.append(isofate_function_escape)
from escape import *
from planets_parameters import *
import mors
from constants import *
from XUV_flux import *
from Fxuv_functions import *
from isofunks import phi_E

##### Path to the directory #####
data_fwlmors_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/'

##### Extract the time, luminosity and flux from stellar evolution files #####
Lxuv_filename = 'Lxuv_vs_age_Star_1.0Msun_Omega_1.0.txt'
Lxuv_age, Lxuv_Lxuv, Lxuv_Fxuv = LoadLxuvfluxes(data_fwlmors_dir,Lxuv_filename, a_earth*au2cm)
Star_filename = 'Star_1.0Msun_Myr_Omega_1.0.pickle'
Star_age, Star_Lxuv, Star_Fxuv = LoadStarfluxes(data_fwlmors_dir,Star_filename, a_earth*au2cm)

##### Compute the escape with different formula using the data from stellar evolution files #####
mara_escape_Lxuv = []
zephyrus_escape_Lxuv = []
isofate_escape_Lxuv = []
Vpot = (G*Me)/(Re**3)
for i in range(len(Lxuv_age)) : 
    mara_escape_Lxuv.append((dMdt_EL_Attia2021('no',a_earth*au2cm,e_earth,Me*1e3,Ms*1e3,Re*1e2,Re*1e2,Lxuv_Lxuv[i],Lxuv_Fxuv[i],0.15,'no'))/1e3)
    zephyrus_escape_Lxuv.append(dMdt_EL_Lopez2012('no',a_earth*au2m,e_earth,Me,Ms,0.15,Re,Lxuv_Fxuv[i]*ergcm2stoWm2))
    isofate_escape_Lxuv.append(0.15*Lxuv_Fxuv[i]*ergcm2stoWm2/(4*Vpot))

mara_escape_Star = []
zephyrus_escape_Star = []
isofate_escape_Star = []
for i in range(len(Star_age)) :
    mara_escape_Star.append((dMdt_EL_Attia2021('no',a_earth*au2cm,e_earth,Me*1e3,Ms*1e3,Re*1e2,Re*1e2,Star_Lxuv[i],Star_Fxuv[i],0.15,'no'))/1e3)
    zephyrus_escape_Star.append(dMdt_EL_Lopez2012('no',a_earth*au2m,e_earth,Me,Ms,0.15,Re,Star_Fxuv[i]*ergcm2stoWm2))
    isofate_escape_Star.append(0.15*Star_Fxuv[i]*ergcm2stoWm2/(4*Vpot))

# Combine time and flux values into a single array
data_mara_Lxuv = np.column_stack((Lxuv_age, mara_escape_Lxuv))
data_mara_Star = np.column_stack((Star_age, mara_escape_Star))
data_zephyrus_Lxuv = np.column_stack((Lxuv_age, zephyrus_escape_Lxuv))
data_zephyrus_Star = np.column_stack((Star_age, zephyrus_escape_Star))  
data_isofate_Lxuv = np.column_stack((Lxuv_age, isofate_escape_Lxuv))
data_isofate_Star = np.column_stack((Star_age, isofate_escape_Star))  


# # Save the data to a .dat file
np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Attia+2021_with_Lxuv.txt', data_mara_Lxuv, header='Time [Myr]                 MLR (Lxuv) [kg s-1]', comments='')
np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Attia+2021_with_Star.txt', data_mara_Star, header='Time [Myr]                 MLR (Star) [kg s-1]', comments='')
print(f"Data saved : Attia+2021 Lxuv")

np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Lopez+2012_with_Lxuv.txt', data_zephyrus_Lxuv, header='Time [Myr]                 MLR (Lxuv) [kg s-1]', comments='')
np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Lopez+2012_with_Star.txt', data_zephyrus_Star, header='Time [Myr]                 MLR (Star) [kg s-1]', comments='')
print(f"Data saved : Lopez+2012 (Zephyrus)")

np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Cherubim+2024_with_Lxuv.txt', data_isofate_Lxuv, header='Time [Myr]                 MLR (Lxuv) [kg s-1]', comments='')
np.savetxt('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/MLR_computation_Cherubim+2024_with_Star.txt', data_isofate_Star, header='Time [Myr]                 MLR (Star) [kg s-1]', comments='')
print(f"Data saved : Cherubim+2024 (IsoFate)")