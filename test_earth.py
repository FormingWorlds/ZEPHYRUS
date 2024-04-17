''''
test_earth.py
testing EL escape on Generating plots for testing atmospheric escape
'''

import numpy as np
import matplotlib.pyplot as plt
from constants import *
import XUV_flux
import escape

##### Initial parameters ##### 
time_simulation_array = np.logspace(5,10, 1000)/s2yr                # Simulation time [s]
time_simulation_year = time_simulation_array*s2yr                   # in years (not Gyr)
initial_incident_flux = 14.67                                       # Fxuv received on Earth at t = 0.01 Gyr = 10e6 years = 10 Myr -> see Fig 9. Wordsworth+18 [W.m-2]
present_day_earth_XUV_flux = 4.64e-3                                # Fxuv received on Earth today [W.m-2]
age_earth_year = 4.543e9                                            # Age of the Earth [years]

##### 1. Compute the incident flux Fxuv received by the Earth from the Sun ##### 
incident_xuv_flux_FXUV = []
incident_xuv_flux_Johnstone = []
incident_xuv_flux_Baraffe = []

for i in range(len(time_simulation_array)) : 
    FXUV = XUV_flux.Fxuv(time_simulation_array[i], initial_incident_flux, t_sat = 50e6, beta = -1.23)
    Johnstone, L_EUV = XUV_flux.Fxuv_Johnstone(time_simulation_array[i],au2m, 'G5')
    Baraffe,age,Lstar = XUV_flux.Fxuv_Baraffe_Sun(time_simulation_array[i], au2m)
    incident_xuv_flux_FXUV.append(FXUV)
    incident_xuv_flux_Johnstone.append(Johnstone)
    incident_xuv_flux_Baraffe.append(Baraffe)

print('Incident XUV flux computation [W/m2] : ')
print('Fxuv (adapted from Ribas+2005) : ',incident_xuv_flux_FXUV[0])
print('Fxuv (Johnstone+2021)          : ',incident_xuv_flux_Johnstone[0])
print('Fxuv (Baraffe+2015)          : ',incident_xuv_flux_Baraffe[0])

print(L_EUV)

# # Plots Fxuv tracks
# plt.rcParams['figure.figsize'] = [8,6]
# plt.loglog(time_simulation_year,incident_xuv_flux_FXUV, label = 'F$_{XUV}$ (adapted from Ribas+2005)')
# plt.loglog(time_simulation_year,incident_xuv_flux_Johnstone, label = 'Johnstone+2021')
# plt.loglog(time_simulation_year,incident_xuv_flux_Baraffe, label = 'Baraffe+2015')
# plt.loglog(time_simulation_year[0],initial_incident_flux, 'o', color = 'darkgrey', label = 'Earth t = 0.01 Gyr')
# plt.loglog(age_earth_year,present_day_earth_XUV_flux, 'o', color = 'black', label = 'Earth today')
# plt.grid(alpha=0.15)
# #plt.ylim([0.000000000000000000000001,1000])
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Planetary incident flux F$_{XUV}$ [W.m$^{2}$]', fontsize=15)
# plt.legend(loc='upper right')
# plt.savefig('plots/plot_xuv_flux_comparison_models.png',dpi=180) 

plt.plot(age,Lstar*Ls_ergs,label = 'Baraffe+2015')
#plt.plot(time_simulation_year,(incident_xuv_flux_Johnstone*(4*np.pi*au2m**2))/Ls_ergs,label = 'Johnstone+2021')
#plt.plot(time_simulation_year,L_EUV,label = 'Johnstone+2021')
plt.legend()
plt.grid(alpha=0.15)
plt.xlabel('log(Time) [years]', fontsize=15)
plt.ylabel('L$_{XUV}$ [erg.s-1]', fontsize=15)
plt.savefig('plots/plot_test_baraffe_lin.png',dpi=180)