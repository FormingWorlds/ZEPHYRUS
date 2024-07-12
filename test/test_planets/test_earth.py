''''
Emma Postolec
test_earth.py
testing EL atmospheric escape on Earth : plots
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
eps = 0.15                                                          # efficiency factor for EL escape

##### 1. Compute the incident flux Fxuv received by the Earth from the Sun ##### 

# Define vectorized versions of the functions
vectorized_Fxuv = np.vectorize(XUV_flux.Fxuv)
vectorized_Fxuv_Johnstone = np.vectorize(XUV_flux.Fxuv_Johnstone)
vectorized_Fxuv_Baraffe_Sun = np.vectorize(XUV_flux.Fxuv_Baraffe_Sun)
# Calculate all values at once
incident_xuv_flux_Ribas = vectorized_Fxuv(time_simulation_array, initial_incident_flux, t_sat=50e6, beta=-1.23)
incident_xuv_flux_Johnstone, luminosity_Johnstone = vectorized_Fxuv_Johnstone(time_simulation_array, au2m, 'G5')
incident_xuv_flux_Baraffe, luminosity_Baraffe = vectorized_Fxuv_Baraffe_Sun(time_simulation_array, au2m)

##### 2. Compute hydrodynamic escape EL ##### 

# Define vectorized versions of the functions
vectorized_Baumeister2023 = np.vectorize(escape.dMdt_EL_Baumeister_2023)
vectorized_LehmerCatling2017 = np.vectorize(escape.dMdt_EL_Lehmer_Catling_2017)
vectorized_LopezFortneyMiller2012 = np.vectorize(escape.dMdt_EL_Lopez_Fortney_Miller_2012)

# Calculate all values at once
EL_Baumeister2023_flux_Baraffe2015 = vectorized_Baumeister2023(eps, Re, 0.1*Re, incident_xuv_flux_Baraffe, Me)
EL_Baumeister2023_flux_Johnstone2021 = vectorized_Baumeister2023(eps, Re, 0.1*Re, incident_xuv_flux_Johnstone, Me)
EL_Baumeister2023_flux_Ribas2005 = vectorized_Baumeister2023(eps, Re, 0.1*Re, incident_xuv_flux_Ribas, Me)
EL_LehmerCatling2017_flux_Baraffe2015 = vectorized_LehmerCatling2017(eps, 0.1*Re, incident_xuv_flux_Baraffe, Me)
EL_LehmerCatling2017_flux_Johnstone2021 = vectorized_LehmerCatling2017(eps, 0.1*Re, incident_xuv_flux_Johnstone, Me)
EL_LehmerCatling2017_flux_Ribas2005 = vectorized_LehmerCatling2017(eps, 0.1*Re, incident_xuv_flux_Ribas, Me)
EL_LopezFortneyMiller2012_flux_Baraffe2015 = vectorized_LopezFortneyMiller2012('yes', a_earth * au2m, e_earth, Me, Ms, eps, 0.1*Re, incident_xuv_flux_Baraffe)
EL_LopezFortneyMiller2012_flux_Johnstone2021 = vectorized_LopezFortneyMiller2012('yes', a_earth * au2m, e_earth, Me, Ms, eps, 0.1*Re, incident_xuv_flux_Johnstone)
EL_LopezFortneyMiller2012_flux_Ribas2005 = vectorized_LopezFortneyMiller2012('yes', a_earth * au2m, e_earth, Me, Ms, eps, 0.1*Re, incident_xuv_flux_Ribas)

##### Plot : Stellar luminosity vs time #####

# # Plot : Baraffe+2015 (Lbol)
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,luminosity_Baraffe, color = 'steelblue', label = 'Baraffe+2015')
# plt.legend()
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('L$_{bol}$ [Watts]', fontsize=15)
# plt.savefig('plots/test_earth/plot_test_luminosity_baraffe',dpi=180)

# # Plot : Johnstone+2021 (Leuv)
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,luminosity_Johnstone,color = 'orange', label = 'Johnstone+2021')
# plt.legend()
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('L$_{EUV}$ [Watts]', fontsize=15)
# plt.savefig('plots/test_earth/plot_test_luminosity_johnstone',dpi=180)

# # Plot : Comparison Baraffe+2015 vs Johnstone+2021
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,luminosity_Baraffe,color = 'steelblue',label = 'Baraffe+2015 (L$_{bol}$)')
# plt.loglog(time_simulation_year,luminosity_Johnstone, color = 'orange',label = 'Johnstone+2021 (L$_{EUV}$)')
# plt.legend()
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Luminosity [Watts]', fontsize=15)
# plt.savefig('plots/test_earth/plot_test_luminosity_baraffe_and_johnstone.png',dpi=180)


# ##### Plot : Stellar flux vs time #####
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,incident_xuv_flux_Baraffe, color = 'steelblue', label = 'Baraffe+2015 (F$_{bol}$)')
# plt.loglog(time_simulation_year,incident_xuv_flux_Johnstone,color = 'orange', label = 'Johnstone+2021')
# plt.loglog(time_simulation_year,incident_xuv_flux_Ribas,color = 'green', label = 'Ribas+2005 (adapted)')
# plt.loglog(time_simulation_year[0],initial_incident_flux, 'o', color = 'darkgrey', label = 'Earth t = 0.01 Gyr')
# plt.loglog(age_earth_year,present_day_earth_XUV_flux, 'o', color = 'black', label = 'Earth today')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Planetary incident flux F$_{XUV}$ [W.m$^{2}$]', fontsize=15)
# plt.legend()
# plt.savefig('plots/test_earth/plot_xuv_flux_comparison_models.png',dpi=180) 

# ##### Plot : Mass loss rate vs time 

# ## Plot : Mass loss rate vs time Baumeister+2023
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Baraffe2015, color = 'steelblue', label = 'Baraffe+2015 (F$_{bol}$)')
# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Johnstone2021,color = 'orange', label = 'Johnstone+2021')
# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Ribas2005,color = 'green', label = 'Ribas+2005 (adapted)')
# plt.axvline(5e7, ls = '--', lw = 0.75, color = 'darksalmon', label = '50 Myr')
# plt.axvline(x = age_earth_year, ls='--', lw=0.7, color = 'gray', label = 'Earth today')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=15)
# plt.title('Baumeister+2023')
# plt.legend()
# plt.savefig('plots/test_earth/plot_mass_loss_rate_vs_time_EL_BAUMEISTER_2023.png',dpi=180)

# ## Plot : Mass loss rate vs time Lehmer,Catling+2017
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Baraffe2015, color = 'steelblue', label = 'Baraffe+2015 (F$_{bol}$)')
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Johnstone2021,color = 'orange', label = 'Johnstone+2021')
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Ribas2005,color = 'green', label = 'Ribas+2005 (adapted)')
# plt.axvline(5e7, ls = '--', lw = 0.75, color = 'darksalmon', label = '50 Myr')
# plt.axvline(x = age_earth_year, ls='--', lw=0.7, color = 'gray', label = 'Earth today')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=15)
# plt.title('Lehmer,Catling+2017')
# plt.legend()
# plt.savefig('plots/test_earth/plot_mass_loss_rate_vs_time_EL_lehmercatling_2017.png',dpi=180)

# ## Plot : Mass loss rate vs time Lopez,Fortney,Miller+2012
plt.figure(figsize=(8,6))
plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Baraffe2015, color = 'steelblue', label = 'Baraffe+2015 (F$_{bol}$)')
plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Johnstone2021,color = 'orange', label = 'Johnstone+2021')
plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Ribas2005,color = 'green', label = 'Ribas+2005 (adapted)')
plt.axvline(5e7, ls = '--', lw = 0.75, color = 'darksalmon', label = '50 Myr')
plt.axvline(x = age_earth_year, ls='--', lw=0.7, color = 'gray', label = 'Earth today')
plt.grid(alpha=0.15)
plt.xlabel('Time [years]', fontsize=15)
plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=15)
#plt.title('Lopez,Fortney,Miller+2012')
plt.legend()
plt.savefig('plots/test_earth/plot_mass_loss_rate_vs_time_EL_LopezFortneyMiller2012.png',dpi=180)

# ## Plot : Mass loss rate vs time - Comparison 3 formula

# plt.figure(figsize=(12,8))

# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Baraffe2015, ls='solid', color = 'black', label = 'EL : Baumeister+2023')
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Baraffe2015, ls='dotted', color = 'black', label = 'EL : Lehmer,Catling+2017')
# plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Baraffe2015, ls='dashed', color = 'black', label = 'EL : Lopez,Fortney,Miller+2012')

# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Baraffe2015, ls='solid', color = 'steelblue', label = 'F$_{XUV}$ : Baraffe+2015 (F$_{bol}$)')
# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Johnstone2021, ls='solid', color = 'orange', label = 'F$_{XUV}$ : Johnstone+2021')
# plt.loglog(time_simulation_year,EL_Baumeister2023_flux_Ribas2005, ls='solid', color = 'green', label = 'F$_{XUV}$ : Ribas+2005 (adapted)')

# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Baraffe2015, ls='dotted', color = 'steelblue')
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Johnstone2021, ls='dotted', color = 'orange')
# plt.loglog(time_simulation_year,EL_LehmerCatling2017_flux_Ribas2005, ls='dotted', color = 'green')

# plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Baraffe2015, ls='dashed', color = 'steelblue')
# plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Johnstone2021, ls='dashed', color = 'orange')
# plt.loglog(time_simulation_year,EL_LopezFortneyMiller2012_flux_Ribas2005,  ls='dashed', color = 'green')

# plt.axvline(5e7, ls = '--', lw = 0.75, color = 'darksalmon', label = 't$_{sat}$ = 50 Myr')
# plt.axvline(x = age_earth_year, ls='--', lw=0.7, color = 'gray', label = 'Earth today')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=20)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=20)
# plt.legend(ncol=3)
# plt.savefig('plots/test_earth/plot_mass_loss_rate_vs_time_EL_comparison.png',dpi=180)


# # Define labels and functions
# labels = ['Baraffe+2015', 'Johnstone+2021', 'Ribas+2005 (adapted)']
# functions = [incident_xuv_flux_Baraffe, incident_xuv_flux_Johnstone, incident_xuv_flux_Ribas]

# # Define plot types
# plot_types = [
#     ('Luminosity vs time', luminosity_Baraffe, luminosity_Johnstone, luminosity_Johnstone, 'Luminosity [Watts]', 'plot_test_luminosity_baraffe_and_johnstone.png'),
#     ('Stellar flux vs time', incident_xuv_flux_Baraffe, incident_xuv_flux_Johnstone, incident_xuv_flux_Ribas, 'Planetary incident flux F$_{XUV}$ [W.m$^{2}$]', 'plot_xuv_flux_comparison_models.png'),
#     ('Mass loss rate vs time', EL_Baumeister2023_flux_Baraffe2015, EL_LehmerCatling2017_flux_Baraffe2015, EL_LopezFortneyMiller2012_flux_Baraffe2015, 'Mass loss rate [kg.s$^{-1}$]', 'plot_mass_loss_rate_vs_time_EL_comparison.png')
# ]

# # Plot
# for plot_type, data1, data2, data3, ylabel, filename in plot_types:
#     plt.figure(figsize=(12, 8))
#     for label, function in zip(labels, functions):
#         plt.loglog(time_simulation_year, function, label=label)
#     plt.grid(alpha=0.15)
#     plt.xlabel('Time [years]', fontsize=15)
#     plt.ylabel(ylabel, fontsize=15)
#     plt.title(plot_type, fontsize=15)
#     plt.legend()
#     plt.savefig(f'plots/test_earth/{filename}', dpi=180)
#     plt.close()
