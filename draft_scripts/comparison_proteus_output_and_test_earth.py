
import numpy as np
import matplotlib.pyplot as plt

import os
import sys
import re
import pandas as pd

constants_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/constants.py')
sys.path.append(constants_dir)
from constants import *


planets_parameters_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/planets_parameters.py')
sys.path.append(planets_parameters_dir)
from planets_parameters import *

escape_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/escape.py')
sys.path.append(escape_dir)
from escape import *

xuv_flux_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/XUV_flux.py')
sys.path.append(xuv_flux_dir)
from XUV_flux import *

# ##### Initial parameters ##### 
# time_simulation_array = np.logspace(5,10, 1000)/s2yr                # Simulation time [s]
# time_simulation_year = time_simulation_array*s2yr                   # in years (not Gyr)
# initial_incident_flux = 14.67                                       # Fxuv received on Earth at t = 0.01 Gyr = 10e6 years = 10 Myr -> see Fig 9. Wordsworth+18 [W.m-2]
# present_day_earth_XUV_flux = 4.64e-3                                # Fxuv received on Earth today [W.m-2]
# age_earth_year = 4.543e9                                            # Age of the Earth [years]
eps = 0.15                                                          # efficiency factor for EL escape

# ##### 1. Compute the incident flux Fxuv received by the Earth from the Sun ##### 

# # Define vectorized versions of the functions
# vectorized_Fxuv = np.vectorize(Fxuv_Ribas2005)
# vectorized_Fxuv_Johnstone = np.vectorize(Fxuv_Johnstone)
# vectorized_Fxuv_Baraffe_Sun = np.vectorize(Fxuv_Baraffe_Sun)
# # Calculate all values at once
# incident_xuv_flux_Ribas = vectorized_Fxuv(time_simulation_array, initial_incident_flux, t_sat=50e6, beta=-1.23)
# incident_xuv_flux_Johnstone, luminosity_Johnstone = vectorized_Fxuv_Johnstone(time_simulation_array, au2m, 'G2')
# incident_xuv_flux_Baraffe, luminosity_Baraffe = vectorized_Fxuv_Baraffe_Sun(time_simulation_array, au2m)

# ##### 2. Compute hydrodynamic escape EL ##### 

# # Define vectorized versions of the functions
# vectorized_Lopez2012 = np.vectorize(dMdt_EL_Lopez2012)

# # Calculate all values at once
# EL_Lopez2012_flux_Baraffe2015 = vectorized_Lopez2012('yes', a_earth * au2m, e_earth, Me, Ms, eps,Re, incident_xuv_flux_Baraffe)
# EL_Lopez2012_flux_Johnstone2021 = vectorized_Lopez2012('yes', a_earth * au2m, e_earth, Me, Ms, eps,Re, incident_xuv_flux_Johnstone)
# EL_Lopez2012_flux_Ribas2005 = vectorized_Lopez2012('yes', a_earth * au2m, e_earth, Me, Ms, eps,Re, incident_xuv_flux_Ribas)




# __________________________ Functions __________________________ #                 # I can put this function in a document_handler.py file later ?

def open_flux_files(path) :                                                         
    data_line1 = open(path)                                                         # Load the stellar spectrum from PROTEUS simulations   
    first_line = data_line1.readline()                                              # Read the first line of the file 
    numbers = re.findall(r'\d+\.\d+', first_line)                                   # Only select the float numbers on the string
    time_step = float(numbers[0])                                                   # Convert the numbers form the string to a float value
    data = np.loadtxt(path,skiprows=1)                                              # Load the stellar spectrum from PROTEUS simulations without the first line
    all_wavelength = data[:,0]                                                      # Extract the wavelenght and flux column  for all wavelength range -> to get the bolometric flux 
    all_flux = data[:,1]                                                        
    integrated_all_flux = (np.trapz(all_flux, x=all_wavelength))                    # Integrate the flux over all wavelenght (bolometric flux)
    selected_data = data[(data[:,0] >= 0.517e+00) & (data[:,0] <= 9.200e+01)]       # Select the data in the X-ray + XUV wavelenght range : 0.517 nm < data < 92 nm (EUV range Mors code)
    xuv_wavelength = selected_data[:,0]                                             # Extract the wavelenght [nm] and flux [ergs/cm**2/s/nm] column
    xuv_flux = selected_data[:,1]                                                   
    integrated_xuv_flux = (np.trapz(xuv_flux, x=xuv_wavelength))                    # Integrate the flux over all wavelenght in the XUV range [ergs/cm**2/s]
    return integrated_xuv_flux,time_step, integrated_all_flux, all_wavelength, all_flux, xuv_wavelength, xuv_flux                                                # Return the mean XUV flux value and the corresponding time_step

# __________________________ Extracting XUV flux from PROTEUS simulations + compute EL escape __________________________ #  

files = os.listdir('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/XUV_tracks_from_proteus/')                                 # Open all the files in this directory
sflux_files = [file for file in files if file.endswith('.sflux')]                   # Select only the '.sflux' files
def extract_number(file_name):                                                      # Function to sort the files names in the increasing order -> I can put this function in a document_handler.py file later ?
    number_part = file_name.split('.')[0]                                           # Split the filename at the dot and take the first part (the number)
    return int(number_part)
sorted_files = sorted(sflux_files, key=extract_number)                              # Sort the files using the custom key

time_flux_escape = []    
time_flux_escape_all = []                                                             # Create a list with 3 columns : time_step, XUV_flux, Mass_loss_rate_EL
for i in sorted_files:
    integrated_xuv_flux, time_step, integrated_all_flux, all_wavelength, all_flux, xuv_wavelength, xuv_flux = open_flux_files('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/XUV_tracks_from_proteus/'+i) 
    mass_loss_rate = dMdt_EL_Lopez2012('yes', a_earth * au2m, e_earth, Me, Ms, eps,Re, integrated_xuv_flux*ergcm2stoWm2)
    mass_loss_rate_all = dMdt_EL_Lopez2012('yes', a_earth * au2m, e_earth, Me, Ms, eps,Re, integrated_all_flux*ergcm2stoWm2)
    time_flux_escape.append([time_step, integrated_xuv_flux*ergcm2stoWm2, mass_loss_rate])
    time_flux_escape_all.append([time_step, integrated_all_flux*ergcm2stoWm2, mass_loss_rate_all])

df = pd.DataFrame(time_flux_escape, columns = ['Time step [Myr]',                   # Put the time step, XUV flux and Escape mass loss rate in a dataframe -> easier to manipulate
                                               'XUV Flux [W.m-2]',
                                               'Mass loss rate (EL) [kg.s-1]'])
#print(df)

df_all = pd.DataFrame(time_flux_escape_all, columns = ['Time step [Myr]',                   # Put the time step, XUV flux and Escape mass loss rate in a dataframe -> easier to manipulate
                                               'Bolometric XUV Flux [W.m-2]',
                                               'Mass loss rate (EL) [kg.s-1]'])

# Save DataFrames to text files
save_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations_from_PROTEUS_simulation_for_Earth/'
df.to_csv(save_dir+'EL_escape_only_xuv_flux_from_PROTEUS_simulation.txt', sep='\t', index=False)
df_all.to_csv(save_dir+'EL_escape_all_wavelength_flux_from_PROTEUS_simulation.txt', sep='\t', index=False)

print("DataFrames saved to text files successfully.")

# ##### Plot : Stellar luminosity vs time #####

# ##### Plot : Stellar flux vs time #####
# plt.figure(figsize=(12,8))
# plt.loglog(time_simulation_year,incident_xuv_flux_Baraffe, color = 'steelblue', label = 'Baraffe+2015 (F$_{bol}$)')
# plt.loglog(time_simulation_year,incident_xuv_flux_Johnstone,color = 'orange', label = 'Johnstone+2021')
# plt.loglog(time_simulation_year,incident_xuv_flux_Ribas,color = 'green', label = 'Ribas+2005 (adapted)')
# plt.loglog(df['Time step [Myr]']*1e6,df['XUV Flux [W.m-2]'], ls='-', color = 'red', label = 'PROTEUS - Mors')
# plt.loglog(df_all['Time step [Myr]']*1e6,df_all['Bolometric XUV Flux [W.m-2]'], ls='-', color = 'purple', label = 'PROTEUS - Mors (F$_{bol}$)')
# plt.loglog(time_simulation_year[0],initial_incident_flux, 'o', color = 'darkgrey', label = 'Earth t = 0.01 Gyr')
# plt.loglog(age_earth_year,present_day_earth_XUV_flux, 'o', color = 'black', label = 'Earth today')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Planetary incident flux [W.m$^{2}$]', fontsize=15)
# plt.legend()
# plt.savefig('plots/test_comparison_proteus/plot_comparison_proteus_xuv_flux_comparison_models.png',dpi=180) 

# # ##### Plot : Mass loss rate vs time 

# ## Plot : Mass loss rate vs time Lopez,Fortney,Miller+2012
# plt.figure(figsize=(8,6))
# plt.loglog(time_simulation_year,EL_Lopez2012_flux_Baraffe2015, color = 'steelblue', label = 'Bolometric  model : Baraffe+2015')
# plt.loglog(time_simulation_year,EL_Lopez2012_flux_Ribas2005,color = 'green', label = 'XUV model 1 : Ribas+2005')
# plt.loglog(time_simulation_year,EL_Lopez2012_flux_Johnstone2021,color = 'orange', label = 'XUV model 2 : Johnstone+2021')

# plt.loglog(df['Time step [Myr]']*1e6,df['Mass loss rate (EL) [kg.s-1]'], ls='-', color = 'red', label = 'PROTEUS - Mors (XUV model 2)')
# plt.loglog(df_all['Time step [Myr]']*1e6,df_all['Mass loss rate (EL) [kg.s-1]'], ls='-', color = 'purple', label = 'PROTEUS - Mors (Bolometric model)')
# #plt.axvline(5e7, ls = '--', lw = 0.75, color = 'darksalmon', label = '50 Myr')
# plt.axvline(x = age_earth_year, ls='--', lw=0.7, color = 'gray', label = 't = 4.568 Gyr (Earth today)')

# plt.grid(alpha=0.15)
# plt.xlabel('Time [years]', fontsize=15)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=15)
# plt.legend()

# # Constants
# seconds_per_year = 365.25 * 24 * 3600
# # Conversion function
# def kg_s_to_earth_mass_per_year(mass_loss_rate_kg_s):
#     return mass_loss_rate_kg_s * (seconds_per_year / Me)

# # Create second y-axis

# ax2 = plt.gca().twinx()
# ax2.set_yscale('log')
# ax2.set_ylabel('Mass loss rate [M$_{\oplus}$.year$^{ -1}$]', fontsize=15)

# # Set ticks for the second y-axis
# yticks = plt.gca().get_yticks()
# yticks_earth_mass = kg_s_to_earth_mass_per_year(yticks)
# ax2.set_yticks(yticks)
# ax2.set_yticklabels([f'$10^{{{int(np.log10(y))}}}$' for y in yticks_earth_mass])



# plt.savefig('plots/test_comparison_proteus/plot_comparison_proteus_MLR_comparison_models.png',dpi=180)

# ## Extracted spectra 


# ## plot 1 
# plt.figure(figsize=(8,6))

# # All the stellar spectrum at each time step
# for i in sorted_files:
#     plt.loglog(all_wavelength, all_flux, color='darkblue')# 'Mors : F$_{bol}$')
#     #plt.axhline(y=integrated_all_flux, ls='-', lw=1.5)
#     plt.loglog(xuv_wavelength,xuv_flux, color='salmon')# 'Mors : F$_{XUV}$')
#     #plt.hlines(y=integrated_xuv_flux,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength),ls='-', lw=1.5)
# plt.loglog(all_wavelength, all_flux, color='darkblue', label= 'Mors : F$_{bol}$')
# plt.loglog(xuv_wavelength,xuv_flux, color='salmon', label= 'Mors : F$_{XUV}$')
# # Integrated flux
# #plt.axhline(y=integrated_all_flux * flux_conversion, ls='-', lw=1.5, label= 'Integrated F$_{bol}$')
# #plt.hlines(y=integrated_xuv_flux * flux_conversion,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength), color='orange', ls='-', lw=1.5, label= 'Integrated F$_{XUV}$')
# # XUV Range
# plt.axvline(x = min(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
# plt.axvline(x = max(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
# plt.text(2e0,1e-8,'X-rays + XUV range', fontsize= 10, color='gray')

# plt.grid(alpha=0.15)
# plt.xlabel('Wavelength [nm]', fontsize=15)
# plt.ylabel('Stellar flux at 1 AU [erg/c$^{2}$/s]', fontsize=15)

# plt.legend(loc='upper right')
# plt.savefig('plots/test_comparison_proteus/plot_comparison_raw_output_spectra.png',dpi=180)

# ## plot 2 
# plt.figure(figsize=(8,6))

# # All the stellar spectrum at each time step
# for i in sorted_files:
#     plt.loglog(all_wavelength, all_flux, color='darkblue')
#     plt.loglog(xuv_wavelength,xuv_flux, color='salmon')
#     plt.axhline(y=integrated_all_flux, ls='-', lw=1.5)
#     plt.hlines(y=integrated_xuv_flux,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength),ls='-', lw=1.5)

# plt.loglog(all_wavelength, all_flux, color='darkblue', label= 'Mors : F$_{bol}$')
# plt.loglog(xuv_wavelength,xuv_flux, color='salmon', label= 'Mors : F$_{XUV}$')
# # Integrated flux
# plt.axhline(y=integrated_all_flux , color= 'steelblue', ls='-', lw=1.5, label= 'Integrated F$_{bol}$')
# plt.hlines(y=integrated_xuv_flux ,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength), color='orange', ls='-', lw=1.5, label= 'Integrated F$_{XUV}$')
# # XUV Range
# plt.axvline(x = min(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
# plt.axvline(x = max(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
# plt.text(2e0,1e-8,'X-rays + XUV range', fontsize= 10, color='gray')

# plt.grid(alpha=0.15)
# plt.xlabel('Wavelength [nm]', fontsize=15)
# plt.ylabel('Stellar flux at 1 AU [erg/cm$^{2}$/s]', fontsize=15)

# # Create a secondary y-axis
# ax1 = plt.gca()
# ax2 = ax1.twinx()
# ax2.set_ylabel('Stellar flux at 1 AU [W/m$^2$]', fontsize=15)
# ax2.set_yscale('log')
# #ax2.set_ylim(ax1.get_ylim())
# #ax2.set_yticks(ax1.get_yticks())
# #ax2.set_yticklabels([f'{x * 0.0001:.1e}' for x in ax1.get_yticks()])

# ax1.legend(loc='upper right')
# plt.savefig('plots/test_comparison_proteus/plot_comparison_integrated_value.png',dpi=180)

# mlr_earth_today =  (0.15* np.pi * (6.378e6)**3 * 14.67)/(6.6743e-11*5.9722e24)
# print(mlr_earth_today, "kg.s-1")
# print(mlr_earth_today * (seconds_per_year / Me) , "Me.yr")

# def Fxuv_Ribas2005(t, F0, t_sat = 5e8, beta = -1.23):
#     '''
#     Function taken from IsoFate code for tests
#     Calculates incident XUV flux
#     Adapted from Ribas et al 2005
#     Consistent with empirical data from MUSCLES spectra for early M dwarfs

#     Inputs:
#         - t: time/age [s]
#         - F0: initial incident XUV flux [W/m2]
#         - t_sat: saturation time [yr]; change this for different stellar types (M1:500Myr, G:50Myr)
#         - beta: exponential term [ndim]

#     Output: incident XUV flux [W/m2]
#     '''
#     if t*s2yr < t_sat:
#         return F0
#     else:
#         return F0*(t*s2yr/t_sat)**beta


# test = Fxuv_Ribas2005(4.568e9/s2yr, initial_incident_flux, t_sat = 5e8, beta = -1.23)
# print(test)