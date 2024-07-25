''''
Emma Postolec
main.py
Main script : import XUV flux from PROTEUS simulation
'''
import numpy as np
from constants import *
from escape import *
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# __________________________ Functions __________________________ #                 # I can put this function in a document_handler.py file later ?

def open_flux_files(path) :                                                         
    data_line1 = open(path)                                                         # Load the stellar spectrum from PROTEUS simulations   
    first_line = data_line1.readline()                                              # Read the first line of the file 
    numbers = re.findall(r'\d+\.\d+', first_line)                                   # Only select the float numbers on the string
    time_step = float(numbers[0])                                                   # Convert the numbers form the string to a float value
    data = np.loadtxt(path,skiprows=1)                                              # Load the stellar spectrum from PROTEUS simulations without the first line
    selected_data = data[(data[:,0] >= 1.000e+01) & (data[:,0] <= 1.000e+02)]       # Select the data in the XUV wavelenght range : 10 nm < data < 100 nm
    xuv_wavelength = selected_data[:,0]                                             # Extract the wavelenght [nm] and flux [ergs/cm**2/s/nm] column
    xuv_flux = selected_data[:,1]                                                   
    #mean_xuv_flux = np.mean((xuv_wavelength * xuv_flux) * ergcm2stoWm2)             # Compute the mean value of XUV flux received by the planet at 1 AU in W/m2
    integrated_xuv_flux = np.trapz(xuv_flux * ergcm2stoWm2, x=xuv_wavelength)        # Integrate the flux over all wavelenght in the XUV range
    return xuv_flux,time_step,integrated_xuv_flux                                                        # Return the mean XUV flux value and the corresponding time_step

# __________________________ Extracting XUV flux from PROTEUS simulations + compute EL escape __________________________ #  

files = os.listdir('data/XUV_tracks_from_proteus/')                                 # Open all the files in this directory
sflux_files = [file for file in files if file.endswith('.sflux')]                   # Select only the '.sflux' files
def extract_number(file_name):                                                      # Function to sort the files names in the increasing order -> I can put this function in a document_handler.py file later ?
    number_part = file_name.split('.')[0]                                           # Split the filename at the dot and take the first part (the number)
    return int(number_part)
sorted_files = sorted(sflux_files, key=extract_number)                              # Sort the files using the custom key

epsilon = 0.15                                                                      # Efficiency parameter for escape
time_flux_escape = []                                                               # Create a list with 3 columns : time_step, XUV_flux, Mass_loss_rate_EL
for i in sorted_files:
    XUV_flux, time_step, integrated_xuv_flux = open_flux_files('data/XUV_tracks_from_proteus/'+i) 
    print(i)
    print(XUV_flux)
    print(integrated_xuv_flux)
    print(time_step)
    mass_loss_rate = dMdt_EL_Lopez2012('yes',a_earth*au2m,e_earth,Me,Ms,epsilon,Re,XUV_flux)
    time_flux_escape.append([time_step, XUV_flux, mass_loss_rate])
    print(integrated_xuv_flux)

print('NO ERROR UNTIL HERE')

# df = pd.DataFrame(time_flux_escape, columns = ['Time step [Myr]',                   # Put the time step, XUV flux and Escape mass loss rate in a dataframe -> easier to manipulate
#                                                'XUV Flux [W.m-2]',
#                                                'Mass loss rate (EL) [kg.s-1]'])

# # __________________________ Plots __________________________ #  

# plt.figure(figsize=(12,8))                                                          # Plot : XUV flux vs time 
# plt.plot(df['Time step [Myr]'],df['XUV Flux [W.m-2]'], ls='-', color = 'green', label = 'Mors')
# plt.grid(alpha=0.15)
# plt.legend()
# plt.xlabel('Time [Myr]', fontsize=20)
# plt.ylabel('Planetary incident flux F$_{XUV}$ [W.m$^{2}$] at 1 AU', fontsize=20)
# plt.savefig('plots/test_building_main/test_flux_from_proteus.png',dpi=180)

# plt.figure(figsize=(12,8))                                                          # Plot : XUV flux vs time - loglog scale
# plt.loglog(df['Time step [Myr]'],df['XUV Flux [W.m-2]'], ls='-', color = 'green', label = 'Mors')
# plt.grid(alpha=0.15)
# plt.legend()
# plt.xlabel('Time [Myr]', fontsize=20)
# plt.ylabel('Planetary incident flux F$_{XUV}$ [W.m$^{2}$] at 1 AU', fontsize=20)
# plt.savefig('plots/test_building_main/test_flux_from_proteus_loglog.png',dpi=180)

# plt.figure(figsize=(12,8))                                                          # Plot : Mass loss rate vs time - loglog scale
# plt.plot(df['Time step [Myr]'],df['Mass loss rate (EL) [kg.s-1]'], ls='-', color = 'blue')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [Myr]', fontsize=20)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=20)
# plt.savefig('plots/test_building_main/test_MLR_computation_from_proteus.png',dpi=180)
 
# plt.figure(figsize=(12,8))                                                          # Plot : Mass loss rate vs time - loglog scale
# plt.loglog(df['Time step [Myr]'],df['Mass loss rate (EL) [kg.s-1]'], ls='-', color = 'blue')
# plt.grid(alpha=0.15)
# plt.xlabel('Time [Myr]', fontsize=20)
# plt.ylabel('Mass loss rate [kg.s$^{-1}$]', fontsize=20)
# plt.savefig('plots/test_building_main/test_MLR_computation_from_proteus_loglog.png',dpi=180)
