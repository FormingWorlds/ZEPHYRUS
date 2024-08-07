
import numpy as np
import matplotlib.pyplot as plt

import mors 

import os
import sys
import re
zephyrus_dir = os.path.dirname('../src/zephyrus/')
sys.path.extend([zephyrus_dir])
from constants import *
from planets_parameters import *
from escape import *
from XUV_flux import *

########################### Functions ###############################

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

def extract_number(file_name):                                                      # Function to sort the files names in the increasing order -> I can put this function in a document_handler.py file later ?
    number_part = file_name.split('.')[0]                                           # Split the filename at the dot and take the first part (the number)
    return int(number_part)


########################### Path to directories ###############################

proteus_data                    = '../data/PROTEUS_simulation_Fxuv_sun_earth/'
path_plot                       = '../plots/test_comparison_proteus/'


########################### Initialization #####################################

simulation_time         = np.logspace(6,10, 100)/s2yr               # Simulation time [s]
simulation_time_Myr     = np.logspace(6,10, 100)/1e6                # Simulation time [Myr]

escape                  = []                                        # Mass Loss Rate              [kg s-1]
epsilon                 = 0.15                                      # Escape efficiency factor    [dimensionless]

########################### Extract Fxuv from different models #####################################

# Baraffe+2015
Fbol_baraffe                = Fbol_Baraffe_Sun(simulation_time,au2m)
Fxuv_Baraffe                = Fbol_baraffe/1e3                                                        # [W m-2]

# Johnstone+2021
Fxuv_johnstone2021          = Fxuv_Johnstone_Sun(simulation_time, a_earth*au2cm)                     # [W m-2]

# MORS
Sun                         = mors.Star(Mstar=1.0, Omega=1.0)                         # Extract luminosities using the mors.Star() function
Sun_age                     = Sun.Tracks['Age']
Sun_Lxuv                    = Sun.Tracks['Lx'] + Sun.Tracks['Leuv']                         # XUV luminosity            [erg s-1]
Sun_Fxuv                    = (Sun_Lxuv/(4 * np.pi * a_earth*au2cm **2)) * ergcm2stoWm2     # XUV flux                  [W m-2]

# Ribas+2005
vectorized_Fxuv             = np.vectorize(Fxuv)
Fxuv_ribas                  = vectorized_Fxuv(simulation_time, Fxuv_earth_10Myr, t_sat=50e6, beta=-1.23)

########################### Escape computations ####################################

baraffe_escape          = [EL_escape(0,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_Baraffe]
baraffe_escape_bol      = [EL_escape(0,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv*1e3) for fxuv in Fxuv_Baraffe]
johnstone_escape        = [EL_escape(0,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_johnstone2021]
mors_escape             = [EL_escape(0,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Sun_Fxuv]
ribas_escape            = [EL_escape(0,a_earth*au2m,e_earth,Me,Ms,epsilon,Re,Re,fxuv) for fxuv in Fxuv_ribas]



####################################################################################
# Proteus simulation 
########################### Organize the files #####################################
files = os.listdir(proteus_data)                        # Open all the files in this directory
sflux_files = [file for file in files if file.endswith('.sflux')]                   # Select only the '.sflux' files
sorted_files = sorted(sflux_files, key=extract_number)                              # Sort the files using the custom key
########################### Compute the escape #####################################
for i in sorted_files:
    integrated_xuv_flux, time_step, integrated_all_flux, all_wavelength, all_flux, xuv_wavelength, xuv_flux = open_flux_files(proteus_data+i)  
    mass_loss_rate      = EL_escape(0, a_earth * au2m, e_earth, Me, Ms, epsilon, Re, Re, integrated_xuv_flux*ergcm2stoWm2)
    mass_loss_rate_bol  = EL_escape(0, a_earth * au2m, e_earth, Me, Ms, epsilon, Re, Re, integrated_all_flux*ergcm2stoWm2)
    escape.append({
        'Time': time_step,
        'Fxuv': integrated_xuv_flux,
        'MLR': mass_loss_rate,
        'MLR_bol': mass_loss_rate_bol
    })

times = [entry['Time'] for entry in escape]
mlr = [entry['MLR'] for entry in escape]
mlr_bol = [entry['MLR_bol'] for entry in escape]



########################### Plot #####################################

# Escape vs time
fig, ax1 = plt.subplots(figsize=(10, 7))
ax1.loglog(simulation_time_Myr, baraffe_escape, color='steelblue', label='Baraffe+2015 ($F_{XUV}$ = $F_{bol}$/$10^3$)')
ax1.loglog(simulation_time_Myr, johnstone_escape, color='orange', label='Johnstone+2021 ($F_{XUV}$)')
ax1.loglog(Sun_age, mors_escape, color='gold', label='MORS ($F_{XUV}$)')
ax1.loglog(simulation_time_Myr, ribas_escape, color='green', label='Ribas+2005 ($F_{XUV}$)')
ax1.loglog(times, mlr, color='red', label='Proteus simulation using old MORS ($F_{XUV}$)')
ax1.loglog(simulation_time_Myr, baraffe_escape_bol, color='violet', label='Baraffe+2015 ($F_{bol}$)')
ax1.loglog(times, mlr_bol, color='purple', label='Proteus simulation using old MORS ($F_{bol}$)')
ax1.axvline(x=age_earth/1e6, color='dimgrey', linestyle='--', linewidth=1)
ax1.text(3600, 1.3e5, 'Earth today', color='dimgrey', rotation=90, verticalalignment='bottom')
ax1.set_xlabel('Time [Myr]', fontsize=15)
ax1.set_ylabel(r'Mass loss rate [kg $s^{-1}$]', fontsize=15)
ax1.set_title('Comparison of EL escape for Sun-Earth system', fontsize=15)
ax1.grid(alpha=0.4)
ax1.legend(loc='upper right')
ax2 = ax1.twinx()
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0] / s2yr) / Me, (ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'Mass loss rate [M$_{\oplus}$ $yr^{-1}$]', fontsize=15)
textstr = (r'$\epsilon$ = 0.15' '\n' r'$R_p = R_{\mathrm{XUV}} = R_{\oplus}$' '\n' r'$M_p = M_{\oplus}$' '\n' r'a = a$_{\mathrm{Earth}}$' '\n' r'e = e$_{\mathrm{Earth}}$')
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
ax1.text(1.2, 3e4, textstr, fontsize=14, verticalalignment='top', bbox=props)
plt.savefig(path_plot + 'escape_vs_time_proteus_sim.pdf', dpi=180)


# Extracted spectra 
plt.figure(figsize=(8,6))
for i in sorted_files:
    plt.loglog(all_wavelength, all_flux, color='darkblue')
    plt.loglog(xuv_wavelength,xuv_flux, color='salmon')
    plt.axhline(y=integrated_all_flux, ls='-', lw=1.5)
    plt.hlines(y=integrated_xuv_flux,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength),ls='-', lw=1.5)
plt.loglog(all_wavelength, all_flux, color='darkblue', label= 'Mors : F$_{bol}$')
plt.loglog(xuv_wavelength,xuv_flux, color='salmon', label= 'Mors : F$_{XUV}$')
# Integrated flux
plt.axhline(y=integrated_all_flux , color= 'steelblue', ls='-', lw=1.5, label= 'Integrated F$_{bol}$')
plt.hlines(y=integrated_xuv_flux ,xmin=min(xuv_wavelength),xmax=max(xuv_wavelength), color='orange', ls='-', lw=1.5, label= 'Integrated F$_{XUV}$')
# XUV Range
plt.axvline(x = min(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
plt.axvline(x = max(xuv_wavelength), ls='--', lw=0.7, color = 'gray')
plt.text(2e0,1e-8,'X-rays + XUV range', fontsize= 10, color='gray')
plt.grid(alpha=0.15)
plt.xlabel('Wavelength [nm]', fontsize=15)
plt.ylabel('Stellar flux at 1 AU [erg/cm$^{2}$/s]', fontsize=15)
ax1 = plt.gca()
ax2 = ax1.twinx()
ax2.set_ylabel('Stellar flux at 1 AU [W/m$^2$]', fontsize=15)
ax2.set_yscale('log')
ax1.legend(loc='upper right')
plt.savefig(path_plot + 'plot_comparison_integrated_value.pdf',dpi=180)
plt.show()
