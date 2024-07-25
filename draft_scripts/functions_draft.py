''''
Emma Postolec
functions_draft.py
This file contains functions used during the development stage of Zephyrus.
'''

import re
import pandas as pd

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
