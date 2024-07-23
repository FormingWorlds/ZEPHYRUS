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
import mors
from constants import *
from planets_parameters import *
from XUV_flux import *


def compute_and_store_fluxes(filename, time_simulation, Fxuv_function, *args, **kwargs):
    """
    Compute flux using the provided XUV flux function and store the results in a .dat file.

    Parameters:
    - filename          : Name of the .dat file to store the results.
    - time_simulation   : Array of time points for the simulation.
    - Fxuv_function     : The XUV flux function to vectorize and compute.
    - args              : Positional arguments to pass to the XUV flux function.
    - kwargs            : Keyword arguments to pass to the XUV flux function.
    """
    # Vectorize the XUV flux function
    vectorized_function = np.vectorize(Fxuv_function)
    
    # Calculate flux values
    results = vectorized_function(time_simulation, *args, **kwargs)
    
    # Handle multiple outputs from the function
    if isinstance(results, tuple):
        incident_flux = results[0]  # Assuming the first item is the flux
    else:
        incident_flux = results
    
    # Ensure incident_flux is a one-dimensional array
    if incident_flux.ndim > 1:
        incident_flux = incident_flux.flatten()
    
    # Combine time and flux values into a single array
    data = np.column_stack((time_simulation, incident_flux))
    
    # Save the data to a .dat file
    np.savetxt(filename, data, header='Time [s]                 Flux [W m-2]', comments='')
    
    print(f"Data saved to {filename}")

def read_flux_data(filename):
    """
    Read time and flux data from a .dat file.

    Parameters:
    - filename: Name of the .dat file to read the data from.

    Returns:
    - time: Array of time values.
    - flux: Array of flux values.
    """
    # Read the data from the file
    data = np.loadtxt(filename, skiprows=1)
    
    # Separate the time and flux columns
    time = data[:, 0]
    flux = data[:, 1]
    
    return time, flux

def LoadLxuvfluxes(data_dir, filepath, d):    
    """
    Load data extracted with the Lxuv() function from FWL-Mors

    Parameters:
    - datadir  : Path of the directory containing the file with data
    - filename : Name of the file to read the data from
    - d        : Distance of the orbiting planet                        [cm]

    Returns:
    - ages           : Array of time values             [Myr]
    - xuv_luminosity : Array of XUV luminosity values   [erg s-1]
    - xuv_flux       : Array of XUV flux values         [erg s-1 cm-2]
    """
    data = np.loadtxt(data_dir+filepath, skiprows=1)

    ages = data[:, 0]
    xuv_luminosity = data[:, 1]
    xuv_flux = xuv_luminosity/(4*np.pi*(d)**2)

    return ages, xuv_luminosity, xuv_flux

def LoadStarfluxes(data_dir, filepath, d):    
    """
    Load data extracted with the Lxuv() function from FWL-Mors

    Parameters:
    - datadir  : Path of the directory containing the file with data
    - filename : Name of the file to read the data from
    - d        : Distance of the orbiting planet                        [cm]

    Returns:
    - Star_age      : Array of time values             [Myr]
    - Star_Lxuv     : Array of XUV luminosity values   [erg s-1]
    - Star_Fxuv     : Array of XUV flux values         [erg s-1 cm-2]
    """
    Star = mors.Load(data_dir+filepath)

    print(Star.Units['Lx'])

    Star_age = Star.Tracks['Age']
    Star_Lxuv = Star.Tracks['Lx']+Star.Tracks['Leuv']
    Star_Fxuv = (Star.Tracks['Lx']+Star.Tracks['Leuv'])/(4*np.pi*(d)**2)

    return Star_age, Star_Lxuv, Star_Fxuv