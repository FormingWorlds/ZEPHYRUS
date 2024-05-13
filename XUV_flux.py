''''
Emma Postolec
XUV_flux.py
Compute the XUV flux evolution of a star
'''

from constants import *
import numpy as np

def Fxuv(t, F0, t_sat = 5e8, beta = -1.23):
    '''
    Function taken from IsoFate code for tests
    Calculates incident XUV flux
    Adapted from Ribas et al 2005
    Consistent with empirical data from MUSCLES spectra for early M dwarfs

    Inputs:
        - t: time/age [s]
        - F0: initial incident XUV flux [W/m2]
        - t_sat: saturation time [yr]; change this for different stellar types (M1:500Myr, G:50Myr)
        - beta: exponential term [ndim]

    Output: incident XUV flux [W/m2]
    '''
    if t*s2yr < t_sat:
        return F0
    else:
        return F0*(t*s2yr/t_sat)**beta

def Fxuv_Johnstone(t, d, stellar_type):
    '''
    Function taken from IsoFate code for tests
    Calculates incident XUV flux
    Adapted from Johnstone et al 2021 semi-empirical XUV tracks

    Inputs:
       - t: time/age [s]
       - d: orbital distance [m]
       - stellar_type: 'M1', 'K5', 'G5' [str]
    Output: incident XUV flux [W/m2]
    '''

    if stellar_type == 'M1':
        path = 'data/RotationXUVTracks/TrackGrid_MstarPercentile/0p5Msun_50percentile_basic.dat'
    elif stellar_type == 'K5':
        path = 'data/RotationXUVTracks/TrackGrid_MstarPercentile/0p7Msun_50percentile_basic.dat'
    elif stellar_type == 'G5':
        path = 'data/RotationXUVTracks/TrackGrid_MstarPercentile/1p0Msun_50percentile_basic.dat'
    else:
        print('Stellar type not supported. Please input "M1", "K5" or "G5"')

    data = np.loadtxt(path, unpack = True)
    age = data[0]*1e6/s2yr # [s]
    L_EUV = (data[4] + data[5] + data[6])*ergpersecondtowatt # [W]
    F_EUV = L_EUV/(4*np.pi*d**2) # [W/m2]

    return np.interp(t, age, F_EUV), np.interp(t, age, L_EUV)

def Fxuv_Baraffe_Sun(t,d):
    '''
    Computes incident XUV flux for the Sun only 
    Using stellar evolution files from Baraffe et al. 2015

    Inputs:
        - t: time/age [s]
        - d: orbital distance [m]

    Output: incident XUV flux [W/m2]
    '''
    
    path = 'data/Baraffe_stellar_evolution/baraffe_XUV_Sun.dat'
    data = np.loadtxt(path, usecols=(1,3), skiprows=3)
    logt  = data[:,0] # time in years
    logLstar = data[:,1] # luminosity of the star in solar unit
    age = (10**logt)/s2yr # time in years
    Lstar = ((10**logLstar)*Ls_ergs)*ergpersecondtowatt # W
    Flux= (Lstar)/(4*np.pi*d**2) # W/m2

    return np.interp(t,age,Flux),np.interp(t,age,Lstar)


