''''
Emma Postolec
XUV_flux.py
Compute the XUV flux evolution of a star
'''
import numpy as np
import pathlib
from zephyrus.constants import *
from zephyrus.planets_parameters import *
import mors

########################### IsoFATE Fxuv functions ###########################
#### The following functions are taken from the IsoFATE code written by Colin Cherubim 
#### (Cherubim+2024) to test different XUV flux models in this code. 

def Fxuv(t, F0, t_sat = 5e8, beta = -1.23):
    '''
    Function taken from IsoFate code for tests (Cherubim+2024)
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
    Function taken from IsoFate code for tests (Cherubim+2024)
    Calculates incident XUV flux
    Adapted from Johnstone et al 2021 semi-empirical XUV tracks

    Inputs:
       - t: time/age [s]
       - d: orbital distance [m]
       - stellar_type: 'M1', 'K5', 'G5' [str]
    Output: incident XUV flux [W/m2]
    '''

    data_path = str(pathlib.Path(__file__).parent.absolute())+'/data/Johnstone_2021/'

    if stellar_type == 'M1':
        path = data_path + '0p5Msun_50percentile_basic.dat'
    elif stellar_type == 'K5':
        path = data_path + '0p7Msun_50percentile_basic.dat'
    elif stellar_type == 'G5':
        path = data_path + '1p0Msun_50percentile_basic.dat'
    elif stellar_type == 'Sun':
        path = data_path + '1p0Msun_1p0OmegaSun_basic.dat'
    else:
        print('Stellar type not supported. Please input "M1", "K5","G5" or "Sun"')

    data = np.loadtxt(path, unpack = True)
    age = data[0]*1e6/s2yr # [s]
    L_EUV = (data[4] + data[5] + data[6])*ergpersecondtowatt # [W]
    F_EUV = L_EUV/(4*np.pi*d**2) # [W/m2]

    return np.interp(t, age, F_EUV), np.interp(t, age, L_EUV)

def Fxuv_SF(t,d):
    '''
    Function taken from IsoFate code for tests (Cherubim+2024)
    Calculates incident XUV flux
    Adapted from Sanz-Forcada et al 2011 semi-empirical study for M to F stars

    Inputs:
       - t: time/age [s]
       - d: orbital distance [m]
    Output: incident XUV flux [W/m2]
    '''
    L_EUV = 10**(22.12 - 1.24*np.log10(t*s2yr/1e9))
    F_EUV = L_EUV/(4*np.pi*d**2) # [W/m2]
    return F_EUV

def Fxuv_Ribas(t):
    '''
    Function taken from IsoFate code for tests (Cherubim+2024)
    Calculates XUV luminosity
    Adapted from Ribas et al 2005 Sun in time program (eq 1)

    Inputs:
       - t: time/age [s]
    Output: incident XUV flux [W/m2]
    '''
    tau = t*s2yr/1e9
    F = 29.7*tau**-1.23
    F = F*ergcm2stoWm2 # [W/m2]
    L = F*4*np.pi*au2m**2
    return F # [W/m2

def Fxuv_hazmat(t, d, activity):
    '''
    Function taken from IsoFate code for tests (Cherubim+2024)
    Semi-empirical XUV flux estimates based on HAZMAT and MUSCLES
    programs. Consistent with Fxuv power law approximation.

    Inputs: t (time, s); d (orbital distance, m); activity ('high', 'medium', 'low')
    Outputs: incident planetary XUV flux [W/m2]
    '''
    flux_10myr_uq = 45400 # [W/m2]
    flux_10myr_med = 41681
    flux_10myr_lq = 30355
    flux_45myr_uq = 36661
    flux_45myr_med = 3661
    flux_45myr_lq = 30288
    flux_120myr_uq = 47629
    flux_120myr_med = 21601
    flux_120myr_lq = 21601
    flux_650myr_uq = 48673
    flux_650myr_med = 20213
    flux_650myr_lq = 4144
    flux_5000myr_uq = 3964
    flux_5000myr_med = 1146
    flux_5000myr_lq = 1100

    if activity == 'high':
        if t*s2yr/1e6 < 10:
            return flux_10myr_uq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 45:
            return flux_45myr_uq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 120:
            return flux_120myr_uq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 650:
            return flux_650myr_uq*(0.515*Rs)**2/(d**2)
        else:
            return flux_5000myr_uq*(0.515*Rs)**2/(d**2)
    elif activity == 'medium':
        if t*s2yr/1e6 < 10:
            return flux_10myr_med*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 45:
            return flux_45myr_med*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 120:
            return flux_120myr_med*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 650:
            return flux_650myr_med*(0.515*Rs)**2/(d**2)
        else:
            return flux_5000myr_med*(0.515*Rs)**2/(d**2)
    elif activity == 'low':
        if t*s2yr/1e6 < 10:
            return flux_10myr_lq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 45:
            return flux_45myr_lq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 120:
            return flux_120myr_lq*(0.515*Rs)**2/(d**2)
        elif t*s2yr/1e6 < 650:
            return flux_650myr_lq*(0.515*Rs)**2/(d**2)
        else:
            return flux_5000myr_lq*(0.515*Rs)**2/(d**2)
        
###############################################################################
        
########################### Zephyrus Fxuv functions ###########################

def Fxuv_Johnstone_Sun(t,d):
    '''
    Computes the incident XUV flux received by a planet at a distance d for a Sun-like star only (1.0 Msun, 1.0 OmegaSun)
    Using stellar evolution files from Johnstone et al. 2021

    Inputs:
        - t : Time/age                           [s]
        - d : Orbital distance of the planet     [cm]

    Output: 
        - Incident XUV flux                      [W m-2]
    '''

    path_to_file    = str(pathlib.Path(__file__).parent.absolute())+'/data/Johnstone_2021/1p0Msun_1p0OmegaSun_basic.dat'

    data        = np.loadtxt(path_to_file, unpack = True)
    age         = data[0]*1e6/s2yr                                          # [s]
    Lxuv        = (data[3] + data[4] + data[5] + data[6])                   # [erg s-1]
    Fxuv        = Lxuv/(4*np.pi*d**2)                                       # [erg s-1 cm-2]
    Fxuv_SI     = Fxuv * ergcm2stoWm2                                       # [W m-2]
    
    return np.interp(t,age,Fxuv_SI)


def Fxuv_mors(M_star,Omega_star,d):
    '''
    Computes the incident XUV flux received by a planet at a distance d for a given star mass and rotation rate
    Using stellar evolution files from MORS (Star() function)

    Inputs:
        - M_star : Mass of the star               [Msun]   
        - Omega_star : Rotation rate of the star  [Omega sun]
        - t : Time/age                            [s]
        - d : Orbital distance of the planet      [AU]

    Output: 
        - Age of the star                         [Myr]
        - Incident XUV flux                       [W m-2]
    '''
    star_data   = mors.Star(Mstar=M_star, Omega=Omega_star)                                                     # Extract luminosities using the mors.Star() function
    Star_age    = star_data.Tracks['Age']                                                                       # Age of the Star            [Myr]
    Star_Fxuv   = ((star_data.Tracks['Lx']+star_data.Tracks['Leuv'])/(4*np.pi*(d*au2cm)**2)) * ergcm2stoWm2     # XUV flux                   [W m-2]

    return Star_age,Star_Fxuv
