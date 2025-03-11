''''
Emma Postolec
escape.py
Main functions to compute atmospheric escape.
'''
import numpy as np
from zephyrus.constants import *
from zephyrus.planets_parameters import *

########################################################### Energy-Limited escape (EL) ###########################################################

def EL_escape(tidal_contribution:bool, a:float,e:float,
                Mp:float,Ms:float,epsilon:float,
                Rp:float,Rxuv:float,Fxuv:float, scaling:int=2):

    ''''
    Compute the mass-loss rate for Energy-Limited (EL) escape.
    Based on the formula from Lopez, Fortney & Miller, 2012 (Equation 2,3,4).
    Alternatively based on the scaling from Lehmer & Catling (2017), Equation 1.

    Inputs :
        - tidal_contribution : bool
            True K_tilde=1, False K_tide = tidal correction factor (0 < K_tide < 1)
        - a : float
            planetary semi-major axis [m]
        - e : float
            planetary eccentricty [dimensionless]
        - Mp : float
            planetary mass [kg]
        - Ms : float
            Stellar mass [kg]
        - epsilon : float
            efficiency factor (varies typically 0.1 < epsilon < 0.6) [dimensionless]
        - Rp : float
            planetary radius [m]
        - Rxuv : float
            planetary radius at which XUV radiation are opticaly thick (defined at 20 mbar in Baumeister et  al. 2023) [m]
        - Fxuv : float
            XUV incident flux received on the planet from the host star [W m-2]

        - scaling : int
            Planet radius scaling exponent (2: R_xuv^2 * R_int, 3: R_xuv^3)

    Output :
        - escape_EL : float
            Mass-loss rate for EL escape [kg s-1]
    '''
    # Tidal contribution
    if tidal_contribution:                 # Take into account tidal contributions : Ktide
        Rhill = a * (1-e) * (Mp/(3*Ms))**(1/3)
        ksi = Rhill/Rxuv
        K_tide = 1 - (3/(2*ksi)) + (1/(2*(ksi**3)))
    else :                                           # No tidal contributions : Ktide = 1
        K_tide = 1

    # Radius
    match scaling:
        case 2:
            R_cubed = Rp * Rxuv**2
        case 3:
            R_cubed = Rxuv**3
        case _:
            raise ValueError(f"Invalid radius exponent: {scaling}")

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * np.pi * R_cubed * Fxuv) / (G * Mp * K_tide)

    return escape_EL

def dMdt_Cherubim2024(Mp:float,Rp:float,epsilon:float,Feuv:float):

    ''''
    Compute the mass-loss rate for Energy-Limited (EL) escape.
    Formula from Cherubim et al., 2024 (Equation 2).

    Inputs :
        - Mp : float
            planetary mass [kg]
        - Rp : float
            planetary radius [m]
        - epsilon : float
            efficiency factor (varies typically 0.1 < epsilon < 0.6) [dimensionless]
        - Feuv : float
            EUV incident flux received on the planet from the host star [W m-2]

    Output : float
        Mass-flux for EL escape [kg m-2 s-1]
    '''

    # Mass-flux for EL escape

    Vpot = (G * Mp) / Rp
    escape_EL = (epsilon * Feuv) / (4 * Vpot)

    return escape_EL

def dMdt_EL_Attia2021(tidal_contribution:bool,a:float,e:float,Mp:float,Ms:float,
                      Rp:float,Lxuv:float,Fxuv:float,
                      Rxuv:float=None,epsilon:float=None,orbit_effect:float=1.0):

    ''''
    Compute the mass-loss rate for Energy-Limited (EL) escape.
    Formula from Attia et al. 2021 (Equation 25,26,27,28,29).
    (ATTENTION : All units in CGS)

    Inputs :
        - tidal_contribution : bool
            True K_tilde=1, False K_tide = tidal correction factor (0 < K_tide < 1)
        - a : float
            planetary semi-major axis [m]
        - e : float
            planetary eccentricty [dimensionless]
        - Mp : float
            planetary mass [kg]
        - Ms : float
            Stellar mass [kg]
        - Rp : float
            planetary radius [m]
        - Lxuv : float
            XUV luminosity received on the planet from the host star [erg s-1]
        - Fxuv : float
            XUV incident flux received on the planet from the host star [W m-2]
        - Rxuv (optional) : float
            planetary radius at which XUV radiation are opticaly thick computed from  [m]
        - epsilon (optional) : float
            efficiency factor (varies typically 0.1 < epsilon < 0.6) [dimensionless]
        - orbit_effect (optional) : float

    Output : Mass-loss rate for EL escape [g s-1]
    '''

    # Take into account tidal contributions : Ktide
    if tidal_contribution == 'yes':
        ksi = (Mp/(3*Ms))**(1/3) * (a/Rp)  * (1+((e**2)/2))
        K_tide = 1 - (3/(2*ksi)) + (1/(2*(ksi**3)))
    # No tidal contributions : Ktide = 1
    else :
        K_tide = 1

    if Rxuv == None:
        Rxuv = Rp * 10 ** (-0.185 * np.log10(G_cgs * Mp / Rp)+ 0.021 * np.log10(Fxuv)+2.42)
        if Rxuv > 0 :
            Rxuv = Rxuv
        else :
            Rxuv = Rp

    if epsilon == None:
        v = np.log10(G_cgs * Mp / Rp)
        if v <= 13.11:
            epsilon = 10**(-0.50 - 0.44 * (v - 12.00))
        else:
            epsilon = 10**(-0.98 - 7.29 * (v - 13.11))

    # Take into account orbital contributions : orbit_effect
    if orbit_effect != 1.0:
        orbit_effect = np.sqrt(1 - e**2)

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * Lxuv * Rp * Rxuv**2) / (4 * G_cgs * Mp * K_tide * orbit_effect * a**2)

    return escape_EL

########################################################### Diffusion-Limited escape (DL) ###########################################################

## COMMENT THIS SECTION LATER
def dMdt_DL_Baumeister_2023(Ratm) :
    baj = bCO2 * chiCO2 + bCO * chiCO + bHO2 * chiHO2
    escape_DL = 4 * np.pi * (Ratm**2) *(mH2/Na) * baj * chiH2 * ((1/Ha)-(1/HH2))
    return escape_DL

def dMdt_ELtoDL_Baumeister_2023(Rp,Mp,Ratm) :
    EL = dMdt_EL_baumeister(Rp,Mp,Ratm)
    DL = dMdt_DL_baumeister(Ratm)
    fel = 1/((1+ np.exp(-(chiH2-chi0)/(w))))
    transition = fel * EL + (1-fel) * DL
    return transition
