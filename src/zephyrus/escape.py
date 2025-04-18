''''
Emma Postolec, Harrison Nichols
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
