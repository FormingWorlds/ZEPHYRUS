''''
Emma Postolec
escape.py
Main functions to compute atmospheric escape.
'''
import numpy as np
from constants import *
from planets_parameters import *

########################################################### Energy-Limited escape (EL) ###########################################################

def dMdt_EL_Lopez2012(tidal_contribution,a,e,Mp,Ms,epsilon,Rxuv,Fxuv):

    ''''
    Compute the mass-loss rate for Energy-Limited (EL) escape.
    Formula from Lopez, Fortney & Miller, 2012 (Equation 2,3,4). Based on the formulation in Erkaev et al. 2007.
    
    Inputs :
        - tidal_contribution : 'yes'or 'no' -> K_tide = tidal correction factor (0 < K_tide < 1)                                           [dimensionless]
        - a                  : planetary semi-major axis                                                                                   [m]
        - e                  : planetary eccentricty                                                                                       [dimensionless]
        - Mp                 : planetary mass                                                                                              [kg]
        - Ms                 : Stellar mass                                                                                                [kg]
        - epsilon            : efficiency factor (varies typically 0.1 < epsilon < 0.6)                                                    [dimensionless]
        - Rxuv               : planetary radius at which XUV radiation are opticaly thick (defined at 20 mbar in Baumeister et  al. 2023)  [m]
        - Fxuv               : XUV incident flux received on the planet from the host star                                                 [W.m-2]
    Constants :
        - G                  : gravitational constant                                                                                      [m3.kg-1.s-2] 

    Output : Mass-loss rate for EL escape [kg.s-1]
    '''

    # Take into account tidal contributions : Ktide
    if tidal_contribution == 'yes':
        Rhill = a * (1-e) * (Mp/(3*Ms))**(1/3)
        ksi = Rhill/Rxuv
        K_tide = 1 - (3/(2*ksi)) + (1/(2*(ksi**3)))
    # No tidal contributions : Ktide = 1
    else :
        K_tide = 1

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * np.pi * (Rxuv**3) * Fxuv) / (G * Mp * K_tide)
    
    return escape_EL

def dMdt_EL_Attia2021(tidal_contribution,a,e,Mp,Ms,Rp,Rxuv,Lxuv,Fxuv,epsilon,orbit_effect):

    ''''
    Compute the mass-loss rate for Energy-Limited (EL) escape.
    Formula from Attia et al. 2021 (Equation 25,26,27,28,29). 
    (ATTENTION : All units in CGS)
    
    Inputs :
        - tidal_contribution : 'yes'or 'no' -> K_tide = tidal correction factor (0 < K_tide < 1)  [dimensionless]
        - a                  : planetary semi-major axis                                          [cm]
        - e                  : planetary eccentricty                                              [dimensionless]
        - Mp                 : planetary mass                                                     [g]
        - Ms                 : stellar mass                                                       [g]
        - Rp                 : planetary radius                                                   [cm]
        - Rxuv               : 'computation'or user value (usually planet radius)                 [cm] 
        - Lxuv               : XUV luminosity received on the planet from the host star           [erg.s-1]
        - Fxuv               : XUV incident flux received by the planet from the host star        [erg.cm-2.s-1]
        - epsilon            : 'computation'or user value : 0.1 < epsilon < 0.6                   [dimensionless]
        - orbit_effect       : 'yes'or 'no'                                                       [dimensionless]
    Constants :
        - G_cgs              : gravitational constant in cgs units                                [cm3.g-1.s-2] 

    Output : Mass-loss rate for EL escape [kg.s-1]
    '''

    # Take into account tidal contributions : Ktide
    if tidal_contribution == 'yes':
        ksi = (Mp/(3*Ms))**(1/3) * (a/Rp)  * (1+((e**2)/2))  
        K_tide = 1 - (3/(2*ksi)) + (1/(2*(ksi**3)))
    # No tidal contributions : Ktide = 1
    else :
        K_tide = 1
    #print('K_tide = ', K_tide)

    if Rxuv == 'computation' :
        # Compute Rxuv 
        Rxuv = Rp * 10 ** (-0.185 * np.log10(G_cgs * Mp / Rp)+ 0.021 * np.log10(Fxuv)+2.42)
        if Rxuv > 0 :
            Rxuv = Rxuv
            #print('Rxuv = ', Rxuv/1e2, 'Rearth [m]', '(Rxuv =', (Rxuv*100)/Re ,'Rearth)')     
        else :
            Rxuv = Rp
            #print('Rxuv = Rp =', Rp, 'Rearth')
    else :
        Rxuv = Rxuv
    #print('Rxuv = ',Rxuv,' [cm]')
          
    if epsilon == 'computation' :
        # Compute epsilon (escape efficiency factor)
        v = np.log10(G_cgs * Mp / Rp)
        if v <= 13.11: 
            epsilon = 10**(-0.50 - 0.44 * (v - 12.00))
            print(r'$\epsilon$ = ', epsilon)
        else:
            epsilon = 10**(-0.98 - 7.29 * (v - 13.11))
            print(r'$\epsilon$ = ', epsilon)
    else :
        epsilon = epsilon
    #print('Epsilon =',epsilon)

    # Take into account orbital contributions : orbit_effect
    if orbit_effect=='yes':
        orbit_effect = np.sqrt(1 - e**2)
    # No orbital contributions : orbit_effect = 1
    else :
        orbit_effect = 1
    #print('orbit_effect =',orbit_effect)

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * Lxuv * Rp * Rxuv**2) / (4 * G_cgs * Mp * K_tide * orbit_effect * a**2)

    return escape_EL

########################################################### Diffusion-Limited escape (DL) ###########################################################

## COMMENT THIS SECTION 
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