''''
escape.py
Main functions for atmospheric escape.
'''

from constants import *

########################################################### Energy-Limited escape (EL) ###########################################################

def dMdt_EL_Baumeister_2023(Rp,Mp,Ratm) :

    ''''
    Compute the mass-loss rate for EL escape.
    Formula from Baumeister et  al. 2023 (Equation 15)
    
    Inputs :
        - epsilon : efficiency factor (varies typically 0.1 < epsilon < 0.6)                                    [dimensionless]
        - Rp      : planetary radius                                                                            [m]
        - Ratm    : planet radius at the top of the atmosphere (defined at 20 mbar in Baumeister et  al. 2023)  [m]
        - Fxuv    : XUV incident flux received on the planet from the host star                                 [W.m-2]
        - G       : gravitational constant                                                                      [m3.kg-1.s-2] 
        - Mp      : planetary mass                                                                              [kg]

    Output : Mass-loss rate for EL escape [kg.s-1]
    '''

    escape_EL = (epsilon * np.pi * Rp * (Ratm**2) *Fxuv) / (G*Mp)
    return escape_EL

def dMdt_EL_Lehmer_Catling_2017(Rxuv,Mp) :

    ''''
    Compute the mass-loss rate for EL escape.
    Formula from Lehmer & Catling, 2017 (Equation 1)
    
    Inputs :
        - epsilon : efficiency factor (varies typically 0.1 < epsilon < 0.6)                                                    [dimensionless]
        - Rxuv    : planetary radius at which XUV radiation are opticaly thick (defined at 20 mbar in Baumeister et  al. 2023)  [m]
        - Fxuv    : XUV incident flux received on the planet from the host star                                                 [W.m-2]
        - G       : gravitational constant                                                                                      [m3.kg-1.s-2] 
        - Mp      : planetary mass                                                                                              [kg]

    Output : Mass-loss rate for EL escape [kg.s-1]
    '''

    escape_EL = (epsilon * np.pi * (Rxuv**3) *Fxuv) / (G*Mp)
    return escape_EL

def dMdt_EL_Lopez_Fortney_Miller_2012(tidal_contribution,a,e,Mp,Ms,Rxuv,epsilon,Fxuv):

    ''''
    Compute the mass-loss rate for EL escape.
    Formula from Lopez, Fortney & Miller, 2012 (Equation 2,3,4). Based on the formulation in Erkaev et al. 2007.
    
    Inputs :
        - tidal_contribution : 'yes'or 'no'
        - a                  : planetary semi-major axis                                                                                   [dimensionless]
        - e                  : planetary eccentricty                                                                                       [dimensionless]
        - Mp                 : planetary mass                                                                                              [kg]
        - Ms                 : Stellar mass                                                                                                [kg]
        - Rxuv               : planetary radius at which XUV radiation are opticaly thick (defined at 20 mbar in Baumeister et  al. 2023)  [m]
        - epsilon            : efficiency factor (varies typically 0.1 < epsilon < 0.6)                                                    [dimensionless]
        - Fxuv               : XUV incident flux received on the planet from the host star                                                 [W.m-2]
        - G                  : gravitational constant                                                                                      [m3.kg-1.s-2] 
        - K_tide             : tidal correction factor (0 < K_tide < 1)                                                                    [dimensionless]

    Output : Mass-loss rate for EL escape [kg.s-1]
    '''

    # Take into account tidal contributions : Ktide
    if tidal_contribution == 'yes':
        Rhill = a * (1-e) * (Mp/(3*Ms))**(1.3)
        ksi = Rhill/Rxuv
        K_tide = 1 - (3/(2*ksi)) + (1/(2*(ksi**3)))
    # No tidal contributions : Ktide = 1
    else :
        K_tide = 1

    # Mass-loss rate for EL escape
    escape_EL = (epsilon * np.pi * (Rxuv**3) *Fxuv) / (G*Mp*K_tide)
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