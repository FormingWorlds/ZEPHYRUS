''''
Emma Postolec
constants.py
This file contains some often used physical constants and planets parameters.
'''

##### Physical constants #####
kb  = 1.38e-23                    # Boltzmann constant        [m2.kg.s-2.K-1] = [J.K-1]
G   = 6.6743e-11                  # Gravitational constant    [m3.kg-1.s-2]
c   = 2.99792458e8                # Speed of light            [m.s-1]

##### Earth parameters #####
Re                = 6.378e6       # Earth radius                          [m]
Me                = 5.9722e24     # Earth mass                            [kg]
Me_atm            = 5.15e18       # Mass of the Earth atmopshere          [kg]
Fxuv_earth_today  = 4.64e-3       # Stellar flux received on Earth today  [W.m-2]
age_earth         = 4.543e9       # Age of the Earth                      [years]
e_earth           = 0.017         # Earth eccentricity                    [dimensionless]
a_earth           = 1             # Earth semi-major axis                 [au]

##### Sun parameters #####
Rs      = 6.957e8                     # Solar radius        [m]
Ms      = 1.98847e30                  # Solar mass          [kg]
Ls      = 3.828e26                    # Solar luminosity    [W]
Ls_ergs = 3.839e33                    # Solar luminosity    [erg/s]
age_sun = 4.603e9                     # Age of the Sun      [years]

##### Convertion of units #####
s2yr                = 1/(3600*24*365)       # convert [seconds]   to [years]
erg2joule           = 1e-7                  # convert [ergs]      to [Joules]
au2m                = 1.496e11              # convert [au]        to [m]
ergpersecondtowatt  = 1e-7                  # convert [erg/s]     to [W]
ergcm2stoWm2        = 1e-3                  # convert [erg/cm2/s] to [W/m2]

##### TOI-561 parameters #####
## Star (Weiss+2021)
R_TOI561               = 0.832*Rs        # TOI-561 radius                          [m]
R_TOI561_errorbar      = 0.019*Rs        # Errorbars on TOI-561 radius             [m]
M_TOI561               = 0.805*Ms        # TOI-561 mass                            [kg]
M_TOI561_errorbar      = 0.030*Ms        # Errorbars on TOI-561 mass               [kg]
L_TOI561               = 0.522*Ls        # TOI-561 luminosity (bolometric ?)       [W]
L_TOI561_errorbar      = 0.017*Ls        # Errorbars on TOI-561 luminosity         [W]
age_TOI561             = 10e9            # TOI-561 age                             [years]
age_TOI561_errorbar    = 3e9             # Errorbars on TOI-561 age                [years]
## Planet b (Brinkman+2023)
R_TOI561b               = 1.37*Re        # TOI-561b radius                          [m]
R_TOI561b_errorbar      = 0.04*Re        # Errorbars on TOI-561b radius             [m]
M_TOI561b               = 2.24*Me        # TOI-561b mass                            [kg]
M_TOI561b_errorbar      = 0.20*Me        # Errorbars on TOI-561b mass               [kg]
#Fxuv_TOI561b_today      =                # Stellar flux received on TOI-561b today  [W.m-2]
#age_TOI561b             =                # Age of the TOI-561b                      [years]
e_TOI561b               = 0              # TOI-561b eccentricity                    [dimensionless]
a_TOI561b               = 0.0106         # TOI-561b semi-major axis                 [au]
a_TOI561b_errorbar      = 0.0004         # Errorbars onTOI-561b semi-major axis     [au]







