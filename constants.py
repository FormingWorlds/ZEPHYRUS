''''
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
Fxuv_earth_today  = 4.64e-3       # Stellar flux received on Earth today  [W.m-2]
age_earth         = 4.543e9       # Age of the Earth                      [years]
e_earth           = 0.017         # Earth eccentricity                    [dimensionless]
a_earth           = 1             # Earth semi-major axis                 [au]

##### Sun parameters #####
Rs      = 6.957e8                     # Solar radius        [m]
Ms      = 1.98847e30                  # Solar mass          [kg]
Ls      = 3.828e26                    # Solar luminosity    [W]
Ls_ergs = 3.839e33                    # Solar luminosity    [erg/s]

##### Convertion of units #####
s2yr                = 1/(3600*24*365)       # convert [seconds] to [years]
erg2joule           = 1e-7                  # convert [ergs] to [Joules]
au2m                = 1.496e11              # convert [au] to [m]
ergpersecondtowatt  = 1e-7                  # convert [erg/s] to [W]












