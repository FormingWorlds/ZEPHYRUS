''''
Emma Postolec
constants.py
This file contains star-planet system parameters for running escape simulations.
'''

######################################### Sun-Earth system #########################################

# Sun parameters 
Rs      = 6.957e8                      # Solar radius                          [m]
Ms      = 1.98847e30                   # Solar mass                            [kg]
Ls      = 3.828e26                     # Solar luminosity                      [W]
Ls_ergs = 3.839e33                     # Solar luminosity                      [erg s-1]
age_sun = 4.603e9                      # Age of the Sun                        [yr]

# Earth parameters 
Re                = 6.378e6            # Earth radius                          [m]
Me                = 5.9722e24          # Earth mass                            [kg]
Me_atm            = 5.15e18            # Mass of the Earth atmopshere          [kg]
Fxuv_earth_10Myr  = 14.67              # Fxuv received on Earth at t = 10 Myr -> see Fig 9. Wordsworth+18 [W m-2]
Fxuv_earth_today  = 4.64e-3            # Stellar flux received on Earth today  [W m-2]
age_earth         = 4.543e9            # Age of the Earth                      [yr]
e_earth           = 0.017              # Earth eccentricity                    [dimensionless]
a_earth           = 1                  # Earth semi-major axis                 [au]

