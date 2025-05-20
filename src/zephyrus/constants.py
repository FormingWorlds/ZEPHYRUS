''''
Emma Postolec
constants.py
This file contains some often used physical constants and unit conversions.
'''

######################################### Physical constants #########################################
kb      = 1.38e-23                    # Boltzmann constant                     [m2 kg s-2 K-1] = [J K-1]
G       = 6.6743e-11                  # Gravitational constant                 [m3 kg-1 s-2]
G_cgs   = 6.6743e-8                   # Gravitational constant in cgs units    [cm3 g-1 s-2]
c       = 2.99792458e8                # Speed of light                         [m s-1]
M_H     = 0.001008                    # Molar mass of H                        [kg/mol]
M_O     = 0.015999                    # Molar mass of O                        [kg/mol] 
R       = 8.31446                     # Ideal gas constant                     [m2 kg s-2 K-1 mol-1] = [J K-1 mol-1]
amu     = 1.66053906660e-27           # Atomic mass unit                       [kg]

######################################### Units conversions #########################################
s2yr                = 1/(3600*24*365)       # convert [seconds]      to [years]
erg2joule           = 1e-7                  # convert [ergs]         to [Joules]
au2m                = 1.496e11              # convert [au]           to [m]
au2cm               = 1.496e13              # convert [au]           to [cm]
ergpersecondtowatt  = 1e-7                  # convert [erg s-1]      to [W]
ergcm2stoWm2        = 1e-3                  # convert [erg s-1 cm-2] to [W m-2]
