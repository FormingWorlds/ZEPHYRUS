'''
!!! info "`constants.py`"
    Often used physical constants and unit conversions.<br>
    Authors: Emma Postolec, Mara Attia
'''

######################################### Physical constants #########################################
kb      = 1.380649e-23                # Boltzmann constant (exact SI value)    [m2 kg s-2 K-1] = [J K-1]
kb_cgs  = 1.380649e-16                # Boltzmann constant in cgs units        [erg K-1]
G       = 6.6743e-11                  # Gravitational constant                 [m3 kg-1 s-2]
G_cgs   = 6.6743e-8                   # Gravitational constant in cgs units    [cm3 g-1 s-2]
c       = 2.99792458e8                # Speed of light (exact SI value)        [m s-1]
h_planck = 6.62607015e-34             # Planck constant (exact SI value)       [J s]
m_p     = 1.67262192369e-27           # Proton mass (CODATA 2018)              [kg]
amu     = 1.66053906660e-27           # Atomic mass constant (CODATA 2018)     [kg]

# One proton crossing the planetary surface per Julian year: the smallest
# mass-loss rate with physical content on any planetary reservoir. Reported
# beside every escape verdict and never applied, since whether a rate is
# negligible is the caller's decision.
rate_floor = m_p / 3.15576e7          # Rate floor                             [kg s-1]

######################################### Units conversions #########################################
s2yr                = 1/(3600*24*365)       # convert [seconds]      to [years]
erg2joule           = 1e-7                  # convert [ergs]         to [Joules]
ev2joule            = 1.602176634e-19       # convert [eV]           to [Joules] (exact SI value)
au2m                = 1.496e11              # convert [au]           to [m]
au2cm               = 1.496e13              # convert [au]           to [cm]
ergpersecondtowatt  = 1e-7                  # convert [erg s-1]      to [W]
ergcm2stoWm2        = 1e-3                  # convert [erg s-1 cm-2] to [W m-2]
