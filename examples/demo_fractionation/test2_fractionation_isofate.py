import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from zephyrus.constants import kb, G, M_H, M_O
from zephyrus.planets_parameters import Me, Re, Ms, Rs
from zephyrus.fractionation import Fractionation_binary_mixture
from zephyrus.binary_diffusion_coefficients import b_H_O
from isofate_binary import isocalc
from isofunks import Insolation, Luminosity

f_atm = 0.1
Mp = 1.0 * Me
Mstar = 1.0 * Ms
F0 = 1.0e-3
R_star = 1.0*Rs # [m]
T_star = 5600 # [K]
L = Luminosity(R_star, T_star) # [W]
Fp = Insolation(L, a) # [W/m^2]

a = SemiMajor(M_star, P) # [m]
d = a

solutions = (f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', species = 'H/O', rad_evol = True,
mu = mu_solar, eps = 0.15, activity = 'medium', flux_model = 'power law', stellar_type = 'M1',
Rp_override = False, t_sat = 5e8, f_atm_final = 'null', n_TO_final = 'null', 
n_steps = int(1e6), t0 = 1e6, rho_rcb = 1.0, Johnson = False, RR = True, f_pred = False, thermal = True, 
beta = -1.23)