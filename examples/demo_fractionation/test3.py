from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.fractionation import *
from zephyrus.binary_diffusion_coefficients import b_H_O

## Initialization of the parameters
Te = 2880  # Earth temperature [K]
time = [1,2,3,4,5]
n_H = [1e30, 1e20, 1e18, 1e15, 1e10]  # Number of H particules in the atmosphere [mol]
n_O = [1e30, 1e20, 1e18, 1e15, 1e10]  # Number of O particules in the atmosphere [mol]
b_HO = b_H_O(Te)
Phi = 1.0e35  # Total flux of particles in the atmosphere [particles/m2/s]

for t, i in enumerate(time):
    ge = Acceleration_of_gravity(Mp=Me, Rp=Re)  # Gravitational field strength [m/s2]
    H_light_H = Scale_height_single_species(T=Te, g=ge, M_i=M_H)
    H_heavy_O = Scale_height_single_species(T=Te, g=ge, M_i=M_O)
    X_light[i], X_heavy[i] = Molar_concentration_binary_mixture(n_light=n_H[i], n_heavy=n_O[i])
    H_atm = Scale_height_binary_mixture(T=Te, g=ge, X_light=X_light[i], M_light=M_H, X_heavy=X_heavy[i], M_heavy=M_O) 

    # Test individual flux functions
    Phi_diffusion_light = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_H)
    Phi_diffusion_heavy = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_O)
    Phi_crit = Phi_critical(T=Te, g=ge, M_light=M_H, M_heavy=M_O, b=b_HO, x_light=X_light[i], n_light=n_H[i], n_heavy=n_O[i])
    Phi_light, Phi_heavy = Number_flux(Phi, Phi_crit, n_light=n_H[i], n_heavy=n_O[i], M_light=M_H, M_heavy=M_O, X_light=X_light[i], X_heavy=X_heavy[i], Phi_diffusion_light=Phi_diffusion_light, Phi_diffusion_heavy=Phi_diffusion_heavy)

    # Test fractionation functions
    Phi_dl, Phi_df, Phi_c, Phi_l, Phi_h, Phi_tot =  Fractionation_binary_mixture(n_light=n_H[i], n_heavy=n_O[i], Mp=Me, Rp=Re, b=b_HO, T=Te, M_light=M_H, M_heavy=M_O, Phi=Phi)
