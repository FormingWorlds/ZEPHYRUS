from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.fractionation import *
from zephyrus.binary_diffusion_coefficients import b_H_O

## Initialization of the parameters
Te = 288  # Earth temperature [K]
n_H = 1.0e20  # Number of H particules in the atmosphere [mol]
n_O = 1.0e20  # Number of O particules in the atmosphere [mol]
b_HO = b_H_O(Te)
Phi = 1.0e33  # Total flux of particles in the atmosphere [particles/m2/s]

# Test the general function 
ge = Acceleration_of_gravity(Mp=Me, Rp=Re)  # Gravitational field strength [m/s2]
H_light_H = Scale_height_single_species(T=Te, g=ge, M_i=M_H)
H_heavy_O = Scale_height_single_species(T=Te, g=ge, M_i=M_O)
X_light, X_heavy = Molar_concentration_binary_mixture(n_light=n_H, n_heavy=n_O)
H_atm = Scale_height_binary_mixture(T=Te, g=ge, X_light=X_light, M_light=M_H, X_heavy=X_heavy, M_heavy=M_O)

# Test individual flux functions
Phi_diffusion_light = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_H)
Phi_diffusion_heavy = Diffusion_flux(b=b_HO, T=Te, g=ge, M_i=M_O)
Phi_crit = Phi_critical(T=Te, g=ge, M_light=M_H, M_heavy=M_O, b=b_HO, x_light=X_light, n_light=n_H, n_heavy=n_O)
Phi_light, Phi_heavy = Number_flux(Phi, Phi_crit, n_light=n_H, n_heavy=n_O, M_light=M_H, M_heavy=M_O, X_light=X_light, X_heavy=X_heavy, Phi_diffusion_light=Phi_diffusion_light, Phi_diffusion_heavy=Phi_diffusion_heavy)

# Test fractionation functions
Phi_dl, Phi_df, Phi_c, Phi_l, Phi_h, Phi_tot =  Fractionation_binary_mixture(n_light=n_H, n_heavy=n_O, Mp=Me, Rp=Re, b=b_HO, T=Te, M_light=M_H, M_heavy=M_O, Phi=Phi)

# Print the results
print('------------------------------------------------------------')
print('Test of the fractionation functions for a binary mixture')
print('------------------------------------------------------------')
print(f'Binary diffusion coefficient for H and O: {b_HO:.2e} m-1 s-1')
print(f'Acceleration of gravity: {ge:.2f} m/s2')
print(f"Scale height of light species (H): {H_light_H:.2f} m")
print(f"Scale height of heavy species (O): {H_heavy_O:.2f} m")
print(f'Molar concentration of light species (H): {X_light:.2f} [1]')
print(f'Molar concentration of heavy species (O): {X_heavy:.2f} [1]')
print(f"Scale height of the atmosphere: {H_atm:.2f} m")
print('------------------------------------------------------------')
print(f'Diffusion flux of light species (H): {Phi_diffusion_light:.2e} particles/m2/s')
print(f'Diffusion flux of heavy species (O): {Phi_diffusion_heavy:.2e} particles/m2/s')
print(f'Critical flux: {Phi_crit:.2e} particles/m2/s')
print(f'Flux : {Phi:.2e} particles/m2/s')
print('------------------------------------------------------------')
print(f'Flux of light species (H): {Phi_light:.2e} particles/m2/s')
print(f'Flux of heavy species (O): {Phi_heavy:.2e} particles/m2/s')
print('------------------------------------------------------------')
print(f'Fractionation total flux: {Phi_tot:.2e} particles/m2/s')
print('------------------------------------------------------------')