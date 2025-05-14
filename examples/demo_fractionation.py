import numpy as np
import matplotlib.pyplot as plt
import mors

from zephyrus.constants import *
from zephyrus.planets_parameters import *
from zephyrus.fractionation import *
from zephyrus.binary_diffusion_coefficients import b_H_O

## Initialization of the parameters
Te = 288  # Earth temperature [K]
b_HO = b_H_O(Te)
m_O = 16.0 
m_H = 1.0
n_O = 2.0
n_H = 1.0  
Phi = 1.0e2  # Diffusion flux [particles/m2/s]

# Test the general function 
H = Scale_height(T=288, g=9.81, m_i=m_O)  # Scale height of a gas mixture [m]
x_H, x_O = Molar_concentration(n_1=n_H,n_2=n_O)  # Molar concentration of specie i in a gas mixture
g = Acceleration_of_gravity(Mp=Me, Rp=Re)  # Gravitational field strength [m/s2]
A = Planetary_surface_area(Rp=Re)  # Surface area of the planet [m2]

# Test the fractionation implementation
Phi_d_O = Diffusion_flux_i(b=b_HO, T=Te, g=g, m_i=m_O)  # Diffusion flux 
Phi_d_H = Diffusion_flux_i(b=b_HO, T=Te, g=g, m_i=m_H)  # Diffusion flux
Phi_crit = Phi_critical(b=b_HO, T=Te, g=g, m1=m_H, m2=m_O, x1=x_H)  # Critical diffusion flux
Phi_H, Phi_O = Number_flux(Phi, Phi_crit,  m1=m_H, m2=m_O, x1=x_H, x2=x_O, Phi_diffusion_1=Phi_d_H, Phi_diffusion_2=Phi_d_O)
MLR_H = Mass_loss_rate(A=A, Phi=Phi_H)  # Mass loss rate of H
MLR_O = Mass_loss_rate(A=A, Phi=Phi_O)  # Mass loss rate of O   
MLR_total = Mass_loss_rate(A=A, Phi=Phi)  # Total mass loss rate

a,b,c = Fractionation_binary(n1=n_H, n2=n_O, Mp=Me, Rp=Re, b=b_HO, T=Te, g=g, m1=m_H, m2=m_O, Phi=Phi)

print("------------------------------------------")
print("Results of general functions : ")
print("Scale height H =", H, "m")
print("Molar concentration x_O =", x_O)
print("Molar concentration x_H =", x_H)
print("Planet's gravity g =", g, "m/s2")
print("Surface area of the planet A =", A, "m2")
print("------------------------------------------")

print("Results of the fractionation sub-functions : ")
print("Diffusion flux Phi_d =", Phi_d_O, "particles/m2/s")
print("Critical diffusion flux Phi_crit =", Phi_crit, "particles/m2/s")
print("Number flux of H =", Phi_H, "particles/m2/s")
print("Number flux of O =", Phi_O, "particles/m2/s")
print("Mass loss rate of H MLR_H =", MLR_H, "kg/s")
print("Mass loss rate of O MLR_O =", MLR_O, "kg/s")
print("Total mass loss rate MLR_total =", MLR_total, "kg/s")
print("------------------------------------------")

print("Results of the fractionation function : ")
print("Mass loss rate of light species MLR_light =", a, "kg/s")
print("Mass loss rate of heavy species MLR_heavy =", b, "kg/s")     
print("Total mass loss rate MLR_total =", c, "kg/s")
print("------------------------------------------")
