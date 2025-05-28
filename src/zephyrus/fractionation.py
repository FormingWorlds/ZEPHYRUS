
''''
Emma Postolec, Collin Cherubim
fractionation.py
Functions to compute fractionation during atmospheric escape based on IsoFATE functions written by Collin Cherubim (Cherubim et al. 2024).
'''
import numpy as np
from zephyrus.constants import G, R, kb, amu

##################### General functions #####################

def Acceleration_of_gravity(Mp,Rp):
    '''
    Compute the gravitational field strength of a planet.
    
    Inputs:
        - Mp: float
            Mass of the planet [kg]
        - Rp: float
            Radius of the planet [m]
            
    Output: 
        - g : float 
            Gravitational field strength [m/s2]
    '''
    g = G * Mp / (Rp**2)  # Gravitational field strength [m/s2]
    return g

def Planetary_surface_area(Rp):
    '''
    Compute the surface area of a planet.
    
    Inputs:
        - Rp: float
            Radius of the planet [m]
            
    Output: 
        - A : float 
            Surface area of the planet [m2]
    '''
    A = 4 * np.pi * (Rp**2)  # Surface area of the planet [m2]
    return A

##################### Binary mixture #####################

def Effective_scale_height(T, g, M_i):
    '''
    Compute the effective scale height of a single species i at the base of the escaping region.
    
    Inputs:
        - T: float
            Temperature of the atmosphere at the base of the escaping region [K]
        - g: float
            Acceleration of gravity [m/s2]
        - M_i: float
            Molar mass of species i [kg/mol]
            
    Output: 
        - H : float 
            Effective scale height of species i at the base of the escaping region [m]
    '''
    H = R * T / (g * M_i)  # Scale height of species i     [m]
    return H

def Molar_concentration_binary_mixture(n_light, n_heavy):
    '''
    Compute the molar concentration of specie i in a binary gas mixture.
    
    Inputs:
        - n_light,n_heavy: float
            Total number of moles of light/heavy species [mol]
            
    Output: 
        - x_light,x_heavy : float 
            Molar concentration of light/heavy species  [1]
    '''
    N_tot = n_light + n_heavy
    x_light = n_light/(N_tot)
    x_heavy = n_heavy/(N_tot)
    return x_light,x_heavy

def Diffusion_flux(b,H_i):
    '''
    Compute the equivalent diffusion flux of species i relative to the primary escaping light species (usually H) in a binary gas mixture.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018

    Inputs:
        - b: float
            Binary diffusion coefficient for the considered binary mixture [particles/m/s]
        - H_i: float
            Effective scale height of species i at the base of the escaping region [m]

    Output:
        - Phi_diffusion : float
            Equivalent diffusion flux for species i [particles/m2/s]
    '''
    Phi_diffusion_i = b/H_i     # Diffusion flux of species i   [particles/m2/s]
    return Phi_diffusion_i

def Phi_critical(b, x_light, H_light, mu_light, mu_heavy):
    '''
    Compute the critical mass flux for a binary gas mixture.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - b: float
            Binary diffusion coefficient for the considered binary mixture [particles/m/s]
        - x_light: float
            Molar concentration of light species [1]
        - H_light: float
            Effective scale height of light species at the base of the escaping region [m]
        - mu_light/mu_heavy: float
            Molecular mass of light/heavy species [kg/particles]

    Output: 
        - Phi_crit : float 
            Critical mass flux [kg/m2/s] 
    '''
    Phi_crit = (b*x_light*(mu_heavy - mu_light))/H_light         # Critical mass flux    [kg/m2/s]
    return Phi_crit 

def Number_flux(Phi, Phi_crit, mu_light, mu_heavy, X_light, X_heavy, Phi_diffusion_light, Phi_diffusion_heavy):
    '''
    Compute the number flux of light and heavy species in a binary gas mixture.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018

    Inputs:
        - Phi: float
            Mass flux for the escape [kg/m2/s]
        - Phi_crit: float
            Critical mass flux [kg/m2/s]
        - mu_light/mu_heavy: float
            Molecular mass of light/heavy species [kg/particles]
        - X_light/X_heavy: float
            Molar concentration of light/heavy species [1]          
        - Phi_diffusion_light/Phi_diffusion_heavy: float
            Diffusion flux of light/heavy species [particles/m2/s]
    Output:
        - Phi_light/Phi_heavy: float
            Number flux of light/heavy species [particles/m2/s]
    '''
    mean_particle_mass = mu_light*X_light + mu_heavy*X_heavy  # Mean particle mass of the escaping outflow [kg/particle]  

    # Compute Phi_light and Phi_heavy
    if Phi < Phi_crit:
        Phi_light = Phi/mu_light
        Phi_heavy = 0
    else:    
       Phi_light = (X_light*Phi + X_light*X_heavy*(mu_heavy - mu_light)*Phi_diffusion_heavy)/mean_particle_mass
       Phi_heavy = (X_heavy*Phi + X_light*X_heavy*(mu_light - mu_heavy)*Phi_diffusion_light)/mean_particle_mass

    return Phi_light, Phi_heavy

def Fractionation_binary_mixture(Mp, Rp, T, M_light, M_heavy, b, X_light, X_heavy, mu_light, mu_heavy, Phi):
    '''
    Compute the fractionation of light and heavy species in a binary gas mixture undergoing atmospheric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - Mp: float
            Mass of the planet [kg]
        - Rp: float
            Radius of the planet [m]
        - T: float
            Temperature of the atmosphere at the base of the escaping region [K]
        - M_light/M_heavy: float
            Molar mass of light/heavy species [kg/mol]
        - b: float
            Binary diffusion coefficient for the considered binary mixture [particles/m/s]
        - X_light/X_heavy: float     
            Molar concentration of species light or heavy  [1]
        - mu_light/mu_heavy: float
            Molecular mass of light/heavy species [kg/particles]
        - Phi: float
            Mass flux for the escape [kg/m2/s]
    
    Output:
        - Phi_crit: float
            Critical mass flux [kg/m2/s]
        - Phi_light/Phi_heavy: float
            Number flux of light/heavy species [particles/m2/s]
        - Phi_total: float
            Total number flux [particles/m2/s]
    '''
    # Compute the gravity and scale height of the planet
    g = Acceleration_of_gravity(Mp, Rp)
    H_light = Effective_scale_height(T, g, M_light)
    H_heavy = Effective_scale_height(T, g, M_heavy)

    # Compute the diffusion flux of light and heavy species
    Phi_diffusion_light = Diffusion_flux(b,H_light)
    Phi_diffusion_heavy = Diffusion_flux(b,H_heavy)

    # Compute the critical mass flux
    Phi_crit = Phi_critical(H_light, mu_light, mu_heavy, b, X_light)

    # Compute the number flux of light and heavy species    
    Phi_light, Phi_heavy = Number_flux(Phi, Phi_crit, mu_light, mu_heavy, X_light, X_heavy, Phi_diffusion_light, Phi_diffusion_heavy)
    Phi_total = Phi_light + Phi_heavy

    return Phi_crit, Phi_light, Phi_heavy, Phi_total

##################### Ternary mixture #####################
