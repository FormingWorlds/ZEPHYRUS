
''''
Emma Postolec, Collin Cherubim
fractionation.py
Functions to compute fractionation during atmospheric escape based on IsoFATE functions written by Collin Cherubim (Cherubim et al. 2024).
'''
import numpy as np
from zephyrus.constants import kb, G, amu, R

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

def Mass_species(n, M):
    '''
    Compute the mass of a species in the atmosphere.
    
    Inputs:
        - n: float
            Number of moles of the species [mol]
        - M: float
            Molar mass of the species [kg/mol]
            
    Output: 
        - m : float 
            Mass of the species [kg]
    '''
    m = n * M  # Mass of the species [kg]
    return m

##################### Binary mixture #####################

def Scale_height_single_species(T, g, M_i):
    '''
    Compute the scale height of a single species.
    
    Inputs:
        - T: float
            Temperature of the atmosphere [K]
        - g: float
            Acceleration of gravity [m/s2]
        - M_i: float
            Molar mass of species i [kg/mol]
            
    Output: 
        - H : float 
            Scale height of the atmosphere [m]
    '''
    H = R * T / (M_i * g)  # Scale height of species i     [m]
    return H

def Scale_height_binary_mixture(T, g, X_light, M_light, X_heavy, M_heavy):
    '''
    Compute the scale height of a gas mixture.
    
    Inputs:
        - T: float
            Temperature of the atmosphere [K]
        - g: float
            Acceleration of gravity [m/s2]
        - X_light: float
            Molar fraction of light species [1]
        - M_light: float
            Molar mass of light species [kg/mol]
        - X_heavy: float
            Molar fraction of heavy species [1]                                             
        - M_heavy: float
            Molar mass of heavy species [kg/mol]       
        
    Output: 
        - H : float 
            Scale height of the atmosphere [m]
    '''
    mu_mix = X_light*M_light + X_heavy*M_heavy  # Mean molecular weight (MMW) of the atmosphere [kg/mol]
    H = kb * T / (mu_mix * amu *  g)            # Scale height of the atmosphere     [m]
    return H

def Molar_concentration_binary_mixture(n_light, n_heavy):
    '''
    Compute the molar concentration of specie i in a gas mixture.
    
    Inputs:
        - n_light,n_heavy: float
            Total number of moles of species light/heavy [mol]
            
    Output: 
        - x_light,x_heavy : float 
            Molar concentration of species light or heavy  [1]
    '''
    N_tot = n_light + n_heavy
    x_light = n_light/(N_tot)
    x_heavy = n_heavy/(N_tot)
    return x_light,x_heavy

def Diffusion_flux(b, T, g, M_i):
    '''
    Compute the equivalent diffusion flux of species i relative to  primary escaping species the in a binary gas mixture.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    Inputs:
        - b: float
            Binary diffusion coefficient for the considered binary mixture[particles/m/s]
        - T: float
            Temperature of the gas mixture [K]
        - g: float
            Acceleration of gravity [m/s2]       
        - m_i: float
            Mass of species i [kg/particle]
    Output:
        - Phi_diffusion : float
            Diffusion flux of species i [particles/m2/s]
    '''
    H_i = Scale_height_single_species(T, g, M_i)     # Scale height of light species [m]
    Phi_diffusion_i = b/H_i     # Diffusion flux of species i   [?]
    return Phi_diffusion_i

def Phi_critical(T, g, M_light, M_heavy, b, x_light, n_light, n_heavy):
    '''
    Compute the critical mass flux for a binary gas mixture undergoing atmospehric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - T: float
            Temperature of the gas mixture [K]
        - g: float
            Acceleration of gravity [m/s2]
        - M_light/M_heavy: float
            Molar mass of light/heavy species [kg/mol]
        - b: float
            Binary diffusion coefficient for the considered binary mixture[particles/m/s]
        - x_light: float
            Molar fraction of light species [1]
        - n_light/n_heavy: float
            Number of moles of light/heavy species [mol]

    Output: 
        - Phi_crit : float 
            Critical mass flux [particles/m2/s] 
    '''
    H_light = Scale_height_single_species(T, g, M_light)     # Scale height of light species [m]
    m_light = Mass_species(n_light, M_light)                 # Mass of light species [kg]
    m_heavy = Mass_species(n_heavy, M_heavy)                 # Mass of heavy species [kg]
    Phi_crit = b*x_light*(m_heavy - m_light)/H_light         # Critical mass flux    [particles/m2/s]
    return Phi_crit 

def Number_flux(Phi, Phi_crit, n_light, n_heavy, M_light, M_heavy, X_light, X_heavy, Phi_diffusion_light, Phi_diffusion_heavy):
    '''
    Compute the number flux of light and heavy species in a binary gas mixture undergoing atmospheric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    Inputs:
        - Phi: float
            Number flux [particles/m2/s]
        - Phi_crit: float
            Critical mass flux [particles/m2/s]
        - n_light,n_heavy: float
            Total number of moles of species light/heavy [mol]
        - M_light/M_heavy: float
            Molar mass of light/heavy species [kg/mol]
        - X_light/X_heavy: float
            Molar concentration of light/heavy species [1]          
        - Phi_diffusion_light/Phi_diffusion_heavy: float
            Diffusion flux of light/heavy species [particles/m2/s]
    Output:
        - Phi_light: float
            Number flux of light species [particles/m2/s]   
        - Phi_heavy: float
            Number flux of heavy species [particles/m2/s]
    '''
    m_light = Mass_species(n_light, M_light)                 # Mass of light species [kg]
    m_heavy = Mass_species(n_heavy, M_heavy)                 # Mass of heavy species [kg]
    # Compute the mean particle mass of the escaping outflow
    mean_particle_mass = m_light*X_light + m_heavy*X_heavy  # Mean particle mass of the escaping outflow [kg/particle]  

    # Compute Phi_light
    if Phi < Phi_crit:
        Phi_light = Phi/m_light
    else:    
       Phi_light = (X_light*Phi + X_light*X_heavy*(m_heavy - m_light)*Phi_diffusion_heavy)/mean_particle_mass

    # Compute Phi_heavy
    if Phi < Phi_crit:
          Phi_heavy = 0
    else:
        Phi_heavy = (X_heavy*Phi + X_light*X_heavy*(m_light - m_heavy)*Phi_diffusion_light)/mean_particle_mass
   
    return Phi_light, Phi_heavy

def Mass_loss(A, Phi):
    '''
    Compute the mass loss of the atmosphere undergoing atmospheric escape.
    
    Inputs:
        - A: float
            Surface area of the planet [m2]
        - Phi: float
            Number flux [particles/m2/s]
            
    Output: 
        - Mass_loss : float 
            Mass loss [particle/s]
    '''
    Mass_loss = A * Phi  # Mass loss rate [kg/s]
    return Mass_loss

def Fractionation_binary_mixture(n_light, n_heavy, Mp, Rp, b, T, M_light, M_heavy, Phi):
    '''
    Compute the fractionation of light and heavy species in a binary gas mixture undergoing atmospheric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - n_light/n_heavy: float     
            Total number of moles of light/heavy species [mol]
        - Mp: float
            Mass of the planet [kg]
        - Rp: float
            Radius of the planet [m]
        - b: float
            Binary diffusion coefficient for the considered binary mixture [particles/m/s]
        - T: float
            Temperature of the gas mixture [K]
        - M_light/M_heavy: float
            Molar mass of light/heavy species [kg/mol]
        - Phi: float
            Number flux [particle/m2/s]
    
    Output:
        - Mass_loss_light: float  
            Mass loss of light species [kg]
        - Mass_loss_heavy: float
            Mass loss of heavy species [kg]
        - Mass_loss_total: float
            Total mass loss [kg]
    '''
    # Compute the molar concentration of light (1) and heavy (2) species
    X_light, X_heavy = Molar_concentration_binary_mixture(n_light, n_heavy)

    # Compute the acceleration of gravity
    g = Acceleration_of_gravity(Mp, Rp)

    # Compute the diffusion flux of light (1) and heavy (2) species
    Phi_diffusion_light = Diffusion_flux(b, T, g, M_light)
    Phi_diffusion_heavy = Diffusion_flux(b, T, g, M_heavy)

    # Compute the critical mass flux
    Phi_crit = Phi_critical(T, g, M_light, M_heavy, b, X_light, n_light, n_heavy)

    # Compute the number flux of light (1) and heavy (2) species    
    Phi_light, Phi_heavy = Number_flux(Phi, Phi_crit, n_light, n_heavy, M_light, M_heavy, X_light, X_heavy, Phi_diffusion_light,Phi_diffusion_heavy)
    Phi_total = Phi_light + Phi_heavy

    return Phi_diffusion_light, Phi_diffusion_heavy, Phi_crit, Phi_light, Phi_heavy, Phi_total

##################### Ternary mixture #####################