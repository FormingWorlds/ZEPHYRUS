''''
Emma Postolec, Collin Cherubim
fractionation.py
Functions to compute fractionation during atmospheric escape based on IsoFATE functions written by Collin Cherubim (Cherubim et al. 2024).
'''
import numpy as np
from zephyrus.constants import kb, G

def Molar_concentration(N_i, N_tot):
    '''
    Compute the molar concentration of specie i in a gas mixture.
    
    Inputs:
        - N_i: float
            Total number of moles of species i [kg/particle]
        - N_tot: float
            Total number of moles in the entire gas mixture [kg/particle]
            
    Output: 
        - x_i : float 
            Molar concentration of species i [1]
    '''
    x_i = N_i/(N_tot)
    return x_i

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

def Diffusion_flux_i(b, T, g, m_i):
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
            Diffusion flux of species i [?]
    '''
    Scale_height_i = kb * T / (g * m_i)     # Scale height of species i     [m]
    Phi_diffusion_i = b /Scale_height_i     # Diffusion flux of species i   [?]
    return Phi_diffusion_i

def Phi_critical(b, H1, m1, m2, x1):
    '''
    Compute the critical mass flux for a binary gas mixture undergoing atmospehric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - b: float
            Binary diffusion coefficient for the considered binary mixture[particles/m/s]
        - H1: float
            Scale height of light species [m]
        - m1/m2: float
            Molecular mass of light(1) and heavy(2) species [kg/particle]
        - x1: float
            Molar concentration of light/heavy species [1]
    
    Output: 
        - Phi_crit : float 
            Critical mass flux [?] 
    '''

    Phi_crit = b*x1*(m2 - m1)/H1    # Critical mass flux    [?]
    return Phi_crit 

def Number_flux(Phi, Phi_crit, m1, m2, x1, x2, Phi_diffusion_1, Phi_diffusion_2):
    '''
    Compute the number flux of light and heavy species in a binary gas mixture undergoing atmospheric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    Inputs:
        - Phi: float
            Mass flux [kg/m2/s]
        - Phi_crit: float
            Critical mass flux [kg/m2/s]
        - m1/m2: float
            Molecular mass of light(1) and heavy(2) species [kg/particle]
        - x1/x2: float
            Molar concentration of light/heavy species [1]          
        - Phi_diffusion_1/Phi_diffusion_2: float
            Diffusion flux of light/heavy species [kg/m2/s]
    Output:
        - Phi_light: float
            Number flux of light species [particles/m2/s]   
        - Phi_heavy: float
            Number flux of heavy species [particles/m2/s]
    '''
    # Compute the mean particle mass of the escaping outflow
    mean_particle_mass = m1*x1 + m2*x2   

    # Compute Phi_light
    if Phi < Phi_crit:
        Phi_light = Phi/m1
    else:    
       Phi_light = (x1*Phi + x1*x2*(m2 - m1)*Phi_diffusion_2)/mean_particle_mass

    # Compute Phi_heavy
    if Phi < Phi_crit:
          Phi_heavy = 0
    else:
        Phi_heavy = (x2*Phi + x1*x2*(m1 - m2)*Phi_diffusion_1)/mean_particle_mass
   
    return Phi_light, Phi_heavy

def Mass_loss_rate(A, Phi):
    '''
    Compute the mass loss rate of a gas undergoing atmospheric escape.
    
    Inputs:
        - A: float
            Surface area of the planet [m2]
        - Phi: float
            Mass flux [kg/m2/s]
            
    Output: 
        - MLR : float 
            Mass loss rate [kg/s]
    '''
    MLR = - A * Phi  # Mass loss rate [kg/s]
    return MLR

def Fractionation_binary(n1, n2, Mp, Rp, b, T, m1, m2, Phi):
    '''
    Compute the fractionation of light and heavy species in a binary gas mixture undergoing atmospheric escape.
    Adapted from IsoFATE (Cherubim et al.2024) and based on formulation in Cherubim et al.2024 and Wordsworth et al. 2018
    
    Inputs:
        - n1/n2: float     
            Total number of moles of light(1) and heavy(2) species [kg/particle]
        - Mp: float
            Mass of the planet [kg]
        - Rp: float
            Radius of the planet [m]
        - b: float
            Binary diffusion coefficient for the considered binary mixture[particles/m/s]
        - T: float      
            Temperature of the gas mixture [K]
        - g: float
            Gravitational field strength [m/s2]
        - m1/m2: float
            Molecular mass of light(1) and heavy(2) species [kg/particle]
        - Phi: float
            Mass flux [kg/m2/s]
    
    Output:
        - MLR_light: float  
            Mass loss rate of light species [kg/s]
        - MLR_heavy: float
            Mass loss rate of heavy species [kg/s]
        - MLR_total: float
            Total mass loss rate [kg/s]
    '''
    # Compute the molar concentration of light (1) and heavy (2) species
    N_tot = n1 + n2
    x1 = Molar_concentration(n1, N_tot)
    x2 = Molar_concentration(n2, N_tot)

    # Compute the acceleration of gravity
    g = Acceleration_of_gravity(Mp, Rp)

    # Compute the diffusion flux of light (1) and heavy (2) species
    Phi_diffusion_1 = Diffusion_flux_i(b, T, g, m1)
    Phi_diffusion_2 = Diffusion_flux_i(b, T, g, m2)

    # Compute the critical mass flux
    Phi_crit = Phi_critical(b, Phi_diffusion_1, m1, m2, x1)

    # Compute the number flux of light (1) and heavy (2) species    
    Phi_light, Phi_heavy = Number_flux(Phi, Phi_crit, m1, m2, x1, x2, Phi_diffusion_1, Phi_diffusion_2)

    # Compute the mass loss rate of light (1) and heavy (2) species and the total mass loss rate
    Planetary_surface_area = Planetary_surface_area(Rp)
    MLR_light = Mass_loss_rate(Planetary_surface_area, Phi_light)  # Mass loss rate of light species [kg/s]
    MLR_heavy = Mass_loss_rate(Planetary_surface_area, Phi_heavy)  # Mass loss rate of heavy species [kg/s]
    MLR_total = Mass_loss_rate(Planetary_surface_area, Phi)        # Total mass loss rate [kg/s]

    # Check the total mass loss rate
    if MLR_total == MLR_light + MLR_heavy:
        pass
    else:
        raise ValueError(f"Total mass loss rate does not equal the sum of light and heavy species mass loss rates: {MLR_total} != {MLR_light + MLR_heavy}")

    return MLR_light, MLR_heavy, MLR_total

##################### Ternary mixture #####################