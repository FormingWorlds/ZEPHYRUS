"""
!!! info "`collision.py`"
    Function to compute fractional atmospheric loss during giant impacts.<br>
    Author(s): Anna Grace Ulses
"""

from zephyrus.constants import G 
import numpy as np 

def mass_loss(v_c:float, M_i:float, M_T:float, rho_i:float, rho_t:float, R_i:float , R_t:float, b:float):
    '''
    Fractional atmospheric mass loss due to giant impacts given by:

    $X \approx 0.64 [(\frac{v_c}{v_{esc}})^2 (\frac{M_i}{M_{tot}})^{\frac{1}{2}} (\frac{\rho_i}{\rho_t})^{\frac{1}{2}} f_M(b)]^{0.65}$

    where subscripts (i) indicate the impactor, subscripts (t) are for the target, and (T) indicates a combined value
    Source: Kegerreis et al 2020b

    Parameters 
    ----------
    v_c : float 
        Collision velocity between impactor and target [m/s]
    M_ : float 
        Mass of impactor (i) or total (T) [kg]
    rho_ : float 
        Density of impactor (i) or target (t) [kg/m^3]
    R_ : float 
        Radius of impactor (i) or target (t) [m]
    b : impact parameter 
        Defined as np.sin(beta) where beta is the angle of impact between the impactor and the target. b is unitless, beta is [deg]
    
    Returns 
    -------
    mass_loss : float
        Fractional mass loss [0,1] of target body's atmosphere [unitless]

    References
    ----------
    1. Kegerreis J.A, Eke V.R, Catling D.C, Massey R.J, Teodoro L.F.A, Zahnle K.J (2020). 
       Atmospheric Erosion by Giant Impacts onto Terrestrial Planets: A Scaling Law for any Speed, Angle, Mass, and Density
    '''
    #interacting volume, here we assume it is the same as the interacting mass of the colliding pair 
    fm = 1/4*((R_t + R_i)**3/(R_t**3 + R_i**3)) * (1-b)**2 * (1 + 2*b)

    v_esc = np.sqrt((2*G * (M_t + M_i))/(R_t + R_i)) #escape velocity
    M_tot = M_i + M_t

    bracket = (v_c/v_esc)**2 * (M_i/M_tot)**(1/2) * (rho_i/rho_t)**(1/2) * fm 
    mass_loss = 0.64 * bracket**0.65 #fractional atmospheric mass loss of target body, empirical (Kerregeis et al 2020)
    return mass_loss