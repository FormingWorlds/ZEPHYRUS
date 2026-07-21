#mass loss function during impacts

def mass_loss(v_c, M_i, M_T, rho_i, rho_t, R_i, R_t, b):
    '''
    v_c = collision velocity
    b = impact parameter 
    _i = impactor
    _t = target
    _T = total
    '''
    #interacting volume, here we assume it is the same as the interacting mass of the colliding pair 
    fm = 1/4*((R_t + R_i)**3/(R_t**3 + R_i**3)) * (1-b)**2 * (1 + 2*b)

    v_esc = np.sqrt((2*G * (M_T + M_i))/(R_t + R_i)) #escape velocity

    bracket = (v_c/v_esc)**2 * (M_i/M_T)**(1/2) * (rho_i/rho_t)**(1/2) * fm 
    mass_loss = 0.64 * bracket**0.65 #fractional atmospheric mass loss of target body, empirical (Kerregeis et al 2020)
    return mass_loss