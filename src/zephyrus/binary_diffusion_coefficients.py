''''
Emma Postolec, Collin Cherubim
binary_diffusion_coefficients.py
Functions to get binary diffusion coefficients for different gas mixtures.
''''

def b_H_O(T):
    '''
    Binary diffusion coefficient for H in O [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HO : float 
            Binary diffusion coefficient for H in O [molecules/m/s] from Wordsworth et al. 2018
    '''
    b_HO = 4.8e19*T**0.75
    return b_HO

def b_H_H2O(T):
    '''
    Binary diffusion coefficient for H in H2O [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HH2O : float 
            Binary diffusion coefficient for H in H2O [molecules/m/s] from Wordsworth et al. 2018
    '''
    b_HH2O = 6.6e19*T**0.7
    return b_HH2O

def b_H_CO2(T):
    '''
    Binary diffusion coefficient for H in CO2 [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HCO2 : float 
            Binary diffusion coefficient for H in CO2 [molecules/m/s] from Odert et al. 2018 (Zahnle and Kasting, 1986)
    '''
    b_HCO2 = 8.4e19*T**0.6
    return b_HCO2

def b_H_O2(T):
    '''
    Binary diffusion coefficient for H in O2 [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HO2 : float 
            Binary diffusion coefficient for H in O2 [molecules/m/s] from Zahnle and Kasting, 2023 
    '''
    b_HO2 = 6.5e21*(T/1000)**0.7
    return b_HO2

def b_H_N2(T):
    '''
    Binary diffusion coefficient for H in N2 [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HN2 : float 
            Binary diffusion coefficient for H in N2 [molecules/m/s] from Zahnle and Kasting, 2023 
    '''
    b_HN2 = 6.5e21*(T/1000)**0.75
    return b_HN2