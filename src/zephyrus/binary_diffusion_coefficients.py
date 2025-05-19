''''
Emma Postolec, Collin Cherubim
binary_diffusion_coefficients.py
Functions to get binary diffusion coefficients for different gas mixtures.
'''

##################### H as a background gas ######################

# Fractionaiton with atomic species

def b_H_H(T):
    '''
    Binary diffusion coefficient for H in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output:     
        - b_HH : float 
            Binary diffusion coefficient for H in H [molecules/m/s] from Weissman and Masson 1962
    '''
    b_HH = 8.295e19*T**0.728
    return b_HH

def b_H_O(T):
    '''
    Binary diffusion coefficient for O in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HO : float 
            Binary diffusion coefficient for O in H [molecules/m/s] from Wordsworth et al. 2018
    '''
    b_HO = 4.8e19*T**0.75
    return b_HO

def b_H_C(T):
    '''
    Binary diffusion coefficient for C in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K] 
    Output:
        - b_HC : float 
            Binary diffusion coefficient for C in H [molecules/m/s] approximated from b_H_O using Genda, Ikoma 2008 formulation of Equation C3 in Appendix C
    '''
    b_HC = 4.85e19*T**0.75
    return b_HC

def b_H_N(T):
    '''
    Binary diffusion coefficient for N in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output:
        - b_HN : float 
            Binary diffusion coefficient for N in H [molecules/m/s] approximated from b_H_O using Genda, Ikoma 2008 formulation of Equation C3 in Appendix C
    '''
    b_HN = 4.85e19*T**0.75
    return b_HN

def b_H_S(T):   
    '''
    Binary diffusion coefficient for S in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output:     
        - b_HS : float      
            Binary diffusion coefficient for S in H [molecules/m/s] approximated from b_H_O using Genda, Ikoma 2008 formulation of Equation C3 in Appendix C
    '''
    b_HS = 4.73e19*T**0.75
    return b_HS

# Fractionation with molecular species

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
    b_HO2 = 6.5e21*(T/1000)**0.75
    return b_HO2

def b_H_H2 (T):
    '''
    Binary diffusion coefficient for H2 in H [molecules/m/s]

    Input: 
        - T : float
            Temperature [K] 
    Output:
        - b_HH2 : float 
            Binary diffusion coefficient for H2 in H [molecules/m/s] from Masson and Marrero 1970 (computed using Table VII and Eq. 136)
    '''
    b_HH2 = 8.29e19*T**0.728
    return b_HH2

def b_H_H2O(T):
    '''
    Binary diffusion coefficient for H in H2O [molecules/m/s]

    Input: 
        - T : float
            Temperature [K]
    Output: 
        - b_HH2O : float 
            Binary diffusion coefficient for H in H2O [molecules/m/s] from Wordsworth et al. 2018 (Zahnle and Kasting, 1986 / Mason and Weissman, 1972)
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


# To compute to take into account all species available in PROTEUS simulations 

def b_H_N2(T):
def b_H_S2(T):
def b_H_S02(T):
def b_H_H2S(T):
def b_H_NH3(T):
def b_H_CH4(T):
def b_H_CO(T):
