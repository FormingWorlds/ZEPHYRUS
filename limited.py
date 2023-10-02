#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Converts stellar XUV luminosity [J⋅s^-1] to a hydrodynamic luminosity of escape [kg⋅s^-1] with some efficiency eta. 
beta^3 dependence as in Owen R. Lehmer and David C. Catling 2017 ApJ 845 130,taking into account that material 
is lifted from XUV photosphere rather than surface.
"""
## Physical constants
G         = 6.674e-11 # [m^3⋅kg^−1⋅s^−2]
eV        = 1.602e-19 # [J]
## Example parameters (later to be explicitly calculated)
eta_cons  = 0.1       # evap efficiency, conservative estimate, Ly alpha cooling 
eta_UB    = 1         # 0.5 * 2 : delta(h*nu) * (Teq boost)
beta_cons = 1         # Rxuv =  beta * Rpl
beta      = 2         # reasonable beta: 2 - 10
## Example inputs 
Mearth    = 5.972e27  # [kg]
Lxuv      = 135.7e21  # [J⋅s^-1] XUV luminosity e.g. active solar mass star (Kubyshinka 2018)
semiaxis  = 1.496e11  # [m] 1 AU Earth semiaxis

import numpy as np

    def dMdtEL(Mpl, Rpl, Lxuv, semiaxis):
        '''Making a good energy-limited escape calculation is a game of finding the right efficiency. Here outputs are a conservative
        estimate and an upper bound estimate. The efficiencies will be explicitly calculated. E.g. the excess energy is a function of
        the stellar spectrum. Beta is a function of the scale height of the atmosphere and can deviate significantly from 1.
        '''
        cons = 0.25 * eta_cons * (Lxuv * Rpl^2 /semiaxis^2) / (G * Mpl /Rpl)
        UB   = (eta_UB/eta_cons) * beta^3 * cons
        return cons, UB
