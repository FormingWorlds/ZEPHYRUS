import numpy as np

import mors
from constants import *
from escape import *

initial_incident_flux = 14.67                                       # Fxuv received on Earth at t = 0.01 Gyr = 10e6 years = 10 Myr -> see Fig 9. Wordsworth+18 [W.m-2]
present_day_earth_XUV_flux = 4.64e-3                                # Fxuv received on Earth today [W.m-2] --> value from Luger&Barnes+2015 (from Ribas+2005)
age_earth_year = 4.543e9                                            # Age of the Earth [years]
eps = 0.15                                                          # efficiency factor for EL escape


mlr_hand_10_Myr_kg_s = (eps * np.pi * (Re**3) * initial_incident_flux) / (G * Me * 1)
mlr_hand_10_Myr_kg_yr = mlr_hand_10_Myr_kg_s * (s2yr)
mlr_hand_10_Myr_Me_yr = mlr_hand_10_Myr_kg_yr / Me

mlr_hand_today_kg_s = (eps * np.pi * (Re**3) * present_day_earth_XUV_flux) / (G * Me * 1)
mlr_hand_today_kg_yr = mlr_hand_today_kg_s * (s2yr)
mlr_hand_today_Me_yr = mlr_hand_today_kg_yr / Me

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Computation by hand for Fxuv from literature : ')
print('MLR on Earth at t=10 Myr  : ',mlr_hand_10_Myr_kg_s, 'kg/s', '|', mlr_hand_10_Myr_kg_yr, 'kg/yr', '|', mlr_hand_10_Myr_Me_yr, 'Me/yr')
print('MLR on Earth at today     : ',mlr_hand_today_kg_s, 'kg/s', '|', mlr_hand_today_kg_yr, 'kg/yr', '|', mlr_hand_today_Me_yr, 'Me/yr')

MLR_earth_10_Myr_kg_s = dMdt_EL_Lopez2012('no', a_earth * au2m, e_earth, Me, Ms, eps,Re, initial_incident_flux)
MLR_earth_10_Myr_kg_yr = MLR_earth_10_Myr_kg_s * (s2yr)
MLR_earth_10_Myr_Me_yr = MLR_earth_10_Myr_kg_yr / Me

MLR_earth_today_kg_s = dMdt_EL_Lopez2012('no', a_earth * au2m, e_earth, Me, Ms, eps,Re, present_day_earth_XUV_flux)
MLR_earth_today_kg_yr = MLR_earth_today_kg_s * (s2yr)
MLR_earth_today_Me_yr = MLR_earth_today_kg_yr / Me

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Computation with function for Fxuv from literature : ')
print('MLR on Earth at t=10 Myr  : ',MLR_earth_10_Myr_kg_s, 'kg/s', '|', MLR_earth_10_Myr_kg_yr, 'kg/yr', '|', MLR_earth_10_Myr_Me_yr, 'Me/yr')
print('MLR on Earth today        : ',MLR_earth_today_kg_s, 'kg/s', '|', MLR_earth_today_kg_yr, 'kg/yr', '|', MLR_earth_today_Me_yr, 'Me/yr')
print('--------------------------------------------------------------------------------------------------------------------------------------')

Fxuv_mors_earth_today = 2e-2 
MLR_Mors_earth_today_kg_s = dMdt_EL_Lopez2012('no', a_earth * au2m, e_earth, Me, Ms, eps,Re, Fxuv_mors_earth_today)
MLR_Mors_earth_today_Me_yr = MLR_Mors_earth_today_kg_s * (s2yr/Me)

print('--------------------------------------------------------------------------------------------------------------------------------------')
print('Computation with Fxuv from Mors (by eye) : ')
print('MLR on Earth today        : ',MLR_Mors_earth_today_kg_s, 'kg/s', '|', MLR_Mors_earth_today_Me_yr, 'Me/yr')
print('--------------------------------------------------------------------------------------------------------------------------------------')


Lxuv = mors.Lxuv(Mstar=1.0, Age=age_earth_year, Omega=10.0)
print(Lxuv)