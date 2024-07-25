import numpy as np
import matplotlib.pyplot as plt

import sys
import os

# To import functions from mors
mors_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/MORS-master/src/mors')
sys.path.append(mors_dir)

import mors

age = 40
omega = 1.0
print('Sun-like star : Mstar=1.0, Age = ', age, '[Myr], Omega = ', omega)
star = mors.Star(Mstar=1.0, Age=age, Omega=omega)
print('Lxuv from Star() = ', star.Lx(age)+star.Leuv(age), '[erg s-1]')
Lxuv = mors.Lxuv(Mstar=1.0, Age=age, Omega=omega)
print('Lxuv from Lxuv() = ', Lxuv['Lxuv'], '[erg s-1]')

#star.Save('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/Sun_1Msun.pickle')
#star.PrintAvailableTracks() # print all the units of all quantities
#print(star.Units['Lx'])
#print(star.Lx(150.0))
#print(star.Tracks['Age'])

#star = mors.Load("/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/Sun_1Msun.pickle")

# plt.loglog(star.Tracks['Age'], star.Tracks['Lx'], color='orange',label='Sun')
# plt.ylabel(r'Lx [erg $s^{-1}$]')
# plt.xlabel('Age [Myr]')
# plt.legend()
# plt.grid(alpha=0.5)
# plt.show()

#Lxuv = mors.Lxuv(Mstar=1.0, Age=1.0, Omega=1.0)
#print(Lxuv)
