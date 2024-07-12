''''
Emma Postolec
new.py
This file is a test to verify and validate escape computation.
'''

import numpy as np
import matplotlib.pyplot as plt

from escape import dMdt_EL_Lopez2012
from constants import *
from planets_parameters import *

# Initialization 
masses = [0.5,1.0,2.0,3.0,4.0,5.0]                # [Me]
MLR_today = []                                # Mass Loss Rate (MLR) [W.m-2]
MLR_10myr = []
epsilon = 0.15

# Computation escape for Earth-like planet at 
for i in masses : 
    MLR_today.append(dMdt_EL_Lopez2012('no',a_earth*au2m,e_earth,i*Me,Ms,epsilon,Re,Fxuv_earth_today))
    MLR_10myr.append(dMdt_EL_Lopez2012('no',a_earth*au2m,e_earth,i*Me,Ms,epsilon,Re,Fxuv_earth_10Myr))

MLR_10myr_Meyr = [(x / s2yr) / Me for x in MLR_10myr]
MLR_today_Meyr = [(x / s2yr) / Me for x in MLR_today]

# Save the data to a text file
np.savetxt('MLR_values_10Myr_and_today_my_computation.txt', np.column_stack((masses, MLR_10myr, MLR_10myr_Meyr, MLR_today, MLR_today_Meyr)), 
           header='masses   MLR_10myr [kg/s]    MLR_10myr [Me/yr]    MLR_today [kg/s]   MLR_today [Me/yr]')


# Create the plot
fig, ax1 = plt.subplots()

# Plotting the data on the primary y-axis
ax1.plot(masses, MLR_10myr, 'o', color='blue', label='t=0.01 Gyr')
ax1.plot(masses, MLR_10myr, color='blue')
ax1.plot(masses, MLR_today, 'o', color='orange', label='t=4.568 Gyr')
ax1.plot(masses, MLR_today, color='orange')
# Adding labels and title for primary y-axis
ax1.set_xlabel('Mp [Me]')
ax1.set_ylabel(r'MLR [kg $s^{-1}$]')
ax1.set_title('EL escape for Sun-Earth system')
# Adding grid and legend
ax1.grid(alpha=0.4)
ax1.legend()
# Setting logarithmic scale for primary y-axis
ax1.set_yscale('log')
# Adding a second y-axis with the converted scale and labels
ax2 = ax1.twinx()
# Apply the conversion factor to the y-axis values
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0] / s2yr) / Me,(ylims[1] / s2yr) / Me)
print(ylims)
ax2.set_yscale('log')
ax2.set_ylabel(r'MLR [Me $yr^{-1}$]')

# Saving the plot
plt.savefig('plots/test_new/MLR_vs_Me_earth.png', dpi=180)

# Show the plot
plt.show()