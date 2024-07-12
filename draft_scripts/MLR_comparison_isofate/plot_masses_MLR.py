import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# To import functions from src/zephyrus
main_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/')
src_zephyrus_folder_path = os.path.join(main_dir, 'src/zephyrus/')
sys.path.append(src_zephyrus_folder_path)

from constants import *
from planets_parameters import *

data_dir_MLR = '/data/MLR_computations/'
MLR_me = np.loadtxt(main_dir + data_dir_MLR + "MLR_values_10Myr_and_today_zephyrus.txt")
MLR_isofate = np.loadtxt(main_dir + data_dir_MLR +'MLR_values_10Myr_and_today_isofate.txt')

masses = MLR_me[:,0]
MLR_10myr_kgs_me = MLR_me[:,1]
MLR_10myr_Meyr_me = MLR_me[:,2]
MLR_today_kgs_me = MLR_me[:,3]
MLR_today_Meyr_me = MLR_me[:,4]

MLR_10myr_kgs_isofate = MLR_isofate[:,1]
MLR_10myr_Meyr_isofate = MLR_isofate[:,2]
MLR_today_kgs_isofate = MLR_isofate[:,3]
MLR_today_Meyr_isofate = MLR_isofate[:,4]


# Create the plot
fig, ax1 = plt.subplots()

# Plotting the data on the primary y-axis
ax1.plot(masses, MLR_10myr_kgs_me, 'o', color='blue', label=' Zephyrus, t = 0.01 Gyr')
ax1.plot(masses, MLR_10myr_kgs_me, color='blue')
ax1.plot(masses, MLR_10myr_kgs_isofate, 'o', color='green', label=' IsoFate, t = 0.01 Gyr')
ax1.plot(masses, MLR_10myr_kgs_isofate, color='green')
ax1.plot(masses, MLR_today_kgs_me, 'o', color='orange', label=' Zephyrus, t = 4.568 Gyr')
ax1.plot(masses, MLR_today_kgs_me, color='orange')
ax1.plot(masses, MLR_today_kgs_isofate, 'o', color='red', label=' IsoFate, t = 4.568 Gyr')
ax1.plot(masses, MLR_today_kgs_isofate, color='red')

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
plot_dir_MLR = '/plots/test_MLR_comparison/'
plt.savefig(main_dir + plot_dir_MLR + 'plot_compare_MLR_vs_Me_earth_isofate_vs_me.png', dpi=180)

# Show the plot
plt.show()