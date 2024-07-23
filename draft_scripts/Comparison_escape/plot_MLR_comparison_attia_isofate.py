import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
# To import functions from mors
mors_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/MORS-master/src/mors')
sys.path.append(mors_dir)
constants_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/constants.py')
sys.path.append(constants_dir)
planets_parameters_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/planets_parameters.py')
sys.path.append(planets_parameters_dir)
escape_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/escape.py')
sys.path.append(escape_dir)
from escape import *
from planets_parameters import *
import mors
from constants import *
xuv_flux_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/src/zephyrus/XUV_flux.py')
sys.path.append(xuv_flux_dir)
from XUV_flux import *
fxuv_functions_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/tools/Fxuv_functions.py')
sys.path.append(fxuv_functions_dir)
from Fxuv_functions import *

data_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/MLR_computations/'

data_Lxuv_mara = np.loadtxt(data_dir + 'MLR_computation_Attia+2021_with_Lxuv.txt', skiprows=1)
Time_Lxuv_mara = data_Lxuv_mara[:, 0] 
MLR_Lxuv_mara = data_Lxuv_mara[:, 1]  
data_Star_mara = np.loadtxt(data_dir + 'MLR_computation_Attia+2021_with_Star.txt', skiprows=1)
Time_Star_mara = data_Star_mara[:, 0] 
MLR_Star_mara = data_Star_mara[:, 1]  

data_Lxuv_zephyrus = np.loadtxt(data_dir + 'MLR_computation_Lopez+2012_with_Lxuv.txt', skiprows=1)
Time_Lxuv_zephyrus = data_Lxuv_zephyrus[:, 0] 
MLR_Lxuv_zephyrus = data_Lxuv_zephyrus[:, 1]  
data_Star_zephyrus = np.loadtxt(data_dir + 'MLR_computation_Lopez+2012_with_Star.txt', skiprows=1)
Time_Star_zephyrus = data_Star_zephyrus[:, 0] 
MLR_Star_zephyrus = data_Star_zephyrus[:, 1]  

data_Lxuv_isofate = np.loadtxt(data_dir + 'MLR_computation_Cherubim+2024_with_Lxuv.txt', skiprows=1)
Time_Lxuv_isofate = data_Lxuv_isofate[:, 0] 
MLR_Lxuv_isofate = data_Lxuv_isofate[:, 1]  
data_Star_isofate = np.loadtxt(data_dir + 'MLR_computation_Cherubim+2024_with_Star.txt', skiprows=1)
Time_Star_isofate = data_Star_isofate[:, 0] 
MLR_Star_isofate = data_Star_isofate[:, 1] 

MLR_earth_today_zephyrus = dMdt_EL_Lopez2012('no',a_earth*au2m,e_earth,Me,Ms,0.15,Re,Fxuv_earth_today)
MLR_earth_today_mara = dMdt_EL_Attia2021('no',a_earth*au2cm,e_earth,Me*1e3,Ms*1e3,Re*1e2,Re*1e2,((Fxuv_earth_today/ergcm2stoWm2)*4*np.pi*au2cm**2),Fxuv_earth_today/ergcm2stoWm2,0.15,'no')
Vpot = (G*Me)/(Re**3)
MLR_earth_today_isofate = 0.15*Fxuv_earth_today/(4*Vpot)

# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 8))

# Plotting the data on the primary y-axis
ax1.loglog(Time_Lxuv_mara, MLR_Lxuv_mara, ':', color='black', label='Flux from mors.Lxuv()')
ax1.loglog(Time_Star_mara, MLR_Star_mara, color='black', label='Flux from mors.Star()',)
ax1.axvline(x=age_earth/1e6, color='lightgrey', linestyle='--', linewidth=1)
ax1.loglog(age_earth/1e6, MLR_earth_today_zephyrus,'o', color='black', label='Earth today (t = 4.543 Gyr)')
ax1.loglog(age_earth/1e6, MLR_earth_today_mara/1e3,'o', color='red')
ax1.loglog(age_earth/1e6, MLR_earth_today_zephyrus,'o', color='steelblue')
ax1.loglog(age_earth/1e6, MLR_earth_today_isofate,'o', color='gold')
#ax1.loglog(age_earth/1e6, 1423.221809,'o', color='green')
ax1.loglog(Time_Lxuv_mara, MLR_Lxuv_mara, ':', color='red')
ax1.loglog(Time_Lxuv_zephyrus, MLR_Lxuv_zephyrus,':', color='steelblue')
ax1.loglog(Time_Lxuv_isofate, MLR_Lxuv_isofate,':', color='gold')
ax1.loglog(Time_Star_mara, MLR_Star_mara, color='red', label='Mara : Attia+2021')
ax1.loglog(Time_Star_zephyrus, MLR_Star_zephyrus, color='steelblue', label='Zephyrus : Lopez+2012')
ax1.loglog(Time_Star_isofate, MLR_Star_isofate, color='gold', label='IsoFATE : Cherubim+2024')



# Adding labels and title for primary y-axis
ax1.set_xlabel('Time [Myr]')
ax1.set_ylabel(r'MLR [kg $s^{-1}$]')
ax1.set_title('Comparison for EL escape for Sun-Earth system (Zephyrus vs IsoFATE vs Attia+2021)')
# Adding grid and legend
ax1.grid(alpha=0.4)
ax1.legend(ncol=2, loc='upper right')
# Setting logarithmic scale for primary y-axis
ax1.set_yscale('log')
# Adding a second y-axis with the converted scale and labels
ax2 = ax1.twinx()
# Apply the conversion factor to the y-axis values
ylims = ax1.get_ylim()
ax2.set_ylim((ylims[0] / s2yr) / Me,(ylims[1] / s2yr) / Me)
ax2.set_yscale('log')
ax2.set_ylabel(r'MLR [$M_{\oplus}$ $yr^{-1}$]')

# Add a text box with LaTeX writing
textstr = r'$\epsilon$ = 0.15' '\n' r'$R_p = R_{\mathrm{XUV}} = R_{\oplus}$' '\n' r'$M_p = M_{\oplus}$'
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
# Place the text box in the upper left corner of the plot
ax1.text(1.2, 6e2, textstr, fontsize=14,
         verticalalignment='top', bbox=props)

# Saving the plot
plot_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/plots/test_MLR_comparison/'
plt.savefig(plot_dir+'comparison_MLR_zephyrus_vs_attia2021.png', dpi=180)

# Show the plot
plt.show()
