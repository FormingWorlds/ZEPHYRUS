import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# To import functions from mors
mors_dir = os.path.dirname('/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/MORS-master/src/mors')
sys.path.append(mors_dir)

import mors

# List of different masses
masses = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
colors = ['brown', 'red', 'orange', 'gold', 'yellow', 'palegoldenrod']
labels = [f'm={mass} M$_{{\\odot}}$' for mass in masses]

# Omega parameter
omega = 1.0

# Define the save directory and output file name
save_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/'
plot_save_path = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/plots/test_mors/'

# Function to extract Lxuv luminosity from the dictionary
def get_lxuv(mstar, age, omega):
    result = mors.Lxuv(Mstar=mstar, Age=age, Omega=omega)
    print(result)
    return result['Lxuv']

# Function to create and save the data to a text file
def create_and_save_data(filepath, mstar, omega):
    # Generate a range of ages from 10 Myr to 4900 Myr with a spacing of 100 Myr
    ages = np.logspace(np.log10(1e6),np.log10(5e9), 1000)/1e6
    print(ages)
    
    # Calculate the Lxuv flux for each age
    xuv_luminosity = []
    for age in ages:
        lxuv = get_lxuv(mstar, age, omega)
        print(f"Mass: {mstar}, Age: {age}, Lxuv: {lxuv}")  # Debugging info
        xuv_luminosity.append(lxuv)
    
    # Save the ages and their corresponding Lxuv fluxes to a txt file
    np.savetxt(filepath, np.column_stack((ages, xuv_luminosity)), header="Age[Myr]   Lxuv[erg.s-1]")
    print(f"Data saved to {filepath}")

# Function to read data from a text file
def read_data(filepath):
    data = np.loadtxt(filepath, skiprows=1)
    ages = data[:, 0]
    xuv_luminosity = data[:, 1]
    return ages, xuv_luminosity

# Function to get filename based on parameters
def get_filename(mass, omega):
    return f'Lxuv_vs_age_Star_{mass}Msun_Omega_{omega}.txt'

plt.figure(figsize=(10, 8))

# Loop over each mass to create/read data and plot
for mass, color, label in zip(masses, colors, labels):
    filename = get_filename(mass, omega)
    filepath = os.path.join(save_dir, filename)
    
    if not os.path.isfile(filepath):
        # If the file does not exist, create and save the data
        create_and_save_data(filepath, mass, omega)
    else:
        print(f"File {filepath} already exists. Reading data from the file.")
    
    # Read the data from the file
    ages, xuv_luminosity = read_data(filepath)
    
    # Plot the data
    plt.loglog(ages, xuv_luminosity, color=color, label=label)

# Add labels, legend, and grid
plt.ylabel(r'XUV luminosity L$_{{\mathrm{XUV}}}$ [erg $s^{-1}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title(f'Stellar Evolution Tracks from MORS (Omega = {omega})', fontsize=15)

# Annotate specific points and lines (like the ZAMS line)
plt.axvline(x=4603, color='black', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e29, 'Today', rotation=90, verticalalignment='bottom')

# Save and show the plot
plt.savefig(plot_save_path+'MORS_function_Lxuv_vs_age_masses.png', dpi=180)
plt.show()
print(f"Plot saved to {plot_save_path}")
