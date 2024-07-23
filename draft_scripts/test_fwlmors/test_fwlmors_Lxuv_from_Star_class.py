import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import mors

# Directory to save the files
save_dir = '/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/data/Mors_stellar_evolution_tracks/'

# List of different masses
masses = [0.2, 0.4, 0.6, 0.8 ,1.0, 1.2]
#masses = [1.0, 1.2]
colors = ['brown', 'red', 'orange', 'gold', 'yellow', 'palegoldenrod']
labels = [f'm={mass} M$_{{\odot}}$' for mass in masses]

# Age and Omega parameters
omega = 1.0

# Function to get filename based on parameters
def get_filename(mass, omega):
    return f'Star_{mass}Msun_Myr_Omega_{omega}.pickle'

# Generate and save tracks if not already saved
for mass in masses:
    print(mass)
    filename = get_filename(mass, omega)
    filepath = os.path.join(save_dir, filename)
    
    if not os.path.exists(filepath):
        star = mors.Star(Mstar=mass, Omega=omega)
        star.Save(filepath)

plt.figure(figsize=(10, 8))

# Plot each mass track from saved files
# Plot each mass track from saved files
for mass, color, label in zip(masses, colors, labels):
    filename = get_filename(mass, omega)
    filepath = os.path.join(save_dir, filename)
       
    star = mors.Load(filepath)
    plt.loglog(star.Tracks['Age'], star.Tracks['Lx']+star.Tracks['Leuv'], color=color, label=label)

# Add labels, legend, and grid
plt.ylabel(r'XUV luminosity L$_{{\mathrm{XUV}}}$ [erg $s^{-1}$]', fontsize=15)
plt.xlabel('Age [Myr]', fontsize=15)
plt.legend()
plt.grid(alpha=0.5)
plt.title(f'Stellar Evolution Tracks from MORS (Omega = {omega})', fontsize=15)

# Annotate specific points and lines (like the ZAMS line)
#plt.axhline(y=1, color='black', linestyle='--', linewidth=0.7)
plt.axvline(x=4603, color='black', linestyle='--', linewidth=0.7)
plt.text(3600, 1.1e29, 'Today', rotation=90, verticalalignment='bottom')

plt.savefig(f'/Users/emmapostolec/Documents/PHD/SCIENCE/CODES/ZEPHYRUS/plots/test_mors/Mors_tracks_Lxuv_Omega_{omega}.png',dpi=180)
plt.show()
