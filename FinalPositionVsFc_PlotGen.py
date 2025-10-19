# Sami Alawadhi, David M. Rutkowski (2025) 
# Vavylonis Group
# Department of Physics, Lehigh University
# If you use any part of the scripts in this package, please cite:
# TBD

import numpy as np 
from matplotlib import pyplot as plt

# Final displacements were obtained manually from individual displacement curves and gathered in an array corresponding to each contractile force in the simulations.

# Data: Forces and displacements in two conditions - to get displacement vs. Fc
forces = np.arange(0.208, 0.509, 0.025)  # Generate array of force values used in the simulations.
septin_mediated = np.array([
    [30.29, 96.9, 29.7, 29.6, 96.8, 60, 30.77, 86.4, 89.1, 90.9, 29.74, 29.69],
    [87.2, 113.3, 116.2, 97.5, 110, 89.8, 98.4, 114.8, 109.9, 113.7, 69, 93.8],
    [104.6, 120.9, 120.5, 105.5, 123.7, 110.9, 116.1, 116.4, 100.7, 122.3, 98.9, 130.1],
    [128.2, 121.7, 135.9, 119.2, 120.9, 111.7, 125.4, 125.7, 119.4, 121.8, 117.1, 130.1],
    [136.4, 131.3, 115.6, 114.2, 131.5, 121.9, 133.4, 126.6, 116.5, 137.1, 113.5, 126.5],
    [138.7, 150.8, 142.6, 136.5, 122.5, 121.1, 136.5, 147.3, 140.1, 141.3, 139.4, 134.9],
    [147.7, 146.4, 144, 139.6, 141.8, 141.5, 132.1, 149.6, 142.4, 148.9, 145.4, 136.3],
    [139.8, 144.6, 149, 142.1, 134.8, 153.8, 141.2, 142.5, 155.8, 149.9, 150.6, 146],
    [144.2, 147.1, 148.3, 143, 139.4, 141.9, 140.6, 139.7, 141, 138, 148.9, 146.8],
    [155.4, 157.5, 157.8, 153.9, 131, 109.9, 150.2, 164.8, 147.9, 161.1, 156.5, 163.2],
    [162.8, 161.2, 162.4, 121.6, 109.4, 162.9, 165.4, 115.4, 151.6, 161, 167.9, 147.8],
    [174, 161.7, 157.4, 162, 132, 167, 171, 157.5, 141.8, 110.2, 169.5, 158.7],
    [151.7, 170.5, 153.6, 171.5, 130.8, 171.9, 170.3, 163, 153.7, 126.3, 153.7, 124],
]) # measurement of CB cell max displacement.
septin_inhibited = np.array([
    [30.57, 113, 109.2, 120.7, 106.3, 31, 30.3, 116.6, 81.9, 99.5, 30.66, 30.32],
    [99, 104.8, 111, 119, 131, 105.5, 118.9, 133.9, 111.3, 94.6, 114.2, 93.8],
    [131.6, 130.2, 140.2, 110.7, 134, 105.3, 124.6, 133.9, 88.1, 100, 110.4, 105.4],
    [138.1, 137.6, 116.1, 111.1, 99.9, 105, 116.9, 134.5, 82, 126.7, 139, 105.2],
    [138.9, 152.6, 105.3, 106, 105.8, 93.9, 140.5, 146.1, 81.9, 82.1, 129.3, 105.8],
    [93.9, 131.1, 144.3, 134.9, 134.1, 94.1, 87.9, 132.4, 88.5, 88.4, 139.4, 105.7],
    [154.1, 137.7, 116.6, 111.4, 82.2, 93.9, 116.4, 158.1, 111.5, 99.8, 122.5, 94],
    [145.5, 105.6, 105.3, 150.2, 127.6, 154.7, 105.4, 140.4, 104.7, 82.2, 121.9, 105.6],
    [121.5, 122.6, 93.3, 144.3, 105.6, 111, 99.1, 150, 105.1, 83, 118.7, 145.5],
    [154.4, 121.8, 87.5, 116.8, 105.7, 121.9, 104.5, 155, 99.8, 152.5, 139.3, 139.5],
    [150.7, 110.9, 144.1, 121.6, 111.1, 110.6, 133, 127.9, 87.8, 111.4, 139, 127.6],
    [154.4, 145.8, 105.6, 119.1, 105.7, 88.1, 122.1, 158.9, 81.3, 111.2, 122.2, 133.8],
    [145.4, 138.1, 105.8, 106.4, 127.7, 132.6, 104.3, 145.5, 94.3, 111.8, 116.7, 145.2],
]) #measurement of noCB cell max displacement

# Define forces and equilibrium bond length
forces = np.arange(0.208, 0.509, 0.025)
equilibrium_bond_length = 0.273043478
forces_per_unit_length = forces / equilibrium_bond_length

# Function to calculate means and standard errors of the mean (SEM)
def compute_stats(data):
    means = np.mean(data, axis=1)
    sems = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])
    return means, sems

# Calculate means and SEMs
CB_means, CB_sems = compute_stats(septin_mediated)
noCB_means, noCB_sems = compute_stats(septin_inhibited)

# Function to plot error bars
def plot_error_bars(x, means, sems, color, label, marker):
    plt.errorbar(x, means, yerr=sems, fmt=marker, color=color, label=label, markersize=8)
    plt.plot(x, means, linestyle='-', color=color)

# Define color-blind friendly colors
colors = {
    "With compartmentalization": "#0006B2",  # Blue
    "Without compartmentalization": "#D50000"  # Red
}

# Plotting with SEM
plot_error_bars(forces_per_unit_length, CB_means, CB_sems, colors["With compartmentalization"], 'With compartmentalization', 'o')
plot_error_bars(forces_per_unit_length - 0.015, noCB_means, noCB_sems, colors["Without compartmentalization"], 'Without compartmentalization', 'o')

plt.gca().tick_params(axis='both', labelsize=18)

plt.xlabel('Contractile force per unit length (nN/\u03BCm)', fontsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.ylabel('Displacement (\u03BCm)', fontsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.legend(loc='lower right', fontsize=18)

plt.tight_layout()  # Adjusts subplot parameters; helps ensure that everything fits well within the display area.

# Set the aspect of the plot to be equal
plt.gca().set_aspect('auto')

plt.show()
