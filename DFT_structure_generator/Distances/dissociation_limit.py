##########################################################################
#                          DISSOCIATION LIMIT                            #
#      DISTANCE SCANNER FOR MOLECULAR FRAGMENTS IN DFT CALCULATIONS      #
#                          ------------------                            #
#                          v.1.0.1 / 01.07.23                            #
#                          ETTORE BARTALUCCI                             #
##########################################################################


# SECTION X: EXTRACT AND COMPUTE NUMERICAL DISSOCIATION LIMIT FOR DISTANCE SCANNING LOOP

import numpy as np
from scipy.optimize import curve_fit

def extrapolate_dissociation_energy(distances, energies):
    def exponential_decay(x, a, b, c):
        return a * np.exp(-b * x) + c

    # Fit the energy data to an exponential decay function
    popt, _ = curve_fit(exponential_decay, distances, energies)

    # Extrapolate the energy to large distances
    extrapolated_energy = exponential_decay(np.inf, *popt)
    return extrapolated_energy

# Load energy data from a text file
data = np.loadtxt('energy_data.txt')
distances = data[:, 0]  # Column 1 contains distances
energies = data[:, 1]  # Column 2 contains energies

# Perform energy extrapolation
extrapolated_energy = extrapolate_dissociation_energy(distances, energies)

# Print the extrapolated dissociation energy
print("Extrapolated dissociation energy:", extrapolated_energy)

# Find the dissociation limit (distance where energy becomes negligibly small)
threshold_energy = 1e-6  # Define the threshold energy value
dissociation_limit_index = np.argmax(np.abs(energies) < threshold_energy)
dissociation_limit = distances[dissociation_limit_index]

# Print the dissociation limit
print("Dissociation limit:", dissociation_limit, "Å")
