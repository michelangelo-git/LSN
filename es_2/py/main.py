import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load the data from files
n_steps = 100
walk_blocks = 100
x, y, z= np.loadtxt("../LSN/es_2/OUTPUT/latticeRW.txt", unpack = True)
average_lattice = y
x, y, z = np.loadtxt("../LSN/es_2/OUTPUT/contRW.txt", unpack = True)
average_solid = y 
steps = np.arange(1, n_steps + 1)

# Define the fitting function
def fit_function(N, k):
    return k * np.sqrt(N)

# Perform the fit
params_lattice, _ = curve_fit(fit_function, steps, average_lattice)
params_solid, _ = curve_fit(fit_function, steps, average_solid)

# Extract the fitting parameter k
k_lattice = params_lattice[0]
k_solid = params_solid[0]

# Plot the results
plt.figure(figsize=(12, 6))

# Plot for lattice walker
plt.subplot(1, 2, 1)
plt.plot(steps, average_lattice, 'o', label='Lattice Walker Data')
plt.plot(steps, fit_function(steps, k_lattice), '-', label=f'Fit: k={k_lattice:.3f}')
plt.xlabel('Number of Steps (N)')
plt.ylabel('Average Distance')
plt.title('Lattice Walker')
plt.legend()

# Plot for solid angle walker
plt.subplot(1, 2, 2)
plt.plot(steps, average_solid, 'o', label='Solid Angle Walker Data')
plt.plot(steps, fit_function(steps, k_solid), '-', label=f'Fit: k={k_solid:.3f}')
plt.xlabel('Number of Steps (N)')
plt.ylabel('Average Distance')
plt.title('Solid Angle Walker')
plt.legend()

plt.tight_layout()
plt.show()

