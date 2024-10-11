import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the histogram data
file_path = '../OUTPUT/psi_squared_histogram.csv'
data = pd.read_csv(file_path)

# Extract data from the DataFrame
x_bins = data['x_bin']
counts = data['Count']

total_samples = counts.sum()
bin_width = x_bins[1] - x_bins[0]
counts /= (total_samples*bin_width)

# Plot the histogram
plt.figure(figsize=(10, 6))
plt.bar(x_bins, counts, width=(x_bins[1] - x_bins[0]), color='blue', alpha=0.7, label='Sampled $|\Psi_T(x)|^2$')


def Vpot(x):
    return (x**2 - 2.5)*x**2
    #return 0.5*x**2

hbar = 1
m = 1
a = 10
N = 1000 # number of iterations

# Step sizes
x = np.linspace(-a/2, a/2, N)
dx = x[1] - x[0] # the step size
V = Vpot(x)

# The central differences method: f" = (f_1 - 2*f_0 + f_-1)/dx^2

CDiff = np.diag(np.ones(N-1),-1)-2*np.diag(np.ones(N),0)+np.diag(np.ones(N-1),1)
# np.diag(np.array,k) construct a "diagonal" matrix using the np.array
# The default is k=0. Use k>0 for diagonals above the main diagonal, 
# and k<0 for diagonals below the main diagonal

# Hamiltonian matrix
H = (-(hbar**2)*CDiff)/(2*m*dx**2) + np.diag(V)

# Compute eigenvectors and their eigenvalues
E,psi = np.linalg.eigh(H)

# Take the transpose & normalize
psi = np.transpose(psi)
psi = psi/np.sqrt(dx)

print("Ground state energy: ", E[0])

# Plot a few things
plt.plot(x,(psi[0])**2)

# Add labels and title
plt.xlabel('x')
plt.ylabel('Count')
plt.title('Histogram of Sampled $|\Psi_T(x)|^2$')
plt.legend()

# Show grid
plt.grid()

# Save the plot as an image file
plt.savefig('psi_squared_histogram.png', dpi=300)

plt.xlim((-3,3))
plt.ylim((-0.6,0.6))
# Display the plot
plt.show()

