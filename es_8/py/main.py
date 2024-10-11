import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    blocks = []
    avg = []
    sigma = []
    with open(filename, 'r') as f:
        next(f)  # skip header
        for line in f:
            cols = line.split()
            blocks.append(int(cols[0]))
            avg.append(float(cols[1]))
            sigma.append(float(cols[2]))
    return np.array(blocks), np.array(avg), np.array(sigma)

def plot_data(blocks, avg, sigma, title, ylabel, filename):
    plt.figure(figsize=(10, 6))
    plt.errorbar(blocks, avg, yerr=sigma, fmt='o-', label='Average with uncertainty')
    plt.xlabel('Block number')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename)
    plt.show()

# Read and plot call options (direct sampling)
blocks, c_avg, c_sigma = read_data('../OUTPUT/BS_call_direct.dat')
plot_data(blocks, c_avg, c_sigma, 'Call Option Price (Direct Sampling)', 'Option Price', 'call_direct.png')

# Read and plot put options (direct sampling)
blocks, p_avg, p_sigma = read_data('../OUTPUT/BS_put_direct.dat')
plot_data(blocks, p_avg, p_sigma, 'Put Option Price (Direct Sampling)', 'Option Price', 'put_direct.png')

# Read and plot call options (discretized sampling)
blocks, c_avg, c_sigma = read_data('../OUTPUT/BS_call_discrete.dat')
plot_data(blocks, c_avg, c_sigma, 'Call Option Price (Discretized Sampling)', 'Option Price', 'call_discrete.png')

# Read and plot put options (discretized sampling)
blocks, p_avg, p_sigma = read_data('../OUTPUT/BS_put_discrete.dat')
plot_data(blocks, p_avg, p_sigma, 'Put Option Price (Discretized Sampling)', 'Option Price', 'put_discrete.png')

