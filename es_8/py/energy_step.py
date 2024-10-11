import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the histogram data
file_path = '../LSN/es_8/OUTPUT/energy_history.csv'
data = pd.read_csv(file_path)

# Extract data from the DataFrame
step = data['SA_step']
energy = data['Energy']
err = data['Energy_Error']

# Plot the histogram
plt.figure(figsize=(10, 6))
plt.error(step, energy,yerr=sigma, fmt='o-', color='blue', alpha=0.7, label='Energy-SAstep')



# Add labels and title
plt.xlabel('SA_step')
plt.ylabel('Energy')
plt.title('Energy')
plt.legend()

# Show grid
plt.grid()

#plt.xlim((-3,3))
#plt.ylim((-0.6,0.6))
# Display the plot
plt.show()

