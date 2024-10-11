#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include "random_setup.h"

using namespace std;

// ========================================
// Function Definitions
// ========================================

// Function to calculate the exponential part of the wavefunction
double exp_psi(double x, double mu, double sigma) {
    return exp(-pow(x + mu, 2) / (2 * sigma * sigma));
}

// Function to calculate the wavefunction ψ(x) as the sum of two exponentials
double psi(double x, double mu, double sigma) {
    return exp_psi(x, mu, sigma) + exp_psi(x, -mu, sigma);
}

// Function to calculate the potential energy V(x)
double pot(double x) {
    return pow(x, 4) - 2.5 * x * x;
}

// Function to calculate the kinetic energy component using finite differences
double kin(double x, double mu, double sigma) {
    double dx = 0.1; // Step size for finite difference
    return (psi(x + dx, mu, sigma) - 2 * psi(x, mu, sigma) + psi(x - dx, mu, sigma)) / (dx * dx);
}

// Function to calculate the total energy E(x) = T + V
double energy(double x, double mu, double sigma) {
    double psi_val = psi(x, mu, sigma);
    if (psi_val == 0) return 0; // Prevent division by zero
    return -0.5 * kin(x, mu, sigma) / psi_val + pot(x);
}

// Metropolis algorithm for sampling x
bool metro(double &x_current, double mu, double sigma, Random &rnd, double max_step) {
    // Propose a new position 'move' uniformly in [x_current - max_step, x_current + max_step]
    double move = x_current + max_step - 2 * rnd.Rannyu() * max_step;

    // Calculate the ratio of ψ(move) to ψ(x_current)
    double psi_ratio = psi(move, mu, sigma) / psi(x_current, mu, sigma);

    // Acceptance probability is the square of the ratio
    double ratio = pow(psi_ratio, 2);

    // Decide whether to accept the move
    if (rnd.Rannyu() < ratio) {
        x_current = move; // Accept the move
        return true;
    }
    else {
        return false; // Reject the move
    }
}

// Function to calculate the statistical error using block averaging
double error(double av, double av2, int block) {
    if (block < 2)
        return 0;
    return sqrt(abs(pow(av / block, 2) - av2 / block) / block);
}

// ========================================
// Function to Compute State Energy with Block Averaging
// ========================================

/**
 * @brief Computes the state energy using block averaging and collects x samples for histogram.
 * 
 * @param mu Variational parameter mu.
 * @param sigma Variational parameter sigma.
 * @param rnd Random number generator object.
 * @param NBlocks Number of blocks for block averaging.
 * @param steps Number of Metropolis steps per block.
 * @param x_samples Vector to store sampled x values for histogram.
 * @return pair<double, double> Pair containing the average energy and its statistical error.
 */
pair<double, double> state_energy(double mu, double sigma, Random &rnd, int NBlocks, int steps, vector<double> &x_samples) {
    double globE = 0;    // Global sum of energies
    double glob2E = 0;   // Global sum of squared energies
    double x = 0.0;      // Current position, initialized at 0 for better sampling
    int equi_steps = 1000; // Number of equilibration steps to stabilize the system
    double accepted = 0; // Counter for accepted moves during equilibration
    double max_step = 2.5; // Initial maximum step size for Metropolis moves

    // ========================================
    // Equilibration Phase
    // ========================================
    // Perform equilibration steps to reach a stable distribution
    for (int i = 0; i < equi_steps; i++) {
        accepted += metro(x, mu, sigma, rnd, max_step);
    }
    // Reset acceptance counter after equilibration
    accepted = 0;

    // ========================================
    // Block Averaging Phase
    // ========================================
    // Vector to store average energy per block (not used further, but kept for potential extensions)
    vector<double> block_energies(NBlocks, 0.0);

    // Loop over each block
    for (int b = 0; b < NBlocks; b++) {
        double ave = 0; // Sum of energies in the current block

        // Loop over each step within the block
        for (int i = 0; i < steps; i++) {
            // Perform a Metropolis move and count accepted moves
            accepted += metro(x, mu, sigma, rnd, max_step);

            // Calculate energy at the current position
            double e = energy(x, mu, sigma);
            ave += e; // Accumulate energy

            // Store the current x for histogram sampling
            x_samples.push_back(x);
        }

        // Adjust the step size based on acceptance rate to maintain efficiency
        double acceptance_rate = accepted / static_cast<double>(steps);
        if (acceptance_rate < 0.4)
            max_step *= 1.15; // Increase step size if acceptance rate is too low
        else if (acceptance_rate > 0.6)
            max_step /= 1.15; // Decrease step size if acceptance rate is too high

        // Reset acceptance counter for the next block
        accepted = 0;

        // Calculate average energy for the current block
        ave /= steps;
        block_energies[b] = ave;

        // Accumulate global sums for averaging and error calculation
        globE += ave;
        glob2E += ave * ave;
    }

    // Calculate overall average energy and its squared average
    double avg = globE / NBlocks;
    double avg2 = glob2E / NBlocks;

    // Calculate the statistical error
    double err = error(avg, avg2, NBlocks);

    return make_pair(avg, err); // Return average energy and error
}

// ========================================
// Function for Final Energy Sampling after Optimization
// ========================================

/**
 * @brief Performs final energy sampling with optimized parameters and writes results to a file.
 * 
 * @param rnd Random number generator object.
 * @param final_mu Optimized variational parameter mu.
 * @param final_sigma Optimized variational parameter sigma.
 * @param new_blocks Number of blocks for final energy sampling.
 * @param new_steps_per_block Number of Metropolis steps per block.
 * @param filename Output filename for final energy data.
 */
void final_energy_sampling(Random &rnd, double final_mu, double final_sigma, int new_blocks, int new_steps_per_block, const string &filename) {
    // Open the output file and write the header
    ofstream energy_final_file(filename);
    energy_final_file << "Block,Average_Energy,Energy_Error\n";

    double globE_final = 0.0;   // Global sum of energies
    double glob2E_final = 0.0;  // Global sum of squared energies
    double x = 0.0;              // Current position, initialized at 0

    // Loop over each block
    for (int b = 0; b < new_blocks; b++) {
        double aveE = 0.0; // Sum of energies in the current block

        // Loop over each step within the block
        for (int i = 0; i < new_steps_per_block; i++) {
            // Perform a Metropolis move with a fixed step size of 2.0
            metro(x, final_mu, final_sigma, rnd, 2.0);

            // Calculate energy at the current position
            double e = energy(x, final_mu, final_sigma);
            aveE += e; // Accumulate energy
        }

        // Calculate average energy for the current block
        aveE /= new_steps_per_block;

        // Accumulate global sums for averaging and error calculation
        globE_final += aveE;
        glob2E_final += aveE * aveE;

        // Calculate cumulative average and squared average up to the current block
        double avg_final = globE_final / (b + 1);
        double avg2_final = glob2E_final / (b + 1);

        // Calculate the statistical error
        double err_final = error(avg_final, avg2_final, b + 1);

        // Write the block data to the output file
        energy_final_file << (b + 1) << "," << avg_final << "," << err_final << "\n";
    }

    // Close the output file
    energy_final_file.close();
}

// ========================================
// Function to Sample ψ² and Store in a Histogram
// ========================================

/**
 * @brief Samples the square of the wavefunction ψ²(x) and stores the results in a histogram.
 * 
 * @param rnd Random number generator object.
 * @param final_mu Optimized variational parameter mu.
 * @param final_sigma Optimized variational parameter sigma.
 * @param num_bins Number of bins in the histogram.
 * @param x_min Minimum x value for the histogram.
 * @param x_max Maximum x value for the histogram.
 * @param bin_measures Number of samples to collect for the histogram.
 * @param filename Output filename for the histogram data.
 */
void sample_psi_squared_histogram(Random &rnd, double final_mu, double final_sigma, int num_bins, double x_min, double x_max, int bin_measures, const string &filename) {
    // Initialize histogram bins to zero
    vector<int> histogram(num_bins, 0);
    double bin_width = (x_max - x_min) / num_bins; // Width of each bin
    double x = 0.0; // Current position, initialized at 0

    // Sampling loop to fill the histogram
    for (int b = 0; b < bin_measures; b++) {
        // Perform a Metropolis move with a fixed step size of 2.0
        metro(x, final_mu, final_sigma, rnd, 2.0);

        // Determine which bin the current x falls into
        int bin = static_cast<int>((x - x_min) / bin_width);

        // Ensure the bin index is within the valid range
        if (bin >= 0 && bin < num_bins) {
            histogram[bin]++; // Increment the count for the corresponding bin
        }
    }

    // Open the histogram output file and write the header
    ofstream hist_file(filename);
    hist_file << "x_bin,Count\n";

    // Loop over each bin to calculate the center x value and write the count
    for (int i = 0; i < num_bins; i++) {
        double x_center = x_min + (i + 0.5) * bin_width; // Center of the bin
        hist_file << x_center << "," << histogram[i] << "\n";
    }

    // Close the histogram output file
    hist_file.close();
}

// ========================================
// Main Function
// ========================================

int main() {
    // ========================================
    // Initialization of Random Number Generator
    // ========================================
    Random rnd;         // Create a Random object
    rnd_setup(rnd);    // Initialize the Random object (assumed to be defined in "random_setup.h")

    // ========================================
    // File Outputs Setup
    // ========================================
    // Open the energy history file and write the header
    ofstream energy_file("../OUTPUT/energy_history.csv");
    energy_file << "SA_step,Energy,Energy_Error,Mu,Sigma\n";

    // Open the parameter trajectory file and write the header
    ofstream traj_file("../OUTPUT/parameter_trajectory.csv");
    traj_file << "SA_step,Mu,Sigma\n";

    // ========================================
    // Simulated Annealing Parameters
    // ========================================
    int blocks = 10;                // Number of blocks for energy calculation during SA
    int steps_per_block = 500;      // Number of Metropolis steps per block during SA
    double T_initial = 2.0;         // Initial temperature for simulated annealing
    double cooling_rate = 0.8;      // Cooling rate (exponential cooling factor)
    int total_SA_steps = 1000;      // Total number of simulated annealing steps
    double exploration_width = T_initial; // Initial exploration width based on temperature

    // ========================================
    // Variational Parameters Initialization
    // ========================================
    double mu = 1.0;         // Initial variational parameter mu
    double sigma = 1.0;      // Initial variational parameter sigma
    double currentEnergy = 0.0; // Current energy value

    // ========================================
    // History Tracking Vectors (Optional for Analysis)
    // ========================================
    vector<double> energy_history;   // Vector to store energy history
    vector<double> energy_errors;    // Vector to store energy errors
    vector<double> mu_history;       // Vector to store mu history
    vector<double> sigma_history;    // Vector to store sigma history

    // ========================================
    // Initial Energy Calculation Before SA
    // ========================================
    vector<double> x_samples_initial; // Vector to store initial x samples (optional)
    pair<double, double> initial = state_energy(mu, sigma, rnd, blocks, steps_per_block, x_samples_initial);
    currentEnergy = initial.first;    // Set the current energy
    double currentError = initial.second; // Set the current energy error

    // ========================================
    // Logging Initial Step
    // ========================================
    energy_history.push_back(currentEnergy);
    energy_errors.push_back(currentError);
    mu_history.push_back(mu);
    sigma_history.push_back(sigma);

    // Write initial energy and parameters to output files
    energy_file << "0," << currentEnergy << "," << currentError << "," << abs(mu) << "," << sigma << "\n";
    traj_file << "0," << abs(mu) << "," << sigma << "\n";

    // ========================================
    // Simulated Annealing Loop
    // ========================================
    double T = T_initial; // Current temperature

    // Loop over each simulated annealing step
    for (int sa = 1; sa <= total_SA_steps; sa++) {
        // Update exploration width based on current temperature
        exploration_width = sqrt(T);

        // ========================================
        // Parameter Proposal
        // ========================================
        // Propose a new mu by adding a random perturbation
        double new_mu = mu + 2 * exploration_width * (rnd.Rannyu() - 0.5);
        if (new_mu > 1.5) new_mu = mu; // Prevent mu from exceeding 1.5

        // Propose a new sigma by adding a random perturbation
        double new_sigma = sigma + 2 * exploration_width * (rnd.Rannyu() - 0.5);
        if (new_sigma <= 0.2 || new_sigma > 3) new_sigma = sigma; // Ensure sigma remains within physical bounds

        // ========================================
        // Energy Calculation for New Parameters
        // ========================================
        vector<double> x_samples_new; // Vector to store new x samples (optional)
        pair<double, double> new_energy_pair = state_energy(new_mu, new_sigma, rnd, blocks, steps_per_block, x_samples_new);
        double newEnergy = new_energy_pair.first;   // New energy
        double newError = new_energy_pair.second;   // New energy error

        // ========================================
        // Metropolis Criterion for Parameter Acceptance
        // ========================================
        double deltaE = newEnergy - currentEnergy; // Change in energy
        bool accept = false; // Flag to determine if new parameters are accepted

        if (deltaE < 0) {
            // If the new energy is lower, accept the new parameters
            accept = true;
        }
        else {
            // If the new energy is higher, accept with probability exp(-deltaE / T)
            double probability = exp(-deltaE / T);
            if (rnd.Rannyu() < probability)
                accept = true;
        }

        // ========================================
        // Update Parameters if Accepted
        // ========================================
        if (accept) {
            mu = new_mu;               // Update mu
            sigma = new_sigma;         // Update sigma
            currentEnergy = newEnergy; // Update current energy
            currentError = newError;   // Update current energy error
        }

        // ========================================
        // Record Data
        // ========================================
        energy_history.push_back(currentEnergy);
        energy_errors.push_back(currentError);
        mu_history.push_back(mu);
        sigma_history.push_back(sigma);

        // Write the current step's data to the energy and trajectory files
        energy_file << sa << "," << currentEnergy << "," << currentError 
                    << "," << abs(mu) << "," << sigma << endl;
        traj_file << sa << "," << abs(mu) << "," << sigma << endl;

        // ========================================
        // Temperature Update and Progress Logging
        // ========================================
        if (sa % 10 == 0) {
            // Apply cooling rate every 10 steps
            T *= cooling_rate;
        }

        if (sa % 100 == 0) {
            // Print progress every 100 steps
            cout << "SA Step: " << sa 
                 << " | Energy: " << currentEnergy
                 << " | Mu: " << mu 
                 << " | Sigma: " << sigma 
                 << " | Temperature: " << T << endl;
        }
    }

    // ========================================
    // Close Output Files after SA Loop
    // ========================================
    energy_file.close();
    traj_file.close();

    // ========================================
    // Final Output after Simulated Annealing
    // ========================================
    cout << "\nFinal Energy: " << currentEnergy 
         << " | Mu: " << abs(mu) 
         << " | Sigma: " << abs(sigma) << endl;

    // ========================================
    // Final Sampling after Optimization
    // ========================================
    // Define parameters for final energy sampling
    int final_blocks = 40;            // Number of blocks for final energy sampling
    int final_steps = 200;            // Number of Metropolis steps per block

    // Perform final energy sampling and write results to file
    final_energy_sampling(rnd, mu, sigma, final_blocks, final_steps, "../OUTPUT/final_energy.csv");

    // ========================================
    // Histogram Sampling of ψ²(x)
    // ========================================
    const int num_bins = 400;           // Number of histogram bins
    const double x_min = -5.0;          // Minimum x value for histogram
    const double x_max = 5.0;           // Maximum x value for histogram
    int bin_measures = num_bins * 5000; // Total number of samples for histogram

    // Perform histogram sampling and write results to file
    sample_psi_squared_histogram(rnd, mu, sigma, num_bins, x_min, x_max, bin_measures, "../OUTPUT/psi_squared_histogram.csv");

    return 0; // End of program
}

