#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "random.h"
#include "random_setup.h"

using namespace std;

// Function to calculate the statistical uncertainty (sigma)
double calculateSigma(double average, double averageSquared, int blockCount) {
    if (blockCount <= 1) return 0;
    return sqrt(abs(averageSquared - average * average) / blockCount);
}

// Integrand function F(x) = π * cos(πx / 2) / 2
double integrandFunction(double x) {
    return M_PI * cos(M_PI * x / 2) / 2;
}

// Cumulative Distribution Function (CDF): CDF(x) = 2x - x^2
double cumulativeDistributionFunction(double x) {
    return 2 * x - x * x;
}

// Find the inverse CDF using the bisection method
double inverseCDF(double probability, double tolerance = 1e-8) {
    double lowerBound = 0.0;
    double upperBound = 1.0;
    double midpoint;

    while (upperBound - lowerBound > tolerance) {
        midpoint = (lowerBound + upperBound) / 2.0;
        if (cumulativeDistributionFunction(midpoint) < probability) {
            lowerBound = midpoint;
        } else {
            upperBound = midpoint;
        }
    }

    return midpoint;
}

//////////////////////////////////////////////////////
// File Handling Functions

// Create a file with a header for storing lattice or continuum walker data
int createFileHeader(const string& filename) {
    ofstream outFile(filename);
    if (outFile.is_open()) {
        outFile << "#STEP\t#AVG\tSIGMA" << endl;
        outFile.close();
        return 1;
    } else {
        cerr << "Problem creating " << filename << "." << endl;
        return -1;
    }
}

// Save averaged data to a file
int saveData(double average, double averageSquared, int blockCount, int step, const string& filename) {
    ofstream outFile(filename, ios::app);
    if (outFile.is_open()) {
        outFile << step << "\t" << average << "\t" << calculateSigma(average, averageSquared, blockCount) << endl;
        outFile.close();
        return 1;
    } else {
        cerr << "Problem writing to " << filename << "." << endl;
        return -1;
    }
}

//////////////////////////////////////////////////////
// Walker Movement Functions

// Determine a random direction (0, 1, or 2 corresponding to x, y, z)
int getRandomDirection(Random& rnd) {
    return int(3.0 * rnd.Rannyu());
}

// Determine a random step forward or backward
int getRandomStep(Random& rnd) {
    return (rnd.Rannyu() < 0.5) ? 1 : -1;
}

// Calculate the squared radius from the origin
float calculateRadiusSquared(const float* position) {
    float radiusSquared = 0.0f;
    for (int i = 0; i < 3; i++)
        radiusSquared += pow(position[i], 2);
    return radiusSquared;
}

// Update the position based on random angles theta and phi
void updatePosition(float* position, Random& rnd) {
    float theta = acos(1 - 2 * rnd.Rannyu());
    float phi = 2 * M_PI * rnd.Rannyu();
    position[0] += sin(theta) * cos(phi);
    position[1] += sin(theta) * sin(phi);
    position[2] += cos(theta);
}

//////////////////////////////////////////////////////
// Main Function

int main() {
    // Initialize random number generator
    Random rnd;
    rnd_setup(rnd);

    // Simulation parameters
    const int numberOfBlocks = 40;
    const int stepsPerBlock = 1000;
    double average = 0.0, averageSquared = 0.0, sigmaValue = 0.0;

    // =========================
    // Uniform Sampling
    // =========================
    ofstream uniformDataFile("../OUTPUT/uniform.txt");
    for (int block = 1; block <= numberOfBlocks; block++) {
        double blockAverage = 0.0;
        for (int step = 0; step < stepsPerBlock; step++) {
            blockAverage += integrandFunction(rnd.Rannyu());
        }
        blockAverage /= stepsPerBlock;
        average += blockAverage;
        averageSquared += blockAverage * blockAverage;
        uniformDataFile << average / block << " " << calculateSigma(average / block, averageSquared / block, block) << endl;
    }
    uniformDataFile.close();
    average /= numberOfBlocks;
    sigmaValue = calculateSigma(average, averageSquared / numberOfBlocks, numberOfBlocks);

    cout << "== Integral with uniform sampling: I = " << average << " ± " << sigmaValue << endl;

    // Reset averages for importance sampling
    average = 0.0;
    averageSquared = 0.0;

    // =========================
    // Importance Sampling
    // =========================
    ofstream importanceDataFile("../OUTPUT/importance.txt");
    for (int block = 1; block <= numberOfBlocks; block++) {
        double blockAverage = 0.0;
        for (int step = 0; step < stepsPerBlock; step++) {
            double x = inverseCDF(rnd.Rannyu());
            blockAverage += integrandFunction(x) / (2 - 2 * x);
        }
        blockAverage /= stepsPerBlock;
        average += blockAverage;
        averageSquared += blockAverage * blockAverage;
        importanceDataFile << average / block << " " << calculateSigma(average / block, averageSquared / block, block) << endl;
    }
    importanceDataFile.close();
    average /= numberOfBlocks;
    sigmaValue = calculateSigma(average, averageSquared / numberOfBlocks, numberOfBlocks);

    cout << "== Integral with importance sampling: I = " << average << " ± " << sigmaValue << endl;

    // =========================
    // Lattice Walker Simulation
    // =========================
    string latticeFile = "../OUTPUT/latticeRW.txt";
    createFileHeader(latticeFile);

    const int latticeBlockSize = 100;
    const int latticeBlocks = 20;
    const int walkSteps = 100; // Number of steps the walker takes

    double averageRadiusLattice[walkSteps][latticeBlocks] = {0};
    float positionLattice[3] = {0.0f, 0.0f, 0.0f};

    cout << "Lattice walker simulation started... ";
    for (int block = 0; block < latticeBlocks; block++) {
        for (int bstep = 0; bstep < latticeBlockSize; bstep++) {
            // Reset position to origin for each walk
            positionLattice[0] = positionLattice[1] = positionLattice[2] = 0.0f;
            for (int step = 0; step < walkSteps; step++) {
                // Move in a random direction with a random step
                int direction = getRandomDirection(rnd);
                int stepDirection = getRandomStep(rnd);
                positionLattice[direction] += stepDirection;
                // Accumulate squared radius
                averageRadiusLattice[step][block] += calculateRadiusSquared(positionLattice);
            }
        }
        // Compute average radius for each step in the current block
        for (int step = 0; step < walkSteps; step++) {
            averageRadiusLattice[step][block] = sqrt(averageRadiusLattice[step][block] / latticeBlockSize);
        }
    }

    // Compute cumulative averages and save to file
    for (int step = 0; step < walkSteps; step++) {
        double cumulativeAverage = 0.0;
        double cumulativeAverageSquared = 0.0;
        for (int block = 0; block < latticeBlocks; block++) {
            cumulativeAverage += averageRadiusLattice[step][block];
            cumulativeAverageSquared += pow(averageRadiusLattice[step][block], 2);
        }
        saveData(cumulativeAverage / latticeBlocks, cumulativeAverageSquared / latticeBlocks, latticeBlocks, step + 1, latticeFile);
    }
    cout << "Done." << endl;

    // =========================
    // Continuum Walker Simulation
    // =========================
    // Reset averageRadius array
    for (int step = 0; step < walkSteps; step++)
        for (int block = 0; block < latticeBlocks; block++)
            averageRadiusLattice[step][block] = 0.0;

    string continuumFile = "../OUTPUT/contRW.txt";
    createFileHeader(continuumFile);

    cout << "Continuum walker simulation started... ";
    for (int block = 0; block < latticeBlocks; block++) {
        for (int bstep = 0; bstep < latticeBlockSize; bstep++) {
            // Reset position to origin for each walk
            positionLattice[0] = positionLattice[1] = positionLattice[2] = 0.0f;
            for (int step = 0; step < walkSteps; step++) {
                // Update position based on random angles
                updatePosition(positionLattice, rnd);
                // Accumulate squared radius
                averageRadiusLattice[step][block] += calculateRadiusSquared(positionLattice);
            }
        }
        // Compute average radius for each step in the current block
        for (int step = 0; step < walkSteps; step++) {
            averageRadiusLattice[step][block] = sqrt(averageRadiusLattice[step][block] / latticeBlockSize);
        }
    }

    // Compute cumulative averages and save to file
    for (int step = 0; step < walkSteps; step++) {
        double cumulativeAverage = 0.0;
        double cumulativeAverageSquared = 0.0;
        for (int block = 0; block < latticeBlocks; block++) {
            cumulativeAverage += averageRadiusLattice[step][block];
            cumulativeAverageSquared += pow(averageRadiusLattice[step][block], 2);
        }
        saveData(cumulativeAverage / latticeBlocks, cumulativeAverageSquared / latticeBlocks, latticeBlocks, step + 1, continuumFile);
    }	
    cout << "Done." << endl;

    return 0;
}

