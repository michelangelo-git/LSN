#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

double sigma(double integral, double integral2, int n); // computes std dev for data blocking

//1
void _integral(string integrand, Random& rnd); // computes the integral of two functions, which one is decided by string integrand
void chi_test(Random& rnd); //performs the chi test on the given random number generator

//2
double which_dist(char dist, Random& rnd); // returns a random variable from the distribution that matches the char dist
void generate_dist(string distribution, Random& rnd, int* N_values, int dice_throws); // generates histograms for different distributions

//3
bool throw_needle(Random& rnd, double d, double L); // simulates a buffon needle throw, sampling the angle with a rejection technique
double pi_estimate(double probability, double d, double L); // converts the probability found in a block into a pi estimate
void buffon(int pi_throws, int n_points, Random& rnd); // simulates the buffon experiment



