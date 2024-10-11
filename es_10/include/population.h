#ifndef __population__
#define __population__

#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <unordered_set>
#include <mpi.h>
#include "random.h"
#include "random_setup.h"
#include "Salesman.h"

using namespace std;
using namespace arma;

class population {
public:
	population(int numSalesmen, int numCities, int rank);
	~population();	
	Random rnd;
	Random migControl;
	int numSalesmen;
	int numCities;	
	field <Salesman> salesmen;
	void _sort();
	void selection(double _pow);
	void migrate(int size, int rank);
	void crossover(double crossoverAcceptance);
	void mutation(double mutationAcceptance);
	double _sigma_dist();
	void printUpdate(double progress);
	int olympics(int size, int rank);
	void bestHalf_clear(int rank);
	void bestHalf(int rank);
};

#endif
