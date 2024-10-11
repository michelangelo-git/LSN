#ifndef __population__
#define __population__

#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <unordered_set>
#include "random.h"
#include "random_setup.h"
#include "Salesman.h"

using namespace std;
using namespace arma;

class population {
public:
	population(int numSalesmen, int numCities);
	~population();	
	Random rnd;
	int numSalesmen;
	int numCities;	
	field <Salesman> salesmen;
	void _sort();
	void selection(double _pow);
	void crossover(double crossoverAcceptance);
	void mutation(double mutationAcceptance);
	double _sigma_dist();
	void printUpdate(double progress);
	void bestHalf_clear(string type);
	void bestHalf(string type);
};

#endif
