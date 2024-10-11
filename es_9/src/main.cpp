#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <unordered_set>
#include "random.h"
#include "random_setup.h"
#include "Salesman.h"
#include "population.h"

using namespace std;
using namespace arma;

int main(int argc, char *argv[]) {
	//setting up map
	Random rnd;
	rnd_setup(rnd);
	string output_name = "square";
	bool circle = 0;
	if (argc > 1 && strcmp(argv[1], "circle") == 0){
		circle = 1;
		output_name = "circle";
	}
	
	int numCities = 34;
    vec x(numCities), y(numCities);
    for (int i = 0; i < numCities; i++) {
        if (circle){
		     double angle = 2.* 3.141592653589793 * rnd.Rannyu();
		     x(i) = cos(angle);
		     y(i) = sin(angle); // Ensure y is within the circumference
        }
        else{
		     x(i) = rnd.Rannyu(); //cities in the square
		     y(i) = rnd.Rannyu();
        }
    }
    vector<city> italy(numCities);
    for (int i = 0; i < numCities; i++) {
        italy[i] = { x(i), y(i)};
    }
	// Setting up the army of capitalism
	const int numSalesmen = 50;
	//population constructor 
	population pop(numSalesmen, numCities);
	pop.bestHalf_clear(output_name);
	// GA loop
	int generation = 1;
	double selectionPow = 2.5;
	double crossoverAcceptance = 0.5;
	double mutationAcceptance = 0.05;
	int maxExe = 1.5*pow(10,3);

	while (generation <= maxExe) {
		// Update the distance for all salesmen
		for (int i = 0; i < numSalesmen; i++)
			 pop.salesmen(i).updateDistance(italy);
		
		pop._sort();
		pop.bestHalf(output_name);
		
		pop.selection(selectionPow);
		
		if(generation%1000 == 0){
			pop.printUpdate(double(generation)/maxExe);
			if(pop._sigma_dist() < 1)
				mutationAcceptance*=1.1;
			}
				
		pop.crossover(crossoverAcceptance);
		pop.mutation(mutationAcceptance);
		
		generation++;
   }
   pop._sort();
	pop.salesmen(0).save(italy, "../OUTPUT/best_" + output_name + ".dat");
	cout << "Best distance :" << pop.salesmen(0).distance << endl;

	return 0;
}

