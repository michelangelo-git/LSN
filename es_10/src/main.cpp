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
#include "population.h"

using namespace std;
using namespace arma;

int mapSize(string filename){
	int numCities = 0;
	ifstream cityData;
		 cityData.open(filename);
		 if (!cityData) {
			  cout << "Problems reading " << filename << endl;
			  return -1;
		 }
	 string line;
	 while (getline(cityData, line)) {
	     numCities++;
	 }
	 cityData.close();
	 return numCities;
}

int initMap(vector<city>& map, int numCities, string filename){
	ifstream cityData;
	cityData.open(filename);
	if (!cityData) {
		cout << "Problems reading " << filename << endl;
		return 0;
	}
	vec x(numCities), y(numCities); 
	for (int i = 0; i < numCities; i++) {
		cityData >> x(i) >> y(i);
		map[i] = {x(i), y(i)};
	}
	cityData.close();
	return 1;
}

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Open file to count the number of cities
	int numCities = mapSize("../INPUT/cap_prov_ita.dat");
	if(numCities < 1){
		MPI_Finalize();
		return -1;
	}

	// Open the file to read city coordinates
	vector<city> italy(numCities);
	initMap(italy, numCities, "../INPUT/cap_prov_ita.dat");

	// Setting up the army of capitalism
	const int numSalesmen = 300;
	//population constructor 
	population pop(numSalesmen, numCities, rank);
	pop.bestHalf_clear(rank);
	// GA loop
	int generation = 1;
	double selectionPow = 3.2;
	double crossoverAcceptance = 0.4;
	double mutationAcceptance = 0.1;
	int maxExe = pow(10,4);
	int migrationTime = 1000;

	while (generation <= maxExe) {
		// Update the distance for all salesmen
		for (int i = 0; i < numSalesmen; i++)
			 pop.salesmen(i).updateDistance(italy);

		// Sort the salesmen based on their distance
		pop._sort();
		pop.bestHalf(rank);
		// Migration
		if (generation % migrationTime == 0){
			pop.migrate(size, rank);
			for(int i = 1; i < 5; i++)
			pop.salesmen(i) = pop.salesmen(0); 
			}
		//selection
		pop.selection(selectionPow);
		
		if(generation%500 == 0 && rank == 0)
			pop.printUpdate(double(generation)/maxExe);
				

		pop.mutation(mutationAcceptance);
    	pop.crossover(crossoverAcceptance);
		
		generation++;
   }
   pop._sort();

	MPI_Barrier(MPI_COMM_WORLD);
	int best_rank = pop.olympics(size, rank);
	if(rank == best_rank){
		pop.salesmen(0).save(italy, rank);
		cout << "Process " << rank << " with distance " << pop.salesmen(0).distance << endl;
	}

	MPI_Finalize();
	return 0;
}

