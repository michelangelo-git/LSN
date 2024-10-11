#include "population.h"

population :: population(int n_Salesmen, int n_cities, int rank)
    : numSalesmen(n_Salesmen), numCities(n_cities) // Initialize member variables
{
	// Set up random seeds
	rnd_setup_multi(rnd, rank);
	rnd_setup(migControl);

	// Resize population pool and initialize it
	salesmen.set_size(numSalesmen);
	for(int i = 0; i < numSalesmen; ++i)
		salesmen(i).initialize(numCities, rnd);
}

population :: ~population(){}

// Custom sort function for arma::field<Salesman>
void population :: _sort(){
	 // Simple selection sort for arma::field<Salesman>
	 for (int i = 0; i < numSalesmen - 1; i++) {
	     int minIndex = i;
	     for (int j = i + 1; j < numSalesmen; j++) {
	         if (salesmen(j).distance < salesmen(minIndex).distance) {
	             minIndex = j;
	         }
	     }
	     // Swap the salesmen at i and minIndex
	     if (minIndex != i) {
	         Salesman temp = salesmen(i);
	         salesmen(i) = salesmen(minIndex);
	         salesmen(minIndex) = temp;
	     }
	 }
}

double population :: _sigma_dist() {
	 int dim = numSalesmen;
	 if (dim <= 1) return 0;
	 double avg = 0;
	 for (int i = 0; i < dim; i++) {
	     avg += salesmen(i).distance;
	 }
	 avg /= dim;
	 double diff = 0;
	 for (int i = 0; i < dim; i++) {
	     diff += pow(salesmen(i).distance - avg, 2);
	 }
	 return sqrt(diff / (dim-1));
}

void population :: printUpdate(double progress){
	vec distances(numSalesmen);
	for (int i = 0; i < numSalesmen; ++i)
		distances(i) = salesmen(i).distance;
	double sigma = _sigma_dist();
	cout << progress << " sigma: " << sigma;
		cout << " best: " << salesmen(0).distance << endl;
}

void population :: selection(double _pow){
	// Selection
	field <Salesman> clones = salesmen;
	for (int i = 1; i < numSalesmen; i++) {
		 int son = int( (numSalesmen) * pow(rnd.Rannyu(), _pow));
		 clones(i) = salesmen(son);
	}
	
	salesmen = clones;
};

//migration
void population :: migrate(int size, int rank){
	if(size < 2) return;
	
	int* couples = new int[size];
	for(int i = 0; i < size; i++)
		couples[i] = i;
	for(int i = 0; i < 22; i++){
		int r1 = migControl.Rannyu() * size;
		int r2 = migControl.Rannyu() * size;
		swap(couples[r1], couples[r2]);
	}
	int _vec_size = int(salesmen(0).visitedCities.n_elem);
	vec _local(_vec_size);
	vec _send(_vec_size);

	for(int i = 0; i < _vec_size; i++){
		_local(i) = salesmen(0).visitedCities(i);
		_send(i) = salesmen(0).visitedCities(i);
	}
	for(int m = 0; m < size - 1; m += 2) {
		if(rank == couples[m]) {           
			MPI_Sendrecv(_send.memptr(), _vec_size, MPI_DOUBLE, couples[m+1], 0, 
					           _local.memptr(), _vec_size, MPI_DOUBLE, couples[m+1], 1, 
					           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		}else if(rank == couples[m+1]) {
			MPI_Sendrecv(_send.memptr(), _vec_size, MPI_DOUBLE, couples[m], 1, 
					           _local.memptr(), _vec_size, MPI_DOUBLE, couples[m], 0, 
					           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}

	for(int i = 0; i < _vec_size; i++)
		salesmen(0).visitedCities(i) = _local(i);
		
	delete[] couples;
}

void population :: crossover(double crossoverAcceptance) {
    int crossed = 0;
    while (crossed < numSalesmen) {
        if (rnd.Rannyu() < crossoverAcceptance) {
            int cut = 1 + int(rnd.Rannyu() * (numCities - 2));
            int dad = 1 + int(rnd.Rannyu() * (numSalesmen - 1));
            int mom;
            do {
                mom = 1 + int(rnd.Rannyu() * (numSalesmen - 1));
            } while (mom == dad);

            // Perform crossover only once and assign offspring
            Salesman offspring1 = salesmen(dad);
            Salesman offspring2 = salesmen(mom);

            offspring1.crossover(salesmen(dad), salesmen(mom), cut);
            offspring2.crossover(salesmen(mom), salesmen(dad), cut);

            // Assign offspring back to population
            salesmen(dad) = offspring1;
            salesmen(mom) = offspring2;
        }
        crossed++;
    }
}


//mutation
void population :: mutation(double mutationAcceptance){
	for (int mutant = 1; mutant < numSalesmen; mutant++) {
		//shift by n positions	
		if (rnd.Rannyu() < mutationAcceptance) {
			int start = 1 + int((numCities-1) * rnd.Rannyu());
			salesmen(mutant).quickSwap(start);
		}
		//shift by n positionss
		if (rnd.Rannyu() < mutationAcceptance) {
			int shiftBox = int(rnd.Rannyu() * (numCities-1));
			int shift = 1 + int(rnd.Rannyu() * shiftBox);
			int start = 1 + int((numCities - shiftBox) * rnd.Rannyu());
			salesmen(mutant).shift(start, shift, shiftBox);
		}
		//swap blocks of the path
		if (rnd.Rannyu() < mutationAcceptance) {
			int n =  + int(rnd.Rannyu() * (numCities-1)/2 );
			int start = 1 + int((numCities - 2*n) * rnd.Rannyu());
			salesmen(mutant).swapCities(start, n);
		}
		//order mutation
		if (rnd.Rannyu() < mutationAcceptance) {
			int length = 1 + int(rnd.Rannyu() * (numCities-2));
			int start = 1 + int((numCities - length - 1) * rnd.Rannyu());
			salesmen(mutant).reversePath(start, length);
		}
	} 
}

int population :: olympics(int size, int rank){
	
	cout << "Process " << rank << " distance " <<  salesmen(0).distance << endl;

	double dist;
	double *all_distances = new double[size];
	// Gather all distances by broadcasting
	for (int i = 0; i < size; i++) {
		dist = salesmen(0).distance;
		MPI_Bcast(&dist, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
		all_distances[i] = dist;
	}
	// Now, each process has the distances of all processes
	int best_rank = rank;
	double best_distance = numCities * 1000;

	for (int i = 0; i < size; i++) {
		if (all_distances[i] < best_distance) {
			best_distance = all_distances[i];
			best_rank = i;
		}
	}
	
	return best_rank;
}


void population :: bestHalf_clear(int rank){
	string spec = to_string(rank);
	ofstream best("../OUTPUT/"+ spec +"_best_individual.dat");
	best.close();
}

void population :: bestHalf(int rank){
	string spec = to_string(rank);
	
	fstream best("../OUTPUT/"+ spec +"_best_individual.dat", ios::app);
	best << salesmen(0).distance << endl;
	best.close();
}

