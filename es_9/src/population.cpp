#include "population.h"

population :: population(int n_Salesmen, int n_cities)
    : numSalesmen(n_Salesmen), numCities(n_cities) // Initialize member variables
{
	// Set up random seeds
	rnd_setup(rnd);
	// Resize population pool and initialize it
	salesmen.set_size(numSalesmen);
	for(int i = 0; i < numSalesmen; ++i)
		salesmen(i).initialize(numCities, rnd);
}

population :: ~population(){}

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

void population :: bestHalf_clear(string type){
	ofstream besthalf("../OUTPUT/"+ type +"_bestHalf.dat");
	besthalf.close();
	ofstream best("../OUTPUT/"+ type +"_best_individual.dat");
	best.close();
}

void population :: bestHalf(string type){
	ofstream besthalf("../OUTPUT/"+ type +"_bestHalf.dat", ios::app);
	double half_best = 0.;
	for(int i = 0; i < 0.5*numSalesmen; i++)
		half_best += salesmen(i).distance;
	half_best /= (0.5*numSalesmen);
	besthalf << half_best << endl;	
	besthalf.close();
	
	fstream best("../OUTPUT/"+ type +"_best_individual.dat", ios::app);
	best << salesmen(0).distance << endl;
	best.close();
}


void population :: selection(double _pow){
	
	field <Salesman> clones = salesmen;
	
	for (int i = 1; i < numSalesmen; i++) {
		 int son = int( (numSalesmen) * pow(rnd.Rannyu(), _pow));
		 clones(i) = salesmen(son);
	}
	
	salesmen = clones;
};
//crossover
void population :: crossover(double crossoverAcceptance){

int crossed = 0;
	while(crossed < numSalesmen){
		if(rnd.Rannyu() < crossoverAcceptance){ 
			int cut = 1 + int(rnd.Rannyu()*(numCities-2));
			int dad = 1+int(rnd.Rannyu()*(numSalesmen-1));
			int mom = 1+int(rnd.Rannyu()*(numSalesmen-1));
			while(mom == dad) mom = int(rnd.Rannyu()*numSalesmen); 
			salesmen(mom).crossover(salesmen(dad), salesmen(mom), cut);
			salesmen(dad).crossover(salesmen(mom), salesmen(dad), cut);
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

