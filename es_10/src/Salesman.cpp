#include "Salesman.h"



Salesman::Salesman() {}
Salesman::~Salesman() {}

void Salesman :: initialize(int numCities, Random &rnd) {
     visitedCities.set_size(numCities);
     for(int i = 0; i < numCities; i++)
     		visitedCities(i) = i;
     int shuffles = numCities*4;
     for(int i = 0; i < shuffles; i++){
     		int city1 = 1+int( rnd.Rannyu()*(numCities-1) );
     		int city2 = 1+int( rnd.Rannyu()*(numCities-1) );
     		swap(visitedCities(city1), visitedCities(city2));
     	}
 }

 /*void Salesman :: visitCity(Random &rnd, int c, int numCities) {
     int index = c + int(rnd.Rannyu()*double(numCities -c));
     //cout << "\nindex " << index << " - VC(c) " <<  visitedCities(c) << " - VC(I) " <<  visitedCities(index) << endl  ;
     int temp = visitedCities(c);
     visitedCities(c) = visitedCities(index);
     visitedCities(index) = temp;       
 }*/
 
 
 double Salesman :: cost(const city& last, const city& next) {
 	return sqrt(pow(last.x - next.x,2) + pow(last.y - next.y,2)); // the ideal path is slightly shorter than 2PI
 }

 void Salesman :: updateDistance(const vector<city>& map) {
distance = 0;
int numCities = int(visitedCities.n_elem);
for (int i = 0; i < numCities-1; i++) {
    distance += cost(map[visitedCities(i)], map[visitedCities(i + 1)]);
}
// Add the distance back to the starting city to complete the cycle
distance += cost(map[visitedCities(numCities-1)], map[visitedCities(0)]);
 }

 
 void Salesman::crossover(Salesman& parent1, Salesman& parent2, int cut) {
 int numCities = int(visitedCities.n_elem);

 // Ensure cut is within bounds
 if (cut < 0 || cut >= numCities) {
     throw std::out_of_range("Cut index is out of range.");
 }
 // Copy the first part from parent1
 for (int i = 0; i < cut; i++) {
     visitedCities(i) = parent1.visitedCities(i);
 }

 // Copy the remaining part from parent2, excluding cities in parent1
 int index = cut;
 for (int i = 0; i < numCities; i++) {
 	if(index == numCities)
 		break;
     bool found = false;
     for (int vc = 0; vc <= cut; vc++) {
         if (parent2.visitedCities(i) == visitedCities(vc)) {
             found = true;
             break;
         }
     }
     if (!found) {
         visitedCities(index++) = parent2.visitedCities(i);
     }
 }
 
 if(index != numCities){
 	//cout << "problems with CO: " << index << endl;
 	visitedCities = parent1.visitedCities; 
 }

}

 
 void Salesman :: quickSwap(int start){
 	if(start >= int(visitedCities.n_elem -1)) return;
 	swap (visitedCities(start), visitedCities(start+1));
 }
 
 void Salesman :: shift(int start, int shift, int shiftBox){
    int numCities = int(visitedCities.n_elem);
    if (start + shiftBox >= numCities) return;
    
    vec temp(shiftBox);
    for(int i = 0; i < shiftBox; i++)
    	temp(i) = visitedCities(start+i);
    
    for(int i = 0; i < shiftBox; i++)
    	visitedCities(start + i) = temp((i + shift) % shiftBox);  
    }
 
 
 void Salesman :: swapCities(int start, int pathLenght){
 	if(start + 2*pathLenght >= visitedCities.n_elem) return;
 	
 	for(int i = 0; i < pathLenght; i++)
 		swap(visitedCities(start+i), visitedCities(start + pathLenght +i));
 }
 
 void Salesman :: reversePath(int start, int n){
if(start + n > int(visitedCities.n_elem)) return;
for(int i = 0; i < n; i++)
	swap( visitedCities(start+i), visitedCities(start+n--));
 }
 
 void Salesman :: printPath(){
     visitedCities.print("Visited Cities:");
 }
 
 void Salesman :: save(const vector<city>& map, int rank){
 
 ofstream best;
 best.open("../OUTPUT/bestSalesman.dat");
 best << "#Rank: " << rank << endl;
 best << "#Distance: " << distance << endl;
 best << "#Path coordinates: \n# x:\ty:" << endl;
 int numCities = int(visitedCities.n_elem);
 for(int i = 0; i < numCities; i++){
 	best << " " << map[visitedCities(i)].x << "\t" << map[visitedCities(i)].y << endl;
 }
 best.close();
 
 }

