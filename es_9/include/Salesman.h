#ifndef __Salesman__
#define __Salesman__

#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <algorithm>
#include <unordered_set>
#include "random.h"

using namespace std;
using namespace arma;

class city {
public:
    double x, y;
};

class Salesman {
public:
    Salesman();
    ~Salesman();
    double distance;
    vec visitedCities;
    void initialize(int numCities, Random& rnd); 
    double cost(const city& last, const city& next);
    void updateDistance(const vector<city>& map);   //calculate the distance of the current path  
    //crossover    
    void crossover(Salesman& parent1, Salesman& parent2, int cut);    
    //mutations    
    void quickSwap(int start); // swap two contiguous cities
    void swapCities(int start, int pathLenght); // swap a n-long path with a contiguous m-long path
    void shift(int start, int n, int m); //shift by n an m-long path
    void reversePath(int start, int n); // reverse a n-long path
    void printPath();
    void save(const vector<city>& map, string output_name);
};

#endif
