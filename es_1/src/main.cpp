#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "lib.h"
#include "random.h"
#include "random_setup.h"

using namespace std;

int main(){

Random rnd;
rnd_setup(rnd);

// integral in [0,1] of f(r) = r
_integral("straight", rnd);
// integral in [0,1] of f(r) = (r^2 - 0.5)
_integral("parabolic", rnd);

// 2 X^2 test
chi_test(rnd);

// histograms
int N_values[4] = {1, 2, 10, 100};
int dice_throws = pow(10,4);
generate_dist("uniform", rnd, N_values, dice_throws);
generate_dist("exponential", rnd, N_values, dice_throws);
generate_dist("cauchy", rnd, N_values, dice_throws);


//Buffon experiment
int pi_throws = pow(10,4);
int n_blocks = 40;
buffon(pi_throws, n_blocks, rnd);


return 0;

}

