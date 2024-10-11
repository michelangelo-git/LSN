
#ifndef __Random_setup__
#define __Random_setup__

// This class contains a function to setup the seeds of the random number generator

int rnd_setup(Random& rnd);
int rnd_setup_multi(Random& rnd, int rank);

#endif
