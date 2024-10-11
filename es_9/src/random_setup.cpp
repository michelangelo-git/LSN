#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "random_setup.h"

using namespace std;
 
int rnd_setup(Random& rnd){

   int seed[4];
   int p1, p2;
   ifstream Primes("../INPUT/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../INPUT/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();
   return 0;
}

int rnd_setup_multi(Random& rnd, int rank){

   int seed[4];
   int p1, p2;
   ifstream Primes("../INPUT/Primes");
   if (Primes.is_open()){
      for(int i = 0; i < rank; i++)
         Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../INPUT/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            // change the seed depending on rank
            for(int i = 0; i < 4; i++){
               seed[i] += rank;
            }
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();
   return 0;
}