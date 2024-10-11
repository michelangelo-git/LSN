#include "lib.h"

using namespace std;

double sigma(double integral, double integral2, int n){
   if(n < 2)
      return 0;
  
   return sqrt( (integral2- integral*integral) / (n-1) );
}


////////////////////////////////////////////////////////////////////
//

void _integral(string integrand_func, Random& rnd){
	string filename = "../OUTPUT/" + integrand_func + "_integral_01.txt";
	char which_func = integrand_func[0];
	const int throws = 100000;
	const int n_blocks = 50;
	const int block_size = throws / n_blocks;
	double integral = 0;
	double integral2 = 0;

	ofstream data(filename);
	for(int blk = 1; blk <= n_blocks; blk++){
		double temp_integral = 0;
	   for(int k = 0; k < block_size; k++){
	      if(which_func == 's')
	      	temp_integral += rnd.Rannyu();
	      else if(which_func == 'p')
				temp_integral += pow(rnd.Rannyu() - 0.5, 2);
	   }
	   temp_integral /= block_size;
	   integral += temp_integral;
	   integral2 += temp_integral*temp_integral;   

	   data << integral/blk << " " << sigma(integral/blk, integral2/blk, blk) << endl;
	}
	
	cout << integrand_func << ":\t" << integral/n_blocks 
	     << " +/- " << sigma(integral/n_blocks, integral2/n_blocks, n_blocks) << endl;
	
	data.close();
}

void chi_test(Random& rnd){

	int Xi_throws = pow(10,4);
	int bins = 100;
	int n_blocks = Xi_throws / 100;
	double Xi_bins[bins] = {0};

	ofstream chi_data("../OUTPUT/chi_data.txt");

	for(int blk = 0; blk < n_blocks; blk++){
		double expected = Xi_throws * (blk + 1) / bins;
		for(int i = 0; i < Xi_throws; i++){
			int bin_index = int(rnd.Rannyu() * 100);
			Xi_bins[bin_index]++;
		}
		double Xi = 0;
		for(int i = 0; i < bins; i++){
			Xi += pow(Xi_bins[i] - expected, 2) / expected;
		}
		chi_data << Xi << endl;
	}
	chi_data.close();
}

/////////////////////////////////////////////////////////////
//2: histogram with different distributions

double which_dist(char dist, Random& rnd){
	if(dist == 'u')
		return rnd.Rannyu();
	if(dist == 'e')
		return rnd.Exponential(1.);
	if(dist == 'c')
		return rnd.CauchyLorentz(0., 1.);
	return 0;
}

void generate_dist(string distribution, Random& rnd, int* N_values, int dice_throws){

	char dist = distribution[0];
	string filename = "../OUTPUT/" + distribution + ".txt";
	ofstream data;
	data.open(filename);

	for(int n = 0; n < dice_throws; n++){
		for(int i = 0; i < 4; i++){
			double out = 0;
			for(int k = 0; k < N_values[i]; k++){
				out += which_dist(dist, rnd);
			}
			out /= N_values[i];
			data << out << " ";
		}
		data << endl;
	}
	data.close();
}

//////////////////////////////////////////////////////////////
//3: buffon simulation

bool throw_needle(Random& rnd, double d, double L){
	// Throwing a needle with a rejection technique 
	double x, y;
	double y_center = rnd.Rannyu(-d/2, d/2);
	bool inside = false;
	while(!inside){
		x = rnd.Rannyu(-L, L);
		y = rnd.Rannyu(-L, L);
		if((x*x + y*y) > 0 && (x*x + y*y) <= L*L)
			inside = true;
	}
	// finding the y intersection with the circle of radius = L
	y = 0.5*y*L / sqrt(x*x + y*y);
	
	// checking intersection with lines
	if(y_center > 0) y = y_center + abs(y);
	else y = y_center - abs(y);
	
	if(abs(y) >= d/2) return 1;
	return 0;
}

double pi_estimate(double probability, double L, double d){
	//This function evaluates PI, given a probability estimate
	return 2 * L / (d * probability);
}


void buffon(int pi_throws, int n_blocks, Random& rnd){
	double d = 1.;
	double L = 1.;
	
	const int pi_block_size = pi_throws / n_blocks;
	double pi_estimates = 0;
	double pi_estimates2 = 0;
	ofstream pi_data("../OUTPUT/PI_data.txt");

	for(int k = 1; k <= n_blocks; k++){
		double _temp_pi = 0;
		for (int step = 0; step < pi_block_size; step++)
			_temp_pi += throw_needle(rnd, d, L); //checking if each needle intersects a line
		_temp_pi /= pi_block_size; // Probability = successes / attempts
		_temp_pi = pi_estimate(_temp_pi, d, L); // converting Probability to Pi
		pi_estimates += _temp_pi;
		pi_estimates2 += _temp_pi*_temp_pi;
		pi_data << pi_estimates/k << " " << sigma(pi_estimates/k, pi_estimates2/k, k) << endl;
	}

	cout << "Buffon estimate of PI:\t" << pi_estimates/n_blocks
		  << "\t+/- " << sigma(pi_estimates/n_blocks, pi_estimates2/n_blocks, n_blocks) << endl; 

	pi_data.close();
	
}




