#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "random.h"
#include "random_setup.h"

using namespace std;

double s_sample(float S0, float r, float sigma, float t0, float t1, Random &rnd) {
    float z = rnd.Gauss(0, 1);
    return S0 * exp((r - sigma * sigma / 2) * (t1 - t0) + sigma * z * sqrt(t1 - t0));
}

double compute_sigma(double avg, double avg2, int blk) {
    if (blk <= 1) return 0;
    return sqrt(abs(avg2 - avg * avg) / blk);
}

int fileHeader(string filename) {
    ofstream out;
    out.open(filename);
    if (out.is_open()) {
        out << "#BLOCK\t#AVG\tSIGMA" << endl;
        out.close();
        return 1;
    } else {
        cerr << "Problem reading " << filename << "." << endl;
        return -1;
    }
}

int save(double avg, double avg2, int block, string filename) {
    ofstream out;
    out.open(filename, ios::app);
    if (out.is_open()) {
        out << block << "\t" << avg << "\t" << compute_sigma(avg, avg2, block) << endl;
        out.close();
        return 1;
    } else {
        cerr << "Problem reading " << filename << "." << endl;
        return -1;
    }
}

int main() {
    Random rnd;
    rnd_setup(rnd);

    // Market parameters
    const float r = 0.1;
    const float sigma = 0.25;
    float S0 = 100.;
    float t0 = 0.;
    float t_delivery = 1.0;
    float k = 100.;

    // Block parameters
    const int blocks = 100;
    const int steps = 2000;

    // Variables for storing averages
    double c_avg = 0, c_avg2 = 0;
    double p_avg = 0, p_avg2 = 0;

    // Direct sampling
    fileHeader("../OUTPUT/BS_call_direct.dat");
    fileHeader("../OUTPUT/BS_put_direct.dat");

    for (int blk = 1; blk <= blocks; blk++) {
        double call = 0, put = 0;

        for (int i = 0; i < steps; i++) {
            double St = s_sample(S0, r, sigma, t0, t_delivery, rnd);
            call += exp(-r) * max(0., St - k);
            put += exp(-r) * max(0., k - St);
        }

        call /= steps;
        c_avg += call;
        c_avg2 += call * call;
        save(c_avg / blk, c_avg2 / blk, blk, "../OUTPUT/BS_call_direct.dat");

        put /= steps;
        p_avg += put;
        p_avg2 += put * put;
        save(p_avg / blk, p_avg2 / blk, blk, "../OUTPUT/BS_put_direct.dat");
    }

    cout << "call: " << c_avg / blocks << endl;
    cout << "put: " << p_avg / blocks << endl;

    // Reset averages for the discretized sampling
    c_avg = 0;
    c_avg2 = 0;
    p_avg = 0;
    p_avg2 = 0;

    // Discretized sampling
    fileHeader("../OUTPUT/BS_call_discrete.dat");
    fileHeader("../OUTPUT/BS_put_discrete.dat");

    for (int blk = 1; blk <= blocks; blk++) {
        double call = 0, put = 0;

        for (int i = 0; i < steps; i++) {
            double St = S0;
            for (double t = 0.; t < t_delivery; t += 0.01) {
                St = s_sample(St, r, sigma, t, t + 0.01, rnd);
            }
            call += exp(-r) * max(0., St - k);
            put += exp(-r) * max(0., k - St);
        }

        call /= steps;
        c_avg += call;
        c_avg2 += call * call;
        save(c_avg / blk, c_avg2 / blk, blk, "../OUTPUT/BS_call_discrete.dat");

        put /= steps;
        p_avg += put;
        p_avg2 += put * put;
        save(p_avg / blk, p_avg2 / blk, blk, "../OUTPUT/BS_put_discrete.dat");
    }

    cout << "\n\ncall: " << c_avg / blocks << endl;
    cout << "put: " << p_avg / blocks << endl;

    return 0;
}
