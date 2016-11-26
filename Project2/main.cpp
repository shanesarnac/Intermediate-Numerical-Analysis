#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <math.h>
#include <vector>

#include "data_point.h"
#include "integrate.h"
#include "interpolation.h"
#include "ode_solver.h"

using namespace std;

vector<vector<double> > f_and_g_values;


double f1(double nu, double v1_i, double v2_i, double u1_i, double u2_i, double u3_i) {
	return u2_i;
}

double f2(double nu, double v1_i, double v2_i, double u1_i, double u2_i, double u3_i) {
	return u3_i;
}

double f3(double nu, double v1_i, double v2_i, double u1_i, double u2_i, double u3_i) {
	return -0.5*u1_i*u3_i;
}

double g1(double nu, double v1_i, double v2_i, double u1_i, double u2_i, double u3_i) {
	return v2_i;
}

double g2(double nu, double v1_i, double v2_i, double u1_i, double u2_i, double u3_i) {
	double pr = 4.0; // can change
	return -(pr / 2.0)*u1_i*v2_i;
}

//double F(double nu, double v1_i, double v2_i, double u1_i) {
	//cout << "nu = " << nu << endl;
	//int f_size = f_values[0].size();
	//for (int i = 0; i < f_size; i++) {
		//cout << "nu, f_values[0][i] = " << nu << "," << f_values[0][i] << endl;
		//if (nu == f_values[0][i]) {
			//return f_values[1][i];
		//}
	//}
	//cout << "Error: nu value not present" << endl;
	//return 0.0;
//}

int main() {
	double h = 0.1;
	double max = 10.0;
	double t0 = 0;
	double x10 = 1; // G(0) = 1
	double x20 = -0.534787; // Guess for G'(0)
	double x30 = 0; // F(0) = 0
	double x40 = 0; // F'(0) = 0
	double x50 = 0.332057; // Guess for F''(0)
	
	/* F and G value indexing:
	 * 		0: mu values
	 * 		1: G values
	 * 		2: G' values
	 * 		3: F values
	 * 		4: F' values
	 * 		5: F'' values*/
	f_and_g_values = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	//cout << "mu,G, G',F,F',F''" << endl;
	//for (unsigned int i = 0; i < f_and_g_values[0].size(); i++) {
		//cout << f_and_g_values[0][i] << "," << f_and_g_values[1][i] << "," << f_and_g_values[2][i] << "," << f_and_g_values[3][i] <<",";
		//cout << f_and_g_values[4][i] << "," << f_and_g_values[5][i] << endl;
	//}
	
	//rk4_3_coupled(h, f1, f2, f3, max, t0, x0, y0, z0);
	//cout << "mu,F,F',F''" << endl;
	//for (unsigned int i = 0; i < f_values[0].size(); i++) {
		//cout << f_values[0][i] << "," << f_values[1][i] << "," << f_values[2][i] << "," << f_values[3][i] << endl;
	//}
	
	//h = 0.1;
	//z0 = f_values[1][0]; // Guess for G'(0)
	//g_values = rk4_3_coupled(h, g1, g2, F, max, t0, x0, y0, z0);
	//cout << "mu,G,G',F?" << endl;
	//for (unsigned int i = 0; i < g_values[0].size(); i++) {
		//cout << g_values[0][i] << "," << g_values[1][i] << "," << g_values[2][i] << "," << g_values[3][i] << endl;
	//}
	
}

#endif
