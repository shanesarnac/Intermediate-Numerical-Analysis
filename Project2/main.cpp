#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <math.h>
#include <vector>

#include "data_point.h"
#include "integrate.h"
#include "interpolation.h"
#include "ode_solver.h"
#include "root_finding.h"

using namespace std;

vector<vector<double> > f_and_g_values;
double prandtl;


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
	double pr = prandtl; // can change
	return -(pr / 2.0)*u1_i*v2_i;
}

double f_prime_estimate(double neta) {
	Interpolation interp;
	vector<Data_Point> f_prime_values;
	
	for (unsigned int i = 40; i < 55; i++) {
		f_prime_values.push_back(Data_Point(f_and_g_values[0][i], f_and_g_values[4][i]));		
	}
	
	/*for (unsigned int i = 0; i < f_prime_values.size(); i++) {
		cout << "neta, F' => " << f_prime_values[i].getX() << "," << f_prime_values[i].getY() << endl;
	}*/
	//return 0;
	return interp.laGrange(neta, f_prime_values);	
}

double g_estimate(double neta) {
	Interpolation interp;
	double h = 0.1;
	double max = 10.0;
	double t0 = 0;
	double x10 = 1; // G(0) = 1
	
	double x20 = -0.184096; // Guess for G'(0) for prandtl = 0.2
	//double x20 = -0.332057; // Guess for G'(0) for prandtl = 1.0
	//double x20 = -0.534787; // Guess for G'(0) for prandtl = 4.0
	//double x20 = -0.675580; // Guess for G'(0) for prandtl = 8.0
	
	double x30 = 0; // F(0) = 0
	double x40 = 0; // F'(0) = 0
	double x50 = 0.332057; // Guess for F''(0)
	
	prandtl = 0.2;
	//prandtl = 1.0;
	//prandtl= 4.0;
	//prandtl= 8.0;
	
	vector<vector<double> > estimates = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	vector<Data_Point> g_values;
	for (unsigned int i = 75; i < 95; i++) {
		g_values.push_back(Data_Point(estimates[0][i], estimates[1][i]));		
	}
	
	/*for (unsigned int i = 0; i < g_values.size(); i++) {
		cout << "neta, G => " << g_values[i].getX() << "," << g_values[i].getY() << endl;
	}*/
	//return 0;
	return interp.laGrange(neta, g_values);	
}

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
	prandtl = 4.0;
	f_and_g_values = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	
	// Finding non-dimensionalized velocities
	// u component and v component
	/*cout << "mu, horizontal velocity, vertical velocity" << endl;
	for (int i = 0; i < (int) f_and_g_values[0].size(); i++) {
		double mu = f_and_g_values[0][i];
		double f_prime = f_and_g_values[4][i];
		double f = f_and_g_values[3][i];
		double u_value = f_prime;
		double v_value = 0.5*(mu*f_prime - f);
		cout << mu <<"," << u_value << "," << v_value << endl;
	}*/
	
	// Finding non-dimensionalized temperature for varying Pr values
	prandtl= 0.2;
	x20 = -0.184096; // Need G to approach 0, so make guess like before
	vector<vector<double> > g_pr_fifth = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	
	prandtl= 1.0;
	x20 = -0.332057;
	vector<vector<double> > g_pr_1 = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	
	prandtl= 4.0;
	x20 = -0.534787;
	vector<vector<double> > g_pr_4 = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	
	prandtl= 8.0;
	x20 = -0.675580;
	vector<vector<double> > g_pr_8 = rk4_5_coupled(h, g1, g2, f1, f2, f3, max, t0, x10, x20, x30, x40, x50);
	
	/**cout << "mu, G with Pr = 0.2, G with Pr = 1, G with Pr = 4, G with Pr = 8" << endl;
	for (int i = 0; i < (int) g_pr_fifth[0].size(); i++) {
		double mu = g_pr_fifth[0][i];
		double pr_fifth = g_pr_fifth[1][i];
		double pr_1 = g_pr_1[1][i];
		double pr_4 = g_pr_4[1][i];
		double pr_8 = g_pr_8[1][i];
		cout << mu << "," << pr_fifth << "," << pr_1 << "," << pr_4 << "," << pr_8 << endl;
	}**/
	
	// Finding neta_m Using Bisection and laGrange polynomials

	//double neta_m_estimate = bisection(f_prime_estimate, 0.98, 4.5, 4.6);
	//cout << "neta_m = " << neta_m_estimate << endl;
	double neta_t_estimate = bisection(g_estimate, 0.02, 9.0, 7.0);
	cout << "neta_t = " << neta_t_estimate << endl; 
	
	
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
