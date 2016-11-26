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

double f1(double nu, double u1_i, double u2_i, double u3_i) {
	return u2_i;
}

double f2(double nu, double u1_i, double u2_i, double u3_i) {
	return u3_i;
}

double f3(double nu, double u1_i, double u2_i, double u3_i) {
	return -0.5*u1_i*u3_i;
}

double g1(double nu, double 

int main() {
	double h = 0.1;
	double max = 20.0;
	double t0 = 0;
	double x0 = 0; 
	double y0 = 0;
	double z0 = 0.332057; // Guess
	
	rk4_3_coupled(h, f1, f2, f3, max, t0, x0,y0,z0);
}

#endif
