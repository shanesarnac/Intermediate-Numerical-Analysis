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

double dtheta_dsigma(double sigma, double theta) {
	double delta = 1.0/3.0;
	return delta*exp(theta) - theta;
}

int main() {
	double h = 0.01;
	Data_Point initial_data = Data_Point(0.0, 0.0);
	double x_max = 20.0;
	vector<Data_Point> theta_estimates = rk4(h, dtheta_dsigma, initial_data, x_max);
	
	int num_points = theta_estimates.size();
	cout << "sigma, delta" << endl;
	for (int i = 0; i < num_points; i++) {
		cout << theta_estimates[i].getX() << ", " << theta_estimates[i].getY() << endl;
	}

	
}

#endif
