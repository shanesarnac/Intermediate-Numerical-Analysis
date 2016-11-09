#ifndef ODE_SOLVER_CPP
#define ODE_SOLVER_CPP

#include "ode_solver.h"

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0) {
	double acc = 0.0;
	double y_ip1 = y_0;
	double y_i;
	cout << "y_{i + 1}(" << x_0 << ") = " << y_ip1 << endl;
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_i = y_ip1;
		double k1 = h*f(x_i, y_i);
		double k2 = h*f(x_i + h/2.0, y_i + k1/2.0);
		double k3 = h*f(x_i + h/2.0, y_i + k2/2.0);
		double k4 = h*f(x_i + h, y_i + k3);
		y_ip1 = y_i + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
		cout << "y_{i + 1}(" << x_i + h << ") = " << y_ip1 << endl;
	}
	
	return acc;
}



#endif
