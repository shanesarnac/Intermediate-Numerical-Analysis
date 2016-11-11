#ifndef ODE_SOLVER_CPP
#define ODE_SOLVER_CPP

#include "ode_solver.h"

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0) {
	cout << "RK4" << endl;
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
		cout << "k1 = " << k1 << endl;
		cout << "k2 = " << k2 << endl;
		cout << "k3 = " << k3 << endl;
		cout << "k4 = " << k4 << endl;
		cout << "y_{i + 1}(" << x_i + h << ") = " << y_ip1 << endl << endl;
	}
	
	return acc;
}

double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1){
	cout << "Adams-Bashforth Two-Step" << endl;
	double acc = 0.0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		/*cout << "x_i = " << x_i << endl;
		cout << "y_i = " << y_i << endl;
		cout << "y_{i  -1} = " << y_im1 << endl;*/
		
		y_ip1 = y_i + (h/2.0)*(3*f(x_i, y_i) - f(x_i - h, y_im1));
		y_im1 = y_i;
		y_i = y_ip1;

		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

// Note: Only works for y' = y + x (had to solve for y_ip1)
double adams_moulton_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1) {
	cout << "Adams-Moulton Two Step" << endl;
	double acc = 0.0;
	double x_im1 = x_0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		y_ip1 = (y_i + (h/12.0) * (5.0*(x_i + h) + 8.0*f(x_i, y_i) - f(x_im1, y_im1)))*pow(1.0 - (5.0*h)/12.0, -1.0);
		y_im1 = y_i;
		x_im1 = x_i;
		y_i = y_ip1;
		
		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

double improved_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0) {
	cout << "Improved Euler's Method" << endl;
	double y_i = y_0;
	double y_bar, y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_bar = y_i + h*f(x_i, y_i);
		y_ip1 = y_i + (h/2.0)*(f(x_i,y_i) + f(x_i + h, y_bar));
		y_i = y_ip1;
		cout << "y bar = " << y_bar << endl;
		cout << "y(" << x_i + h<< ") = " << y_ip1 << endl << endl;
	}
	cout << endl;
	
	return 0.0;
}



#endif
