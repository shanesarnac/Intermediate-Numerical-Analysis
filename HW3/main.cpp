
#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <math.h>

#include "differentiation.h"
#include "integrate.h"

using namespace std;

double function1(double x) {
	return pow(x,2)*exp(-pow(x,2));
}

double function2(double x) {
	return pow(x, 2.0) * exp(-x);
}


int main() {
	double f_data[] = {9.025013, 11.02318, 13.46374, 16.44465};
	cout << "f'(x0 = " << 1.1 << ") = " << three_point_foward(0.1, f_data[0], f_data[1], f_data[2]) << " using three point forward" << endl;
	cout << "f'(x0 = " << 1.2 << ") = " << three_point_middle(0.1, f_data[0], f_data[2]) << " using three point middle" << endl;
	cout << "f'(x0 = " << 1.3 << ") = " << three_point_middle(0.1, f_data[1], f_data[3]) << " using three point middle" << endl;
	cout << "f'(x0 = " << 1.4 << ") = " << three_point_backward(0.1, f_data[1], f_data[2], f_data[3]) << " using three point backward" << endl;
	
	cout << "integrate from 0 to 2 of x^2 e^(-x^2) = " << midpoint_rule(0.25, function1, 0, 2) << " using midpoint rule.";
	cout << endl;
	
	cout << "integrate from 0 to 2 of x^2 e^(-x^2) = " << trapezoidal_integration(0.25, function1, 0, 2) << " using trapezoidal rule.";
	cout << endl;
	
	cout << "integrate from 0 to 2 of x^2 e^(-x^2) = " << simpsons_one_third_integration(0.25, function1, 0, 2) << " using simpson's 1/3 rule.";
	cout << endl;
	
	cout << "integrate from 0 to 1 of x^2 e^(-x) = " << gaussian_quadrature_n2(function2, 0, 1) << " using gaussian quadrature with n = 2";
	cout << endl;
	
	cout << "integrate from 0 to 1 of x^2 e^(-x) = " << gaussian_quadrature_n3(function2, 0, 1) << " using gaussian quadrature with n = 3";
	cout << endl;
	
	
	return 0;
}

#endif
