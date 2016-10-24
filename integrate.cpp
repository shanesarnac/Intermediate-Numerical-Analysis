
#ifndef INTEGRATE_CPP
#define INTEGRATE_CPP

#include "integrate.h"

double trapezoidal_integration(double h, double (*f)(double), double x0, double x1) {
	double area = 0.0;
	
	for (double i = 0; i < x1; i+=h) {
		area += (h/2)*(f(i) + f(i + h));
	}
	
	return area;
}

double simpsons_one_third_integration(double h, double (*f)(double), double x0, double x2) {
	double sum = 0.0;
	int index;
	int n = (int) ((x2 - x0)/h);
	double x[n+1];
	double y[n+1];
	
	x[0] = x0;
	y[0] = f(x0);
	x[n] = x2;
	y[n] = f(x2);
	
	index = 1;
	for (double i = x0 + h; i < x2; i+=h) {
		x[index] = i;
		y[index] = f(x[index]);
		index++;
	}

	for (int i = 1; i < n; i+=2) {
		sum += 4*y[i];
	}

	for (int i = 2; i < n-1; i+=2) {
		sum += 2*y[i];
	}
	
	sum = (f(x0) + f(x2) + sum)*(h/3.0);
	
	return sum;
}

double midpoint_rule(double h, double (*f)(double), double x0, double x2) {
	double sum = 0.0;
	double mid;
	
	for (double i = x0; i < x2; i+=h) {
		mid = (2*i + h)/2.0;
		sum += f(mid);
	}
	
	return sum*h;
}

double gaussian_quadrature_n2(double (*f)(double), double a, double b) {
	double sum = 0.0;
	double t1 = -1/sqrt(3);
	double t2 = 1/sqrt(3);
	double x1 = ((b-a)*t1 + a + b)/2.0;
	double x2 = ((b-a)*t2 + a + b)/2.0;
	sum = f(x1) + f(x2);
	cout << "sum = " << sum << endl;
	return sum*(b-a)/2.0;
}

double gaussian_quadrature_n3(double (*f)(double), double a, double b) {
	double sum = 0.0;
	double t1 = -sqrt(3.0/5.0);
	double t2 = 0;
	double t3 = sqrt(3.0/5.0);
	double x1 = ((b-a)*t1 + a + b)/2.0;
	double x2 = ((b-a)*t2 + a + b)/2.0;
	double x3 = ((b-a)*t3 + a + b)/2.0;
	sum = (5.0/9.0)*f(x1) + (8.0/9.0)*f(x2) + (5.0/9.0)*f(x3);
	cout << "sum = " << sum << endl;
	return sum*(b - a)/2.0;
}

#endif
