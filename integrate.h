
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <iostream>
#include <math.h>

using namespace std;

double trapezoidal_integration(double h, double (*f)(double), double x0, double x1);
double simpsons_one_third_integration(double h, double (*f)(double), double x0, double x2);
double midpoint_rule(double h, double (*f)(double), double x0, double x2);

double gaussian_quadrature_n2(double (*f)(double), double a, double b);
double gaussian_quadrature_n3(double (*f)(double), double a, double b);

#endif
