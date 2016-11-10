#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <iostream>
#include <math.h>

using namespace std;

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0);

double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1);



#endif 
