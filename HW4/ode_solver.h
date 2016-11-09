#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <iostream>
#include <math.h>

using namespace std;

double rk4(double h, double (*f)(double,double), double x_0, double x_max, double y_0);



#endif 
