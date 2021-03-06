#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

#include <iostream>
#include <math.h>
#include "data_point.h"
#include <vector>

using namespace std;

vector<Data_Point> rk4(double h, double (*f)(double,double), Data_Point initial_data, double x_max);

double rk4_2_coupled(double h, double (*f1)(double,double,double), double (*f2)(double, double, double), double t_max, double t_0, double x_0, double y_0);

vector<vector<double> > rk4_3_coupled(double h, double (*f1)(double,double,double,double), double (*f2)(double, double, double, double),
						double (*f3)(double, double, double, double),double max, double t_0, double x_0, double y_0, double z_0);
						
vector<vector<double> > rk4_5_coupled(
		double h, 
		double (*f1)(double,double,double,double, double, double), 
		double (*f2)(double,double,double,double, double, double),
		double (*f3)(double,double,double,double, double, double), 
		double (*f4)(double,double,double,double, double, double), 
		double (*f5)(double,double,double,double, double, double),
		double max, 
		double t_0, 
		double x1_0, 
		double x2_0, 
		double x3_0, 
		double x4_0, 
		double x5_0);
						
vector<vector<double> > rk4_coupled_general(double h, double max, vector<double (*)(double,double,double,double)> functions, vector<double> initial_values);




double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1);

double adams_bashforth_four_step_predictor(double h, double (*f)(double, double), double x_i, double y_i, double x_im1, double y_im1, 
																				double x_im2, double y_im2, double x_im3, double y_im3);

double adams_moulton_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1);

double adams_moulton_three_step_corrector(double h, double (*f)(double, double), double x_max, double x_0, double y_0, 
																double x_1, double y_1, double x_2, double y_2, double x_3, double y_3);



double improved_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0);



#endif 
