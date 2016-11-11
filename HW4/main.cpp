
#include <iostream> 
#include <math.h>
#include "ode_solver.h"

using namespace std;

double f3_real_solution(double x) {
	return exp(x) - x - 1.0;
}

double f3(double x, double y) {
	return x + y;
}

double u1_real_solution(double t) {
	return (1.0/3.0)*exp(5.0*t) - (1.0/3.0)*exp(-t) + exp(2.0*t);
}

double f1(double t, double u1, double u2) {
	return 3.0*u1 + 2.0*u2 - (2.0*pow(t,2.0) + 1)*exp(2.0*t);
}

double u2_real_solution(double t) {
	return (1.0/3.0)*exp(5.0*t) + (2.0/3.0)*exp(-t) + pow(t, 2.0)*exp(2.0*t);
}

double f2(double t, double u1, double u2) {
	return 4.0*u1 + u2 + (pow(t,2.0) + 2.0*t -4.0)*exp(2.0*t);
}



void print_real_solution(double h, double x_0, double x_max, double (*f)(double)) {
	for (double x_i = x_0; x_i <= x_max; x_i+= h) {
		cout << "y(" << x_i << ") = " << f(x_i) << endl;
	}
}

int main() {
	
	double h = 0.1;
	//double x_0 = 0.0;
	//double x_1 = x_0 + h;
	//double x_2 = x_1 + h;
	//double x_3 = x_2 + h;
	//double x_max = 0.5;
	//double y_0 = 0.0;
	//double y_1 = pow(0.1, 2.0)/ 2.0;
	//double y_2 = pow(x_2, 2.0)/ 2.0;
	//double y_3 = pow(x_3, 2.0)/ 2.0;
	
	//rk4(h, f3, x_0, x_max, y_0);
	//adams_bashforth_two_step(h, f3, x_0, x_max, y_0, y_1);
	//adams_moulton_two_step(h, f3, x_0, x_max, y_0, y_1);
	//improved_eulers(h, f3, x_0, x_max, y_0);
	//adams_moulton_three_step_corrector(h, f3, x_max, x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);
	//print_real_solution(h, x_0, x_max, f3_real_solution);
	h = 0.2;
	double t_0 = 0;
	double t_max = 1.0;
	double u1_0 = 1;
	double u2_0 = 1;
	
	rk4_coupled(h, f1, f2, t_max, t_0, u1_0, u2_0);
	print_real_solution(h, t_0, t_max, u1_real_solution);
	
	//rk4_coupled(h, f1, f2, t_max, t_0, u1_0, u2_0);
	//print_real_solution(h, t_0, t_max, u2_real_solution);
	
	
	
	return 0;
}
