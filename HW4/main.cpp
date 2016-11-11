
#include <iostream> 
#include <math.h>
#include "ode_solver.h"

using namespace std;

double real_solution(double x) {
	return exp(x) - x - 1.0;
}

double f1(double x, double y) {
	return x + y;
}

void print_real_solution(double h, double x_0, double x_max) {
	for (double x_i = x_0; x_i <= x_max; x_i+= h) {
		cout << "y(" << x_i << ") = " << real_solution(x_i) << endl;
	}
}

int main() {
	
	double h = 0.1;
	double x_0 = 0.0;
	double x_max = 0.5;
	double y_0 = 0.0;
	double y_1 = pow(0.1, 2.0)/ 2.0;
	
	//rk4(h, f1, x_0, x_max, y_0);
	//adams_bashforth_two_step(h, f1, x_0, x_max, y_0, y_1);
	//adams_moulton_two_step(h, f1, x_0, x_max, y_0, y_1);
	improved_eulers(h, f1, x_0, x_max, y_0);
	print_real_solution(h, x_0, x_max);
	
	return 0;
}
