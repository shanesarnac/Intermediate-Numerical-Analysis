
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
	double x_1 = x_0 + h;
	double x_2 = x_1 + h;
	double x_3 = x_2 + h;
	double x_max = 0.5;
	double y_0 = 0.0;
	double y_1 = pow(0.1, 2.0)/ 2.0;
	double y_2 = pow(x_2, 2.0)/ 2.0;
	double y_3 = pow(x_3, 2.0)/ 2.0;
	
	
	//rk4(h, f1, x_0, x_max, y_0);
	//adams_bashforth_two_step(h, f1, x_0, x_max, y_0, y_1);
	//adams_moulton_two_step(h, f1, x_0, x_max, y_0, y_1);
	//improved_eulers(h, f1, x_0, x_max, y_0);
	adams_moulton_three_step_corrector(h, f1, x_max, x_0, y_0, x_1, y_1, x_2, y_2, x_3, y_3);
	print_real_solution(h, x_0, x_max);
	
	return 0;
}
