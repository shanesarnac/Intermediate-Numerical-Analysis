
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
	
	double h = 0.2;
	double x_0 = 0.0;
	double x_max = 1.0;
	double y_0 = 0.0;
	
	rk4(h, f1, x_0, x_max, y_0);
	print_real_solution(h, x_0, x_max);
	
	return 0;
}
