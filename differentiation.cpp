
#ifndef DIFFERENTIATION_CPP
#define DIFFERENTIATION_CPP

double three_point_foward(double h, double f_x0, double f_x1, double f_x2) {
	// Use the two points immediately after x0 to estimate the slope at x0
	double f_prime = (1/(2*h))*(-3*f_x0 + 4*f_x1 -f_x2);
	return f_prime;	
}

double three_point_backward(double h, double f_x0, double f_x1, double f_x2) {
	// Use the two points immediately before x2 to estimate the slope at x2
	double f_prime = (1/h)*(0.5*f_x0 - 2*f_x1 + 1.5*f_x2);
	return f_prime;	
}

double three_point_middle(double h, double f_mh, double f_ph) {
	// Use the the point immediately before and immediately after x0 to find the slope at x0
	double f_prime = (1/(2*h))*(-f_mh + f_ph);
	return f_prime;
}



#endif
