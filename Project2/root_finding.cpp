#ifndef ROOT_FINDING_CPP
#define ROOT_FINDING_CPP

#include "root_finding.h"

double epsilon = 0.0000000001;

// f(guess1) < target and f(guess2) > target
double bisection(double (*f)(double), double target, double guess1, double guess2) {
	cout << "f(guess1) = " << f(guess1) << endl;
	cout << "f(guess2) = " << f(guess2) << endl;
	if (fabs(target - f(guess1)) < epsilon) {
		return guess1;
	}
	if (fabs(target - f(guess2)) < epsilon) {
		return guess2;
	}
	if((f(guess1) > target && f(guess2) > target) || (f(guess1) < target && f(guess2) < target)) {
		cout << "Error: invalid guesses. Target must be between guesses." << endl;
		return NAN;
	} 
	
	double new_estimate = 0.5*(guess1 + guess2);
	if (fabs(f(new_estimate) - target) < epsilon)  {
		return new_estimate;
	}
	else if (f(new_estimate) < target) {
		return bisection(f, target, new_estimate, guess2);
	}
	else {
		return bisection(f, target, guess1, new_estimate);
	}
	
	return 0.0;
}

#endif
