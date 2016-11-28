#ifndef ROOT_FINDING_CPP
#define ROOT_FINDING_CPP

#include "root_finding.h"

double bisection(double (*f)(double), double target, double guess1, double guess2) {
	if (target == f(guess1)) {
		return guess1;
	}
	if (target == f(guess2)) {
		return guess2;
	}
	if((f(guess1) > target && f(guess2) > target) || (f(guess1) < target && f(guess2) < target)) {
		cout << "Error: invalid guesses. Target must be between guesses." << endl;
		return NAN;
	} 
	
	double new_estimate = 0.5*(guess1 + guess2);
	if (f(new_estimate) == target)  {
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
