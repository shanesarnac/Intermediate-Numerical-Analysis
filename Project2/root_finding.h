#ifndef ROOT_FINDING_H
#define ROOT_FINDING_H

#include <iostream>
#include <stdlib.h> 
#include <math.h>
#include "data_point.h"
#include <vector>

using namespace std;

// f(guess1) < target and f(guess2) > target
double bisection(double (*f)(double), double target, double guess1, double guess2);

#endif
