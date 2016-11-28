#ifndef ROOT_FINDING_H
#define ROOT_FINDING_H

#include <iostream>
#include <math.h>
#include "data_point.h"
#include <vector>

using namespace std;

// guess1 < target and guess2 > target
double bisection(double (*f)(double), double target, double guess1, double guess2);

#endif
