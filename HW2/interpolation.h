#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

#include "data_point.h"

class Interpolation {
public:
	double laGrange(double x, vector<Data_Point> data);
	double nevillesMethod(double alpha, vector<Data_Point> data);
	double newtonsDividedDifferenceMethod(double alpha, vector<Data_Point> data);
	double evaluateSpline(vector<double> coefficients);
	vector<vector<double>> naturalCubicSpline(vector<Data_Point> data);
	vector<vector<double>> clampedCubicSpline(vector<Data_Point> data, double fp0, double fpn);

private:
	
};

#endif 
