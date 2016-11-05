#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <vector>

#include "root_finding.h"
#include "interpolation.h"
#include "data_point.h"

using namespace std;

int main(int argc, const char * argv[]) {
	Interpolation interp;
	vector<Data_Point> data;
	double alpha = 8.4;
	data.push_back(Data_Point(8.1,16.94410));
	data.push_back(Data_Point(8.3, 17.56492));
	data.push_back(Data_Point(8.6, 18.50515));
	data.push_back(Data_Point(8.7, 18.82091));
	
	//interp.laGrange(alpha, data);
	//interp.nevillesMethod(alpha, data);
	//interp.newtonsDividedDifferenceMethod(alpha, data);
	
	vector<Data_Point> data1;
	data1.push_back(Data_Point(-3, 81));
	data1.push_back(Data_Point(-1,1));
	data1.push_back(Data_Point(1,1));
	data1.push_back(Data_Point(3,81));
	
	double fp0 = -81;
	double fpn = 81;
	
	//interp.clampedCubicSpline(data1, fp0, fpn);
	//cout << endl;
	
	vector<Data_Point> data2;
	data2.push_back(Data_Point(-1, 1));
	data2.push_back(Data_Point(0,0));
	data2.push_back(Data_Point(1,1));
	
	fp0 = -3;
	fpn = 3;
	
	//interp.clampedCubicSpline(data2, fp0, fpn);
	
	//cout << endl;
	
	interp.naturalCubicSpline(data2);
	
}


#endif
