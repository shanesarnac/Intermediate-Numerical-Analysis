#ifndef DATA_POINT_H
#define DATA_POINT_H

#include <iostream>

using namespace std;

class Data_Point {
public:
	Data_Point(double X, double Y);
	//Data_Point(const Data_Point* data);
	void setX(double X);
	void setY(double Y);
	double getX();
	double getY();

private:
	double x;
	double y;
	
};

#endif
