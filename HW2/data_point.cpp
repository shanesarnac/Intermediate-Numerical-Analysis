#ifndef DATA_POINT_CPP
#define DATA_POINT_CPP

#include "data_point.h"

using namespace std;

Data_Point::Data_Point(double X, double Y) {
	x = X;
	y = Y;
}

//Data_Point::Data_Point(const Data_Point* data) {
	//x = data->getX();
	//y = data->getY();
//}

void Data_Point::setX(double X) {
	x = X;
}

void Data_Point::setY(double Y) {
	y = Y;
}

double Data_Point::getX() {
	return x;
}

double Data_Point::getY() {
	return y;
}

#endif
