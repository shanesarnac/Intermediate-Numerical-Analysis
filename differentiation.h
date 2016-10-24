
#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H

double three_point_foward(double h, double f_x0, double f_x1, double f_x2);

double three_point_middle(double h, double f_mh, double f_ph);

double three_point_backward(double h, double f_x0, double f_x1, double f_x2);

#endif
