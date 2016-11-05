#ifndef INTERPOLATION_CPP
#define INTERPOLATION_CPP

#include "interpolation.h"

using namespace std;

double Interpolation::laGrange(double x, vector<Data_Point> data) {
	cout << "LaGrange Polynomial method: " << endl;
	
	double estimate = 0;
	int num_points = data.size();
	if (num_points < 2) {
		cout << "Error: Must have more than one point to make an estimation" << endl;
		return NAN;
	}
	
	cout << "estimate = ";
	for (int i = 0; i < num_points; i++) {
		Data_Point temp(data[i].getX(), data[i].getY());
		double product = 1;
		for (int k = 0; k < num_points; k++) {
			Data_Point current(data[k].getX(), data[k].getY());
			if (i != k) {
				product *= (x - current.getX())/(temp.getX() - current.getX());
			}
		}
		product *= temp.getY();
		cout << " + " << product;
		estimate += product;
	}
	cout << endl;
	
	
	cout << "The estimate at x = " << x << " is " << setprecision(8) << estimate << " using a degree " << num_points -1 << " degree polynomial " << endl << endl;
	
	
	return estimate;
}

double Interpolation::nevillesMethod(double alpha, vector<Data_Point> data) {
	
	cout << "Neville's method: " << endl;
	double estimate = 0;
	int num_points = data.size();
	if (num_points < 2) {
		cout << "Error: must have two or more points" << endl;
		return -1;
	}
	
	double q_est[num_points][num_points];
	double x_values[num_points];
	
	
	for (int i = 0; i < num_points; i++) {
		q_est[0][i] = data[i].getY();
		x_values[i] = data[i].getX();
		cout << "x" << i << " = " << x_values[i] << " => ";
		cout << "Q(" << i << ",0)(" << alpha <<") = " << q_est[0][i] << endl;
	}
	
	// degree
	for (int i = 1; i < num_points; i++) {
		cout << endl;
		// last point
		for (int j = i; j < num_points; j++) {
			double xj = x_values[j]; // x_j
			double xjm1 = x_values[j-1]; // x_{j - 1}
			q_est[i][j] = q_est[i-1][j]*(alpha - xjm1)/(xj - xjm1) + q_est[i-1][j-1]*(alpha - xj)/(xjm1 - xj);
			estimate = q_est[i][j];
			cout << "x" << j << " = " << x_values[j] << " => ";
			cout << "Q(" << j << "," << i << ")(" << alpha << ") = " << q_est[i][j] << endl;
		}
		
	}
	
	cout << "Neville's method: " << endl;
	cout << "The estimate at x = " << alpha << " is " << setprecision(8) << estimate << " using a degree " << num_points -1 << " degree polynomial " << endl;
	
	return estimate;
}

double Interpolation::newtonsDividedDifferenceMethod(double alpha, vector<Data_Point> data) {
	double estimate = 0;
	int num_points = data.size();
	if (num_points < 2) {
		cout << "Error: must have two or more data points" << endl;
		return -1;
	}
	
	double divided_differences[num_points][num_points];
	double x_values[num_points];
	
	for (int i = 0; i < num_points; i++) {
		divided_differences[0][i] = data[i].getY();
		x_values[i] = data[i].getX();
		cout << "f[x" << i << ",0] = " << divided_differences[0][i] << endl;
	}
	
	for (int j = 1; j < num_points; j++) {
		cout << endl;
		for (int k = j; k < num_points; k++) {
			divided_differences[j][k] = (divided_differences[j-1][k] - divided_differences[j-1][k-1])/(x_values[k]-x_values[k - j]);
			cout << "f[x" << j << "," << k << "] = " << divided_differences[j][k] << endl;
		}
		
	}
		
	for (int i = 0; i < num_points; i++) {
		double product = 1;
		for (int j = 0; j < i; j++) {
			product *= (alpha - x_values[j]);
		}
		estimate += divided_differences[i][i] * product;
	}	
	
	cout << "f(" << alpha << ") ~ " << estimate << endl;
		
	return estimate;
}

vector<vector<double>> Interpolation::naturalCubicSpline(vector<Data_Point> data) {
	//vector<double[]> coefficients;
	int num_points = data.size();
	double a_n[num_points];
	double h_n[num_points];
	double alpha_n[num_points];
	
	for (int i = 0; i < num_points; i++) {
		a_n[i] = data[i].getY();
	}
	
	for (int i = 0; i < num_points; i++) {
		h_n[i] = data[i+1].getX() - data[i].getX();
	}
	
	for (int i = 1; i < num_points; i++) {
		alpha_n[i] = (3 / h_n[i]) *(a_n[i + 1] - a_n[i]) - (3 / h_n[i - 1])*(a_n[i] - a_n[i - 1]);
	}
	
	double l_n[num_points + 1];
	double mu_n[num_points + 1];
	double z_n[num_points + 1];
	double c_n[num_points + 1];
	double b_n[num_points];
	double d_n[num_points];
	
	l_n[0] = 1;
	mu_n[0] = 0;
	z_n[0] = 0;
	
	for (int i = 1; i < num_points; i++) {
		l_n[i] = 2*(data[i+1].getX() - data[i-1].getX()) - h_n[i - 1] * mu_n[i-1];
		mu_n[i] = h_n[i] / l_n[i];
		z_n[i] = (alpha_n[i] - h_n[i-1] * z_n[i-1]) / l_n[i];
	}
	
	l_n[num_points] = 1;
	z_n[num_points] = 0;
	c_n[num_points] = 0;
	
	
	for (int j = num_points - 1; j >= 0; j--) {
		c_n[j] = z_n[j] - mu_n[j]* c_n[j+1];
		b_n[j] = (a_n[j+1] - a_n[j])/h_n[j] - h_n[j] * (c_n[j+1] + 2*c_n[j]) / 3;
		d_n[j] = (c_n[j+1] - c_n[j])/ (3*h_n[j]);
	}
	
	vector<double> a_list;
	vector<double> b_list;
	vector<double> c_list;
	vector<double> d_list;
	
	for (int i = 0; i < num_points; i++) {
		a_list.push_back(a_n[i]);
		b_list.push_back(b_n[i]);
		c_list.push_back(c_n[i]);
		d_list.push_back(d_n[i]);
		
		cout << "a" << i << " = " << a_n[i] << endl;
		cout << "b" << i << " = " << b_n[i] << endl;
		cout << "c" << i << " = " << c_n[i] << endl;
		cout << "d" << i << " = " << d_n[i] << endl;
	}
	
	vector<vector<double>> coefficients;
	coefficients.push_back(a_list);
	coefficients.push_back(b_list);
	coefficients.push_back(c_list);
	coefficients.push_back(d_list);
	//coefficients.push_back(a_n);
	//coefficients.push_back(b_n);
	//coefficients.push_back(c_n);
	//coefficients.push_back(d_n);
	
	return coefficients;
	
}

vector<vector<double>> Interpolation::clampedCubicSpline(vector<Data_Point> data, double fp0, double fpn) {
	int num_points = data.size();
	double a_n[num_points];
	double h_n[num_points];
	double alpha_n[num_points];
	
	for (int i = 0; i < num_points; i++) {
		a_n[i] = data[i].getY();
	}
	
	for (int i = 0; i < num_points; i++) {
		h_n[i] = data[i+1].getX() - data[i].getX();
	}
	
	alpha_n[0] = 3*(a_n[1] - a_n[0])/ h_n[0] - 3 - fp0;
	alpha_n[num_points-1] = 3*fpn - 3*(a_n[num_points-1] - a_n[num_points - 2])/ h_n[num_points - 2];
	
	for (int i = 1; i < num_points; i++) {
		alpha_n[i] = (3 / h_n[i])*(a_n[i + 1] - a_n[i]) - (3 / h_n[i-1])*(a_n[i] - a_n[i - 1]);
	}
	
	double l_n[num_points + 1];
	double mu_n[num_points + 1];
	double z_n[num_points + 1];
	double c_n[num_points + 1];
	double b_n[num_points];
	double d_n[num_points];
	
	l_n[0] = 2*h_n[0];
	mu_n[0] = 0.5;
	z_n[0] = alpha_n[0] / l_n[0];
	
	for (int i = 1; i < num_points; i++) {
		l_n[i] = 2*(data[i+1].getX() - data[i - 1].getX()) - h_n[i - 1]*mu_n[i-1];
		mu_n[i] = h_n[i] / l_n[i];
		z_n[i] = (alpha_n[i] - h_n[i-1]*z_n[i-1]) / l_n[i];
	}
	
	l_n[num_points] = h_n[num_points-1] * (2 - mu_n[num_points - 1]);
	z_n[num_points] = (alpha_n[num_points] - h_n[num_points - 1] * z_n[num_points - 1]);
	c_n[num_points] = z_n[num_points];
	
	for (int j = num_points - 1; j >= 0; j--) {
		c_n[j] = z_n[j] - mu_n[j]*c_n[j+1];
		b_n[j] = (a_n[j+1] - a_n[j])/h_n[j] - h_n[j] * (c_n[j+1] + 2*c_n[j]) / 3;
		d_n[j] = (c_n[j+1] - c_n[j]) / (3*h_n[j]);
	}
	
	vector<double> a_list;
	vector<double> b_list;
	vector<double> c_list;
	vector<double> d_list;
	
	for (int i = 0; i < num_points; i++) {
		a_list.push_back(a_n[i]);
		b_list.push_back(b_n[i]);
		c_list.push_back(c_n[i]);
		d_list.push_back(d_n[i]);
		
		cout << "a" << i << " = " << a_n[i] << endl;
		cout << "b" << i << " = " << b_n[i] << endl;
		cout << "c" << i << " = " << c_n[i] << endl;
		cout << "d" << i << " = " << d_n[i] << endl;
	}
	
	vector<vector<double>> coefficients;
	coefficients.push_back(a_list);
	coefficients.push_back(b_list);
	coefficients.push_back(c_list);
	coefficients.push_back(d_list);
	
	return coefficients;
}



#endif
