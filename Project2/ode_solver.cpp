#ifndef ODE_SOLVER_CPP
#define ODE_SOLVER_CPP

#include "ode_solver.h"

vector<Data_Point> rk4(double h, double (*f)(double,double), Data_Point initial_data, double x_max) {
	//cout << "RK4" << endl;
	vector<Data_Point> estimates;
	
	double x_0 = initial_data.getX();
	double y_0 = initial_data.getY();
	double y_ip1 = y_0;
	double y_i;
	estimates.push_back(Data_Point(x_0, y_0));
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_i = y_ip1;
		double k1 = h*f(x_i, y_i);
		double k2 = h*f(x_i + h/2.0, y_i + k1/2.0);
		double k3 = h*f(x_i + h/2.0, y_i + k2/2.0);
		double k4 = h*f(x_i + h, y_i + k3);
		y_ip1 = y_i + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
		estimates.push_back(Data_Point(x_i + h, y_ip1));
	}
	
	return estimates;
}

double rk4_2_coupled(double h, double (*f1)(double,double,double), double (*f2)(double, double, double), double t_max, double t_0, double x_0, double y_0) {
	cout << "RK4 2 Coupled" << endl;
	double x_ip1 = x_0;
	double y_ip1 = y_0;
	double x_i, y_i;
	
	cout << "x(" << t_0 << ") = " << x_0 << endl;
	//cout << "y(" << t_0 << ") = " << y_0 << endl;
	
	for(double t_i = t_0; t_i < t_max; t_i += h) {
		x_i = x_ip1;
		y_i = y_ip1;
		
		double k0 = h*f1(t_i, x_i, y_i);
		double l0 = h*f2(t_i, x_i, y_i);
		double k1 = h*f1(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0); 
		double l1 = h*f2(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0); 
		double k2 = h*f1(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1); 
		double l2 = h*f2(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1); 
		double k3 = h*f1(t_i + h, x_i + k2, y_i + l2); 
		double l3 = h*f2(t_i + h, x_i + k2, y_i + l2); 
		
		//cout << "k0 = " << k0 << endl; 
		//cout << "k1 = " << k1 << endl; 
		//cout << "k2 = " << k2 << endl; 
		//cout << "k3 = " << k3 << endl; 
		
		x_ip1 = x_i + (1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3);
		y_ip1 = y_i + (1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3);
		
		cout << "x(" << t_i + h << ") = " << x_ip1 << endl;
		//cout << "y(" << t_i + h << ") = " << y_ip1 << endl;
	} 
	cout << endl;
	
	return 0.0;
}

vector<vector<double> > rk4_3_coupled(double h, double (*f1)(double,double,double,double), double (*f2)(double, double, double, double),
						double (*f3)(double, double, double, double),double max, double t_0, double x_0, double y_0, double z_0) {
	//cout << "RK4 3 Coupled" << endl;
	
	
	vector<vector<double> > function_values;
	for (int i = 0; i < 4; i++) {
		vector<double> initializer;
		function_values.push_back(initializer);
		function_values[i].push_back(0.0);
	}
	function_values[0].push_back(t_0);
	function_values[1].push_back(x_0);
	function_values[2].push_back(y_0);
	function_values[3].push_back(z_0);
	
	
	double x_ip1 = x_0;
	double y_ip1 = y_0;
	double z_ip1 = z_0;
	double x_i, y_i, z_i;
	
	//cout << "mu,F,F',F''" << endl;
	//cout << t_0 << "," << x_0 << "," << y_0 << "," << z_0 << endl;
	//cout << "x(" << t_0 << ") = " << x_0 << endl;
	//cout << "y(" << t_0 << ") = " << y_0 << endl;
	//cout << "z(" << t_0 << ") = " << z_0 << endl;
	
	for(double t_i = t_0; t_i < max; t_i += h) {
		x_i = x_ip1;
		y_i = y_ip1;
		z_i = z_ip1;
		
		double k0 = h*f1(t_i, x_i, y_i, z_i);
		double l0 = h*f2(t_i, x_i, y_i, z_i);
		double m0 = h*f3(t_i, x_i, y_i, z_i);
		
		double k1 = h*f1(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		double l1 = h*f2(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		double m1 = h*f3(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		
		double k2 = h*f1(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		double l2 = h*f2(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		double m2 = h*f3(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		
		double k3 = h*f1(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		double l3 = h*f2(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		double m3 = h*f3(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		
		//cout << "k0 = " << k0 << endl; 
		//cout << "k1 = " << k1 << endl; 
		//cout << "k2 = " << k2 << endl; 
		//cout << "k3 = " << k3 << endl; 
		
		x_ip1 = x_i + (1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3);
		y_ip1 = y_i + (1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3);
		z_ip1 = z_i + (1.0/6.0)*(m0 + 2.0*m1 + 2.0*m2 + m3);
		
		function_values[0].push_back(t_i + h); 
		function_values[1].push_back(x_ip1); 
		function_values[2].push_back(y_ip1); 
		function_values[3].push_back(z_ip1); 
		
		
		//cout << t_i + h << "," << x_ip1 << "," << y_ip1 << "," << z_ip1 << endl;
		
		//cout << "x(" << t_i + h << ") = " << x_ip1 << endl;
		//cout << "y(" << t_i + h << ") = " << y_ip1 << endl;
		//cout << "z(" << t_i + h << ") = " << z_ip1 << endl;
	} 
	cout << endl;
	
	return function_values;			
							
}

vector<vector<double> > rk4_5_coupled(double h, double (*f1)(double,double,double,double, double, double), double (*f2)(double,double,double,double, double, double),
						double (*f3)(double,double,double,double, double, double), double (*f4)(double,double,double,double, double, double), double (*f5)(double,double,double,double, double, double),
						double max, double t_0, double x1_0, double x2_0, double x3_0, double x4_0, double x5_0) {
	//cout << "RK4 3 Coupled" << endl;
	
	
	vector<vector<double> > function_values;
	for (int i = 0; i < 6; i++) {
		vector<double> initializer;
		function_values.push_back(initializer);
		function_values[i].push_back(0.0);
	}
	function_values[0].push_back(t_0);
	function_values[1].push_back(x1_0);
	function_values[2].push_back(x2_0);
	function_values[3].push_back(x3_0);
	function_values[4].push_back(x4_0);
	function_values[5].push_back(x5_0);
	
	
	double x1_ip1 = x1_0;
	double x2_ip1 = x2_0;
	double x3_ip1 = x3_0;
	double x4_ip1 = x4_0;
	double x5_ip1 = x5_0;
	double x1_i, x2_i, x3_i, x4_i, x5_i;
	
	//cout << "mu,F,F',F''" << endl;
	//cout << t_0 << "," << x_0 << "," << y_0 << "," << z_0 << endl;
	//cout << "x(" << t_0 << ") = " << x_0 << endl;
	//cout << "y(" << t_0 << ") = " << y_0 << endl;
	//cout << "z(" << t_0 << ") = " << z_0 << endl;
	
	for(double t_i = t_0; t_i < max; t_i += h) {
		x1_i = x1_ip1;
		x2_i = x2_ip1;
		x3_i = x3_ip1;
		x4_i = x4_ip1;
		x5_i = x5_ip1;

		
		double k0 = h*f1(t_i, x1_i, x2_i, x3_i, x4_i, x5_i);
		double l0 = h*f2(t_i, x1_i, x2_i, x3_i, x4_i, x5_i);
		double m0 = h*f3(t_i, x1_i, x2_i, x3_i, x4_i, x5_i);
		double n0 = h*f4(t_i, x1_i, x2_i, x3_i, x4_i, x5_i);
		double p0 = h*f5(t_i, x1_i, x2_i, x3_i, x4_i, x5_i);

		
		double k1 = h*f1(t_i + 0.5*h, x1_i + 0.5*k0, x2_i + 0.5*l0, x3_i + 0.5*m0, x4_i + 0.5*n0, x5_i + 0.5*p0); 
		double l1 = h*f2(t_i + 0.5*h, x1_i + 0.5*k0, x2_i + 0.5*l0, x3_i + 0.5*m0, x4_i + 0.5*n0, x5_i + 0.5*p0); 
		double m1 = h*f3(t_i + 0.5*h, x1_i + 0.5*k0, x2_i + 0.5*l0, x3_i + 0.5*m0, x4_i + 0.5*n0, x5_i + 0.5*p0); 
		double n1 = h*f4(t_i + 0.5*h, x1_i + 0.5*k0, x2_i + 0.5*l0, x3_i + 0.5*m0, x4_i + 0.5*n0, x5_i + 0.5*p0); 
		double p1 = h*f5(t_i + 0.5*h, x1_i + 0.5*k0, x2_i + 0.5*l0, x3_i + 0.5*m0, x4_i + 0.5*n0, x5_i + 0.5*p0); 
 

		
		double k2 = h*f1(t_i + 0.5*h, x1_i + 0.5*k1, x2_i + 0.5*l1, x3_i + 0.5*m1, x4_i + 0.5*n1, x5_i + 0.5*p1); 
		double l2 = h*f2(t_i + 0.5*h, x1_i + 0.5*k1, x2_i + 0.5*l1, x3_i + 0.5*m1, x4_i + 0.5*n1, x5_i + 0.5*p1); 
		double m2 = h*f3(t_i + 0.5*h, x1_i + 0.5*k1, x2_i + 0.5*l1, x3_i + 0.5*m1, x4_i + 0.5*n1, x5_i + 0.5*p1); 
		double n2 = h*f4(t_i + 0.5*h, x1_i + 0.5*k1, x2_i + 0.5*l1, x3_i + 0.5*m1, x4_i + 0.5*n1, x5_i + 0.5*p1); 
		double p2 = h*f5(t_i + 0.5*h, x1_i + 0.5*k1, x2_i + 0.5*l1, x3_i + 0.5*m1, x4_i + 0.5*n1, x5_i + 0.5*p1); 
		
		
		double k3 = h*f1(t_i + h, x1_i + k2, x2_i + l2, x3_i + m2, x4_i + n2, x5_i + p2); 
		double l3 = h*f2(t_i + h, x1_i + k2, x2_i + l2, x3_i + m2, x4_i + n2, x5_i + p2); 
		double m3 = h*f3(t_i + h, x1_i + k2, x2_i + l2, x3_i + m2, x4_i + n2, x5_i + p2); 
		double n3 = h*f4(t_i + h, x1_i + k2, x2_i + l2, x3_i + m2, x4_i + n2, x5_i + p2); 
		double p3 = h*f5(t_i + h, x1_i + k2, x2_i + l2, x3_i + m2, x4_i + n2, x5_i + p2); 

		
		//cout << "k0 = " << k0 << endl; 
		//cout << "k1 = " << k1 << endl; 
		//cout << "k2 = " << k2 << endl; 
		//cout << "k3 = " << k3 << endl; 
		
		x1_ip1 = x1_i + (1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3);
		x2_ip1 = x2_i + (1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3);
		x3_ip1 = x3_i + (1.0/6.0)*(m0 + 2.0*m1 + 2.0*m2 + m3);
		x4_ip1 = x4_i + (1.0/6.0)*(n0 + 2.0*n1 + 2.0*n2 + n3);
		x5_ip1 = x5_i + (1.0/6.0)*(p0 + 2.0*p1 + 2.0*p2 + p3);
		
		function_values[0].push_back(t_i + h); 
		function_values[1].push_back(x1_ip1); 
		function_values[2].push_back(x2_ip1); 
		function_values[3].push_back(x3_ip1); 
		function_values[4].push_back(x4_ip1); 
		function_values[5].push_back(x5_ip1); 
		
		
		//cout << t_i + h << "," << x_ip1 << "," << y_ip1 << "," << z_ip1 << endl;
		
		//cout << "x(" << t_i + h << ") = " << x_ip1 << endl;
		//cout << "y(" << t_i + h << ") = " << y_ip1 << endl;
		//cout << "z(" << t_i + h << ") = " << z_ip1 << endl;
	} 
	//cout << endl;
	
	return function_values;			
							
}

vector<vector<double> > rk4_coupled_general(double h, double initial_x, double max, vector<double (*)(double,double,double,double)> functions, vector<double> initial_values) {
	//cout << "RK4 Coupled General" << endl;
	
	vector<vector<double> > function_values;
	vector<double> f_i, f_ip1, k0_values, k1_values, k2_values, k3_values;
	
	if (functions.size() != initial_values.size()) {
		cout << "Error: number of functions does not match number of initial values" << endl;
		return function_values;
	}
	
	for (unsigned int i = 0; i < initial_values.size(); i++) {
		vector<double> initializer;
		function_values.push_back(initializer);
		function_values[i].push_back(initial_values[i]);
		
		f_ip1.push_back(initial_values[i]);
		f_i.push_back(0.0); // initialize f_i values to 0
		k0_values.push_back(0.0);
		k1_values.push_back(0.0);
		k2_values.push_back(0.0);
		k3_values.push_back(0.0);
	}
	
	//double x_ip1 = x_0;
	//double y_ip1 = y_0;
	//double z_ip1 = z_0;
	//double x_i, y_i, z_i;
	
	//cout << "mu,F,F',F''" << endl;
	//cout << t_0 << "," << x_0 << "," << y_0 << "," << z_0 << endl;
	////cout << "x(" << t_0 << ") = " << x_0 << endl;
	////cout << "y(" << t_0 << ") = " << y_0 << endl;
	////cout << "z(" << t_0 << ") = " << z_0 << endl;
	
	for(double x_i = initial_x; x_i < max; x_i += h) {
		// Reset values for f_i
		for (unsigned int i = 0; i < initial_values.size(); i++) {
			f_i[i] = f_ip1[i];
		} 
		
		// Reset k0_values
		for (unsigned int i = 0; i < initial_values.size(); i++) {
			//k0_values[i] = h*functions[i](x_i, ;
		} 
		//x_i = x_ip1;
		//y_i = y_ip1;
		//z_i = z_ip1;
		
		//double k0 = h*f1(t_i, x_i, y_i, z_i);
		//double l0 = h*f2(t_i, x_i, y_i, z_i);
		//double m0 = h*f3(t_i, x_i, y_i, z_i);
		
		//double k1 = h*f1(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		//double l1 = h*f2(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		//double m1 = h*f3(t_i + 0.5*h, x_i + 0.5*k0, y_i + 0.5*l0, z_i + 0.5*m0); 
		
		//double k2 = h*f1(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		//double l2 = h*f2(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		//double m2 = h*f3(t_i + 0.5*h, x_i + 0.5*k1, y_i + 0.5*l1, z_i + 0.5*m1); 
		
		//double k3 = h*f1(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		//double l3 = h*f2(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		//double m3 = h*f3(t_i + h, x_i + k2, y_i + l2, z_i + m2); 
		
		////cout << "k0 = " << k0 << endl; 
		////cout << "k1 = " << k1 << endl; 
		////cout << "k2 = " << k2 << endl; 
		////cout << "k3 = " << k3 << endl; 
		
		//x_ip1 = x_i + (1.0/6.0)*(k0 + 2.0*k1 + 2.0*k2 + k3);
		//y_ip1 = y_i + (1.0/6.0)*(l0 + 2.0*l1 + 2.0*l2 + l3);
		//z_ip1 = z_i + (1.0/6.0)*(m0 + 2.0*m1 + 2.0*m2 + m3);
		
		//cout << t_i + h << "," << x_ip1 << "," << y_ip1 << "," << z_ip1 << endl;
		
		////cout << "x(" << t_i + h << ") = " << x_ip1 << endl;
		////cout << "y(" << t_i + h << ") = " << y_ip1 << endl;
		////cout << "z(" << t_i + h << ") = " << z_ip1 << endl;
	} 
	//cout << endl;
	
	return function_values;			
							
}

double adams_bashforth_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1){
	cout << "Adams-Bashforth Two-Step" << endl;
	double acc = 0.0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		/*cout << "x_i = " << x_i << endl;
		cout << "y_i = " << y_i << endl;
		cout << "y_{i  -1} = " << y_im1 << endl;*/
		
		y_ip1 = y_i + (h/2.0)*(3*f(x_i, y_i) - f(x_i - h, y_im1));
		y_im1 = y_i;
		y_i = y_ip1;

		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

// Note: Only works for y' = y + x (had to solve for y_ip1)
double adams_moulton_two_step(double h, double (*f)(double, double), double x_0, double x_max, double y_0, double y_1) {
	cout << "Adams-Moulton Two Step" << endl;
	double acc = 0.0;
	double x_im1 = x_0;
	double y_i = y_1;
	double y_im1 = y_0;
	double y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_0 + h << ") = " << y_1 << endl;
	
	for (double x_i = x_0 + h; x_i + h <= x_max; x_i += h) {
		y_ip1 = (y_i + (h/12.0) * (5.0*(x_i + h) + 8.0*f(x_i, y_i) - f(x_im1, y_im1)))*pow(1.0 - (5.0*h)/12.0, -1.0);
		y_im1 = y_i;
		x_im1 = x_i;
		y_i = y_ip1;
		
		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 2.0) << endl;
	}
	cout << endl;
	
	return acc;
}

double improved_eulers(double h, double (*f)(double, double), double x_0, double x_max, double y_0) {
	cout << "Improved Euler's Method" << endl;
	double y_i = y_0;
	double y_bar, y_ip1;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	
	for (double x_i = x_0; x_i < x_max; x_i += h) {
		y_bar = y_i + h*f(x_i, y_i);
		y_ip1 = y_i + (h/2.0)*(f(x_i,y_i) + f(x_i + h, y_bar));
		y_i = y_ip1;
		//cout << "y bar = " << y_bar << endl;
		cout << "y(" << x_i + h<< ") = " << y_ip1 << endl;
	}
	cout << endl;
	
	return 0.0;
}

double adams_bashforth_four_step_predictor(double h, double (*f)(double, double), double x_i, double y_i, double x_im1, double y_im1, double x_im2, double y_im2, double x_im3, double y_im3) {
	return y_i + (h/24.0)*(55.0*f(x_i, y_i) - 59.0*f(x_im1, y_im1) + 37.0*f(x_im2, y_im2) - 9.0*f(x_im3, y_im3));
}

double adams_moulton_three_step_corrector(double h, double (*f)(double, double), double x_max, double x_0, double y_0, double x_1, double y_1, double x_2, double y_2, double x_3, double y_3) {
	cout << "Adams-Moulton Three Step" << endl;
	double acc = 0.0;
	double x_im3 = x_0;
	double x_im2 = x_1;
	double x_im1 = x_2;
	double y_im3 = y_0;
	double y_im2 = y_1;
	double y_im1 = y_2;
	double y_i = y_3;
	double y_ip1, y_bar;
	
	cout << "y(" << x_0 << ") = " << y_0 << endl;
	cout << "y(" << x_1 << ") = " << y_1 << endl;
	cout << "y(" << x_2 << ") = " << y_2 << endl;
	
	for (double x_i = x_3; x_i + h <= x_max; x_i += h) {
		y_bar = adams_bashforth_four_step_predictor(h, f, x_i, y_i, x_im1, y_im1, x_im2, y_im2, x_im3, y_im3);
		y_ip1 = y_i + (h/24.0)*(9.0*f(x_i + h, y_bar) + 19.0*f(x_i, y_i) -5.0*f(x_im1, y_im1) + f(x_im2, y_im2));
		y_im3 = y_im2;
		y_im2 = y_im1;
		y_im1 = y_i;
		y_i = y_ip1;
		
		x_im3 = x_im2;
		x_im2 = x_im1;
		x_im1 = x_i;
		
		
		cout << "y(" << x_i + h<< ") = " << y_ip1 << " +- " << pow(h, 5.0) << endl;
	}
	cout << endl;
	
	return acc;
	
}



#endif
