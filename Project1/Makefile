OBJ = ode_solver.o main.o
CC = g++
CFLAGS = -Wall -g 
TARGETS = main.cpp ode_solver.cpp data_point.cpp interpolation.cpp integrate.cpp

analysis: $(TARGETS)
	$(CC) $(CFLAGS) $(TARGETS) -lm -o Numerical_Analysis
	
clean: 
	rm *.o


