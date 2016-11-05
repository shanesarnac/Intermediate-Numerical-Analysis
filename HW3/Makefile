OBJ = integrate.o main.o
CC = g++
CFLAGS = -Wall -g 
TARGETS = main.cpp integrate.cpp  differentiation.cpp

analysis: $(TARGETS)
	$(CC) $(CFLAGS) $(TARGETS) -lm -o Numerical_Analysis
	
clean: 
	rm *.o
