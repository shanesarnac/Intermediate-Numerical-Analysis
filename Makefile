OBJ = integrate.o main.o
CC = g++
CFLAGS = -Wall -g 

analysis: main.cpp integrate.cpp 
	$(CC) $(CFLAGS) main.cpp -lm -o Numerical_Analysis
	
clean: 
	rm *.o
