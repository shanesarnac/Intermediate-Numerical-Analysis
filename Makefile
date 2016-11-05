numerical_analysis: main.cpp interpolation.cpp root_finding.cpp data_point.cpp
	g++ -Wall -std=c++11 -g main.cpp interpolation.cpp root_finding.cpp data_point.cpp -lm -o numerical_analysis

clean: 
	rm *.o
