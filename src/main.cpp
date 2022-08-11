#include "../include/linprog.h"
#include "linprog.cpp"
#include <iostream>
#include <string>
#include <unistd.h>

using namespace std;

int main(int argc, char* argv[]){
	Solver solver;
	int o;
	const char* optstring = "i:o:";
	string input_file, output_file;
	while((o = getopt(argc, argv, optstring)) != -1){
		switch(o){
			case 'i':{input_file = optarg;	break;}
			case 'o':{output_file = optarg;	break;}
		}
	}
	solver.parse_input(input_file);
	if (solver.check_b()){
		solver.solve(output_file);
		solver.report();
	}
	else{
		cout << "Requires 2-stage simplex algorithm." << endl;
	}
	return 0;
}
