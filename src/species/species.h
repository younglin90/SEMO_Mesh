#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

class SEMO_Species{
	public:
		string name;
		string phase;
		double rho;
		double mu;
		double sigma;
		
		string rhoType;
		double Pinf;
		double cv;
		double gamma;
		double b;
		double q;
		

};

// class SEMO_Species_Builder{
	// public:
		// // vector<string> name;
		// // vector<string> phase;
		// // vector<double> rho;
		// // vector<double> mu;
		
// };

