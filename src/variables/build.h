#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

#include "../mesh/build.h"

class SEMO_Variables_Builder{
	public:
		void initialization(SEMO_Mesh_Builder& mesh);
		
		double timestep;
		vector<double> pressure;
		vector<vector<double>> velocities;
		vector<double> temperature;
		vector<vector<double>> volFracs;
		vector<vector<double>> massFracs;
		
};

