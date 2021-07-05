#pragma once
#include <iostream>

#include <vector>
#include "../species/species.h"
// #include "../mesh/build.h"
// #include "../controls/build.h"

class SEMO_Mesh_Builder;
class SEMO_Controls_Builder;

class SEMO_Mesh_Save {
private:

public:
    SEMO_Mesh_Save() {
    }
    ~SEMO_Mesh_Save() {
    }
	
	void vtu(SEMO_Mesh_Builder &in);
	
	void vtu(
		string folder, 
		SEMO_Mesh_Builder &in, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species);
		
	void vtu(string folder, int rank, SEMO_Mesh_Builder &in);
	void vtuZlib(SEMO_Mesh_Builder &in, SEMO_Controls_Builder &controls);
	void gnuplot(int iter, double dtime, vector<double>& norm);

};


