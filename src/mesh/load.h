#pragma once
#include <iostream>
#include "../species/species.h" 

// #include "build.h"
class SEMO_Mesh_Builder;
class SEMO_Controls_Builder;

class SEMO_Mesh_Load {
private:

public:
    SEMO_Mesh_Load() {
    }
    ~SEMO_Mesh_Load() {
    }
	
	void OpenFoam(SEMO_Mesh_Builder &in, string folder);
	
	void vtu(SEMO_Mesh_Builder &mesh, string folder);
	
	void vtu(string folder, 
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species);
	
	void boundaryProperties(string file, 
			SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<double>& value);
	void boundaryVelocities(string file, 
			SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<vector<double>>& value);

};


