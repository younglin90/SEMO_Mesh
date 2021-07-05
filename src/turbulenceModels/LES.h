#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <mpi.h>
using namespace std;

#include "../controls/build.h" 
#include "../species/species.h" 
#include "../math/math.h"

// #include "build.h"
class SEMO_Mesh_Builder;
class SEMO_Controls_Builder;

class SEMO_TurbModel_LES {
private:

public:
    SEMO_TurbModel_LES() {
    }
    ~SEMO_TurbModel_LES() {
    }
	
	// void calcEffectiveValues(
		// SEMO_Mesh_Builder &mesh, 
		// SEMO_Controls_Builder &controls,
		// vector<SEMO_Species>& species
		// );
		
	void calcWALE(
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species
		);
		
	void calcSmagorinsky(
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species
		);
		
	void calcDynamicSmagorinsky(
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species
		);

};