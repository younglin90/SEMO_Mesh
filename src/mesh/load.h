#pragma once
#include <iostream>

// #include "build.h"
class SEMO_Mesh_Builder;

class SEMO_Mesh_Load {
private:

public:
    SEMO_Mesh_Load() {
    }
    ~SEMO_Mesh_Load() {
    }
	
	void OpenFoam(SEMO_Mesh_Builder &in);

};


