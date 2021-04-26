#pragma once
#include <iostream>

// #include "build.h"
class SEMO_Mesh_Builder;

class SEMO_Mesh_Save {
private:

public:
    SEMO_Mesh_Save() {
    }
    ~SEMO_Mesh_Save() {
    }
	
	void vtu(SEMO_Mesh_Builder &in);

};


