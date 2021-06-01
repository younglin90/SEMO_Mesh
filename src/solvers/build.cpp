#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::setInitValues(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	
	
	for(auto& cell : mesh.cells){
		// cell.var.resize(5 + controls.nSp,0.0);
		cell.var.resize(controls.nTotalCellVar,0.0);
		cell.var[controls.P] = 101325.0;
		cell.var[controls.U] = 0.0;
		cell.var[controls.V] = 0.0;
		cell.var[controls.W] = 0.0;
		cell.var[controls.T] = 300.0;
		
		if(controls.nSp>1){
			cell.var[controls.VF[0]] = 1.0;
			for(int i=1; i<controls.nSp; ++i){
				cell.var[controls.VF[i]] = 0.0;
			}
		}
		else{
			cell.var[controls.VF[0]] = 0.0;
		}
		
		cell.var[controls.Rho] = 1000.0*cell.var[controls.VF[0]]+1.0*(1.0-cell.var[controls.VF[0]]);
	}
	
	
	for(auto& face : mesh.faces){
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		face.var.resize(controls.nTotalFaceVar,0.0);
	}
	
	
}


