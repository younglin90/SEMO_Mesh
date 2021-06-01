#include "build.h"
#include <cmath>
#include <array>





void SEMO_Solvers_Builder::calcCellTransport(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int nSp = controls.nSp;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		

		vector<double> VF;
		double tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			VF.push_back(cell.var[controls.VF[ns]]);
			tmp1 += cell.var[controls.VF[ns]];
		}
		VF.push_back(1.0 - tmp1);
		
		double mu = 0.0;
		for(int ns=0; ns<nSp; ++ns){
			mu += VF[ns]*species[ns].mu;
		}
		
		cell.var[controls.mu] = mu;
		
	}
	
	
	
	
}



