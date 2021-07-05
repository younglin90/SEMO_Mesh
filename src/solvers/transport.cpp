#include "build.h"
#include <cmath>
#include <array>


#include "../turbulenceModels/LES.h"


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
	
	//========================
	// turbulence models
	if(controls.turbType == "LES"){
			
		SEMO_TurbModel_LES LES;
		if(controls.turbLESModel == "WALE"){
			LES.calcWALE(mesh, controls, species);
		}
		else if(controls.turbLESModel == "smagorinsky"){
			LES.calcSmagorinsky(mesh, controls, species);
		}
		else if(controls.turbLESModel == "dynamicSmagorinsky"){
			LES.calcDynamicSmagorinsky(mesh, controls, species);
		}
		else{
			cerr << "| #Error : not defined LES model" << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			
		}
		
	}
	else if(controls.turbType == "RANS"){
			
		cerr << "| #Error : not yet defined RANS model" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	}
	else if(
	controls.turbType == "laminar" ||
	controls.turbType == "none"
	){
	
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			cell.var[controls.muT] = 0.0;
			cell.var[controls.kSGS] = 0.0;
		}
	}
	else{
		cerr << "| #Error : not defined turb type" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	}
	//========================
		
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		cell.var[controls.mu] += cell.var[controls.muT];
		// cell.var[controls.muEffective] = cell.var[controls.mu] + cell.var[controls.muT];
		
		// cell.var[controls.kappaEffective] = 
			// cell.var[controls.kappa] + cell.var[controls.cp] * cell.var[controls.muT] / controls.PrT;
		
		// cell.var[controls.DEffective] = 
			// cell.var[controls.D] + cell.var[controls.muT] / (cell.var[controls.Rho] * controls.ScT);
	}
	
	
}



