#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcRHS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	
	this->calcFluxes(mesh, controls, residuals);
	
	
	this->sourceTerms(mesh, controls, residuals);
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		SEMO_Cell& cell = mesh.cells[i];
		
		for(int nq=0; nq<controls.nEq; ++nq){
			residuals[i][nq] *= -cell.var[controls.dtPseudo]/cell.volume;
		}
	}
	
	
}


void SEMO_Solvers_Builder::calcSingleRHS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	
	this->calcFluxes(mesh, controls, residuals);
	
	
	// this->sourceTerms(mesh, controls, residuals);
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		SEMO_Cell& cell = mesh.cells[i];
		
		for(int nq=0; nq<controls.nEq; ++nq){
			residuals[i][nq] *= -controls.timeStep/cell.volume;
		}
	}
	
	
}





void SEMO_Solvers_Builder::calcFluxes(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
		


	// calc gradient
	SEMO_Utility_Math math;
	vector<vector<double>> gradU;
	math.calcLeastSquare(mesh, controls.U, controls.fU, gradU);
	math.calcGradientFace(mesh, gradU, controls.U, controls.fU, 
					controls.dUdx, controls.dUdy, controls.dUdz);
	vector<vector<double>> gradV;
	math.calcLeastSquare(mesh, controls.V, controls.fV, gradV);
	math.calcGradientFace(mesh, gradV, controls.V, controls.fV, 
					controls.dVdx, controls.dVdy, controls.dVdz);
	vector<vector<double>> gradW;
	math.calcLeastSquare(mesh, controls.W, controls.fW, gradW);
	math.calcGradientFace(mesh, gradW, controls.W, controls.fW, 
					controls.dWdx, controls.dWdy, controls.dWdz);
	

	
	for(auto& face : mesh.faces){
		
		vector<double> flux(controls.nEq,0.0);
		vector<double> vFlux(controls.nEq,0.0);
		
		vector<double> MFL;
		vector<double> MFR;
		
		double MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			MFL.push_back(face.varL[controls.fMF[ns]]);
			MFnSp += face.varL[controls.fMF[ns]];
		}
		MFL.push_back(1.0 - MFnSp);
		
		MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			MFR.push_back(face.varR[controls.fMF[ns]]);
			MFnSp += face.varR[controls.fMF[ns]];
		}
		MFR.push_back(1.0 - MFnSp);
		
		
		
		this->calcFluxHAUS(
			face.varL[controls.fP],
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW],
			face.varL[controls.fT], MFL,
			face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
			face.varR[controls.fP],
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW],
			face.varR[controls.fT], MFR,
			face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			face.unitNormals,
			flux );


		
		this->calcViscousFlux(
			face.varL[controls.fmu],
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW],
			face.varL[controls.fRho],
			face.varR[controls.fmu],
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW],
			face.varR[controls.fRho],
			face.var[controls.dUdx], face.var[controls.dUdy], face.var[controls.dUdz],
			face.var[controls.dVdx], face.var[controls.dVdy], face.var[controls.dVdz],
			face.var[controls.dWdx], face.var[controls.dWdy], face.var[controls.dWdz],
			face.unitNormals,
			vFlux );
			
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(int nq=0; nq<controls.nEq; ++nq){
				residuals[face.owner][nq] += (flux[nq]-vFlux[nq])*face.area;
				residuals[face.neighbour][nq] -= (flux[nq]-vFlux[nq])*face.area;
			}
		}
		else if(
		face.getType() == SEMO_Types::BOUNDARY_FACE ||
		face.getType() == SEMO_Types::PROCESSOR_FACE
		){
			for(int nq=0; nq<controls.nEq; ++nq){
				residuals[face.owner][nq] += (flux[nq]-vFlux[nq])*face.area;
			}
		}
		
	}
	
	

		
}