#include "build.h"
#include <cmath>
#include <array>



void SEMO_Solvers_Builder::sourceTerms(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls, 
	vector<SEMO_Species>& species,
	vector<vector<double>>& residuals){


	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		//==================
		// real time
		vector<double> Q(controls.nEq,0.0);
		
		Q[0] = cell.var[controls.Rho];
		Q[1] = cell.var[controls.Rho] * cell.var[controls.U];
		Q[2] = cell.var[controls.Rho] * cell.var[controls.V];
		Q[3] = cell.var[controls.Rho] * cell.var[controls.W];
		Q[4] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
		for(int ns=0; ns<controls.nSp-1; ++ns){
			Q[5+ns] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
		}
		
		
		for(int j=0; j<controls.nEq; ++j){
			
			// 1st-order euler
			// double source = 
				// (Q[j] - cell.var[controls.Qn[j]]) 
				// / controls.timeStep * cell.volume;
				
			// 2nd-order
			double source = 
				(1.5*Q[j] - 2.0*cell.var[controls.Qn[j]] + 0.5*cell.var[controls.Qm[j]]) 
				/ controls.timeStep * cell.volume;
			
			residuals[i][j] += source;
		}
		
		//==================
		// gravity
		residuals[i][1] -= cell.var[controls.Rho] * controls.gravityAcceleration[0] * cell.volume;
		residuals[i][2] -= cell.var[controls.Rho] * controls.gravityAcceleration[1] * cell.volume;
		residuals[i][3] -= cell.var[controls.Rho] * controls.gravityAcceleration[2] * cell.volume;
		
		
	}
	
	
	
	//==================
	// surface tension
	
	SEMO_Utility_Math math;
	
	vector<vector<double>> gradAi;
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	
	vector<double> kappa;
	this->calcCurvature(mesh, controls.VF[0], kappa);
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		residuals[i][1] -= cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
		residuals[i][2] -= cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
		residuals[i][3] -= cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
	}
	
	

}