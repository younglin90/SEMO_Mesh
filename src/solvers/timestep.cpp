#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcPseudoTimeStep(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	
	for(auto& cell : mesh.cells){
		
		int i = controls.dtPseudo;
		int Ur = controls.Ur;
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		double C = cell.var[controls.C];
		
		
		cell.var[Ur] = magVel;
		// cell.var[Ur] = max(cell.var[Ur], 0.001);
		cell.var[Ur] = max(cell.var[Ur], 0.1*C);
		cell.var[Ur] = min(cell.var[Ur], C);
		
		double eta = pow(cell.var[Ur]/C,2.0);
		double iegenU = 0.5*magVel*(1.0+eta);
		double iegenC = 0.5*sqrt(pow(magVel*(1.0-eta),2.0)+4.0*pow(cell.var[Ur],2.0));
		
		double maxA=0.0;
		for(auto& j : cell.faces){
			maxA = max(maxA, mesh.faces[j].area);
		}
		
		cell.var[i] = cell.volume / ( iegenU + iegenC ) / maxA;
		
		double viscTime = pow(cell.volume / maxA, 2.0) / 2.0 
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		cell.var[i] = min(cell.var[i], viscTime);
		
		// cell.var[i] = pow(cell.volume,0.3333) / ( iegenU + iegenC );
		// double viscTime = pow(cell.volume,0.6666) / 2.0 
			// / ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		// cell.var[i] = min(cell.var[i], viscTime);
		
		cell.var[i] *= controls.pseudoCo;
		
	}
	

	
}




void SEMO_Solvers_Builder::calcRealTimeStep(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	double minValue=10000.0;
	for(auto& cell : mesh.cells){
		
		int i = controls.dtPseudo;
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		double C = cell.var[controls.C];
		
		cell.var[i] = pow(cell.volume,0.333) / ( magVel + abs(C) );
		
		cell.var[i] *= controls.pseudoCo;
		
		minValue = min(minValue,cell.var[i]);
	}
	

	double minValueReduced;
    MPI_Allreduce(&minValue, &minValueReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	controls.timeStep = minValueReduced;

	
}






void SEMO_Solvers_Builder::calcIncomRealTimeStep(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	double minDt=1000000.0;
	for(auto& cell : mesh.cells){
		
		double maxA=0.0;
		double minA=1000000.0;
		for(auto& j : cell.faces){
			maxA = max(maxA, mesh.faces[j].area);
			minA = min(minA, mesh.faces[j].area);
		}
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		double convTime = controls.maxCo * cell.volume / maxA / (magVel+1.e-200);
		minDt = min(minDt,convTime);
		
		double viscTime = controls.maxCo * 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		minDt = min(minDt,viscTime);
		
		double gravMag = 
			sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 pow(controls.gravityAcceleration[1],2.0) +
                 pow(controls.gravityAcceleration[2],2.0));
		double gravTime = controls.maxCo * sqrt( pow(cell.volume,0.333) / (gravMag+1.e-200) );
		minDt = min(minDt,gravTime);
		
		
		double surfTensCoeff = 1.96;
		double surfTensTime = 
			controls.maxCo * pow(pow(cell.volume,0.333),1.5)*sqrt(1001/(4.0*3.14*surfTensCoeff));
		minDt = min(minDt,surfTensTime);
		
		
		
	}
	

	double minDtReduced;
    MPI_Allreduce(&minDt, &minDtReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	
	controls.timeStep = controls.orgTimeStep;
	if(controls.orgTimeStep > controls.maxTimeStep) controls.timeStep = controls.maxTimeStep;
	if(minDtReduced > controls.maxTimeStep) controls.timeStep = controls.maxTimeStep;
	if(minDtReduced < controls.orgTimeStep) controls.timeStep = minDtReduced;

	
}