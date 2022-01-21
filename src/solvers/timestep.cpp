#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcPseudoTimeStep(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	
	
	
	for(auto& cell : mesh.cells){

		double maxA=0.0;
		double minA=1000000.0;
		for(auto& j : cell.faces){
			maxA = max(maxA, mesh.faces[j].area);
			minA = min(minA, mesh.faces[j].area);
		}
		
		
		int i = controls.dtPseudo;
		int Ur = controls.Ur;
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		double C = cell.var[controls.C];
		
		
		cell.var[Ur] = magVel;
		double Uco = controls.Uco;
		cell.var[Ur] = max(cell.var[Ur], Uco);
		double Uun = controls.Lch/3.141592/controls.timeStep;
		cell.var[Ur] = max(cell.var[Ur], Uun);
		// cell.var[Ur] = max(cell.var[Ur], 0.0001*C);
		cell.var[Ur] = min(cell.var[Ur], C);
		
		
		
		
		
		// cell.var[Ur] = C;
		
		
	
		
		double eta = pow(cell.var[Ur]/C,2.0);
		double iegenU = 0.5*magVel*(1.0+eta);
		double iegenC = 0.5*sqrt(pow(magVel*(1.0-eta),2.0)+4.0*pow(cell.var[Ur],2.0));
		double convTime = controls.allowableCFL * cell.volume / maxA / (iegenU+iegenC);
		
		double viscTime = controls.allowableCFL * 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
			
		double minDt = min(convTime,viscTime);
		
		double gravMag = 
			sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 pow(controls.gravityAcceleration[1],2.0) +
                 pow(controls.gravityAcceleration[2],2.0));
		double gravTime = controls.allowableCFL * sqrt( pow(cell.volume,0.333) / (gravMag+1.e-200) );
		minDt = min(minDt,gravTime);
		
		
		double surfTensCoeff = 1.96;
		double surfTensTime = 
			controls.allowableCFL * pow(pow(cell.volume,0.333),1.5)*sqrt(1001/(4.0*3.14*surfTensCoeff));
		minDt = min(minDt,surfTensTime);
		
		
		cell.var[i] = minDt;
		
	}
	
	
	// double minTime = 1000000.0;
	// for(auto& cell : mesh.cells){

		// minTime = min(minTime,cell.var[controls.dtPseudo]);
		
	// }
	// for(auto& cell : mesh.cells){

		// cell.var[controls.dtPseudo] = minTime;
		
	// }
	
	

	
}




void SEMO_Solvers_Builder::calcRealTimeStep(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	
	double minValue=10000.0;
	for(auto& cell : mesh.cells){

		double maxA=0.0;
		double minA=1000000.0;
		for(auto& j : cell.faces){
			maxA = max(maxA, mesh.faces[j].area);
			minA = min(minA, mesh.faces[j].area);
		}
		
		
		int i = controls.dtPseudo;
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		double C = cell.var[controls.C];
		
		cell.var[i] = controls.pseudoCo * cell.volume / maxA / ( magVel + C );
		
		
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




void SEMO_Solvers_Builder::calcCorantNumberForPrint(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	double& corantNum){
	
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
		
		double convTime = cell.volume / maxA / (magVel+1.e-200);
		minDt = min(minDt,convTime);
		
		double viscTime = 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		minDt = min(minDt,viscTime);
		
		double gravMag = 
			sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 pow(controls.gravityAcceleration[1],2.0) +
                 pow(controls.gravityAcceleration[2],2.0));
		double gravTime = sqrt( sqrt(minA) / (gravMag+1.e-200) );
		minDt = min(minDt,gravTime);
		
		
		double surfTensCoeff = 0.0728;
		double surfTensTime = 
			// pow(pow(cell.volume,0.333),1.5)*sqrt(1001.0/(4.0*3.14*surfTensCoeff));
			sqrt(pow(minA,1.5)*(1000.0+1.2)/(4.0*3.141592*surfTensCoeff));
		minDt = min(minDt,surfTensTime);
		
	}
	

	double minDtReduced;
    MPI_Allreduce(&minDt, &minDtReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	
	corantNum = controls.timeStep / minDtReduced;

	
}






double SEMO_Solvers_Builder::calcTimeStepFromCorantNumber(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
		
		
	double corantNum = controls.maxCo;
	
	
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
		
		double convTime = cell.volume / maxA / (magVel+1.e-200);
		minDt = min(minDt,convTime);
		
		double viscTime = 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		minDt = min(minDt,viscTime);
		
		double gravMag = 
			sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 pow(controls.gravityAcceleration[1],2.0) +
                 pow(controls.gravityAcceleration[2],2.0));
		double gravTime = sqrt( sqrt(minA) / (gravMag+1.e-200) );
		minDt = min(minDt,gravTime);
		
		
		double surfTensCoeff = species[0].sigma;
		double minDX = cell.volume/maxA;
		double C_sigma = sqrt( 6.28*surfTensCoeff/(species[0].rho + species[1].rho)/minDX );
		double surfTensTime = 
			// pow(pow(cell.volume,0.333),1.5)*sqrt(1001.0/(4.0*3.14*surfTensCoeff));
			// sqrt(pow(minA,1.5)*(1000.0+1.2)/(4.0*3.141592*surfTensCoeff));
			// sqrt(pow(minA,1.5)*(cell.var[controls.Rho])/(4.0*3.141592*surfTensCoeff));
			// sqrt(pow(minA,1.5)*(species[0].rho + species[1].rho)/(4.0*3.141592*surfTensCoeff));
			minDX/(C_sigma + magVel);
			
			
		minDt = min(minDt,surfTensTime);
		
	}
	

	double minDtReduced;
    MPI_Allreduce(&minDt, &minDtReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	// cout << minDtReduced << endl;
	
	controls.timeStep = corantNum * minDtReduced;
	
	double corantNum_out = corantNum;
	
	if(controls.timeStep>controls.maxTimeStep) {
		corantNum_out = corantNum * controls.maxTimeStep/controls.timeStep;
		controls.timeStep = controls.maxTimeStep;
	}

	return corantNum_out;
}





void SEMO_Solvers_Builder::calcExplicitCorantNumberForPrint(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	double& corantNum){
	
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
		
		double convTime = cell.volume / maxA / (magVel + cell.var[controls.C] +1.e-200);
		minDt = min(minDt,convTime);
		
		double viscTime = 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		minDt = min(minDt,viscTime);
		
		double gravMag = 
			sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 pow(controls.gravityAcceleration[1],2.0) +
                 pow(controls.gravityAcceleration[2],2.0));
		double gravTime = sqrt( sqrt(minA) / (gravMag+1.e-200) );
		minDt = min(minDt,gravTime);
		
		
		double surfTensCoeff = 0.0728;
		double surfTensTime = 
			// pow(pow(cell.volume,0.333),1.5)*sqrt(1001.0/(4.0*3.14*surfTensCoeff));
			sqrt(pow(minA,1.5)*(1000.0+1.2)/(4.0*3.141592*surfTensCoeff));
		minDt = min(minDt,surfTensTime);
		
	}
	

	double minDtReduced;
    MPI_Allreduce(&minDt, &minDtReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	
	corantNum = controls.timeStep / minDtReduced;

	
}






double SEMO_Solvers_Builder::calcExplicitTimeStepFromCorantNumber(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
		
		
	double corantNum = controls.maxCo;
	
	
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
		
		double convTime = cell.volume / maxA / (magVel + cell.var[controls.C] +1.e-200);
		minDt = min(minDt,convTime);
		
		double viscTime = 0.5 * minA
			/ ( cell.var[controls.mu] + 1.e-200 ) * cell.var[controls.Rho];
		minDt = min(minDt,viscTime);
		
		// double gravMag = 
			// sqrt(pow(controls.gravityAcceleration[0],2.0) +
                 // pow(controls.gravityAcceleration[1],2.0) +
                 // pow(controls.gravityAcceleration[2],2.0));
		// double gravTime = sqrt( sqrt(minA) / (gravMag+1.e-200) );
		// minDt = min(minDt,gravTime);
		
		
		// double surfTensCoeff = species[0].sigma;
		// double minDX = cell.volume/maxA;
		// double C_sigma = sqrt( 6.28*surfTensCoeff/(species[0].rho + species[1].rho)/minDX );
		// double surfTensTime = 
			// // pow(pow(cell.volume,0.333),1.5)*sqrt(1001.0/(4.0*3.14*surfTensCoeff));
			// // sqrt(pow(minA,1.5)*(1000.0+1.2)/(4.0*3.141592*surfTensCoeff));
			// // sqrt(pow(minA,1.5)*(cell.var[controls.Rho])/(4.0*3.141592*surfTensCoeff));
			// // sqrt(pow(minA,1.5)*(species[0].rho + species[1].rho)/(4.0*3.141592*surfTensCoeff));
			// minDX/(C_sigma + magVel);
			
			
		// minDt = min(minDt,surfTensTime);
		
	}
	// cout << "AAAAAAA : " << minDt << endl;
	
	

	double minDtReduced;
    MPI_Allreduce(&minDt, &minDtReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	
	// cout << minDtReduced << endl;
	
	controls.timeStep = corantNum * minDtReduced;
	
	double corantNum_out = corantNum;
	
	if(controls.timeStep>controls.maxTimeStep) {
		corantNum_out = corantNum * controls.maxTimeStep/controls.timeStep;
		controls.timeStep = controls.maxTimeStep;
	}

	return corantNum_out;
}