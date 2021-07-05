#include "build.h"
#include <cmath>
#include <array>
#include <time.h>

void SEMO_Solvers_Builder::calcNormResiduals(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals,
	vector<double>& norm){
	
	norm.clear();
	norm.resize(controls.nEq,0.0);
	
	double totalVol = 0.0;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		for(int j=0; j<controls.nEq; ++j){
			norm[j] += pow(residuals[i][j],2.0)*mesh.cells[i].volume;
			// norm[j] += pow(residuals[i][j],2.0);
		}
		totalVol += mesh.cells[i].volume;
		// ++totalVol;
	}
	
	vector<double> normReduced(controls.nEq,0.0);
	double totalVolReduced;
    MPI_Allreduce(norm.data(), normReduced.data(), controls.nEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&totalVol, &totalVolReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	

	for(int j=0; j<controls.nEq; ++j){
		normReduced[j] /= totalVolReduced;
		norm[j] = sqrt(normReduced[j]);
	}
	
	double dClock = clock() - controls.startClock;
	// dClock /= 1000.0;
	dClock /= CLOCKS_PER_SEC;
 
	SEMO_Mesh_Save save;
	save.gnuplot(controls.iterTotal, dClock, norm);
 
	
}





void SEMO_Solvers_Builder::calcNormResiduals(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<double>& residuals,
	vector<double>& norm){
	
	norm.clear();
	norm.resize(controls.nEq,0.0);
	
	double totalVol = 0.0;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		totalVol += mesh.cells[i].volume;
	}
	
	for(int j=0; j<controls.nEq; ++j){
		norm[j] += residuals[j];
	}
	
	vector<double> normReduced(controls.nEq,0.0);
	double totalVolReduced;
    MPI_Allreduce(norm.data(), normReduced.data(), controls.nEq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&totalVol, &totalVolReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	

	for(int j=0; j<controls.nEq; ++j){
		normReduced[j] /= totalVolReduced;
		norm[j] = sqrt(normReduced[j]);
	}
	
	
	double dClock = clock() - controls.startClock;
	// dClock /= 1000.0;
	dClock /= CLOCKS_PER_SEC;
	
	SEMO_Mesh_Save save;
	save.gnuplot(controls.iterTotal, dClock, norm);
 
	
}