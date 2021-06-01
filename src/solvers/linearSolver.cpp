#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcLinearSolver(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	SEMO_Utility_Math math;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		int matSize = controls.nEq;
		
		
		vector<vector<double>> matrixP(matSize,vector<double>(matSize,0.0));
		vector<vector<double>> matrixT(matSize,vector<double>(matSize,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(matSize,vector<double>(matSize,0.0));
		
		for(int in=0; in<matSize; ++in){
			for(int jn=0; jn<matSize; ++jn){
				matrixA[in][jn] = matrixP[in][jn] + 
					1.5 * cell.var[controls.dtPseudo] / controls.timeStep * matrixT[in][jn];
					
			}
		}
		
		vector<double> vectorB;
		for(auto& j : residuals[i]){
			vectorB.push_back(j);
		}
		
		
		math.GaussSeidelSOR(matrixA);
		
		
		for(int j=0; j<residuals[i].size(); ++j){
			double updateResi = 0.0;
			for(int k=0; k<residuals[i].size(); ++k){
				updateResi += matrixA[j][k]*vectorB[k];
			}
			residuals[i][j] = updateResi;
		}
		
		
	}
	

	
}




void SEMO_Solvers_Builder::calcSingleLinearSolver(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	SEMO_Utility_Math math;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		int matSize = controls.nEq;
		
		
		vector<vector<double>> matrixP(matSize,vector<double>(matSize,0.0));
		vector<vector<double>> matrixT(matSize,vector<double>(matSize,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(matSize,vector<double>(matSize,0.0));
		
		for(int in=0; in<matSize; ++in){
			for(int jn=0; jn<matSize; ++jn){
				matrixA[in][jn] = matrixT[in][jn];
					
			}
		}
		
		vector<double> vectorB;
		for(auto& j : residuals[i]){
			vectorB.push_back(j);
		}
		
		
		math.GaussSeidelSOR(matrixA);
		
		
		for(int j=0; j<residuals[i].size(); ++j){
			double updateResi = 0.0;
			for(int k=0; k<residuals[i].size(); ++k){
				updateResi += matrixA[j][k]*vectorB[k];
			}
			residuals[i][j] = updateResi;
		}
		
		
	}
	

	
}






void SEMO_Solvers_Builder::constructMatrix(
	SEMO_Cell& cell, SEMO_Controls_Builder& controls,
	vector<vector<double>>& matrixP,
	vector<vector<double>>& matrixT){


	int N = matrixP.size();

	double DRDP = cell.var[controls.dRhoDP];
	double DHDP = cell.var[controls.dHtDP];
	double DRDT = cell.var[controls.dRhoDT];
	double DHDT = cell.var[controls.dHtDT];
	vector<double> DRDY, DHDY;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		DRDY.push_back(cell.var[controls.dRhoDMF[ns]]);
		DHDY.push_back(cell.var[controls.dHtDMF[ns]]);
	}
	// cell.var[controls.dRhoDMF[ns]]
	double Rho = cell.var[controls.Rho];
	double U = cell.var[controls.U];
	double V = cell.var[controls.V];
	double W = cell.var[controls.W];
	double T = cell.var[controls.T];
	vector<double> Y;
	for(int ns=0; ns<controls.nSp; ++ns){
		Y.push_back(cell.var[controls.MF[ns]]);
	}
	double Ht = cell.var[controls.Ht];
	
	
	matrixT[0][0] = DRDP;
	matrixT[1][0] = DRDP*U;
	matrixT[2][0] = DRDP*V;
	matrixT[3][0] = DRDP*W;
	matrixT[4][0] = DRDP*Ht + Rho*DHDP - 1.0;
	
	// cout << DRDP << " " << Ht << " " << Rho << " " << DHDP << endl;
	
	matrixT[1][1] = Rho;
	matrixT[4][1] = Rho*U;
	
	matrixT[2][2] = Rho;
	matrixT[4][2] = Rho*V;
	
	matrixT[3][3] = Rho;
	matrixT[4][3] = Rho*W;
	
	matrixT[0][4] = DRDT;
	matrixT[1][4] = DRDT*U;
	matrixT[2][4] = DRDT*V;
	matrixT[3][4] = DRDT*W;
	matrixT[4][4] = DRDT*Ht + Rho*DHDT;
	
	// cout << DRDT*Ht + Rho*DHDT << endl;
	
	for(int ns=0; ns<controls.nSp-1; ++ns){
		int nq = 5 + ns;
		matrixT[nq][0] = DRDP*Y[ns];
		matrixT[nq][4] = DRDT*Y[ns];
		
		// cout << DRDT << " " << Y[ns] << endl;
		
		for(int ns2=0; ns2<controls.nSp-1; ++ns2){
			int nq2 = 5 + ns2;
			matrixT[nq][nq2] = Y[ns]*DRDY[ns2];
		}
		matrixT[nq][nq] += Rho;
		
		matrixT[0][nq] = DRDY[ns];
		matrixT[1][nq] = DRDY[ns]*U;
		matrixT[2][nq] = DRDY[ns]*V;
		matrixT[3][nq] = DRDY[ns]*W;
		matrixT[4][nq] = DRDY[ns]*Ht + Rho*DHDY[ns];
	}
	
	
		// cout << matrixT[4][5] << endl;
	

	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			matrixP[i][j] = matrixT[i][j];
		}
	}
	
	// double beta = 0.5*DRDP;
	// double beta = 1.0/pow(cell.var[controls.Ur],2.0);
	double beta = 1.0/pow(cell.var[controls.Ur],2.0) 
					+ DRDP - 1.0/pow(cell.var[controls.C],2.0);
	
	matrixT[0][0] = beta;
	matrixT[1][0] = beta*U;
	matrixT[2][0] = beta*V;
	matrixT[3][0] = beta*W;
	matrixT[4][0] = beta*Ht + Rho*DHDP - 1.0;
	
	for(int ns=0; ns<controls.nSp-1; ++ns){
		int nq = 5 + ns;
		matrixT[nq][0] = beta*Y[ns];
	}
	




}