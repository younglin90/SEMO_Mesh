#include "LES.h" 

// void SEMO_TurbModel_LES::calcEffectiveValues(
	// SEMO_Mesh_Builder &mesh, 
	// SEMO_Controls_Builder &controls,
	// vector<SEMO_Species>& species
	// ){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// cell.var[controls.muEffective] = cell.var[controls.mu] + cell.var[controls.muT];
		
		// // cell.var[controls.kappaEffective] = 
			// // cell.var[controls.kappa] + cell.var[controls.cp] * cell.var[controls.muT] / controls.PrT;
		
		// // cell.var[controls.DEffective] = 
			// // cell.var[controls.D] + cell.var[controls.muT] / (cell.var[controls.Rho] * controls.ScT);
	// }
	
	
	
// }

void SEMO_TurbModel_LES::calcWALE(
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
	double vkc = 0.05;
	double Ci = 0.007;
	double Cw = 0.325;
	

	// gradient U, V, W
	SEMO_Utility_Math math;
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	vector<vector<double>> S(3,vector<double>(3,0.0));
	vector<vector<double>> Sd(3,vector<double>(3,0.0));
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		double dUdx = gradU[i][0];
		double dUdy = gradU[i][1];
		double dUdz = gradU[i][2];
		double dVdx = gradV[i][0];
		double dVdy = gradV[i][1];
		double dVdz = gradV[i][2];
		double dWdx = gradW[i][0];
		double dWdy = gradW[i][1];
		double dWdz = gradW[i][2];
		
		S[0][0] = dUdx;
		S[0][1] = 0.5*(dUdy+dVdx);
		S[0][2] = 0.5*(dUdz+dWdx);
		S[1][0] = S[0][1];
		S[1][1] = dVdy;
		S[1][2] = 0.5*(dVdz+dWdy);
		S[2][0] = S[0][2];
		S[2][1] = S[1][2];
		S[2][2] = dWdz;

		Sd[0][0] = 0.5*(pow(dUdx,2.0)+pow(dUdx,2.0));
		Sd[0][1] = 0.5*(pow(dUdy,2.0)+pow(dVdx,2.0));
		Sd[0][2] = 0.5*(pow(dUdz,2.0)+pow(dWdx,2.0));
			
		Sd[1][0] = 0.5*(pow(dVdx,2.0)+pow(dUdy,2.0));
		Sd[1][1] = 0.5*(pow(dVdy,2.0)+pow(dVdy,2.0));
		Sd[1][2] = 0.5*(pow(dVdz,2.0)+pow(dWdy,2.0));
			
		Sd[2][0] = 0.5*(pow(dWdx,2.0)+pow(dUdz,2.0));
		Sd[2][1] = 0.5*(pow(dWdy,2.0)+pow(dVdz,2.0));
		Sd[2][2] = 0.5*(pow(dWdz,2.0)+pow(dWdz,2.0));
		
		double tmpDVDx = 1.0/3.0*pow(dUdx+dVdy+dWdz,2.0);
		Sd[0][0] -= tmpDVDx;
		Sd[1][1] -= tmpDVDx;
		Sd[2][2] -= tmpDVDx;
		
        // calculate S_ij * S_ij
		double SS = 
			pow(S[0][0],2.0) + pow(S[1][1],2.0) + pow(S[2][2],2.0) +
			2.0*(pow(S[0][1],2.0) + pow(S[0][2],2.0) + pow(S[1][2],2.0));
		double SSd = 
			pow(S[0][0],2.0) + pow(S[1][1],2.0) + pow(S[2][2],2.0) +
			2.0*(pow(Sd[0][1],2.0) + pow(Sd[0][2],2.0) + pow(Sd[1][2],2.0));
		
		double Sabs = sqrt(2.0*SS);
		
		// calculate Delta_bar
		double Ls = Cw * pow(cell.volume,0.333333333);
		// Ls = min(Ls, vkc*wdist);
		double Ls2 = Ls*Ls;
		
		// calculate muT & kSGS
		double tmpSSSSd = pow(SS,2.5)+pow(SSd,1.25);
		cell.var[controls.muT] = 0.0;
		if( tmpSSSSd != 0.0 ){
			cell.var[controls.muT] = cell.var[controls.Rho] * Ls2 * pow(SSd,1.5) / tmpSSSSd;
		}
		
		cell.var[controls.kSGS] = Ci * pow(cell.volume,0.6666666666) * (0.5*Sabs*Sabs);
		
		
	}
	
	
	
}








void SEMO_TurbModel_LES::calcSmagorinsky(
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
    double C_I = 0.007;
    double C_R = 0.1;
	double C_eta = 0.92;
	

	// gradient U, V, W
	SEMO_Utility_Math math;
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	vector<vector<double>> S(3,vector<double>(3,0.0));
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		double dUdx = gradU[i][0];
		double dUdy = gradU[i][1];
		double dUdz = gradU[i][2];
		double dVdx = gradV[i][0];
		double dVdy = gradV[i][1];
		double dVdz = gradV[i][2];
		double dWdx = gradW[i][0];
		double dWdy = gradW[i][1];
		double dWdz = gradW[i][2];
		
		S[0][0] = dUdx;
		S[0][1] = 0.5*(dUdy+dVdx);
		S[0][2] = 0.5*(dUdz+dWdx);
		S[1][0] = S[0][1];
		S[1][1] = dVdy;
		S[1][2] = 0.5*(dVdz+dWdy);
		S[2][0] = S[0][2];
		S[2][1] = S[1][2];
		S[2][2] = dWdz;

        // calculate S_ij * S_ij
		double SS = 
			pow(S[0][0],2.0) + pow(S[1][1],2.0) + pow(S[2][2],2.0) +
			2.0*(pow(S[0][1],2.0) + pow(S[0][2],2.0) + pow(S[1][2],2.0));
		
		double Sabs = sqrt(2.0*SS);
		
        // calculate Delta_bar
        double DeltaBar2 = pow(cell.volume,0.666666666);
		
		// calculate muT & kSGS
		cell.var[controls.muT] = C_R * cell.var[controls.Rho]*DeltaBar2*Sabs;
		
		cell.var[controls.kSGS] = C_I * DeltaBar2 * (0.5*Sabs*Sabs);
		
		
	}
	
	
	
}



void SEMO_TurbModel_LES::calcDynamicSmagorinsky(
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
    double C_I = 0.007;
	double C_eta = 0.92;
	

	// gradient U, V, W
	SEMO_Utility_Math math;
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	

	vector<vector<double>> S(3,vector<double>(3,0.0));
	vector<vector<double>> L(3,vector<double>(3,0.0));
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		double dUdx = gradU[i][0];
		double dUdy = gradU[i][1];
		double dUdz = gradU[i][2];
		double dVdx = gradV[i][0];
		double dVdy = gradV[i][1];
		double dVdz = gradV[i][2];
		double dWdx = gradW[i][0];
		double dWdy = gradW[i][1];
		double dWdz = gradW[i][2];
		
        // calculate strain rate tensor, Sij
		S[0][0] = dUdx;
		S[0][1] = 0.5*(dUdy+dVdx);
		S[0][2] = 0.5*(dUdz+dWdx);
		S[1][0] = S[0][1];
		S[1][1] = dVdy;
		S[1][2] = 0.5*(dVdz+dWdy);
		S[2][0] = S[0][2];
		S[2][1] = S[1][2];
		S[2][2] = dWdz;

        // calculate S_ij * S_ij
		double SS = 
			pow(S[0][0],2.0) + pow(S[1][1],2.0) + pow(S[2][2],2.0) +
			2.0*(pow(S[0][1],2.0) + pow(S[0][2],2.0) + pow(S[1][2],2.0));
		
		double Sabs = sqrt(2.0*SS);
		
		
        // filtering process f(rho*ui*uj), f(rho), f(ui)
		double Rho = cell.var[controls.Rho];
		double U = cell.var[controls.U];
		double V = cell.var[controls.V];
		double W = cell.var[controls.W];
		double vol = cell.volume;
		L[0][0] = Rho * U * U * vol;
		L[0][1] = Rho * U * V * vol;
		L[0][2] = Rho * U * W * vol;
		L[1][1] = Rho * V * V * vol;
		L[1][2] = Rho * V * W * vol;
		L[2][2] = Rho * W * W * vol;
		
		double sumVol = vol;
		double sumRho = Rho * vol;
		double sumV = U * vol;
		double sumU = V * vol;
		double sumW = W * vol;
		double sumSabs = Sabs * vol;
		
		
		vector<vector<double>> M(3,vector<double>(3,0.0));
		double tmpVol062 = pow(vol,0.66666666);
		for(int j=0; j<3; ++j){
			for(int k=0; k<3; ++k){
				M[j][k] = -tmpVol062*Sabs*S[j][k]*vol;
				S[j][k] *= vol;
			}
		}
		
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
			

			double dUdxSten = gradU[j][0];
			double dUdySten = gradU[j][1];
			double dUdzSten = gradU[j][2];
			double dVdxSten = gradV[j][0];
			double dVdySten = gradV[j][1];
			double dVdzSten = gradV[j][2];
			double dWdxSten = gradW[j][0];
			double dWdySten = gradW[j][1];
			double dWdzSten = gradW[j][2];
			double volSten = cellSten.volume;
			
			vector<vector<double>> tmpS(3,vector<double>(3,0.0));
			tmpS[0][0] = dUdxSten;
			tmpS[0][1] = 0.5*(dUdySten+dVdxSten);
			tmpS[0][2] = 0.5*(dUdzSten+dWdxSten);
			tmpS[1][0] = tmpS[0][1];
			tmpS[1][1] = dVdySten;
			tmpS[1][2] = 0.5*(dVdzSten+dWdySten);
			tmpS[2][0] = tmpS[0][2];
			tmpS[2][1] = tmpS[1][2];
			tmpS[2][2] = dWdzSten;
			
			S[0][0] += tmpS[0][0]*volSten;
			S[0][1] += tmpS[0][1]*volSten;
			S[0][2] += tmpS[0][2]*volSten;
			S[1][0] += tmpS[1][0]*volSten;
			S[1][1] += tmpS[1][1]*volSten;
			S[1][2] += tmpS[1][2]*volSten;
			S[2][0] += tmpS[2][0]*volSten;
			S[2][1] += tmpS[2][1]*volSten;
			S[2][2] += tmpS[2][2]*volSten;
			
			double SSSten = 
				pow(tmpS[0][0],2.0) + pow(tmpS[1][1],2.0) + pow(tmpS[2][2],2.0) +
				2.0*(pow(tmpS[0][1],2.0) + pow(tmpS[0][2],2.0) + pow(tmpS[1][2],2.0));
			
			double SabsSten = sqrt(2.0*SSSten);
			
			
			double tmpVol06 = pow(volSten,0.66666666);
			M[0][0] -= tmpVol06*SabsSten*tmpS[0][0]*volSten;
			M[0][1] -= tmpVol06*SabsSten*tmpS[0][1]*volSten;
			M[0][2] -= tmpVol06*SabsSten*tmpS[0][2]*volSten;
			M[1][0] -= tmpVol06*SabsSten*tmpS[1][0]*volSten;
			M[1][1] -= tmpVol06*SabsSten*tmpS[1][1]*volSten;
			M[1][2] -= tmpVol06*SabsSten*tmpS[1][2]*volSten;
			M[2][0] -= tmpVol06*SabsSten*tmpS[2][0]*volSten;
			M[2][1] -= tmpVol06*SabsSten*tmpS[2][1]*volSten;
			M[2][2] -= tmpVol06*SabsSten*tmpS[2][2]*volSten;
			
			sumSabs += SabsSten*volSten;
			
			
			double RhoSten = cellSten.var[controls.Rho];
			double USten = cellSten.var[controls.U];
			double VSten = cellSten.var[controls.V];
			double WSten = cellSten.var[controls.W];
			
			L[0][0] += RhoSten * USten * USten * volSten;
			L[0][1] += RhoSten * USten * VSten * volSten;
			L[0][2] += RhoSten * USten * WSten * volSten;
			L[1][1] += RhoSten * VSten * VSten * volSten;
			L[1][2] += RhoSten * VSten * WSten * volSten;
			L[2][2] += RhoSten * WSten * WSten * volSten;
			
			sumRho += RhoSten * volSten;
			sumV += USten * volSten;
			sumU += VSten * volSten;
			sumW += WSten * volSten;
			sumVol += volSten;
			
		}
		
		for(int j=0; j<3; ++j){
			for(int k=0; k<3; ++k){
				M[j][k] /= sumVol;
				S[j][k] /= sumVol;
				L[j][k] /= sumVol;
			}
		}
		sumRho /= sumVol;
		sumV /= sumVol;
		sumU /= sumVol;
		sumW /= sumVol;
		sumSabs /= sumVol;
		
		
		L[0][0] -= sumRho * sumU * sumU;
		L[0][1] -= sumRho * sumU * sumV;
		L[0][2] -= sumRho * sumU * sumW;
		L[1][1] -= sumRho * sumV * sumV;
		L[1][2] -= sumRho * sumV * sumW;
		L[2][2] -= sumRho * sumW * sumW;
		L[1][0] = L[0][1];
		L[2][0] = L[0][2];
		L[2][1] = L[1][2];
		
		double sumDelta = pow(sumVol,0.66666666666);
		for(int j=0; j<3; ++j){
			for(int k=0; k<3; ++k){
				M[j][k] = sumDelta * sumSabs * S[j][k] - M[j][k];
			}
		}
		
        // ensemble average for numerator and denominator
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(int j=0; j<3; ++j){
			for(int k=0; k<3; ++k){
				sum1 += L[j][k] * M[j][k];
				sum2 += M[j][k] * M[j][k];
			}
		}
		sum1 /= 9.0;
		sum2 /= 9.0;
		
        // calculate parameter C_R
		double C_R = 0.0;
		if(sum2 != 0.0){
			C_R = -0.5 * sum1 / sum2;
		}
		
        // calculate Delta_bar
        double DeltaBar2 = pow(cell.volume,0.666666666);
		
		// calculate muT & kSGS
		cell.var[controls.muT] = C_R * cell.var[controls.Rho]*DeltaBar2*Sabs;
		
		cell.var[controls.kSGS] = C_I * DeltaBar2 * (0.5*Sabs*Sabs);
		
		
	}
	
	
	
	
}