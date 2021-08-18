#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include <ctime>




void SEMO_Solvers_Builder::calcPressureEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<double>& linAD,
	vector<double>& linAOD,
	vector<double>& linAFL,
	vector<double>& linAFR,
	vector<double>& residuals){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	vector<clock_t> startTime;
	
	SEMO_Utility_Math math;

	// // gradient P
	// vector<vector<double>> gradP;
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
	// vector<double> limGradP;
	// math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] *= limGradP[i]; 
		// }
	// }
	
	
	
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	
	
	
	

	
	// //=========================================
	// vector<double> maxXPF(mesh.cells.size(),0.0);
	// vector<double> maxS(mesh.cells.size(),0.0);	// internal faces
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double xPF = sqrt(
				// pow(face.x - mesh.cells[face.owner].x,2.0)+
				// pow(face.y - mesh.cells[face.owner].y,2.0)+
				// pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// double xNF = sqrt(
				// pow(face.x - mesh.cells[face.neighbour].x,2.0)+
				// pow(face.y - mesh.cells[face.neighbour].y,2.0)+
				// pow(face.z - mesh.cells[face.neighbour].z,2.0));
				
			// maxXPF[face.neighbour] = max(maxXPF[face.neighbour],xNF);
			
			// maxS[face.owner] = max(maxS[face.owner],face.area);
			// maxS[face.neighbour] = max(maxS[face.neighbour],face.area);
			
			
		// }
		// else{
			// double xPF = sqrt(
				// pow(face.x - mesh.cells[face.owner].x,2.0)+
				// pow(face.y - mesh.cells[face.owner].y,2.0)+
				// pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// maxS[face.owner] = max(maxS[face.owner],face.area);
			
		// }
	// }
	// //=========================================
	
	
	// //=========================================
	// vector<vector<double>> gradPLSQ;
	// vector<vector<double>> gradPGG;
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradPLSQ);
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradPGG);
	// vector<vector<double>> gradP;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		// double beta0 = min(1.0, 2.0/AR);
		
		// vector<double> tmpVec;
		// tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][0]);
		// tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][1]);
		// tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][2]);
		
		// gradP.push_back(tmpVec);
		
	// }
	// //=========================================
	
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradPx_send;
		// vector<double> gradPy_send;
		// vector<double> gradPz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradPx_send.push_back(gradP[face.owner][0]);
				// gradPy_send.push_back(gradP[face.owner][1]);
				// gradPz_send.push_back(gradP[face.owner][2]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatasDouble(
					// gradPx_send, gradPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatasDouble(
					// gradPy_send, gradPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatasDouble(
					// gradPz_send, gradPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	// }
	



	// vector<vector<double>> gradU;
	// math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// vector<vector<double>> gradV;
	// math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// vector<vector<double>> gradW;
	// math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);


	

	// // calc gauss-green gradient
	// vector<vector<double>> gradGrav(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double RhoF = 2.0/(1.0/face.varL[controls.fRho]+1.0/face.varR[controls.fRho]);
			// double RhoF = 0.5*(face.varL[controls.fRho]+face.varR[controls.fRho]);
			// gradGrav[face.owner][1] += 
				// 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.owner].volume;
			// gradGrav[face.neighbour][1] -= 
				// 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.neighbour].volume;
		// }
		// else{
			// // double RhoF = 2.0/(1.0/face.varL[controls.fRho]+1.0/face.varR[controls.fRho]);
			// double RhoF = 0.5*(face.varL[controls.fRho]+face.varR[controls.fRho]);
			// gradGrav[face.owner][1] += 
				// 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.owner].volume;
		// }
	// }




	vector<double> cellVolume_recv;
	if(size>1){
		vector<double> cellVolume_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cellVolume_send.push_back(mesh.cells[face.owner].volume);
			}
		}
		mpi.setProcsFaceDatas(
					cellVolume_send, cellVolume_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		cellVolume_send.clear();
		
	}





	// linAD
	vector<double> linAD_send;
	vector<double> linAD_recv;
	if(size>1){
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				linAD_send.push_back(linAD[face.owner]);
			}
		}
		
		mpi.setProcsFaceDatas(
					linAD_send, linAD_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		linAD_send.clear();
		
	}







	vector<double> linA;
	vector<double> linB;
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		
		linA.push_back(0.0);
		linB.push_back(0.0);
			
	}

	vector<double> linAL;
	vector<double> linAR;
	
	int proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		double area = face.area;
		vector<double> nvec(3,0.0);
		nvec[0] = face.unitNormals[0];
		nvec[1] = face.unitNormals[1];
		nvec[2] = face.unitNormals[2];
		
		vector<double> distanceCells(3,0.0);
		distanceCells[0] = face.distCells[0];
		distanceCells[1] = face.distCells[1];
		distanceCells[2] = face.distCells[2];
		
		double wCL = face.wC;
		double wCR = 1.0-face.wC;
		
		double UL = face.varL[controls.fU];
		double VL = face.varL[controls.fV];
		double WL = face.varL[controls.fW];
		double PL = face.varL[controls.fP];
		double RhoL = face.varL[controls.fRho];
		
		double UR = face.varR[controls.fU];
		double VR = face.varR[controls.fV];
		double WR = face.varR[controls.fW];
		double PR = face.varR[controls.fP];
		double RhoR = face.varR[controls.fRho];
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  pow(distanceCells[1],2.0) + 
						  pow(distanceCells[2],2.0));
		// double dPNEff = dPN 
			// /( (nvec[0]*distanceCells[0] + 
			    // nvec[1]*distanceCells[1] + 
				// nvec[2]*distanceCells[2])/dPN );
				
		double nonOrtholimiter = 1.0;
		double cosAlpha = dPN_e/dPN;
		if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			nonOrtholimiter = 0.5;
			// nonOrtholimiter = 0.333;
		}
		else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			nonOrtholimiter = 0.333;
			// nonOrtholimiter = 0.1;
			// nonOrtholimiter = 0.4;
		}
		else if( cosAlpha < 0.342 ){
			nonOrtholimiter = 0.0;
			// nonOrtholimiter = 0.3;
		}
		// nonOrtholimiter = 0.333;
		
		// cout << cosAlpha << " " << nonOrtholimiter << endl;
				
		double Ef[3];
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
				
		// over-relaxed approach
		// double overRelaxCoeff = 1.0;
		// double overRelaxCoeff = Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2];
		// double overRelaxCoeff = 1.0 / (Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2]);
		
		// non-orthogonal, over-relaxed approach
		// double Tf[3];
		// Tf[0] = nvec[0] - Ef[0]*overRelaxCoeff;
		// Tf[1] = nvec[1] - Ef[1]*overRelaxCoeff;
		// Tf[2] = nvec[2] - Ef[2]*overRelaxCoeff;
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		// double dPNEff = dPN / overRelaxCoeff;
		// double dPNEff = dPN;
		
		// double coeffP = area/dPNEff * dPN/dPNEff;
		
		double UnF = wCL*UnL+wCR*UnR;
		// double UnF = 0.5*UnL+0.5*UnR;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			UnF = face.var[controls.Un];
		}
		
		
		// double RhieChowCoeff = 0.8;
		// double RhieChowCoeff = 1.0;
		
		
		double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// tmp1 = (wCL/linAD[face.owner] + wCR/linAD[face.neighbour]);
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// tmp1 = (wCL/linAD[face.owner] + wCR/linAD_recv[proc_num]);
			// ++proc_num;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// tmp1 = 1.0/linAD[face.owner];
		// }
		
		
		// double tmp1 = (wCL/RhoL + wCR/RhoR);
		
		
		// // double tmp1 = (0.5/RhoL + 0.5/RhoR)*controls.timeStep;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double orgPL = mesh.cells[face.owner].var[controls.P];
			// double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			// // tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			
			// // UnF = 
				// // wCL*( 
				// // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (UR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][0])*nvec[0] +
				// // (VR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][1])*nvec[1] +
				// // (WR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][2])*nvec[2] );
				
			// UnF += RhieChowCoeff * (
				// wCL*( 
				// (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// wCR*( 
				// (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // UnF += RhieChowCoeff * (
				// // 0.5*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // 0.5*( 
				// // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // UnF = UnF - tmp1*(orgPR-orgPL)/dPN;
			// // UnF = UnF - tmp1*(orgPR-orgPL)/dPN * overRelaxCoeff;
			// UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			
			// // non-orthogonal, over-relaxed approach
			// double gradPf[3];
			// gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			// // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradP[face.neighbour][0]);
			// // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradP[face.neighbour][1]);
			// // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradP[face.neighbour][2]);
			
			// // gradPf[0] = (newWCL*gradP[face.owner][0]+newWCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (newWCL*gradP[face.owner][1]+newWCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (newWCL*gradP[face.owner][2]+newWCR*gradP[face.neighbour][2]);
			
			// // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // double cellVolL = mesh.cells[face.owner].volume;
			// // double cellVolR = cellVolume_recv[proc_num];
			
			// // // tmp1 = (cellVolL+cellVolR)/
				// // // (cellVolL*RhoL/controls.timeStep + cellVolR*RhoR/controls.timeStep);
				
			// double orgPL = mesh.cells[face.owner].var[controls.P];
			// double orgPR = PR;
			
			// // tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			
			// // UnF = 
				// // wCL*( 
				// // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (UR + 1.0/linAD_recv[proc_num]*gradPx_recv[proc_num])*nvec[0] +
				// // (VR + 1.0/linAD_recv[proc_num]*gradPy_recv[proc_num])*nvec[1] +
				// // (WR + 1.0/linAD_recv[proc_num]*gradPz_recv[proc_num])*nvec[2] );
			
			// UnF += RhieChowCoeff * (
				// wCL*( 
				// (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// wCR*( 
				// (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			// // UnF += RhieChowCoeff * (
				// // 0.5*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // 0.5*( 
				// // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
				
			// // UnF = UnF - tmp1*(orgPR-orgPL)/dPNEff;
			// UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// // non-orthogonal, over-relaxed approach
			// double gradPf[3];
			// gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			// // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradPx_recv[proc_num]);
			// // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradPy_recv[proc_num]);
			// // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradPz_recv[proc_num]);
			// // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			// ++proc_num;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// // tmp1 = 1.0/linAD[face.owner];
			// tmp1 = controls.timeStep/RhoL;
			
		// }
		
		
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double RhoF = RhoL*weightL + RhoR*weightR;
		double RhoF = 1.0;//RhoL*weightL + RhoR*weightR;
		
		// linAL.push_back(-tmp1/dPN * area*overRelaxCoeff);
		// linAR.push_back(-tmp1/dPN * area*overRelaxCoeff);
		
		linAL.push_back(-tmp1/dPN_e * area);
		linAR.push_back(-tmp1/dPN_e * area);
		
		// linAL.back() += (-tmp1/dPN * area*overRelaxCoeff * (nvec[0]*Tf[0]+nvec[1]*Tf[1]+nvec[2]*Tf[2]));
		// linAR.back() += (-tmp1/dPN * area*overRelaxCoeff * (nvec[0]*Tf[0]+nvec[1]*Tf[1]+nvec[2]*Tf[2]));
		
		// linAL.push_back(-0.5*tmp1*coeffP);
		// linAR.push_back(-0.5*tmp1*coeffP);
			
		// convective term
		// linB[face.owner] -= RhoF*UnF*area;
			
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			linA[face.owner]     -= linAL[i];
			linA[face.neighbour] -= linAR[i];
			
			// convective term
			linB[face.owner] -= RhoF*UnF*area;
			linB[face.neighbour] += RhoF*UnF*area;
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			linA[face.owner] -= linAL[i];
			
			// convective term
			linB[face.owner] -= RhoF*UnF*area;
		}
		
	}
	
	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			if(boundary.type[0] == "fixedValue"){
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec(3,0.0);
					nvec[0] = face.unitNormals[0];
					nvec[1] = face.unitNormals[1];
					nvec[2] = face.unitNormals[2];
					
					vector<double> distanceCells(3,0.0);
					distanceCells[0] = face.distCells[0];
					distanceCells[1] = face.distCells[1];
					distanceCells[2] = face.distCells[2];
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
					
					double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
					double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
					
					double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
					double dPN = sqrt(pow(distanceCells[0],2.0) + 
									  pow(distanceCells[1],2.0) + 
									  pow(distanceCells[2],2.0));
								
					double UnF = 0.5*UnL+0.5*UnR;
				
					double RhoL = face.varL[controls.fRho];
					double tmp1 = controls.timeStep/RhoL;
					// double tmp1 = 1.0/linAD[face.owner];
							
					linA[face.owner] += tmp1*area/dPN_e;
					
					double RhoF = 1.0;
					
					// // convective term
					// double orgPL = mesh.cells[face.owner].var[controls.P];
					// double PR = face.varR[controls.fP];
					
					// // UnF = 
						// // ( 
						// // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
						// // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
						// // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] );
					// UnF = 
						// ( 
						// (UL + controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
						// (VL + controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
						// (WL + controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] );
						
					// double corrPL = orgPL;
					// UnF = UnF - tmp1*(PR-corrPL)/dPN_e;
					
					linB[face.owner] -= RhoF*UnF*area;
				}
			}
			else{
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec(3,0.0);
					nvec[0] = face.unitNormals[0];
					nvec[1] = face.unitNormals[1];
					nvec[2] = face.unitNormals[2];
					
					vector<double> distanceCells(3,0.0);
					distanceCells[0] = face.distCells[0];
					distanceCells[1] = face.distCells[1];
					distanceCells[2] = face.distCells[2];
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
					
					double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
					double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
					
					double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
					
					double UnF = 0.5*UnL+0.5*UnR;
				
					double RhoL = face.varL[controls.fRho];
					// double tmp1 = controls.timeStep/RhoL;
					double RhoF = 1.0;
					
					// // convective term
					// double orgPL = mesh.cells[face.owner].var[controls.P];
					// double PR = orgPL 
						// + gradP[face.owner][0]*distanceCells[0]
						// + gradP[face.owner][1]*distanceCells[1]
						// + gradP[face.owner][2]*distanceCells[2];
					
					// // UnF = 
						// // ( 
						// // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
						// // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
						// // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] );
					// UnF = 
						// ( 
						// (UL + controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
						// (VL + controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
						// (WL + controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] );
						
					// double corrPL = orgPL;
					// UnF = UnF - tmp1*(PR-corrPL)/dPN_e;
					
					// convective term
					linB[face.owner] -= RhoF*UnF*area;
				}
			}
			
		}
		
	}
	
	// linear solver : PETSc library
	vector<double> resiVar(mesh.cells.size(),0.0);
	
	if(controls.iterPBs != controls.iterPBsMax-1){
	
		solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			controls.solverP, controls.toleranceP, 
			controls.relTolP, controls.preconditionerP,
			controls.maxIterP);
			
		// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverP, controls.toleranceP, 
			// controls.relTolP, controls.preconditionerP,
			// controls.maxIterP);
		
	}
	else{
	
		solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			controls.solverFinalP, controls.toleranceFinalP, 
			controls.relTolFinalP, controls.preconditionerFinalP,
			controls.maxIterFinalP);
			
		// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalP, controls.toleranceFinalP, 
			// controls.relTolFinalP, controls.preconditionerFinalP,
			// controls.maxIterFinalP);
	}
	

	

	
	vector<vector<double>> gradResiP;
	// math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare(mesh, 0, resiVar, gradResiP);
	math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
	
	// math.gradientGGBoundaryTreatment(mesh, controls, gradResiP);
	
	// //==========================================
	// vector<vector<double>> gradResiPLSQ;
	// vector<vector<double>> gradResiPGG;
	// math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiPLSQ);
	// math.calcGaussGreen(mesh, 0, resiVar, gradResiPGG);
	// vector<vector<double>> gradResiP;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		// double beta0 = min(1.0, 2.0/AR);
		
		// vector<double> tmpVec;
		// tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][0]);
		// tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][1]);
		// tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][2]);
		
		// gradResiP.push_back(tmpVec);
		
	// }
	// //==========================================
	
	
	// update
	double test0 = 0.0;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.P] += controls.prePreURF*resiVar[i];
		
		// double resiU = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][0];
		// double resiV = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][1];
		// double resiW = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][2];
		// double resiU = controls.preVelURF/linAD[i]*gradResiP[i][0];
		// double resiV = controls.preVelURF/linAD[i]*gradResiP[i][1];
		// double resiW = controls.preVelURF/linAD[i]*gradResiP[i][2];
		double resiU = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0];
		double resiV = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1];
		double resiW = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2];
		
		cell.var[controls.U] -= resiU;
		cell.var[controls.V] -= resiV;
		cell.var[controls.W] -= resiW;
		
		

		if( cell.var[controls.P] <= controls.minP ) 
			cell.var[controls.P] = controls.minP;
		if( cell.var[controls.P] >= controls.maxP ) 
			cell.var[controls.P] = controls.maxP;
		
		if( cell.var[controls.U] <= controls.minU ) 
			cell.var[controls.U] = controls.minU;
		if( cell.var[controls.U] >= controls.maxU ) 
			cell.var[controls.U] = controls.maxU;
		
		if( cell.var[controls.V] <= controls.minV ) 
			cell.var[controls.V] = controls.minV;
		if( cell.var[controls.V] >= controls.maxV ) 
			cell.var[controls.V] = controls.maxV;
		
		if( cell.var[controls.W] <= controls.minW ) 
			cell.var[controls.W] = controls.minW;
		if( cell.var[controls.W] >= controls.maxW ) 
			cell.var[controls.W] = controls.maxW;
		
		test0 += controls.prePreURF*resiVar[i];
		// test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0];
		// test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1];
		// test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2];
		
		
		// cell.var[controls.UDV[0]] = resiVar[i];
		// cell.var[controls.UDV[1]] = gradResiP[i][1];
		
		
		// cell.var.at(controls.UDV[0]) += gradP[i][0];
		// cell.var.at(controls.UDV[1]) += gradP[i][1];
		// cell.var.at(controls.UDV[2]) += gradP[i][2];
		
		// cell.var[controls.U] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][0];
		// cell.var[controls.V] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][1];
		// cell.var[controls.W] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][2];
		
		
				
		residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		residuals[1] += pow(resiU,2.0)*cell.volume;
		residuals[2] += pow(resiV,2.0)*cell.volume;
		residuals[3] += pow(resiW,2.0)*cell.volume;
		
	}
	
	
	
	// cout << "A : " << test0 << endl;
	
	


	vector<double> gradResiPx_recv;
	vector<double> gradResiPy_recv;
	vector<double> gradResiPz_recv;
	vector<double> resiVar_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> gradResiPx_send;
		vector<double> gradResiPy_send;
		vector<double> gradResiPz_send;
		vector<double> resiVar_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				gradResiPx_send.push_back(gradResiP[face.owner][0]);
				gradResiPy_send.push_back(gradResiP[face.owner][1]);
				gradResiPz_send.push_back(gradResiP[face.owner][2]);
				resiVar_send.push_back(resiVar[face.owner]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					gradResiPx_send, gradResiPx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradResiPy_send, gradResiPy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradResiPz_send, gradResiPz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		
		mpi.setProcsFaceDatas(
					resiVar_send, resiVar_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	}
	
	
	
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
	
		double area = face.area;
		vector<double> nvec(3,0.0);
		nvec[0] = face.unitNormals[0];
		nvec[1] = face.unitNormals[1];
		nvec[2] = face.unitNormals[2];
		
		vector<double> distanceCells(3,0.0);
		distanceCells[0] = face.distCells[0];
		distanceCells[1] = face.distCells[1];
		distanceCells[2] = face.distCells[2];
		
		double wCL = face.wC;
		double wCR = 1.0-face.wC;
		
		double RhoL = face.varL[controls.fRho];
		
		double RhoR = face.varR[controls.fRho];
		
		double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  pow(distanceCells[1],2.0) + 
						  pow(distanceCells[2],2.0));
				
		double nonOrtholimiter = 1.0;
		double cosAlpha = dPN_e/dPN;
		if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			nonOrtholimiter = 0.5;
		}
		else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			nonOrtholimiter = 0.333;
		}
		else if( cosAlpha < 0.342 ){
			nonOrtholimiter = 0.0;
		}
				
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// tmp1 = (wCL/linAD[face.owner] + wCR/linAD[face.neighbour]);
			
			face.var[controls.Un] -= 
				controls.preVelURF * tmp1*(resiVar[face.neighbour]-resiVar[face.owner])/dPN_e;
				
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradResiP[face.owner][0]+wCR*gradResiP[face.neighbour][0]);
			gradPf[1] = (wCL*gradResiP[face.owner][1]+wCR*gradResiP[face.neighbour][1]);
			gradPf[2] = (wCL*gradResiP[face.owner][2]+wCR*gradResiP[face.neighbour][2]);
			
			face.var[controls.Un] -= controls.preVelURF * 
				nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// tmp1 = (wCL/linAD[face.owner] + wCR/linAD_recv[proc_num]);
			
			face.var[controls.Un] -= 
				controls.preVelURF * tmp1*(resiVar_recv[proc_num]-resiVar[face.owner])/dPN_e;
				
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradResiP[face.owner][0]+wCR*gradResiPx_recv[proc_num]);
			gradPf[1] = (wCL*gradResiP[face.owner][1]+wCR*gradResiPy_recv[proc_num]);
			gradPf[2] = (wCL*gradResiP[face.owner][2]+wCR*gradResiPz_recv[proc_num]);
			
			face.var[controls.Un] -= controls.preVelURF * 
				nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			++proc_num;
			
		}
		
	}
	
	
	
	
	
	
	
	linA.clear();
	linAL.clear();
	linAR.clear();
	linB.clear();
	resiVar.clear();
	
	
}














// // add Rho
// void SEMO_Solvers_Builder::calcPressureEq(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& residuals){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	// vector<clock_t> startTime;
	
	// SEMO_Utility_Math math;

	// // gradient P
	// vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// double PF = 0.5*(face.varL[controls.fP]+face.varR[controls.fP]);
		
		// for(int j=0; j<3; ++j){
			// gradP[face.owner][j] += PF*face.unitNormals[j]*face.area;
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// gradP[face.neighbour][j] -= PF*face.unitNormals[j]*face.area;
			// }
		// }
	// }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] /= mesh.cells[i].volume;
		// }
	// }
	// // processor faces
	// // gradP , 
	// vector<double> gradPx_send;
	// vector<double> gradPy_send;
	// vector<double> gradPz_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// gradPx_send.push_back(gradP[face.owner][0]);
			// gradPy_send.push_back(gradP[face.owner][1]);
			// gradPz_send.push_back(gradP[face.owner][2]);
		// }
	// }
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// mpi.setProcsFaceDatasDouble(
				// gradPx_send, gradPx_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPy_send, gradPy_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPz_send, gradPz_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// gradPx_send.clear();
	// gradPy_send.clear();
	// gradPz_send.clear();
	


	// // linAD
	// vector<double> linAD_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// linAD_send.push_back(linAD[face.owner]);
		// }
	// }
	
	// vector<double> linAD_recv;
	// mpi.setProcsFaceDatasDouble(
				// linAD_send, linAD_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// linAD_send.clear();
	
	
	
	// vector<double> linA;
	// vector<double> linB;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// linA.push_back(0.0);
		// // linB.push_back(0.0);
		// linB.push_back(-
			// (cell.var[controls.Rho]-cell.var[controls.oldRho])
			// *cell.volume/controls.timeStep);
			
	// }


	// vector<double> linAL;
	// vector<double> linAR;
	
	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double area = face.area;
		// vector<double> nvec(3,0.0);
		// nvec[0] = face.unitNormals[0];
		// nvec[1] = face.unitNormals[1];
		// nvec[2] = face.unitNormals[2];
		
		// vector<double> distanceCells(3,0.0);
		// distanceCells[0] = face.distCells[0];
		// distanceCells[1] = face.distCells[1];
		// distanceCells[2] = face.distCells[2];
		
		// double UL = face.varL[controls.fU];
		// double VL = face.varL[controls.fV];
		// double WL = face.varL[controls.fW];
		// double PL = face.varL[controls.fP];
		// double RhoL = face.varL[controls.fRho];
		
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		// double PR = face.varR[controls.fP];
		// double RhoR = face.varR[controls.fRho];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double dPN = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		// double coeffP = area/dPN;
		
		// double UnF = 0.5*(UnL+UnR);
		
		
		// double coeff1 = 0.5*(1.0/RhoL + 1.0/RhoR)*controls.timeStep;
		// double tmp1 = controls.timeStep;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// tmp1 = 0.5*(RhoL/linAD[face.owner] + RhoR/linAD[face.neighbour]);
			// UnF = UnF 
				// - 0.5*(1.0/linAD[face.owner] + 1.0/linAD[face.neighbour])
				// *((PR-PL)/dPN - 0.5*(
				 // gradP[face.owner][0]*nvec[0]+gradP[face.neighbour][0]*nvec[0]
				// +gradP[face.owner][1]*nvec[1]+gradP[face.neighbour][1]*nvec[1]
				// +gradP[face.owner][2]*nvec[2]+gradP[face.neighbour][2]*nvec[2]
				// ));
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// tmp1 = 0.5*(RhoL/linAD[face.owner] + RhoR/linAD_recv[proc_num]);
			// UnF = UnF 
				// - 0.5*(1.0/linAD[face.owner] + 1.0/linAD_recv[proc_num])
				// *((PR-PL)/dPN - 0.5*(
				 // gradP[face.owner][0]*nvec[0]+gradPx_recv[proc_num]*nvec[0]
				// +gradP[face.owner][1]*nvec[1]+gradPy_recv[proc_num]*nvec[1]
				// +gradP[face.owner][2]*nvec[2]+gradPz_recv[proc_num]*nvec[2]
				// ));
			// ++proc_num;
		// }
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double RhoF = RhoL*weightL + RhoR*weightR;
		
		// linAL.push_back(-tmp1*coeffP);
		// linAR.push_back(-tmp1*coeffP);
			
		// // convective term
		// linB[face.owner] -= RhoF*UnF*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA[face.owner]     -= linAL[i];
			// linA[face.neighbour] -= linAR[i];
			
			// // convective term
			// linB[face.neighbour] += RhoF*UnF*area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA[face.owner] -= linAL[i];
			
		// }
		
	// }
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(boundary.type[0] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec(3,0.0);
					// nvec[0] = face.unitNormals[0];
					// nvec[1] = face.unitNormals[1];
					// nvec[2] = face.unitNormals[2];
					
					// vector<double> distanceCells(3,0.0);
					// distanceCells[0] = face.distCells[0];
					// distanceCells[1] = face.distCells[1];
					// distanceCells[2] = face.distCells[2];
				
					// double RhoL = face.varL[controls.fRho];
				
					// double tmp1 = controls.timeStep;
					
					// double dPN = nvec[0]*distanceCells[0] 
							   // + nvec[1]*distanceCells[1] 
							   // + nvec[2]*distanceCells[2];
							   
					// double coeffP = area/dPN;
					
					// linA[face.owner] += tmp1*coeffP;
				// }
			// }
			
		// }
		
	// }
	
	// // linear solver : PETSc library
	// vector<double> resiVar(mesh.cells.size(),0.0);
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverP, controls.toleranceP, 
			// controls.relTolP, controls.preconditionerP,
			// controls.maxIterP);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalP, controls.toleranceFinalP, 
			// controls.relTolFinalP, controls.preconditionerFinalP,
			// controls.maxIterFinalP);
	// }
	

	
	
	// // velocities correction
	// vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double UnL = face.varL[controls.fU]*face.unitNormals[0] +
			             // face.varL[controls.fV]*face.unitNormals[1] +
						 // face.varL[controls.fW]*face.unitNormals[2];
			// double UnR = face.varR[controls.fU]*face.unitNormals[0] +
			             // face.varR[controls.fV]*face.unitNormals[1] +
						 // face.varR[controls.fW]*face.unitNormals[2];
						 
			// double PF = 0.5*(resiVar[face.owner]+resiVar[face.neighbour]);
			
			// for(int j=0; j<3; ++j){
				// gradResiP[face.owner][j] += 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// gradResiP[face.neighbour][j] -= 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			// }
		// }
	// }

	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(boundary.type[0] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double PF = 0.5*(resiVar[face.owner]);
					// for(int j=0; j<3; ++j){
						// gradResiP[face.owner][j] += 
							// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// }
					
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double PF = resiVar[face.owner];
					// for(int j=0; j<3; ++j){
						// gradResiP[face.owner][j] += 
							// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// }
					
				// }
			// }
		// }
	// }
	
	// // processor faces
	// vector<double> sendValues;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// sendValues.push_back(resiVar[face.owner]);
		// }
	// }
	// vector<double> recvValues;
	// mpi.setProcsFaceDatasDouble(
				// sendValues, recvValues,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// int num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
			// double PF = 0.5*(resiVar[face.owner]+recvValues[num]);
			
			// for(int j=0; j<3; ++j){
				// gradResiP[face.owner][j] += 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			// }
			// ++num;
		// }
	// }
	

	
	// // // calc gradient
	// // vector<vector<double>> gradResiP;
	
	// // gradResiP.clear();
	// // gradResiP.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	
	// // for(auto& face : mesh.faces){
		
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double DVar = 
				// // resiVar[face.neighbour] - resiVar[face.owner];
			
			// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
			// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
			// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
			
			// // gradResiP[face.neighbour][0] += face.distCells[0] * DVar;
			// // gradResiP[face.neighbour][1] += face.distCells[1] * DVar;
			// // gradResiP[face.neighbour][2] += face.distCells[2] * DVar;
			
		// // }
		
	// // }
	
	
	// // // boundary
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // if(boundary.type[0] == "fixedValue"){
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double DVar = 
						// // 0.0 - resiVar[face.owner];
					
					// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
					// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
					// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
					
					
				// // }
			// // }
		// // }
	// // }
	
	
	// // // processor faces
	// // vector<double> sendValues;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // sendValues.push_back(resiVar[face.owner]);
		// // }
	// // }
	
	// // vector<double> recvValues;
	// // mpi.setProcsFaceDatasDouble(
				// // sendValues, recvValues,
				// // mesh.countsProcFaces, mesh.countsProcFaces, 
				// // mesh.displsProcFaces, mesh.displsProcFaces);

	// // int num=0;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // double DVar = 
				// // recvValues[num] - resiVar[face.owner];
			
			// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
			// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
			// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
					
			// // ++num;
		// // }
	// // }
	
	
	
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // double tmp0 = 
			// // cell.coeffLeastSquare[0] * gradResiP[i][0] +
			// // cell.coeffLeastSquare[1] * gradResiP[i][1] +
			// // cell.coeffLeastSquare[2] * gradResiP[i][2];
			
		// // double tmp1 = 
			// // cell.coeffLeastSquare[1] * gradResiP[i][0] +
			// // cell.coeffLeastSquare[3] * gradResiP[i][1] +
			// // cell.coeffLeastSquare[4] * gradResiP[i][2];
			
		// // double tmp2 = 
			// // cell.coeffLeastSquare[2] * gradResiP[i][0] +
			// // cell.coeffLeastSquare[4] * gradResiP[i][1] +
			// // cell.coeffLeastSquare[5] * gradResiP[i][2];
			
		// // gradResiP[i][0] = tmp0;
		// // gradResiP[i][1] = tmp1;
		// // gradResiP[i][2] = tmp2;
		
	// // }
	
	
	// // update
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// cell.var[controls.P] += controls.prePreURF*resiVar[i];
		// // cell.var[controls.U] -= controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][0];
		// // cell.var[controls.V] -= controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][1];
		// // cell.var[controls.W] -= controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][2];
		// cell.var[controls.U] -= controls.preVelURF/linAD[i]*gradResiP[i][0];
		// cell.var[controls.V] -= controls.preVelURF/linAD[i]*gradResiP[i][1];
		// cell.var[controls.W] -= controls.preVelURF/linAD[i]*gradResiP[i][2];
		
				
		// residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		// residuals[1] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0],2.0)*cell.volume;
		// residuals[2] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1],2.0)*cell.volume;
		// residuals[3] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2],2.0)*cell.volume;
		
	// }
	
	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	
// }















// void SEMO_Solvers_Builder::calcPressureEq(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& residuals){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	// vector<clock_t> startTime;
	
	// // if(rank==0) startTime.push_back(clock());
	
	
	

	// SEMO_Utility_Math math;
	// vector<vector<double>> gradP;
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	

	// // processor faces
	// // gradP , 
	// vector<double> gradPx_send;
	// vector<double> gradPy_send;
	// vector<double> gradPz_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// gradPx_send.push_back(gradP[face.owner][0]);
			// gradPy_send.push_back(gradP[face.owner][1]);
			// gradPz_send.push_back(gradP[face.owner][2]);
		// }
	// }
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// mpi.setProcsFaceDatasDouble(
				// gradPx_send, gradPx_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPy_send, gradPy_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPz_send, gradPz_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// gradPx_send.clear();
	// gradPy_send.clear();
	// gradPz_send.clear();


	// // volume, linAD
	// vector<double> volume_send;
	// vector<double> linAD_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// volume_send.push_back(mesh.cells[face.owner].volume);
			// linAD_send.push_back(linAD[face.owner]);
		// }
	// }
	
	// vector<double> volume_recv;
	// vector<double> linAD_recv;
	// mpi.setProcsFaceDatasDouble(
				// volume_send, volume_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// linAD_send, linAD_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// volume_send.clear();
	// linAD_send.clear();
	
	
	
	// vector<double> linA(mesh.cells.size(),0.0);
	// vector<double> linB(mesh.cells.size(),0.0);

	// vector<double> linAL;
	// vector<double> linAR;
	
	
	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double area = face.area;
		// vector<double> nvec(3,0.0);
		// nvec[0] = face.unitNormals[0];
		// nvec[1] = face.unitNormals[1];
		// nvec[2] = face.unitNormals[2];
		
		// vector<double> distanceCells(3,0.0);
		// distanceCells[0] = face.distCells[0];
		// distanceCells[1] = face.distCells[1];
		// distanceCells[2] = face.distCells[2];
		
		
		
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // distanceCells[0] = mesh.cells[face.neighbour].x-mesh.cells[face.owner].x;
			// // distanceCells[1] = mesh.cells[face.neighbour].y-mesh.cells[face.owner].y;
			// // distanceCells[2] = mesh.cells[face.neighbour].z-mesh.cells[face.owner].z;
		// // }
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // distanceCells[0] = face.x-mesh.cells[face.owner].x;
			// // distanceCells[1] = face.y-mesh.cells[face.owner].y;
			// // distanceCells[2] = face.z-mesh.cells[face.owner].z;
			// // distanceCells[0] *= 2.0;
			// // distanceCells[1] *= 2.0;
			// // distanceCells[2] *= 2.0;
		// // }
	
		// // if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// // distanceCells[0] = face.x-mesh.cells[face.owner].x;
			// // distanceCells[1] = face.y-mesh.cells[face.owner].y;
			// // distanceCells[2] = face.z-mesh.cells[face.owner].z;
			// // distanceCells[0] *= 2.0;
			// // distanceCells[1] *= 2.0;
			// // distanceCells[2] *= 2.0;
		// // }
		
		
		
		// double UL = face.varL[controls.fU];
		// double VL = face.varL[controls.fV];
		// double WL = face.varL[controls.fW];
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		
		// double PL = face.varL[controls.fP];
		// double PR = face.varR[controls.fP];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double RhoL = face.varL[controls.fRho];
		// double RhoR = face.varR[controls.fRho];
		
	
		
		// // double dPN = 
			// // sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
			                        // // std::begin(distanceCells), 0.0));
									
		// double dPN = sqrt(pow(distanceCells[0],2.0)+pow(distanceCells[1],2.0)+pow(distanceCells[2],2.0));
		// double dPNEff = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		// // dPNEff = area/dPN *(dPNEff/dPN);
		// dPNEff = area/dPN;
		
		

		
		// double UnF = 0.5*(UnL+UnR);
		
		
		

	// // // SLAU
	// // double chat = 10000.0;
	// // double ML = UnL/chat; 
	// // double MR = UnR/chat;
	// // double KLR = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR));
	
    // // double Mbar = ( RhoL*abs(ML)+RhoR*abs(MR) ) / ( RhoL + RhoR );
	// // double Mcy = min(1.0,KLR/chat);
	// // double phi_c = pow((1.0-Mcy),2.0);
	// // double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );

    // // double D_L = ML+(1.0-g_c)*abs(ML);
    // // double D_R = MR-(1.0-g_c)*abs(MR);
    // // double D_rho = Mbar*g_c;
	
	// // double MLP_SLAU = 0.5*(D_L+D_rho);
	// // double MRM_SLAU = 0.5*(D_R-D_rho);
	// // double mdot = RhoL*chat*MLP_SLAU + RhoR*chat*MRM_SLAU
	   // // - 0.5
		 // // *(1.0-0.5*(1.0-cos(3.141592*min(1.0,max(abs(ML),abs(MR))))))
		 // // *(1.0-0.5*(1.0+cos(3.141592*min(abs(PL/PR),abs(PR/PL)))))
		 // // *(PR-PL)
		 // // /chat;
		 
	
	// // // // double MLP  = M_func(ML,1.0,0.0); 
	// // // // double MRM = M_func(MR,-1.0,0.0);
	// // // double preP = pre_func(ML,1.0,0.0); 
	// // // double preM = pre_func(MR,-1.0,0.0);
	// // // double gam = 0.6;
	// // // double PLR = 0.5*(PL+PR) 
				// // // - 1.0*(KLR/chat)*0.5*preP*preM*0.5*(PL+PR)/chat*(UnR-UnL)
				// // // + max(0.2,gam)*(KLR/chat)*0.5*(PL+PR)*(preP+preM-1.0)
				// // // - 0.5*(preP-preM)*(PR-PL);
				
				
	// // UnF = mdot / (0.5*(RhoL+RhoR));
	// // // tmp1 = PLR;
		
		
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		

		// // Rhie-Chow interpolation
		// double tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
		// // double tmp1 = 1.0/(0.5*(RhoL+RhoR))*controls.timeStep;
		// // double tmp1 = 1.0/RhoR*controls.timeStep;
		// // double tmp1 = (weightL/RhoL+weightR/RhoR)*controls.timeStep;
		// // if(abs(UnF) <= 1.e-200) tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
		
		// // if( face.getType() == SEMO_Types::INTERNAL_FACE ){
			// // tmp1 = 0.5*
				// // (mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // +mesh.cells[face.neighbour].volume/linAD[face.neighbour]/RhoR);
			// // // tmp1 = weightL*mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // // +weightR*mesh.cells[face.neighbour].volume/linAD[face.neighbour]/RhoR;
			
			// // // double tmp0;
			// // // tmp0  = gradP[face.owner][0]*nvec[0]+gradP[face.neighbour][0]*nvec[0];
			// // // tmp0 += gradP[face.owner][1]*nvec[1]+gradP[face.neighbour][1]*nvec[1];
			// // // tmp0 += gradP[face.owner][2]*nvec[2]+gradP[face.neighbour][2]*nvec[2];
			// // // UnF = UnF - tmp1*((PR-PL)/dPN - 0.5*tmp0);
		// // }
		// // else if( face.getType() == SEMO_Types::PROCESSOR_FACE ){
			// // tmp1 = 0.5*
				// // (mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // +volume_recv[proc_num]/linAD_recv[proc_num]/RhoR);
			// // // tmp1 = weightL*mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // // +weightR*volume_recv[proc_num]/linAD_recv[proc_num]/RhoR;
				
			// // // double tmp0;
			// // // tmp0  = gradP[face.owner][0]*nvec[0]+gradPx_recv[proc_num]*nvec[0];
			// // // tmp0 += gradP[face.owner][1]*nvec[1]+gradPy_recv[proc_num]*nvec[1];
			// // // tmp0 += gradP[face.owner][2]*nvec[2]+gradPz_recv[proc_num]*nvec[2];
			// // // UnF = UnF - tmp1*((PR-PL)/dPN - 0.5*tmp0);
			
			// // // ++proc_num;
		// // }
		// // else{
			// // tmp1 = mesh.cells[face.owner].volume/linAD[face.owner]/RhoL;
			
		// // }
		
		// // double tmp1 = controls.timeStep;
		
		
		// // // double dPNEff = area/dPN/dPN*
			// // // (std::inner_product(std::begin(nvec), std::end(nvec), 
			                    // // // std::begin(distanceCells), 0.0));

		// // UnF -= tmp1*(PR-PL)/dPN;
		
		// linAL.push_back(-tmp1*dPNEff);
		// linAR.push_back(-tmp1*dPNEff);
		
		// // linA[face.owner] += (-linAL[i]);
		// // linA[face.neighbour] += (-linAR[i]);
			
		// // convective term
		// linB[face.owner] -= UnF*area;
		// // linB[face.neighbour] += UnF*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // linAL[i] = -tmp1*dPNEff;
			// // linAR[i] = -tmp1*dPNEff;
			
			// linA[face.owner] += (-linAL[i]);
			// linA[face.neighbour] += (-linAR[i]);
			
			// // convective term
			// // linB[face.owner] -= UnF*area;
			// linB[face.neighbour] += UnF*area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // linAL[i] = -tmp1*dPNEff;
			// // linAR[i] = -tmp1*dPNEff;
			
			// linA[face.owner] += (-linAL[i]);
			
			// // convective term
			// // linB[face.owner] -= UnF*area;
			
			
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// // convective term
			// // cout << UnF << endl;
			// // linB[face.owner] -= UnF*area;
			
		// }
		
	// }
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(boundary.type[0] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec(3,0.0);
					// nvec[0] = face.unitNormals[0];
					// nvec[1] = face.unitNormals[1];
					// nvec[2] = face.unitNormals[2];
					
					// vector<double> distanceCells(3,0.0);
					// distanceCells[0] = 2.0*face.distCells[0];
					// distanceCells[1] = 2.0*face.distCells[1];
					// distanceCells[2] = 2.0*face.distCells[2];
				
					// double RhoL = face.varL[controls.fRho];
					// // double RhoR = face.varR[controls.fRho];
				
					// double tmp1 = 1.0/RhoL*controls.timeStep;
					// // double tmp1 = mesh.cells[face.owner].volume/linAD[face.owner]/RhoL;
				
					// // double dPN = 
						// // sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
												// // std::begin(distanceCells), 0.0));
					// // double dPNEff = area/dPN/dPN*
						// // (std::inner_product(std::begin(nvec), std::end(nvec), 
											// // std::begin(distanceCells), 0.0));
													
					// double dPN = sqrt(pow(distanceCells[0],2.0)+pow(distanceCells[1],2.0)+pow(distanceCells[2],2.0));
					// double dPNEff = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
					// // dPNEff = area/dPN *(dPNEff/dPN);
					// dPNEff = area/dPN;
					
					
					// linA[face.owner] += tmp1*dPNEff;
				// }
			// }
			
		// }
		
	// }
	
	
	
	// // if(rank==0) cout<< (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	// // if(rank==0) startTime.pop_back();
	
	
	
	// // if(rank==0) startTime.push_back(clock());
	
	// // linear solver : PETSc library
	// vector<double> resiVar(mesh.cells.size(),0.0);
	

	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverP, controls.toleranceP, 
			// controls.relTolP, controls.preconditionerP,
			// controls.maxIterP);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalP, controls.toleranceFinalP, 
			// controls.relTolFinalP, controls.preconditionerFinalP,
			// controls.maxIterFinalP);
	// }
	

	// // if(rank==0) cout<< (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	// // if(rank==0) startTime.pop_back();
	
	// // resiVar.resize(mesh.cells.size(),0.0);

	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // resiVar[i] = linB[i]/linA[i];
	// // }

	// // if(rank==0) startTime.push_back(clock());
	
	
	
	
	
	
	
	// // // velocities correction
	// // vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
	
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double UnL = face.varL[controls.fU]*face.unitNormals[0] +
			             // // face.varL[controls.fV]*face.unitNormals[1] +
						 // // face.varL[controls.fW]*face.unitNormals[2];
			// // double UnR = face.varR[controls.fU]*face.unitNormals[0] +
			             // // face.varR[controls.fV]*face.unitNormals[1] +
						 // // face.varR[controls.fW]*face.unitNormals[2];
			// // // double PF;
			// // // if(UnL+UnR >=0.0){
				// // // PF = resiVar[face.owner];
			// // // }
			// // // else{
				// // // PF = resiVar[face.neighbour];
			// // // }
			
			// // double PF = 0.5*(resiVar[face.owner]+resiVar[face.neighbour]);
			
			// // for(int j=0; j<3; ++j){
				// // gradResiP[face.owner][j] += 
					// // PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// // gradResiP[face.neighbour][j] -= 
					// // PF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			// // }
		// // }
	// // }

	
	// // // boundary
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // if(boundary.type[0] == "fixedValue"){
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double PF = 0.5*(resiVar[face.owner]);
					// // for(int j=0; j<3; ++j){
						// // gradResiP[face.owner][j] += 
							// // PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// // }
					
				// // }
			// // }
			// // else{
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double PF = resiVar[face.owner];
					// // for(int j=0; j<3; ++j){
						// // gradResiP[face.owner][j] += 
							// // PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// // }
					
				// // }
			// // }
		// // }
	// // }
	
	
	// // // processor faces
	// // vector<double> sendValues;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // sendValues.push_back(resiVar[face.owner]);
		// // }
	// // }
	
	
	// // vector<double> recvValues;
	// // mpi.setProcsFaceDatasDouble(
				// // sendValues, recvValues,
				// // mesh.countsProcFaces, mesh.countsProcFaces, 
				// // mesh.displsProcFaces, mesh.displsProcFaces);

	// // int num=0;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){

			// // double UnL = face.varL[controls.fU]*face.unitNormals[0] +
			             // // face.varL[controls.fV]*face.unitNormals[1] +
						 // // face.varL[controls.fW]*face.unitNormals[2];
			// // double UnR = face.varR[controls.fU]*face.unitNormals[0] +
			             // // face.varR[controls.fV]*face.unitNormals[1] +
						 // // face.varR[controls.fW]*face.unitNormals[2];
			// // // double PF;
			// // // if(UnL+UnR >=0.0){
				// // // PF = resiVar[face.owner];
			// // // }
			// // // else{
				// // // PF = recvValues[num];
			// // // }
			
			// // double PF = 0.5*(resiVar[face.owner]+recvValues[num]);
			
			// // for(int j=0; j<3; ++j){
				// // gradResiP[face.owner][j] += 
					// // PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			// // }
			// // ++num;
		// // }
	// // }
	

	
	// // calc gradient
	// vector<vector<double>> gradResiP;
	
	// gradResiP.clear();
	// gradResiP.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	
	// for(auto& face : mesh.faces){
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double DVar = 
				// resiVar[face.neighbour] - resiVar[face.owner];
			
			// gradResiP[face.owner][0] += face.distCells[0] * DVar;
			// gradResiP[face.owner][1] += face.distCells[1] * DVar;
			// gradResiP[face.owner][2] += face.distCells[2] * DVar;
			
			// gradResiP[face.neighbour][0] += face.distCells[0] * DVar;
			// gradResiP[face.neighbour][1] += face.distCells[1] * DVar;
			// gradResiP[face.neighbour][2] += face.distCells[2] * DVar;
			
		// }
		
	// }
	
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(boundary.type[0] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double DVar = 
						// 0.0 - resiVar[face.owner];
					
					// gradResiP[face.owner][0] += face.distCells[0] * DVar;
					// gradResiP[face.owner][1] += face.distCells[1] * DVar;
					// gradResiP[face.owner][2] += face.distCells[2] * DVar;
					
					
				// }
			// }
		// }
	// }
	
	
	// // processor faces
	// vector<double> sendValues;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// sendValues.push_back(resiVar[face.owner]);
		// }
	// }
	
	// vector<double> recvValues;
	// mpi.setProcsFaceDatasDouble(
				// sendValues, recvValues,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);

	// int num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// double DVar = 
				// recvValues[num] - resiVar[face.owner];
			
			// gradResiP[face.owner][0] += face.distCells[0] * DVar;
			// gradResiP[face.owner][1] += face.distCells[1] * DVar;
			// gradResiP[face.owner][2] += face.distCells[2] * DVar;
					
			// ++num;
		// }
	// }
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// double tmp0 = 
			// cell.coeffLeastSquare[0] * gradResiP[i][0] +
			// cell.coeffLeastSquare[1] * gradResiP[i][1] +
			// cell.coeffLeastSquare[2] * gradResiP[i][2];
			
		// double tmp1 = 
			// cell.coeffLeastSquare[1] * gradResiP[i][0] +
			// cell.coeffLeastSquare[3] * gradResiP[i][1] +
			// cell.coeffLeastSquare[4] * gradResiP[i][2];
			
		// double tmp2 = 
			// cell.coeffLeastSquare[2] * gradResiP[i][0] +
			// cell.coeffLeastSquare[4] * gradResiP[i][1] +
			// cell.coeffLeastSquare[5] * gradResiP[i][2];
			
		// gradResiP[i][0] = tmp0;
		// gradResiP[i][1] = tmp1;
		// gradResiP[i][2] = tmp2;
		
	// }
	

	// // if(rank==0) cout<< (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	// // if(rank==0) startTime.pop_back();
	
	// // if(rank==0) startTime.push_back(clock());
	
	// // update
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		// // if(resiVar[i]<-1000.0) {
			// // cout << resiVar[i] << " " << cell.var[controls.U] << " " << cell.var[controls.V] << " " << cell.var[controls.W]<< endl;
		// // }
		
		
		// // resiVar[i] = linB[i]/linA[i];
		
		// cell.var[controls.P] += controls.prePreURF * resiVar[i];
		
		
		// // cell.var[controls.U] -= controls.preVelURF 
				// // *cell.volume/(linAD[i]+linAOD[i])/cell.var[controls.Rho]*gradResiP[i][0];
		// // cell.var[controls.V] -= controls.preVelURF 
				// // *cell.volume/(linAD[i]+linAOD[i])/cell.var[controls.Rho]*gradResiP[i][1];
		// // cell.var[controls.W] -= controls.preVelURF 
				// // *cell.volume/(linAD[i]+linAOD[i])/cell.var[controls.Rho]*gradResiP[i][2];
				
				
		// // cell.var[controls.U] -= controls.preVelURF 
				// // *cell.volume/(linAD[i])/cell.var[controls.Rho]*gradResiP[i][0];
		// // cell.var[controls.V] -= controls.preVelURF 
				// // *cell.volume/(linAD[i])/cell.var[controls.Rho]*gradResiP[i][1];
		// // cell.var[controls.W] -= controls.preVelURF 
				// // *cell.volume/(linAD[i])/cell.var[controls.Rho]*gradResiP[i][2];
		
		// cell.var[controls.U] -= controls.preVelURF 
				// *controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0];
		// cell.var[controls.V] -= controls.preVelURF 
				// *controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1];
		// cell.var[controls.W] -= controls.preVelURF 
				// *controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2];
				
		// // cell.var[controls.U] -= controls.preVelURF 
				// // *controls.timeStep*gradResiP[i][0];
		// // cell.var[controls.V] -= controls.preVelURF 
				// // *controls.timeStep*gradResiP[i][1];
		// // cell.var[controls.W] -= controls.preVelURF 
				// // *controls.timeStep*gradResiP[i][2];
				
				
		// // double avgRho=0.0;
		// // int nRho = 0;
		// // for(auto& j : cell.faces){
			// // if(mesh.faces[j].getType() == SEMO_Types::INTERNAL_FACE){
				// // avgRho += mesh.cells[mesh.faces[j].neighbour].var[controls.Rho];
				// // ++nRho;
			// // }
		// // }
		// // avgRho += cell.var[controls.Rho];
		// // ++nRho;

		// // cell.var[controls.U] -= controls.preVelURF 
				// // *controls.timeStep/(avgRho/(double)nRho)*gradResiP[i][0];
		// // cell.var[controls.V] -= controls.preVelURF 
				// // *controls.timeStep/(avgRho/(double)nRho)*gradResiP[i][1];
		// // cell.var[controls.W] -= controls.preVelURF 
				// // *controls.timeStep/(avgRho/(double)nRho)*gradResiP[i][2];
				
		
		// residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		// residuals[1] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0],2.0)*cell.volume;
		// residuals[2] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1],2.0)*cell.volume;
		// residuals[3] += pow(controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2],2.0)*cell.volume;
		
		
		
		// // cout << resiVar[i] << " " << gradResiP[i][0] << " " << gradResiP[i][1] << " " << gradResiP[i][2] << endl;
		// // if(rank==2) cout << i << " " << resiVar[i] << " " << gradResiP[i][0] << endl;
		
		
		
		// // linA, linAL, linAR, linB
		
		// // if( cell.var[controls.P] <= controls.minP ) 
			// // cell.var[controls.P] = controls.minP;
		// // if( cell.var[controls.P] >= controls.maxP ) 
			// // cell.var[controls.P] = controls.maxP;
		
		// // if( cell.var[controls.U] <= controls.minU ) 
			// // cell.var[controls.U] = controls.minU;
		// // if( cell.var[controls.U] >= controls.maxU ) 
			// // cell.var[controls.U] = controls.maxU;
		
		// // if( cell.var[controls.V] <= controls.minV ) 
			// // cell.var[controls.V] = controls.minV;
		// // if( cell.var[controls.V] >= controls.maxV ) 
			// // cell.var[controls.V] = controls.maxV;
		
		// // if( cell.var[controls.W] <= controls.minW ) 
			// // cell.var[controls.W] = controls.minW;
		// // if( cell.var[controls.W] >= controls.maxW ) 
			// // cell.var[controls.W] = controls.maxW;
		
		
		
	// }
	

	// // if(rank==0) cout<< (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	// // if(rank==0) startTime.pop_back();

// // MPI_Barrier(MPI_COMM_WORLD);
// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
// }