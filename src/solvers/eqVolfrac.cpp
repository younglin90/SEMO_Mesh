#include "build.h"
#include <cmath>
#include <array>
#include <numeric>

void SEMO_Solvers_Builder::calcVolfracEq(
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
	
	SEMO_Utility_Math math;
	
	
	// // gradient P
	// vector<vector<double>> gradP;
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 1, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	// vector<double> limGradP;
	// math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] *= limGradP[i]; 
		// }
	// }
	
	
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	
	
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
		// SEMO_MPI_Builder mpi;
		
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
	





	// vector<double> cellVolume_recv;
	// if(size>1){
		// vector<double> cellVolume_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// cellVolume_send.push_back(mesh.cells[face.owner].volume);
			// }
		// }
		// mpi.setProcsFaceDatasDouble(
					// cellVolume_send, cellVolume_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// cellVolume_send.clear();
		
	// }


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
		
		mpi.setProcsFaceDatasDouble(
					linAD_send, linAD_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		linAD_send.clear();
		
	}






	vector<double> linA(mesh.cells.size(),0.0);
	vector<double> linAL(mesh.faces.size(),0.0);
	vector<double> linAR(mesh.faces.size(),0.0);
	vector<double> linB(mesh.cells.size(),0.0);

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		
		// time marching -> first order euler
		// linA[i] = cell.volume/controls.timeStep;
		
		// linB[i] = -(
			// 1.0*cell.var[controls.VF[0]] - 1.0*cell.var[controls.oldVF[0]]
			// )*cell.volume/controls.timeStep;
			
		// time marching -> second order upwind euler
		linA[i] = 1.5*cell.volume/controls.timeStep;
		
		linB[i] = -(
			1.5*cell.var[controls.VF[0]] - 2.0*cell.var[controls.oldVF[0]] + 0.5*cell.var[controls.Qm[3]]
			)*cell.volume/controls.timeStep;
		
	}
	
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
		double UR = face.varR[controls.fU];
		double VR = face.varR[controls.fV];
		double WR = face.varR[controls.fW];
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		double RhoL = face.varL[controls.fRho];
		double RhoR = face.varR[controls.fRho];
	
		// double tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
	
		// double dPN = 
			// sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
			                        // std::begin(distanceCells), 0.0));
		// double dPNEff = area/dPN/dPN*
			// (std::inner_product(std::begin(nvec), std::end(nvec), 
			                    // std::begin(distanceCells), 0.0));
		
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
		
		
		double PL = face.varL[controls.fP];
		double PR = face.varR[controls.fP];
		
		
		double UnF = wCL*UnL+wCR*UnR;
		// double UnF = 0.5*UnL+0.5*UnR;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			UnF = face.var[controls.Un];
		}
		
		// // double RhieChowCoeff = 0.8;
		// double RhieChowCoeff = 1.0;
		
		// double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
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
		
		double YiL = face.varL[controls.fVF[0]];
		double YiR = face.varR[controls.fVF[0]];
		// double YiF = 0.5*(YiL+YiR);
		double YiF =YiL*weightL+YiR*weightR;
		
		
		
		
		
		// linAL[i] = -0.5*weightL*UnF*area;
		// linAR[i] = +0.5*weightR*UnF*area;
		linAL[i] = -weightL*UnF*area;
		linAR[i] = +weightR*UnF*area;
		
		linB[face.owner] -= YiF*UnF*area;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			linA[face.owner] += (-linAL[i]);
			linA[face.neighbour] += (-linAR[i]);
			
			// convective term
			linB[face.neighbour] += YiF*UnF*area;
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			linA[face.owner] += (-linAL[i]);
			
		}
		
	}
	
	
	// linear solver : PETSc library
	vector<double> resiVar(mesh.cells.size(),0.0);


	if(controls.iterPBs != controls.iterPBsMax-1){
	
		solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			controls.solverVF[0], controls.toleranceVF[0], 
			controls.relTolVF[0], controls.preconditionerVF[0],
			controls.maxIterVF[0]);
	
		// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverVF[0], controls.toleranceVF[0], 
			// controls.relTolVF[0], controls.preconditionerVF[0],
			// controls.maxIterVF[0]);
		
	}
	else{
	
		solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			controls.maxIterFinalVF[0]);
	
		// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			// controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			// controls.maxIterFinalVF[0]);
	}
	
	
	
	// update
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		// resiVar[i] = linB[i]/linA[i];
		
		cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		
		
		residuals[5] += pow(resiVar[i],2.0)*cell.volume;
		
		
		cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
	}
	


	linA.clear();
	linAL.clear();
	linAR.clear();
	linB.clear();
	resiVar.clear();
	
	
}