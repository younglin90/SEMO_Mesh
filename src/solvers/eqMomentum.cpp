

#include "build.h"
#include <cmath>
#include <array>
#include <numeric>


double SEMO_Solvers_Builder::calcMomentumEqs(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	
	// cout << "AAAAAAAAAA" << endl;

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 1;
    int A_n = B_n * B_n;
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	bool boolSkewnessCorrection2 = true;
	// bool boolSkewnessCorrection2 = false;
	bool boolLimiters = true;
	// bool boolLimiters = false;
	int gradIterMax = 1;
	int gradBoundaryIterMax = 0;
	
	// vector<clock_t> startTime;
	
	SEMO_Utility_Math math;
	
	setCellVarMinMax(mesh, controls.P, controls.maximumP, controls.minimumP);
	setCellVarMinMax(mesh, controls.U, controls.maximumU, controls.minimumU);
	setCellVarMinMax(mesh, controls.V, controls.maximumV, controls.minimumV);
	setCellVarMinMax(mesh, controls.W, controls.maximumW, controls.minimumW);
	
	// mesh.cellsGradientVar[controls.P].clear();
	// mesh.cellsGradientVar[controls.U].clear();
	// mesh.cellsGradientVar[controls.V].clear();
	// mesh.cellsGradientVar[controls.W].clear();
	// mesh.cellsGradientVar[controls.T].clear();
	// mesh.cellsGradientVar[controls.VF[0]].clear();
	
	mesh.cellsGradientVar[controls.P].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.U].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.V].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.W].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1){
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				for(int ii=0; ii<3; ++ii){
					mesh.cellsGradientVar[controls.P][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.U][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.V][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.W][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.VF[0]][face.owner][ii] = 0.0;
				}
			}
		}
	}
	
	
	for(int iter=0; iter<gradIterMax; ++iter){
		// math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		// math.calcGaussGreen(mesh, controls.U, controls.fU, mesh.cellsGradientVar[controls.U]);
		// math.calcGaussGreen(mesh, controls.V, controls.fV, mesh.cellsGradientVar[controls.V]);
		// math.calcGaussGreen(mesh, controls.W, controls.fW, mesh.cellsGradientVar[controls.W]);
		// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], mesh.cellsGradientVar[controls.VF[0]]);
			
		// vector<double> dummyVec;
		// math.calcLeastSquare(mesh, "face", "1st", "cell", 
			// controls.P, controls.fP, dummyVec, mesh.cellsGradientVar[controls.P]);
		// math.calcLeastSquare(mesh, "face", "1st", "cell", 
			// controls.U, controls.fU, dummyVec, mesh.cellsGradientVar[controls.U]);
		// math.calcLeastSquare(mesh, "face", "1st", "cell", 
			// controls.V, controls.fV, dummyVec, mesh.cellsGradientVar[controls.V]);
		// math.calcLeastSquare(mesh, "face", "1st", "cell", 
			// controls.W, controls.fW, dummyVec, mesh.cellsGradientVar[controls.W]);
		// math.calcLeastSquare(mesh, "face", "1st", "cell", 
			// controls.VF[0], controls.fVF[0], dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
			
		vector<double> dummyVec;
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.P, controls.fP, dummyVec, mesh.cellsGradientVar[controls.P]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.U, controls.fU, dummyVec, mesh.cellsGradientVar[controls.U]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.V, controls.fV, dummyVec, mesh.cellsGradientVar[controls.V]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.W, controls.fW, dummyVec, mesh.cellsGradientVar[controls.W]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.VF[0], controls.fVF[0], dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
	}
	
	for(int iter=0; iter<gradBoundaryIterMax; ++iter)
	{
		vector<double> dummyVec;
		math.calcLeastSquareOnlyBoundaryCells(mesh, "cell", controls.P, controls.fP, 
			dummyVec, mesh.cellsGradientVar[controls.P]);
		math.calcLeastSquareOnlyBoundaryCells(mesh, "cell", controls.U, controls.fU, 
			dummyVec, mesh.cellsGradientVar[controls.U]);
		math.calcLeastSquareOnlyBoundaryCells(mesh, "cell", controls.V, controls.fV, 
			dummyVec, mesh.cellsGradientVar[controls.V]);
		math.calcLeastSquareOnlyBoundaryCells(mesh, "cell", controls.W, controls.fW, 
			dummyVec, mesh.cellsGradientVar[controls.W]);
		math.calcLeastSquareOnlyBoundaryCells(mesh, "cell", controls.VF[0], controls.fVF[0], 
			dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
	}
	
	// for(int iter=0; iter<gradBoundaryIterMax; ++iter){
		// // boundary face's nodes
		// for(auto& boundary : mesh.boundary){
			// if(boundary.neighbProcNo == -1){
				// vector<int> target(5,0);
				// target[0] = controls.P;
				// target[1] = controls.U;
				// target[2] = controls.V;
				// target[3] = controls.W;
				// target[4] = controls.VF[0];
				// for(auto& tar : target){
					// if(boundary.type[tar] != "zeroGradient") continue;
					// int str = boundary.startFace;
					// int end = str + boundary.nFaces;
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						// auto& cell = mesh.cells[face.owner];
						// double varF = 
							// mesh.cellsGradientVar[tar][face.owner][0]*(face.x-cell.x) +
							// mesh.cellsGradientVar[tar][face.owner][1]*(face.y-cell.y) +
							// mesh.cellsGradientVar[tar][face.owner][2]*(face.z-cell.z);
						// for(int j=0; j<3; ++j){
							// mesh.cellsGradientVar[tar][face.owner][j] += 
								// varF*face.unitNormals[j]*face.area / cell.volume;
						// }
					// }
				// }
			// }
		// }
	// }
	
	
	// vector<vector<double>> gradP_recv;
	// vector<vector<double>> gradU_recv, gradV_recv, gradW_recv;
	// vector<vector<double>> srcGRAV_recv;
	// vector<vector<double>> srcSFT_recv;
	// vector<double> P_recv, U_recv, V_recv, W_recv, Rho_recv, VF_recv, kappa_recv;
	if(size>1){
		// mpi.sendRecvTemporaryCellData(mesh, controls.P, P_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.U, U_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.V, V_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.W, W_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.Rho, Rho_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.VF[0], VF_recv);
		
		// mpi.sendRecvTemporaryVectorData(mesh, gradP, gradP_recv);
		// mpi.sendRecvTemporaryVectorData(mesh, gradU, gradU_recv);
		// mpi.sendRecvTemporaryVectorData(mesh, gradV, gradV_recv);
		// mpi.sendRecvTemporaryVectorData(mesh, gradW, gradW_recv);
		// mpi.sendRecvTemporaryCellVectorData(mesh, controls.sourceGravity, srcGRAV_recv);
		// mpi.sendRecvTemporaryCellVectorData(mesh, controls.sourceSurfaceTension, srcSFT_recv);
		// mpi.sendRecvTemporaryCellData(mesh, controls.kappa, kappa_recv);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumP, mesh.cellsProcVar[controls.maximumP]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumP, mesh.cellsProcVar[controls.minimumP]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumU, mesh.cellsProcVar[controls.maximumU]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumU, mesh.cellsProcVar[controls.minimumU]);
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumV, mesh.cellsProcVar[controls.maximumV]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumV, mesh.cellsProcVar[controls.minimumV]);
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumW, mesh.cellsProcVar[controls.maximumW]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumW, mesh.cellsProcVar[controls.minimumW]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.oldU, mesh.cellsProcVar[controls.oldU]);
		mpi.sendRecvTemporaryCellData(mesh, controls.oldV, mesh.cellsProcVar[controls.oldV]);
		mpi.sendRecvTemporaryCellData(mesh, controls.oldW, mesh.cellsProcVar[controls.oldW]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.P, mesh.cellsProcVar[controls.P]);
		mpi.sendRecvTemporaryCellData(mesh, controls.U, mesh.cellsProcVar[controls.U]);
		mpi.sendRecvTemporaryCellData(mesh, controls.V, mesh.cellsProcVar[controls.V]);
		mpi.sendRecvTemporaryCellData(mesh, controls.W, mesh.cellsProcVar[controls.W]);
		mpi.sendRecvTemporaryCellData(mesh, controls.Rho, mesh.cellsProcVar[controls.Rho]);
		mpi.sendRecvTemporaryCellData(mesh, controls.VF[0], mesh.cellsProcVar[controls.VF[0]]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[0], mesh.cellsProcVar[controls.sourceGravity[0]]);
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[1], mesh.cellsProcVar[controls.sourceGravity[1]]);
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[2], mesh.cellsProcVar[controls.sourceGravity[2]]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[0], mesh.cellsProcVar[controls.sourceSurfaceTension[0]]);
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[1], mesh.cellsProcVar[controls.sourceSurfaceTension[1]]);
		mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[2], mesh.cellsProcVar[controls.sourceSurfaceTension[2]]);
		
		mpi.sendRecvTemporaryCellData(mesh, controls.kappa, mesh.cellsProcVar[controls.kappa]);
		
		mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.P], mesh.cellsProcGradientVar[controls.P]);
		mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.U], mesh.cellsProcGradientVar[controls.U]);
		mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.V], mesh.cellsProcGradientVar[controls.V]);
		mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.W], mesh.cellsProcGradientVar[controls.W]);
		mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.VF[0]], mesh.cellsProcGradientVar[controls.VF[0]]);
		
	}



	// diagonal terms
    vector<int> A_rows(mesh.cells.size()*A_n, 0);
    vector<int> A_cols(mesh.cells.size()*A_n, 0);
    vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
    vector<vector<double>> B(3,vector<double>(mesh.cells.size(), 0.0));
	

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];

        int ijStart_local = B_n*(i) - 1;
        int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
        int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		double Rho = cell.var[controls.Rho];
		double U = cell.var[controls.U];
		double V = cell.var[controls.V];
		double W = cell.var[controls.W];
		double oldRho = cell.var[controls.oldRho];
		double oldU = cell.var[controls.oldU];
		double oldV = cell.var[controls.oldV];
		double oldW = cell.var[controls.oldW];
		double old2U = cell.var[controls.Qm[0]];
		double old2V = cell.var[controls.Qm[1]];
		double old2W = cell.var[controls.Qm[2]];
		double volume = cell.volume;
		double timeStep = controls.timeStep;
		
		int id = 0;
		int i_glo = 0;
		int j_glo = 0;
		int i_loc = i;
		int str_glo = mesh.startCellGlobal*B_n;
		int step_loc = mesh.cells.size();
		int step_glo = mesh.ncellsTotal;
	
		// bounded 2nd order 시간 차분화
		double beta_phi = 1.0;
		for(int ii=0; ii<3; ++ii){
			double beta_phi_tmp = 1.0;
			double phi_SOU = 0.0;
			if(ii==0){
				phi_SOU = ( oldU - old2U ) /
						(abs(U-old2U)+1.e-200)*( U>old2U ? 1.0 : -1.0 );
			}
			if(ii==1){
				phi_SOU = ( oldV - old2V ) /
						(abs(V-old2V)+1.e-200)*( V>old2V ? 1.0 : -1.0 );
			}
			if(ii==2){
				phi_SOU = ( oldW - old2W ) /
						(abs(W-old2W)+1.e-200)*( W>old2W ? 1.0 : -1.0 );
			}
			if(0.0 < phi_SOU && phi_SOU > 1.0){
				beta_phi_tmp = 0.0;
			}
			else{
				beta_phi_tmp = 1.0;
			}
			beta_phi = min(beta_phi,beta_phi_tmp);
		}
		double coeff_old1 = controls.oldTimeStep/(controls.timeStep+controls.oldTimeStep);
		double coeff_old2 = controls.old2TimeStep/(controls.oldTimeStep+controls.old2TimeStep);
		vector<double> cons_np12(3,0.0);
		vector<double> cons_nm12(3,0.0);
		cons_np12[0] = U + beta_phi * coeff_old1 * (U-oldU);
		cons_nm12[0] = oldU + beta_phi * coeff_old2 * (oldU-old2U);
		cons_np12[1] = V + beta_phi * coeff_old1 * (V-oldV);
		cons_nm12[1] = oldV + beta_phi * coeff_old2 * (oldV-old2V);
		cons_np12[2] = W + beta_phi * coeff_old1 * (W-oldW);
		cons_nm12[2] = oldW + beta_phi * coeff_old2 * (oldW-old2W);

        id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
		
        // momentum 방정식
		
		// A_vals[id] =  Rho*volume/timeStep;
		A_vals[id] = Rho*1.5* volume/timeStep;
		// A_vals[id] = Rho*(1.0+beta_phi*coeff_old1)* volume/timeStep;
        
        // B[0][step_loc*0 + i_loc] = -Rho*(U - oldU)*volume / timeStep;
        // B[1][step_loc*0 + i_loc] = -Rho*(V - oldV)*volume / timeStep;
        // B[2][step_loc*0 + i_loc] = -Rho*(W - oldW)*volume / timeStep;
		
        B[0][step_loc*0 + i_loc] = -Rho*(1.5*U - 2.0*oldU + 0.5*old2U)*volume / timeStep;
        B[1][step_loc*0 + i_loc] = -Rho*(1.5*V - 2.0*oldV + 0.5*old2V)*volume / timeStep;
        B[2][step_loc*0 + i_loc] = -Rho*(1.5*W - 2.0*oldW + 0.5*old2W)*volume / timeStep;
		
        // B[0][step_loc*0 + i_loc] = -Rho*(cons_np12[0] - cons_nm12[0])*volume / timeStep;
        // B[1][step_loc*0 + i_loc] = -Rho*(cons_np12[1] - cons_nm12[1])*volume / timeStep;
        // B[2][step_loc*0 + i_loc] = -Rho*(cons_np12[2] - cons_nm12[2])*volume / timeStep;
		
		
		// gravity force terms
        // B[0][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[0]*volume;
        // B[1][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[1]*volume;
        // B[2][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[2]*volume;
        B[0][step_loc*0 + i_loc] += volume*cell.var[controls.sourceGravity[0]];
        B[1][step_loc*0 + i_loc] += volume*cell.var[controls.sourceGravity[1]];
        B[2][step_loc*0 + i_loc] += volume*cell.var[controls.sourceGravity[2]];

		// surface tension force terms
        B[0][step_loc*0 + i_loc] += volume*cell.var[controls.sourceSurfaceTension[0]];
        B[1][step_loc*0 + i_loc] += volume*cell.var[controls.sourceSurfaceTension[1]];
        B[2][step_loc*0 + i_loc] += volume*cell.var[controls.sourceSurfaceTension[2]];
		
		
		
	}
	
	
	
	
	for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		int ijStartL_local = B_n*(face.owner) - 1;
		int ijStartR_local = B_n*(face.neighbour) - 1;
		
        int ijStartL = B_n*mesh.startCellGlobal + ijStartL_local;
        int ijStartR = B_n*mesh.startCellGlobal + ijStartR_local;
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			ijStartR = B_n*mesh.startProcCellGlobal[mesh.neighbProcNo[ip]] + 
					   B_n*(mesh.procNeighbCellNo[ip]) - 1;
		}

		
		double area = face.area;
		vector<double> nvec(3,0.0);
		nvec[0] = face.unitNormals[0];
		nvec[1] = face.unitNormals[1];
		nvec[2] = face.unitNormals[2];
		
		double wCL = face.wC;
		double wCR = 1.0-wCL;
		// double wCL = 0.5;
		// double wCR = 0.5;
		
		double UL_HO = face.varL[controls.fU];
		double VL_HO = face.varL[controls.fV];
		double WL_HO = face.varL[controls.fW];
		double PL_HO = face.varL[controls.fP];
		double RhoL_HO = face.varL[controls.fRho];
		
		double UR_HO = face.varR[controls.fU];
		double VR_HO = face.varR[controls.fV];
		double WR_HO = face.varR[controls.fW];
		double PR_HO = face.varR[controls.fP];
		double RhoR_HO = face.varR[controls.fRho];
		
		double muL = face.varL[controls.fmu];
		double muR = face.varR[controls.fmu];
		// harmonic 평균
		double muF = 1.0/(wCL/muL + wCR/muR);


		double RhoL = mesh.cells[face.owner].var[controls.Rho];
		double PL = mesh.cells[face.owner].var[controls.P];
		double UL = mesh.cells[face.owner].var[controls.U];
		double VL = mesh.cells[face.owner].var[controls.V];
		double WL = mesh.cells[face.owner].var[controls.W];
		double VFL = mesh.cells[face.owner].var[controls.VF[0]];
		double kappaL = mesh.cells[face.owner].var[controls.kappa];
		double UL_old = mesh.cells[face.owner].var[controls.oldU];
		double VL_old = mesh.cells[face.owner].var[controls.oldV];
		double WL_old = mesh.cells[face.owner].var[controls.oldW];
		double maximumPL = mesh.cells[face.owner].var[controls.maximumP];
		double minimumPL = mesh.cells[face.owner].var[controls.minimumP];
		double maximumUL = mesh.cells[face.owner].var[controls.maximumU];
		double minimumUL = mesh.cells[face.owner].var[controls.minimumU];
		double maximumVL = mesh.cells[face.owner].var[controls.maximumV];
		double minimumVL = mesh.cells[face.owner].var[controls.minimumV];
		double maximumWL = mesh.cells[face.owner].var[controls.maximumW];
		double minimumWL = mesh.cells[face.owner].var[controls.minimumW];
		vector<double> gradPL(3,0.0);
		vector<double> gradUL(3,0.0);
		vector<double> gradVL(3,0.0);
		vector<double> gradWL(3,0.0);
		vector<double> gradVFL(3,0.0);
		vector<double> srcGRV_L(3,0.0);
		vector<double> srcSFT_L(3,0.0);
		
		double RhoR = 0.0;
		double PR = 0.0;
		double UR = 0.0;
		double VR = 0.0;
		double WR = 0.0;
		double VFR = 0.0;
		double kappaR = 0.0;
		double UR_old = 0.0;
		double VR_old = 0.0;
		double WR_old = 0.0;
		double maximumPR = 0.0;
		double minimumPR = 0.0;
		double maximumUR = 0.0;
		double minimumUR = 0.0;
		double maximumVR = 0.0;
		double minimumVR = 0.0;
		double maximumWR = 0.0;
		double minimumWR = 0.0;
		vector<double> gradPR(3,0.0);
		vector<double> gradUR(3,0.0);
		vector<double> gradVR(3,0.0);
		vector<double> gradWR(3,0.0);
		vector<double> gradVFR(3,0.0);
		vector<double> srcGRV_R(3,0.0);
		vector<double> srcSFT_R(3,0.0);
		
		for(int ii=0; ii<3; ++ii){
			gradPL[ii] = mesh.cellsGradientVar[controls.P][face.owner][ii];
			gradUL[ii] = mesh.cellsGradientVar[controls.U][face.owner][ii];
			gradVL[ii] = mesh.cellsGradientVar[controls.V][face.owner][ii];
			gradWL[ii] = mesh.cellsGradientVar[controls.W][face.owner][ii];
			gradVFL[ii] = mesh.cellsGradientVar[controls.VF[0]][face.owner][ii];
			srcGRV_L[ii] = mesh.cells[face.owner].var[controls.sourceGravity[ii]];
			srcSFT_L[ii] = mesh.cells[face.owner].var[controls.sourceSurfaceTension[ii]];
		}
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			RhoR = mesh.cells[face.neighbour].var[controls.Rho];
			PR = mesh.cells[face.neighbour].var[controls.P];
			UR = mesh.cells[face.neighbour].var[controls.U];
			VR = mesh.cells[face.neighbour].var[controls.V];
			WR = mesh.cells[face.neighbour].var[controls.W];
			VFR = mesh.cells[face.neighbour].var[controls.VF[0]];
			kappaR = mesh.cells[face.neighbour].var[controls.kappa];
			UR_old = mesh.cells[face.neighbour].var[controls.oldU];
			VR_old = mesh.cells[face.neighbour].var[controls.oldV];
			WR_old = mesh.cells[face.neighbour].var[controls.oldW];
			maximumPR = mesh.cells[face.neighbour].var[controls.maximumP];
			minimumPR = mesh.cells[face.neighbour].var[controls.minimumP];
			maximumUR = mesh.cells[face.neighbour].var[controls.maximumU];
			minimumUR = mesh.cells[face.neighbour].var[controls.minimumU];
			maximumVR = mesh.cells[face.neighbour].var[controls.maximumV];
			minimumVR = mesh.cells[face.neighbour].var[controls.minimumV];
			maximumWR = mesh.cells[face.neighbour].var[controls.maximumW];
			minimumWR = mesh.cells[face.neighbour].var[controls.minimumW];
			for(int ii=0; ii<3; ++ii){
				gradPR[ii] = mesh.cellsGradientVar[controls.P][face.neighbour][ii];
				gradUR[ii] = mesh.cellsGradientVar[controls.U][face.neighbour][ii];
				gradVR[ii] = mesh.cellsGradientVar[controls.V][face.neighbour][ii];
				gradWR[ii] = mesh.cellsGradientVar[controls.W][face.neighbour][ii];
				gradVFR[ii] = mesh.cellsGradientVar[controls.VF[0]][face.neighbour][ii];
				srcGRV_R[ii] = mesh.cells[face.neighbour].var[controls.sourceGravity[ii]];
				srcSFT_R[ii] = mesh.cells[face.neighbour].var[controls.sourceSurfaceTension[ii]];
			}
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			RhoR = mesh.cellsProcVar[controls.Rho][ip];
			PR = mesh.cellsProcVar[controls.P][ip];
			UR = mesh.cellsProcVar[controls.U][ip];
			VR = mesh.cellsProcVar[controls.V][ip];
			WR = mesh.cellsProcVar[controls.W][ip];
			VFR = mesh.cellsProcVar[controls.VF[0]][ip];
			kappaR = mesh.cellsProcVar[controls.kappa][ip];
			UR_old = mesh.cellsProcVar[controls.oldU][ip];
			VR_old = mesh.cellsProcVar[controls.oldV][ip];
			WR_old = mesh.cellsProcVar[controls.oldW][ip];
			maximumPR = mesh.cellsProcVar[controls.maximumP][ip];
			minimumPR = mesh.cellsProcVar[controls.minimumP][ip];
			maximumUR = mesh.cellsProcVar[controls.maximumU][ip];
			minimumUR = mesh.cellsProcVar[controls.minimumU][ip];
			maximumVR = mesh.cellsProcVar[controls.maximumV][ip];
			minimumVR = mesh.cellsProcVar[controls.minimumV][ip];
			maximumWR = mesh.cellsProcVar[controls.maximumW][ip];
			minimumWR = mesh.cellsProcVar[controls.minimumW][ip];
			for(int ii=0; ii<3; ++ii){
				gradPR[ii] = mesh.cellsProcGradientVar[controls.P][ip][ii];
				gradUR[ii] = mesh.cellsProcGradientVar[controls.U][ip][ii];
				gradVR[ii] = mesh.cellsProcGradientVar[controls.V][ip][ii];
				gradWR[ii] = mesh.cellsProcGradientVar[controls.W][ip][ii];
				gradVFR[ii] = mesh.cellsProcGradientVar[controls.VF[0]][ip][ii];
				srcGRV_R[ii] = mesh.cellsProcVar[controls.sourceGravity[ii]][ip];
				srcSFT_R[ii] = mesh.cellsProcVar[controls.sourceSurfaceTension[ii]][ip];
			}
			
		}
		
		double Rho_star = 1.0 / (wCR/RhoL + wCR/RhoR);
		double tmp1 = controls.timeStep / Rho_star;
		
		double dPN = face.magPN;
		double alpha = face.alphaF;
		
		double nonOrtholimiter = 1.0;
		// double cosAlpha = 1.0/alpha;
		// if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			// nonOrtholimiter = 0.5;
		// }
		// else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			// nonOrtholimiter = 0.333;
		// }
		// else if( cosAlpha < 0.342 ){
			// nonOrtholimiter = 0.0;
		// }
		
			// cout << "AAAAAA" << endl;
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		double UnF = wCL*UnL+wCR*UnR;
		
		if(boolSkewnessCorrection){
			// calcInterpolVelSkewness(UnF, gradUL, gradUR, gradVL, gradVR, gradWL, gradWR,
				// wCL, wCR, 
				// nvec,
				// face.vecSkewness);
			double UL_skew = UL;
			double VL_skew = VL;
			double WL_skew = WL;
			double UR_skew = UR;
			double VR_skew = VR;
			double WR_skew = WR;
			for(int ii=0; ii<3; ++ii){
				UL_skew += gradUL[ii]*face.vecSkewness[ii];
				UR_skew += gradUR[ii]*face.vecSkewness[ii];
				VL_skew += gradVL[ii]*face.vecSkewness[ii];
				VR_skew += gradVR[ii]*face.vecSkewness[ii];
				WL_skew += gradWL[ii]*face.vecSkewness[ii];
				WR_skew += gradWR[ii]*face.vecSkewness[ii];
			}
			if(boolLimiters){
				UL_skew = max(minimumUL,min(maximumUL,UL_skew));
				UR_skew = max(minimumUR,min(maximumUR,UR_skew));
				VL_skew = max(minimumVL,min(maximumVL,VL_skew));
				VR_skew = max(minimumVR,min(maximumVR,VR_skew));
				WL_skew = max(minimumWL,min(maximumWL,WL_skew));
				WR_skew = max(minimumWR,min(maximumWR,WR_skew));
			}
			
			UnL = UL_skew*nvec[0] + VL_skew*nvec[1] + WL_skew*nvec[2];
			UnR = UR_skew*nvec[0] + VR_skew*nvec[1] + WR_skew*nvec[2];
			UnF = wCL*UnL+wCR*UnR;
			
			// double UL_skew = UL;
			// double VL_skew = VL;
			// double WL_skew = WL;
			// double UR_skew = UR;
			// double VR_skew = VR;
			// double WR_skew = WR;
			// for(int ii=0; ii<3; ++ii){
				// UL_skew += gradUL[ii]*face.vecPF[ii];
				// VL_skew += gradVL[ii]*face.vecPF[ii];
				// WL_skew += gradWL[ii]*face.vecPF[ii];
				// UR_skew += gradUR[ii]*face.vecNF[ii];
				// VR_skew += gradVR[ii]*face.vecNF[ii];
				// WR_skew += gradWR[ii]*face.vecNF[ii];
			// }
			// if(boolLimiters){
				// UL_skew = max(minimumUL,min(maximumUL,UL_skew));
				// UR_skew = max(minimumUR,min(maximumUR,UR_skew));
				// VL_skew = max(minimumVL,min(maximumVL,VL_skew));
				// VR_skew = max(minimumVR,min(maximumVR,VR_skew));
				// WL_skew = max(minimumWL,min(maximumWL,WL_skew));
				// WR_skew = max(minimumWR,min(maximumWR,WR_skew));
			// }
			// UnL = UL_skew*nvec[0] + VL_skew*nvec[1] + WL_skew*nvec[2];
			// UnR = UR_skew*nvec[0] + VR_skew*nvec[1] + WR_skew*nvec[2];
			// UnF = 0.5*UnL+0.5*UnR;
		}
		calcInterpolVelPressure(UnF, 
			gradPL, gradPR, 
			PL, PR, RhoL, RhoR, wCL, wCR, 
			controls.timeStep,
			alpha, dPN,
			nvec,
			face.unitNomalsPN);
		calcInterpolVelGravity(UnF, srcGRV_L, srcGRV_R,
			RhoL, RhoR, wCL, wCR, 
			controls.timeStep,
			controls.gravityAcceleration,
			alpha,
			nvec,
			face.unitNomalsPN);
		calcInterpolVelSurfTens(UnF, srcSFT_L, srcSFT_R,
			RhoL, RhoR, wCL, wCR, 
			controls.timeStep,
			species[0].sigma, kappaL, kappaR,
			VFL, VFR,
			alpha, dPN,
			nvec,
			face.unitNomalsPN);
		// calcInterpolVelUnsteady(UnF, 
			// wCL, wCR, 
			// face.var[controls.Un],
			// UL_old, UR_old,
			// VL_old, VR_old,
			// WL_old, WR_old,
			// nvec);
			
			
			
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = 1.0 - weightL;
		
		// double RhoF = weightL*RhoL + weightR*RhoR;
		double PF = wCL*PL + wCR*PR;
		double UF_HO = weightL*UL_HO + weightR*UR_HO;
		double VF_HO = weightL*VL_HO + weightR*VR_HO;
		double WF_HO = weightL*WL_HO + weightR*WR_HO;
		
		
		// skewness correction
		if(boolSkewnessCorrection2){
			PF = wCL*PL + wCR*PR;
			for(int ii=0; ii<3; ++ii){
				PF += wCL * gradPL[ii]*face.vecSkewness[ii];
				PF += wCR * gradPR[ii]*face.vecSkewness[ii];
				// PF += 0.5 * gradPL[ii]*face.vecPF[ii];
				// PF += 0.5 * gradPR[ii]*face.vecNF[ii];
				
				// UF_HO += weightL * gradUL[ii]*face.vecSkewness[ii];
				// UF_HO += weightR * gradUR[ii]*face.vecSkewness[ii];
				
				// VF_HO += weightL * gradVL[ii]*face.vecSkewness[ii];
				// VF_HO += weightR * gradVR[ii]*face.vecSkewness[ii];
				
				// WF_HO += weightL * gradWL[ii]*face.vecSkewness[ii];
				// WF_HO += weightR * gradWR[ii]*face.vecSkewness[ii];
				
				// UF_HO += weightL * gradUL[ii]*face.vecPF[ii];
				// UF_HO += weightR * gradUR[ii]*face.vecNF[ii];
				
				// VF_HO += weightL * gradVL[ii]*face.vecPF[ii];
				// VF_HO += weightR * gradVR[ii]*face.vecNF[ii];
				
				// WF_HO += weightL * gradWL[ii]*face.vecPF[ii];
				// WF_HO += weightR * gradWR[ii]*face.vecNF[ii];
			}
			if(boolLimiters){
				PF = max(minimumPL,min(maximumPL,PF));
				PF = max(minimumPR,min(maximumPR,PF));
				// UF_HO = max(minimumUL,min(maximumUL,UF_HO));
				// UF_HO = max(minimumUR,min(maximumUR,UF_HO));
				// VF_HO = max(minimumVL,min(maximumVL,VF_HO));
				// VF_HO = max(minimumVR,min(maximumVR,VF_HO));
				// WF_HO = max(minimumWL,min(maximumWL,WF_HO));
				// WF_HO = max(minimumWR,min(maximumWR,WF_HO));
			}
			

			// PF = 0.5*PL + 0.5*PR;
			// for(int ii=0; ii<3; ++ii){
				// PF += 0.5 * gradPL[ii]*face.vecPF[ii];
				// PF += 0.5 * gradPR[ii]*face.vecNF[ii];
			// }
			// if(boolLimiters){
				// PF = max(minimumPL,min(maximumPL,PF));
				// PF = max(minimumPR,min(maximumPR,PF));
			// }
			
		}


		
		int str_glo_L = mesh.startCellGlobal*B_n;
		int str_glo_R = mesh.startCellGlobal*B_n;
		int step_loc_L = mesh.cells.size();
		int step_loc_R = mesh.cells.size();
		int i_loc_L = face.owner;
		int i_loc_R = face.neighbour;
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			str_glo_R = mesh.startProcCellGlobal[mesh.neighbProcNo[ip]]*B_n;
			step_loc_R = mesh.startProcCellGlobal[mesh.neighbProcNo[ip]+1] - 
			             mesh.startProcCellGlobal[mesh.neighbProcNo[ip]];
			i_loc_R = mesh.procNeighbCellNo[ip];
		}
		
		int id = -1;
		int i_glo = -1;
		int j_glo = -1;

        //------------------------
		id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
		A_vals[id] += ( weightL * RhoL * UnF * area + alpha * muF/dPN*area );
		A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		A_vals.push_back(( weightR * RhoL * UnF * area - alpha * muF/dPN*area ));
			
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( weightR * RhoR * UnF * area - alpha * muF/dPN*area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( weightL * RhoR * UnF * area + alpha * muF/dPN*area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
		vector<double> conv_flux_L(3,0.0);
		vector<double> conv_flux_R(3,0.0);
		vector<double> diff_flux(3,0.0);
		vector<double> diff_flux_L(3,0.0);
		vector<double> diff_flux_R(3,0.0);
		
		// wCL = (float)(round(wCL * 10000000) / 10000000);
		// wCR = 1.0-wCL;
		
		// PF = wCL*PL + wCR*PR;
		// PF = wCL*PL + wCR*PR;
		// for(int ii=0; ii<3; ++ii){
			// PF += wCL * gradPL[ii]*face.vecPF[ii];
			// PF += wCR * gradPR[ii]*face.vecNF[ii];
		// }
		// PF = max(minimumPL,min(maximumPL,PF));
		// PF = max(minimumPR,min(maximumPR,PF));
			// cout << "NONONO " << PF << " " << PL << " " << PR << " " << wCL << endl;
		// PF -= 1.e5;
		// if(PF != 0.0){
			// cout.precision(20);
			// cout << "NONONO " << PF << " " << PL << " " << PR << " " << wCL << endl;
		// }
		
		// convection
		conv_flux_L[0] = (RhoL * UF_HO * UnF + PF * nvec[0]) * area;
		conv_flux_L[1] = (RhoL * VF_HO * UnF + PF * nvec[1]) * area;
		conv_flux_L[2] = (RhoL * WF_HO * UnF + PF * nvec[2]) * area;
		
		conv_flux_R[0] = (RhoR * UF_HO * UnF + PF * nvec[0]) * area;
		conv_flux_R[1] = (RhoR * VF_HO * UnF + PF * nvec[1]) * area;
		conv_flux_R[2] = (RhoR * WF_HO * UnF + PF * nvec[2]) * area;
		
		// viscous
		diff_flux[0] = alpha * muF*(UR-UL)/dPN*area;
		diff_flux[1] = alpha * muF*(VR-VL)/dPN*area;
		diff_flux[2] = alpha * muF*(WR-WL)/dPN*area;
		// viscous, non-orthogonal
		for(int ii=0; ii<3; ++ii){
			diff_flux[0] += ( muF*(wCL*gradUL[ii] + wCR*gradUR[ii])*(nvec[ii] - alpha*face.unitNomalsPN[ii])*area ); 
			diff_flux[1] += ( muF*(wCL*gradVL[ii] + wCR*gradVR[ii])*(nvec[ii] - alpha*face.unitNomalsPN[ii])*area ); 
			diff_flux[2] += ( muF*(wCL*gradWL[ii] + wCR*gradWR[ii])*(nvec[ii] - alpha*face.unitNomalsPN[ii])*area ); 
		}
		
		// 추가 텀2
		diff_flux_L[0] = ( (gradUL[0]*nvec[0]+gradUL[1]*nvec[1]+gradUL[2]*nvec[2])*muF*area );
		diff_flux_L[1] = ( (gradVL[0]*nvec[0]+gradVL[1]*nvec[1]+gradVL[2]*nvec[2])*muF*area );
		diff_flux_L[2] = ( (gradWL[0]*nvec[0]+gradWL[1]*nvec[1]+gradWL[2]*nvec[2])*muF*area );
		
		diff_flux_R[0] = ( (gradUR[0]*nvec[0]+gradUR[1]*nvec[1]+gradUR[2]*nvec[2])*muF*area );
		diff_flux_R[1] = ( (gradVR[0]*nvec[0]+gradVR[1]*nvec[1]+gradVR[2]*nvec[2])*muF*area );
		diff_flux_R[2] = ( (gradWR[0]*nvec[0]+gradWR[1]*nvec[1]+gradWR[2]*nvec[2])*muF*area );
		
		
		
        // ----------------------------
		for(int ii=0; ii<3; ++ii){
			B[ii][step_loc_L*0 + i_loc_L] -= ( conv_flux_L[ii] - diff_flux[ii] - diff_flux_L[ii] );
		}
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(int ii=0; ii<3; ++ii){
				B[ii][step_loc_R*0 + i_loc_R] += ( conv_flux_R[ii] - diff_flux[ii] - diff_flux_R[ii] );
			}
		}
		
		
		
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
		
		
	}
	
	
	
	
	
	
	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
				
				double area = face.area;
				vector<double> nvec(3,0.0);
				nvec[0] = face.unitNormals[0];
				nvec[1] = face.unitNormals[1];
				nvec[2] = face.unitNormals[2];
		
				double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
							 face.varL[controls.fV]*face.unitNormals[1] + 
							 face.varL[controls.fW]*face.unitNormals[2];
				double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
							 face.varR[controls.fV]*face.unitNormals[1] + 
							 face.varR[controls.fW]*face.unitNormals[2];
				
				double UnF = 0.5*UnL+0.5*UnR;
				
				double RhoL = face.varL[controls.fRho];
				double tmp1 = 1.0/RhoL*controls.timeStep;
				
				double UF = 0.5*face.varL[controls.fU] + 0.5*face.varR[controls.fU];
				double VF = 0.5*face.varL[controls.fV] + 0.5*face.varR[controls.fV];
				double WF = 0.5*face.varL[controls.fW] + 0.5*face.varR[controls.fW];
				double PF = 0.5*face.varL[controls.fP] + 0.5*face.varR[controls.fP];
				
				// velocity coefficient
				double coeffUVW = 0.0;
				double coeffUVW_diff = 0.0;
				if( boundary.type[controls.U] == "fixedValue" ){
					coeffUVW = 0.0;
					coeffUVW_diff = 1.0;
				}
				else if( boundary.type[controls.U] == "zeroGradient" ){
					coeffUVW = 1.0;
					coeffUVW_diff = 0.0;
				}
				else if( boundary.type[controls.U] == "slip" ){
					coeffUVW = 0.0;
					coeffUVW_diff = 0.0;
				}
				else if( boundary.type[controls.U] == "noSlip" ){
					coeffUVW = 0.0;
					coeffUVW_diff = 1.0;
				}
				else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
					coeffUVW = 0.0;
					coeffUVW_diff = 1.0;
				}
				else if( boundary.type[controls.U] == "inletOutlet" ){
					double ownNorVel =  
						mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					coeffUVW = (ownNorVel > 0.0) ? 1.0 : 0.0;
					coeffUVW_diff = (ownNorVel > 0.0) ? 0.0 : 1.0;
				}
				else {
					cerr << "| #Error : not defined B.C., var = " << controls.U << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				

				int str_glo = mesh.startCellGlobal*B_n;
				int step_loc = mesh.cells.size();
				int step_glo = mesh.ncellsTotal;
				int i_loc_L = face.owner;
				int id=-1;
				
				
				// properties, viscous term
				double muF = 0.5*face.varL[controls.fmu] + 0.5*face.varR[controls.fmu];
				
				double orgUL = mesh.cells[face.owner].var[controls.U];
				double orgUR = face.varR[controls.fU];
				double orgVL = mesh.cells[face.owner].var[controls.V];
				double orgVR = face.varR[controls.fV];
				double orgWL = mesh.cells[face.owner].var[controls.W];
				double orgWR = face.varR[controls.fW];
				
				
				// // interpolation
				// {
					// double UR = cell.var[controls.U];
					// double VR = cell.var[controls.V];
					// double WR = cell.var[controls.W];
					// double ownNorVel = UR*face.unitNormals[0] +
									   // VR*face.unitNormals[1] +
									   // WR*face.unitNormals[2];
					// if(  boundary.type[controls.U] == "zeroGradient" ||
						// (boundary.type[controls.U] == "inletOutlet" && ownNorVel > 0.0) ){
					
						// UR += mesh.cellsGradientVar[controls.U][face.owner][0]*(face.x-cell.x);
						// UR += mesh.cellsGradientVar[controls.U][face.owner][1]*(face.y-cell.y);
						// UR += mesh.cellsGradientVar[controls.U][face.owner][2]*(face.z-cell.z);
						
						// VR += mesh.cellsGradientVar[controls.V][face.owner][0]*(face.x-cell.x);
						// VR += mesh.cellsGradientVar[controls.V][face.owner][1]*(face.y-cell.y);
						// VR += mesh.cellsGradientVar[controls.V][face.owner][2]*(face.z-cell.z);
						
						// WR += mesh.cellsGradientVar[controls.W][face.owner][0]*(face.x-cell.x);
						// WR += mesh.cellsGradientVar[controls.W][face.owner][1]*(face.y-cell.y);
						// WR += mesh.cellsGradientVar[controls.W][face.owner][2]*(face.z-cell.z);
						
						// UF = UR;
						// VF = VR;
						// WF = WR;
						
						// // UnF = UR*nvec[0]+VR*nvec[1]+WR*nvec[2];
					// }
				// }
				
				
				// face->cell interpolation
				double dPN = 0.5 * face.magPN;
				double alpha = face.alphaF;
				
				
				if( boundary.type[controls.P] == "fixedValue" ){
					double orgPL = mesh.cells[face.owner].var[controls.P];
					double orgPR = boundary.var[controls.P];
					
					UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
					for(int ii=0; ii<3; ++ii){
						UnF += alpha * tmp1 * 
							mesh.cellsGradientVar[controls.P][face.owner][ii]*face.unitNomalsPN[ii];
					}
				}
				
				// {
					
					// if( boundary.type[controls.P] == "fixedValue" ){
						// double orgPL = mesh.cells[face.owner].var[controls.P];
						// double orgPR = boundary.var[controls.P];
						
						// UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
						// for(int ii=0; ii<3; ++ii){
							// UnF += alpha * tmp1 * 
								// mesh.cellsGradientVar[controls.P][face.owner][ii]*face.unitNomalsPN[ii];
						// }
						
						// // PF = boundary.var[controls.P];
					// }
					
					// if( boundary.type[controls.U] == "inletOutlet" ){
						// double ownNorVel =  
							// mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
							// mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
							// mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
						// if(ownNorVel<0.0){
							// double magVel = boundary.var[controls.U];
							// UnF = magVel;
							// UF = magVel*nvec[0];
							// VF = magVel*nvec[1];
							// WF = magVel*nvec[2];	
						// }
						
					// }
					
					
				// }
				
					
				// viscous terms
				double dUdnF = 0.0;
				double dVdnF = 0.0;
				double dWdnF = 0.0;

				dUdnF += alpha * (orgUR-orgUL)/dPN;
				dVdnF += alpha * (orgVR-orgVL)/dPN;
				dWdnF += alpha * (orgWR-orgWL)/dPN;
				for(int ii=0; ii<3; ++ii){
					// face gradients, non-orthogonal
					double tvec = nvec[ii] - alpha * face.unitNomalsPN[ii];
					dUdnF += mesh.cellsGradientVar[controls.U][face.owner][ii] * tvec;
					dVdnF += mesh.cellsGradientVar[controls.V][face.owner][ii] * tvec;
					dWdnF += mesh.cellsGradientVar[controls.W][face.owner][ii] * tvec;
				}
				
				// 추가 텀2
				for(int ii=0; ii<3; ++ii){
					dUdnF += ( mesh.cellsGradientVar[controls.U][face.owner][ii]*nvec[ii] );
					dVdnF += ( mesh.cellsGradientVar[controls.V][face.owner][ii]*nvec[ii] );
					dWdnF += ( mesh.cellsGradientVar[controls.W][face.owner][ii]*nvec[ii] );
				}
				
				
				id = step_loc*(B_n*0+0) + i_loc_L;
				A_vals[id] += coeffUVW * ( RhoL * UnF * area );
				A_vals[id] += coeffUVW_diff * ( muF * alpha/dPN * area );
				
				
				// ----------------------------
				B[0][step_loc*0 + i_loc_L] -= ( RhoL * UF * UnF * area + 
												PF * nvec[0] * area - 
												muF*dUdnF*area );
				B[1][step_loc*0 + i_loc_L] -= ( RhoL * VF * UnF * area + 
												PF * nvec[1] * area - 
												muF*dVdnF*area );
				B[2][step_loc*0 + i_loc_L] -= ( RhoL * WF * UnF * area + 
												PF * nvec[2] * area - 
												muF*dWdnF*area );
			}
		}
	}
	
	
	
	
	
	// linear solver : PETSc library
	vector<vector<double>> resiVar(3,vector<double>(B_n*mesh.cells.size(),0.0));
	
	solveAMGCL("momentum", mesh, B_n, A_rows, A_cols, A_vals, B[0], B[1], B[2], resiVar[0], resiVar[1], resiVar[2]);

	double relaxP = controls.prePreURF;
	double relaxUVW = controls.momVelURF;
	
	double maxP = -1.e9;
	vector<double> maxP_xyz(3,0.0);
	double maxUn = -1.e9;
	vector<double> maxUn_xyz(3,0.0);
	
	double normDelPrim = 0.0;
	double normPrim = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
        int ijStart = B_n*(i) - 1;
        int Astart = A_n*(i) - 1;
        int id = Astart;
		
		int step_loc = mesh.cells.size();
		int i_loc = i;
		
		
		cell.var[controls.U] += controls.momVelURF * resiVar[0][step_loc*0 + i_loc];
		cell.var[controls.V] += controls.momVelURF * resiVar[1][step_loc*0 + i_loc];
		cell.var[controls.W] += controls.momVelURF * resiVar[2][step_loc*0 + i_loc];
		
		// cell.var[controls.U] = cell.var[controls.sourceSurfaceTension[0]];
		// cell.var[controls.V] = cell.var[controls.sourceSurfaceTension[1]];
		// cell.var[controls.W] = cell.var[controls.sourceSurfaceTension[2]];
		
		if(maxP<cell.var[controls.P]){
			maxP = cell.var[controls.P];
			maxP_xyz[0] = cell.x; maxP_xyz[1] = cell.y; maxP_xyz[2] = cell.z;
		}
		double Un = sqrt(
		cell.var[controls.U]*cell.var[controls.U]+
		cell.var[controls.V]*cell.var[controls.V]+
		cell.var[controls.W]*cell.var[controls.W]);
		if(maxUn<Un){
			maxUn = Un;
			maxUn_xyz[0] = cell.x; maxUn_xyz[1] = cell.y; maxUn_xyz[2] = cell.z;
		}
		
		
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
		
		
		normDelPrim += pow(resiVar[0][step_loc*0 + i_loc],2.0);
		normDelPrim += pow(resiVar[1][step_loc*0 + i_loc],2.0);
		normDelPrim += pow(resiVar[2][step_loc*0 + i_loc],2.0);
		
		normPrim += pow(cell.var[controls.U],2.0);
		normPrim += pow(cell.var[controls.V],2.0);
		normPrim += pow(cell.var[controls.W],2.0);
	}
	

	
	
	
	vector<double> maxUn_glo(size,0.0);
	vector<double> maxUn_x_glo(size,0.0);
	vector<double> maxUn_y_glo(size,0.0);
	vector<double> maxUn_z_glo(size,0.0);
	
	MPI_Allgather(&maxUn, 1, MPI_DOUBLE, maxUn_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[0], 1, MPI_DOUBLE, maxUn_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[1], 1, MPI_DOUBLE, maxUn_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[2], 1, MPI_DOUBLE, maxUn_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	
	if(rank==0){
		maxUn = maxUn_glo[0];
		maxUn_xyz[0] = maxUn_x_glo[0]; maxUn_xyz[1] = maxUn_y_glo[0]; maxUn_xyz[2] = maxUn_z_glo[0];
		for(int i=1; i<size; ++i){
			if(maxUn<maxUn_glo[i]){
				maxUn = maxUn_glo[i];
				maxUn_xyz[0] = maxUn_x_glo[i]; maxUn_xyz[1] = maxUn_y_glo[i]; maxUn_xyz[2] = maxUn_z_glo[i];
			}
		}
		
		cout << "-----------------" << endl;
		cout << " maximum Un = " << maxUn << ", xyz = " << maxUn_xyz[0] << ", " << maxUn_xyz[1] << ", " << maxUn_xyz[2] << endl;
		cout << "-----------------" << endl;
	}
	

	
	return sqrt(normDelPrim)/sqrt(normPrim);
	
	
}





















double SEMO_Solvers_Builder::calcMomentumEqs(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	vector<double>& linAD,
	vector<double>& linAOD,
	vector<double>& linAFL,
	vector<double>& linAFR,
	vector<double>& residuals){
		
		
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
		
	// SEMO_MPI_Builder mpi;
	
	
	// int proc_num=0;

	// SEMO_Utility_Math math;
	
	// // gradient P
	// vector<vector<double>> gradP;
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
	// // vector<double> limGradP;
	// // math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // for(int j=0; j<3; ++j){
			// // gradP[i][j] *= limGradP[i]; 
		// // }
	// // }
	
	
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
		
		// mpi.setProcsFaceDatas(
					// gradPx_send, gradPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPy_send, gradPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPz_send, gradPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	// }
	

	
	// // gradient U, V, W
	// vector<vector<double>> gradU;
	// vector<vector<double>> gradV;
	// vector<vector<double>> gradW;
	// // math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// // math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// // math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// // math.calcMGG(mesh, controls.U, controls.fU, 10, 1.e-8, gradU);
	// // math.calcMGG(mesh, controls.V, controls.fV, 10, 1.e-8, gradV);
	// // math.calcMGG(mesh, controls.W, controls.fW, 10, 1.e-8, gradW);
	// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	// vector<double> gradUx_recv, gradUy_recv, gradUz_recv;
	// vector<double> gradVx_recv, gradVy_recv, gradVz_recv;
	// vector<double> gradWx_recv, gradWy_recv, gradWz_recv;
	// if(size>1){
		// // processor faces
		// vector<double> gradUx_send, gradUy_send, gradUz_send;
		// vector<double> gradVx_send, gradVy_send, gradVz_send;
		// vector<double> gradWx_send, gradWy_send, gradWz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradUx_send.push_back(gradU[face.owner][0]);
				// gradUy_send.push_back(gradU[face.owner][1]);
				// gradUz_send.push_back(gradU[face.owner][2]);
				
				// gradVx_send.push_back(gradV[face.owner][0]);
				// gradVy_send.push_back(gradV[face.owner][1]);
				// gradVz_send.push_back(gradV[face.owner][2]);
				
				// gradWx_send.push_back(gradW[face.owner][0]);
				// gradWy_send.push_back(gradW[face.owner][1]);
				// gradWz_send.push_back(gradW[face.owner][2]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// gradUx_send, gradUx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradUy_send, gradUy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradUz_send, gradUz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradUx_send.clear();
		// gradUy_send.clear();
		// gradUz_send.clear();
					
		// mpi.setProcsFaceDatas(
					// gradVx_send, gradVx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradVy_send, gradVy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradVz_send, gradVz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradVx_send.clear();
		// gradVy_send.clear();
		// gradVz_send.clear();
					
		// mpi.setProcsFaceDatas(
					// gradWx_send, gradWx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradWy_send, gradWy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradWz_send, gradWz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradWx_send.clear();
		// gradWy_send.clear();
		// gradWz_send.clear(); 
	// }
	
	
	
	// vector<vector<double>> gradAi;
	// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcMGG(mesh, controls.VF[0], controls.fVF[0], 1, 1.e-8, gradAi);
	// // math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradAi);
	
	// vector<double> kappa;
	// this->calcCurvature(mesh, controls.VF[0], kappa);

	
	
	
	// // diagonal terms
	// vector<double> linA;
	// vector<double> linB0;
	// vector<double> linB1;
	// vector<double> linB2;
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// // time marching -> first order euler
		// linA.push_back(cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		
		// linB0.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.U] - 1.0*cell.var[controls.oldU]
			// )*cell.volume/controls.timeStep);
			
		// linB1.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.V] - 1.0*cell.var[controls.oldV]
			// )*cell.volume/controls.timeStep);
			
		// linB2.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.W] - 1.0*cell.var[controls.oldW]
			// )*cell.volume/controls.timeStep);
		
		// // time marching -> second order upwind euler
		// // linA.push_back(1.5*cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		// // linB0.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.5*cell.var[controls.U] - 2.0*cell.var[controls.oldU] + 0.5*cell.var[controls.Qm[0]]
			// // )*cell.volume/controls.timeStep);
			
		// // linB1.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.5*cell.var[controls.V] - 2.0*cell.var[controls.oldV] + 0.5*cell.var[controls.Qm[1]]
			// // )*cell.volume/controls.timeStep);
			
		// // linB2.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.5*cell.var[controls.W] - 2.0*cell.var[controls.oldW] + 0.5*cell.var[controls.Qm[2]]
			// // )*cell.volume/controls.timeStep);
			
		// // // gravity force terms
		// // linB0.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[0];
		// // linB1.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[1];
		// // linB2.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[2];
			
		// // // surface tension force terms
		// // linB0.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
		// // linB1.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
		// // linB2.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
	// }
	
	
	
	
	




	// // F. Denner, 2013
	// vector<double> A_p_vals(mesh.cells.size(), 0.0);
	// proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double wCL = face.wC;
		// // double wCL = face.wVC;
		// // double wCL = 0.5;
		// double wCR = 1.0-wCL;
		
		// double UnF = 
		    // (wCL*face.varL[controls.fU]*face.unitNormals[0] +
			 // wCL*face.varL[controls.fV]*face.unitNormals[1] +
			 // wCL*face.varL[controls.fW]*face.unitNormals[2] +
			 // wCR*face.varR[controls.fU]*face.unitNormals[0] +
			 // wCR*face.varR[controls.fV]*face.unitNormals[1] +
			 // wCR*face.varR[controls.fW]*face.unitNormals[2]);
			 
		// // if( UnF < std::numeric_limits<double>::min() ){
			// // UnF = std::numeric_limits<double>::min();
		// // }
			 
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = 1.0 - weightL;
		
		// double RhoF = 
			// weightL*face.varL[controls.fRho] +
			// weightR*face.varR[controls.fRho];
		
		// A_p_vals[face.owner] += weightL * RhoF*UnF*face.area / mesh.cells[face.owner].volume;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// A_p_vals[face.neighbour] -=  weightR * RhoF*UnF*face.area / mesh.cells[face.neighbour].volume;
		// }
		
		
	// }

	
	
	
	
	
	// vector<double> linAL(mesh.faces.size(),0.0);
	// vector<double> linAR(mesh.faces.size(),0.0);

	// proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		// double area = face.area;
		// vector<double> nvec;
		// nvec.push_back(face.unitNormals[0]);
		// nvec.push_back(face.unitNormals[1]);
		// nvec.push_back(face.unitNormals[2]);
		
		// vector<double> distanceCells;
		// distanceCells.push_back(face.distCells[0]);
		// distanceCells.push_back(face.distCells[1]);
		// distanceCells.push_back(face.distCells[2]);
		
		// double wCL = face.wC;
		// double wCR = 1.0-face.wC;
	
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
		

		// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		// // double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  // // pow(distanceCells[1],2.0) + 
						  // // pow(distanceCells[2],2.0));
		// double dPN = sqrt(distanceCells[0]*distanceCells[0] + 
						  // distanceCells[1]*distanceCells[1] + 
						  // distanceCells[2]*distanceCells[2]);

		// double nonOrtholimiter = 1.0;
		// // double cosAlpha = dPN_e/dPN;
		// // if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			// // nonOrtholimiter = 0.5;
		// // }
		// // else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			// // nonOrtholimiter = 0.333;
		// // }
		// // else if( cosAlpha < 0.342 ){
			// // nonOrtholimiter = 0.0;
		// // }
				
		// dPN = dPN_e;
				
		// vector<double> Ef(3,0.0);
		// Ef[0] = distanceCells[0]/dPN;
		// Ef[1] = distanceCells[1]/dPN;
		// Ef[2] = distanceCells[2]/dPN;
		// // // over-relaxed
		// // double over_Ef = Ef[0]*nvec[0] + Ef[1]*nvec[1] + Ef[2]*nvec[2];
		// // Ef[0] /= over_Ef; Ef[1] /= over_Ef; Ef[2] /= over_Ef;
		// // double magEf = sqrt(Ef[0]*Ef[0]+Ef[1]*Ef[1]+Ef[2]*Ef[2]);
		// // dPN /= magEf;
		// // // dPN = dPN_e;
			
		// vector<double> Tf(3,0.0);
		// Tf[0] = nvec[0] - Ef[0];
		// Tf[1] = nvec[1] - Ef[1];
		// Tf[2] = nvec[2] - Ef[2];
		
		// double UnF = wCL*UnL+wCR*UnR;
		
		// double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
		// double tmp1 = controls.timeStep / Rho_star;
		// // F.Denner et al., 2013
		// double d_F = 0.0;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// d_F = wCL/(A_p_vals[face.owner]) + wCR/(A_p_vals[face.neighbour]);
		// }
		// else{
			// d_F = 1.0/(A_p_vals[face.owner]);
		// }
		// double d_F_hat = d_F / (2.0 + d_F / tmp1);
		// if(d_F>1.e9) d_F_hat = tmp1;
		// tmp1 = d_F_hat;
		
		// double orgPL = mesh.cells[face.owner].var[controls.P];
		// double orgPR = 0.0;
		// vector<double> gradPL(3,0.0);
		// gradPL[0] = gradP[face.owner][0];
		// gradPL[1] = gradP[face.owner][1];
		// gradPL[2] = gradP[face.owner][2];
		// vector<double> gradPR(3,0.0);
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// orgPR = mesh.cells[face.neighbour].var[controls.P];
			// gradPR[0] = gradP[face.neighbour][0];
			// gradPR[1] = gradP[face.neighbour][1];
			// gradPR[2] = gradP[face.neighbour][2];
			
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// orgPR = PR;
			// gradPR[0] = gradPx_recv[proc_num];
			// gradPR[1] = gradPy_recv[proc_num];
			// gradPR[2] = gradPz_recv[proc_num];
			
		// }
		// for(int ii=0; ii<3; ++ii){
			// UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			// UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
			// // UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*Ef[ii];
			// // UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*Ef[ii];
		// }
			
		// UnF -= tmp1*(orgPR-orgPL)/dPN;
		
		// // // non-orthogonal, over-relaxed approach
		// // for(int ii=0; ii<3; ++ii){
			// // // UnF -= nonOrtholimiter * tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*Tf[ii];
			// // // UnF -= nonOrtholimiter * tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*Tf[ii];
			// // UnF -= nonOrtholimiter * tmp1 * wCL*gradPL[ii]*Tf[ii];
			// // UnF -= nonOrtholimiter * tmp1 * wCR*gradPR[ii]*Tf[ii];
		// // }
		
		
		// // // F.Denner et al., 2018
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // UnF += face.var[controls.Un];
			// // UnF -= 0.5*mesh.cells[face.owner].var[controls.oldU]*nvec[0];
			// // UnF -= 0.5*mesh.cells[face.owner].var[controls.oldV]*nvec[1];
			// // UnF -= 0.5*mesh.cells[face.owner].var[controls.oldW]*nvec[2];
			// // UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldU]*nvec[0];
			// // UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldV]*nvec[1];
			// // UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldW]*nvec[2];
			// // if(controls.iterPBs+1==controls.iterPBsMax){
				// // face.var[controls.Un] = UnF;
			// // }
		// // }
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = 1.0 - weightL;
		
		// double UF = UL*weightL + UR*weightR;
		// double VF = VL*weightL + VR*weightR;
		// double WF = WL*weightL + WR*weightR;
		// double RhoF = RhoL*weightL + RhoR*weightR;
		
		// double PF = wCL*PL + wCR*PR;
		
		
		// double conv_flux_L = weightL*RhoL*UnF*area;
		// double conv_flux_R = weightR*RhoR*UnF*area;
		
		// // properties, viscous term
		// double muL = face.varL[controls.fmu];
		// double muR = face.varR[controls.fmu];
		// double muF = wCL*muL + wCR*muR;
		// double diff_flux = muF/dPN_e*area;
		
		
		// // viscous term, properties
		// double orgUL = mesh.cells[face.owner].var[controls.U];
		// double orgUR = 0.0;
		// double orgVL = mesh.cells[face.owner].var[controls.V];
		// double orgVR = 0.0;
		// double orgWL = mesh.cells[face.owner].var[controls.W];
		// double orgWR = 0.0;
		
		// double gradUf[3];
		// double gradVf[3];
		// double gradWf[3];
		
		// double vNonOrthFlux0 = 0.0;
		// double vNonOrthFlux1 = 0.0;
		// double vNonOrthFlux2 = 0.0;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// orgUR = mesh.cells[face.neighbour].var[controls.U];
			// orgVR = mesh.cells[face.neighbour].var[controls.V];
			// orgWR = mesh.cells[face.neighbour].var[controls.W];
			
			// gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradU[face.neighbour][0]);
			// gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradU[face.neighbour][1]);
			// gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradU[face.neighbour][2]);
			
			// gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradV[face.neighbour][0]);
			// gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradV[face.neighbour][1]);
			// gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradV[face.neighbour][2]);
			
			// gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradW[face.neighbour][0]);
			// gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradW[face.neighbour][1]);
			// gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradW[face.neighbour][2]);
			
			// // // viscous term, non-orthogonal
			// // vNonOrthFlux0 = nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			// // vNonOrthFlux1 = nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			// // vNonOrthFlux2 = nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// orgUR = UR;
			// orgVR = VR;
			// orgWR = WR;
			
			// gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradUx_recv[proc_num]);
			// gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradUy_recv[proc_num]);
			// gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradUz_recv[proc_num]);
			
			// gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradVx_recv[proc_num]);
			// gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradVy_recv[proc_num]);
			// gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradVz_recv[proc_num]);
			
			// gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradWx_recv[proc_num]);
			// gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradWy_recv[proc_num]);
			// gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradWz_recv[proc_num]);
			
			// // // viscous term, non-orthogonal
			// // vNonOrthFlux0 = nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			// // vNonOrthFlux1 = nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			// // vNonOrthFlux2 = nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
			
			
		// }
		
		// // own, ngb
		// linAR[i] = (+conv_flux_R -diff_flux);
		// linA[face.owner] += (+conv_flux_L +diff_flux);
		// // ngb, own
		// linAL[i] = (-conv_flux_L -diff_flux);
		
		// // convective term
		// double conv_flux_x = 0.0;
		// double conv_flux_y = 0.0;
		// double conv_flux_z = 0.0;

		// // convective term	
		// linB0[face.owner] -= UF*RhoL*UnF*area;
		// linB1[face.owner] -= VF*RhoL*UnF*area;
		// linB2[face.owner] -= WF*RhoL*UnF*area;
		
		// // pressure term
		// linB0[face.owner] -= PF*nvec[0]*area;
		// linB1[face.owner] -= PF*nvec[1]*area;
		// linB2[face.owner] -= PF*nvec[2]*area;
		
		// // // viscous term
		// // linB0[face.owner] += muF*(orgUR-orgUL)/dPN_e*area;
		// // linB1[face.owner] += muF*(orgVR-orgVL)/dPN_e*area;
		// // linB2[face.owner] += muF*(orgWR-orgWL)/dPN_e*area;
		
		// // // viscous non-orthogonal term
		// // linB0[face.owner] += vNonOrthFlux0;
		// // linB1[face.owner] += vNonOrthFlux1;
		// // linB2[face.owner] += vNonOrthFlux2;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // ngb, own
			// linA[face.neighbour] += (-conv_flux_R +diff_flux);
			
			// // convective term	
			// linB0[face.neighbour] += UF*RhoR*UnF*area;
			// linB1[face.neighbour] += VF*RhoR*UnF*area;
			// linB2[face.neighbour] += WF*RhoR*UnF*area;
			
			// // pressure term
			// linB0[face.neighbour] += PF*nvec[0]*area;
			// linB1[face.neighbour] += PF*nvec[1]*area;
			// linB2[face.neighbour] += PF*nvec[2]*area;
			
			// // // viscous term
			// // linB0[face.neighbour] -= muF*(orgUR-orgUL)/dPN_e*area;
			// // linB1[face.neighbour] -= muF*(orgVR-orgVL)/dPN_e*area;
			// // linB2[face.neighbour] -= muF*(orgWR-orgWL)/dPN_e*area;
			
			// // // viscous non-orthogonal term
			// // linB0[face.neighbour] -= vNonOrthFlux0;
			// // linB1[face.neighbour] -= vNonOrthFlux1;
			// // linB2[face.neighbour] -= vNonOrthFlux2;
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// ++proc_num;
			
		// }
	// }
	
	
	
	
	
	
	
	
	
	
	

	
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				
				// double area = face.area;
				// vector<double> nvec(3,0.0);
				// nvec[0] = face.unitNormals[0];
				// nvec[1] = face.unitNormals[1];
				// nvec[2] = face.unitNormals[2];
		
				// double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
							 // face.varL[controls.fV]*face.unitNormals[1] + 
							 // face.varL[controls.fW]*face.unitNormals[2];
				// double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
							 // face.varR[controls.fV]*face.unitNormals[1] + 
							 // face.varR[controls.fW]*face.unitNormals[2];
				
				// double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							   // face.unitNormals[1]*face.distCells[1] + 
							   // face.unitNormals[2]*face.distCells[2];
				
				// double dPN = sqrt(pow(face.distCells[0],2.0) + 
								  // pow(face.distCells[1],2.0) + 
								  // pow(face.distCells[2],2.0));
							
				// double UnF = 0.5*UnL+0.5*UnR;
				
				// double RhoL = face.varL[controls.fRho];
				// double tmp1 = 1.0/RhoL*controls.timeStep;
				
				// // for(int ii=0; ii<3; ++ii){
					// // UnF += tmp1 * gradP[face.owner][ii]*nvec[ii];
				// // }
				
				// double UF = 0.5*face.varL[controls.fU] + 0.5*face.varR[controls.fU];
				// double VF = 0.5*face.varL[controls.fV] + 0.5*face.varR[controls.fV];
				// double WF = 0.5*face.varL[controls.fW] + 0.5*face.varR[controls.fW];
				// double PF = 0.5*face.varL[controls.fP] + 0.5*face.varR[controls.fP];
			
				
				// // convective term
				// linB0[face.owner] -= ( RhoL * UF * UnF * area + PF * face.unitNormals[0] * area );
				// linB1[face.owner] -= ( RhoL * VF * UnF * area + PF * face.unitNormals[1] * area );
				// linB2[face.owner] -= ( RhoL * WF * UnF * area + PF * face.unitNormals[2] * area );
			// }
			
			
		// }
		
	// }
	
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar0(mesh.cells.size(),0.0);
	// vector<double> resiVar1(mesh.cells.size(),0.0);
	// vector<double> resiVar2(mesh.cells.size(),0.0);
	
	// string solver;
	// double relTol;
	// double tolerance;
	// string preconditioner;
	// int maxIter;
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
		// solver = controls.solverU;
		// relTol = controls.toleranceU;
		// tolerance = controls.relTolU;
		// preconditioner = controls.preconditionerU;
		// maxIter = controls.maxIterU;
	
	// }
	// else{
		// solver = controls.solverFinalU;
		// relTol = controls.toleranceFinalU;
		// tolerance = controls.relTolFinalU;
		// preconditioner = controls.preconditionerFinalU;
		// maxIter = controls.maxIterFinalU;
	
	// }
	
	// // solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
		// // solver, tolerance, relTol, preconditioner, maxIter);
	// // solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
		// // solver, tolerance, relTol, preconditioner, maxIter);
	// // solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
		// // solver, tolerance, relTol, preconditioner, maxIter);

	// solveHYPRE(mesh, resiVar0, linA, linAL, linAR, linB0,
		// solver, tolerance, relTol, preconditioner, maxIter);
	// solveHYPRE(mesh, resiVar1, linA, linAL, linAR, linB1,
		// solver, tolerance, relTol, preconditioner, maxIter);
	// solveHYPRE(mesh, resiVar2, linA, linAL, linAR, linB2,
		// solver, tolerance, relTol, preconditioner, maxIter);
		
	

	// // update U, V, W
	// double normU = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// cell.var[controls.U] += controls.momVelURF * resiVar0[i];
		// cell.var[controls.V] += controls.momVelURF * resiVar1[i];
		// cell.var[controls.W] += controls.momVelURF * resiVar2[i];
		
		
		// if( cell.var[controls.U] <= controls.minU ) 
			// cell.var[controls.U] = controls.minU;
		// if( cell.var[controls.U] >= controls.maxU ) 
			// cell.var[controls.U] = controls.maxU;
		
		// if( cell.var[controls.V] <= controls.minV ) 
			// cell.var[controls.V] = controls.minV;
		// if( cell.var[controls.V] >= controls.maxV ) 
			// cell.var[controls.V] = controls.maxV;
		
		// if( cell.var[controls.W] <= controls.minW ) 
			// cell.var[controls.W] = controls.minW;
		// if( cell.var[controls.W] >= controls.maxW ) 
			// cell.var[controls.W] = controls.maxW;
		
		// normU += pow(resiVar0[i],2.0);
		// normU += pow(resiVar1[i],2.0);
		// normU += pow(resiVar2[i],2.0);
		
	// }
	
	
	
	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB0.clear();
	// linB1.clear();
	// linB2.clear();
	// resiVar0.clear();
	// resiVar1.clear();
	// resiVar2.clear();
	
	// return sqrt(normU);
	
	
}















// #include "build.h"
// #include <cmath>
// #include <array>
// #include <numeric>




// double SEMO_Solvers_Builder::calcMomentumEqs(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<SEMO_Species>& species,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& linAFL,
	// vector<double>& linAFR,
	// vector<double>& residuals){
		
		
	
	// double coeffPf = 0.05;
	// // double coeffPf = 0.0;
	// double coeffDiss = 0.9;
	
		
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
		
	// SEMO_MPI_Builder mpi;


	// SEMO_Utility_Math math;
	
	// // gradient P
	// vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
	// // vector<double> limGradP;
	// // math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // for(int j=0; j<3; ++j){
			// // gradP[i][j] *= limGradP[i]; 
		// // }
	// // }
	
	
	
	// // // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	
	
	
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
		
		// mpi.setProcsFaceDatas(
					// gradPx_send, gradPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPy_send, gradPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPz_send, gradPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	// }
	

	
	// // gradient U, V, W
	// vector<vector<double>> gradU;
	// vector<vector<double>> gradV;
	// vector<vector<double>> gradW;
	// // math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// // math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// // math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// // math.calcMGG(mesh, controls.U, controls.fU, 10, 1.e-8, gradU);
	// // math.calcMGG(mesh, controls.V, controls.fV, 10, 1.e-8, gradV);
	// // math.calcMGG(mesh, controls.W, controls.fW, 10, 1.e-8, gradW);
	// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	// vector<double> gradUx_recv, gradUy_recv, gradUz_recv;
	// vector<double> gradVx_recv, gradVy_recv, gradVz_recv;
	// vector<double> gradWx_recv, gradWy_recv, gradWz_recv;
	// if(size>1){
		// // processor faces
		// vector<double> gradUx_send, gradUy_send, gradUz_send;
		// vector<double> gradVx_send, gradVy_send, gradVz_send;
		// vector<double> gradWx_send, gradWy_send, gradWz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradUx_send.push_back(gradU[face.owner][0]);
				// gradUy_send.push_back(gradU[face.owner][1]);
				// gradUz_send.push_back(gradU[face.owner][2]);
				
				// gradVx_send.push_back(gradV[face.owner][0]);
				// gradVy_send.push_back(gradV[face.owner][1]);
				// gradVz_send.push_back(gradV[face.owner][2]);
				
				// gradWx_send.push_back(gradW[face.owner][0]);
				// gradWy_send.push_back(gradW[face.owner][1]);
				// gradWz_send.push_back(gradW[face.owner][2]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// gradUx_send, gradUx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradUy_send, gradUy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradUz_send, gradUz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradUx_send.clear();
		// gradUy_send.clear();
		// gradUz_send.clear();
					
		// mpi.setProcsFaceDatas(
					// gradVx_send, gradVx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradVy_send, gradVy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradVz_send, gradVz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradVx_send.clear();
		// gradVy_send.clear();
		// gradVz_send.clear();
					
		// mpi.setProcsFaceDatas(
					// gradWx_send, gradWx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradWy_send, gradWy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradWz_send, gradWz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradWx_send.clear();
		// gradWy_send.clear();
		// gradWz_send.clear();
	// }
	
	
		
	// vector<double> linA;
	// vector<double> linB0;
	// vector<double> linB1;
	// vector<double> linB2;
	
	
	
	// vector<vector<double>> gradAi;
	// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcMGG(mesh, controls.VF[0], controls.fVF[0], 1, 1.e-8, gradAi);
	// // math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradAi);
	
	// vector<double> kappa;
	// this->calcCurvature(mesh, controls.VF[0], kappa);
	
	
	
	
	
	
	
	// //=========================
	// // Un initialization
	// if( controls.iterPBs == 0 && controls.iterMom == 0){
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];

			// vector<double> nvec;
			// nvec.push_back(face.unitNormals[0]);
			// nvec.push_back(face.unitNormals[1]);
			// nvec.push_back(face.unitNormals[2]);
			
			// double wCL = face.wC;
			// double wCR = 1.0-face.wC;
		 
			// double UL = face.varL[controls.fU];
			// double VL = face.varL[controls.fV];
			// double WL = face.varL[controls.fW];
			
			// double UR = face.varR[controls.fU];
			// double VR = face.varR[controls.fV];
			// double WR = face.varR[controls.fW];
			
			// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
			// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
			
			// face.var[controls.Un] = wCL*UnL+wCR*UnR;
			
		// }
	// }
	// //=========================
	
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// // time marching -> first order euler
		// // linA.push_back(cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		// // linB0.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.0*cell.var[controls.U] - 1.0*cell.var[controls.oldU]
			// // )*cell.volume/controls.timeStep);
			
		// // linB1.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.0*cell.var[controls.V] - 1.0*cell.var[controls.oldV]
			// // )*cell.volume/controls.timeStep);
			
		// // linB2.push_back(-
			// // cell.var[controls.Rho]*(
			// // 1.0*cell.var[controls.W] - 1.0*cell.var[controls.oldW]
			// // )*cell.volume/controls.timeStep);
		
		// // time marching -> second order upwind euler
		// linA.push_back(1.5*cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		// linB0.push_back(-
			// cell.var[controls.Rho]*(
			// 1.5*cell.var[controls.U] - 2.0*cell.var[controls.oldU] + 0.5*cell.var[controls.Qm[0]]
			// )*cell.volume/controls.timeStep);
			
		// linB1.push_back(-
			// cell.var[controls.Rho]*(
			// 1.5*cell.var[controls.V] - 2.0*cell.var[controls.oldV] + 0.5*cell.var[controls.Qm[1]]
			// )*cell.volume/controls.timeStep);
			
		// linB2.push_back(-
			// cell.var[controls.Rho]*(
			// 1.5*cell.var[controls.W] - 2.0*cell.var[controls.oldW] + 0.5*cell.var[controls.Qm[2]]
			// )*cell.volume/controls.timeStep);
			
		// // gravity force terms
		// linB0.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[0];
		// linB1.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[1];
		// linB2.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[2];
			
		// // surface tension force terms
		// linB0.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
		// linB1.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
		// linB2.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
	// }
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// vector<double> linAL;
	// vector<double> linAR;

	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// // face.area = 0.00025;
		
		// double area = face.area;
		// vector<double> nvec;
		// nvec.push_back(face.unitNormals[0]);
		// nvec.push_back(face.unitNormals[1]);
		// nvec.push_back(face.unitNormals[2]);
		
		// vector<double> distanceCells;
		// distanceCells.push_back(face.distCells[0]);
		// distanceCells.push_back(face.distCells[1]);
		// distanceCells.push_back(face.distCells[2]);
		
		// double wCL = face.wC;
		// double wCR = 1.0-face.wC;
	
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
		

		// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		// double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  // pow(distanceCells[1],2.0) + 
						  // pow(distanceCells[2],2.0));
		// // double dPNEff = dPN 
			// // /( (nvec[0]*distanceCells[0] + 
			    // // nvec[1]*distanceCells[1] + 
				// // nvec[2]*distanceCells[2])/dPN );

		// double nonOrtholimiter = 1.0;
		// double cosAlpha = dPN_e/dPN;
		// if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			// nonOrtholimiter = 0.5;
		// }
		// else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			// nonOrtholimiter = 0.333;
		// }
		// else if( cosAlpha < 0.342 ){
			// nonOrtholimiter = 0.0;
		// }
		// // nonOrtholimiter = 0.0;
				
				
		// double Ef[3];
		// Ef[0] = distanceCells[0]/dPN;
		// Ef[1] = distanceCells[1]/dPN;
		// Ef[2] = distanceCells[2]/dPN;
				
		// // over-relaxed approach
		// // double overRelaxCoeff = 1.0;
		// // double overRelaxCoeff = Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2];
		// // double overRelaxCoeff = 1.0 / (Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2]);
		
		// // non-orthogonal, over-relaxed approach
		// // double Tf[3];
		// // Tf[0] = nvec[0] - Ef[0]*overRelaxCoeff;
		// // Tf[1] = nvec[1] - Ef[1]*overRelaxCoeff;
		// // Tf[2] = nvec[2] - Ef[2]*overRelaxCoeff;
		// double Tf[3];
		// Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		// Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		// Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		// // double dPNEff = dPN / overRelaxCoeff;
		// // double dPNEff = dPN;
		
		
		
		
		
		

		
		
		// double UnF = wCL*UnL+wCR*UnR;
		// // double UnF = 0.5*UnL+0.5*UnR;
		// // double UnF = face.var[controls.Un];

		// double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double orgPL = mesh.cells[face.owner].var[controls.P];
			// double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			// UnF += (
				// wCL*( 
				// (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// wCR*( 
				// (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			// // UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// double orgPL = mesh.cells[face.owner].var[controls.P];
			// double orgPR = PR;
			
			// UnF += (
				// wCL*( 
				// (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// wCR*( 
				// (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			// UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			// // UnF -=  nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
		// }
		
		

		// // if(
		// // face.getType() == SEMO_Types::INTERNAL_FACE ||
		// // face.getType() == SEMO_Types::PROCESSOR_FACE
		// // ){
			// // UnF = face.var[controls.Un];
		// // }
		
		
		
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double UF = UL*weightL + UR*weightR;
		// double VF = VL*weightL + VR*weightR;
		// double WF = WL*weightL + WR*weightR;
		
		// double PF = wCL*PL + wCR*PR;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			// double barA1 = RhoR*mesh.cells[face.neighbour].volume/controls.timeStep;
			
			// // PF = (barA0*PL + barA1*PR)/(barA0 + barA1);
			
			// PF += coeffPf * 0.5*barA0*barA1/(barA0+barA1)*(UnL-UnR)/area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			// double barA1 = RhoR*mesh.cells[face.owner].volume/controls.timeStep;
			
			// // PF = (barA0*PL + barA1*PR)/(barA0 + barA1);
			
			// PF += coeffPf * 0.5*barA0*barA1/(barA0+barA1)*(UnL-UnR)/area;
		// }
		// // else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// // double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			
			// // PF += coeffPf * 0.5*barA0*(UnL-UnR)/area;
		// // }
		
		
		
		
		
		// linAL.push_back(-weightL*RhoL*UnF*area);
		// linAR.push_back(+weightR*RhoR*UnF*area);
		
		// // properties, viscous term
		// double muL = face.varL[controls.fmu];
		// double muR = face.varR[controls.fmu];
		// double muF = wCL*muL + wCR*muR;
		
		// linAL.back() -= muF/dPN_e*area;
		// linAR.back() -= muF/dPN_e*area;
		
			
		// // pressure term
		// linB0.at(face.owner) -= PF*nvec[0]*area;
		// linB1.at(face.owner) -= PF*nvec[1]*area;
		// linB2.at(face.owner) -= PF*nvec[2]*area;
		
		
		// gradP[face.owner][0] += PF*nvec[0]*area;
		// gradP[face.owner][1] += PF*nvec[1]*area;
		// gradP[face.owner][2] += PF*nvec[2]*area;
		
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA.at(face.owner)     += (-linAL.back());
			// linA.at(face.neighbour) += (-linAR.back());
			
			// // convective term
			// linB0.at(face.owner) -= UF*RhoL*UnF*area;
			// linB1.at(face.owner) -= VF*RhoL*UnF*area;
			// linB2.at(face.owner) -= WF*RhoL*UnF*area;
		
			// linB0.at(face.neighbour) += UF*RhoR*UnF*area;
			// linB1.at(face.neighbour) += VF*RhoR*UnF*area;
			// linB2.at(face.neighbour) += WF*RhoR*UnF*area;
			
			// // pressure term
			// linB0.at(face.neighbour) += PF*nvec[0]*area;
			// linB1.at(face.neighbour) += PF*nvec[1]*area;
			// linB2.at(face.neighbour) += PF*nvec[2]*area;
			
			
			
			// gradP[face.neighbour][0] -= PF*nvec[0]*area;
			// gradP[face.neighbour][1] -= PF*nvec[1]*area;
			// gradP[face.neighbour][2] -= PF*nvec[2]*area;
			
			
			// // viscous term
			// double orgUL = mesh.cells[face.owner].var[controls.U];
			// double orgUR = mesh.cells[face.neighbour].var[controls.U];
			// double orgVL = mesh.cells[face.owner].var[controls.V];
			// double orgVR = mesh.cells[face.neighbour].var[controls.V];
			// double orgWL = mesh.cells[face.owner].var[controls.W];
			// double orgWR = mesh.cells[face.neighbour].var[controls.W];
			
			// linB0.at(face.owner) += muF*(orgUR-orgUL)/dPN_e*area;
			// linB1.at(face.owner) += muF*(orgVR-orgVL)/dPN_e*area;
			// linB2.at(face.owner) += muF*(orgWR-orgWL)/dPN_e*area;
		
			// linB0.at(face.neighbour) -= muF*(orgUR-orgUL)/dPN_e*area;
			// linB1.at(face.neighbour) -= muF*(orgVR-orgVL)/dPN_e*area;
			// linB2.at(face.neighbour) -= muF*(orgWR-orgWL)/dPN_e*area;
			
			// // viscous term, properties
			// double delPhiECF;
			
			// double gradUf[3];
			// gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradU[face.neighbour][0]);
			// gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradU[face.neighbour][1]);
			// gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradU[face.neighbour][2]);
			// // delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
			
			// // gradUf[0] += nvec[0]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// // gradUf[1] += nvec[1]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// // gradUf[2] += nvec[2]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			
			// double gradVf[3];
			// gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradV[face.neighbour][0]);
			// gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradV[face.neighbour][1]);
			// gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradV[face.neighbour][2]);
			// // delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
			
			// // gradVf[0] += nvec[0]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// // gradVf[1] += nvec[1]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// // gradVf[2] += nvec[2]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			
			// double gradWf[3];
			// gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradW[face.neighbour][0]);
			// gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradW[face.neighbour][1]);
			// gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradW[face.neighbour][2]);
			// // delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
			
			// // gradWf[0] += nvec[0]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// // gradWf[1] += nvec[1]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// // gradWf[2] += nvec[2]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			
			
			// // viscous term, non-orthogonal
			// double vNonOrthFlux0 = 
				// nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			// double vNonOrthFlux1 = 
				// nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			// double vNonOrthFlux2 = 
				// nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
			// linB0[face.owner] += vNonOrthFlux0;
			// linB1[face.owner] += vNonOrthFlux1;
			// linB2[face.owner] += vNonOrthFlux2;
				
			// linB0[face.neighbour] -= vNonOrthFlux0;
			// linB1[face.neighbour] -= vNonOrthFlux1;
			// linB2[face.neighbour] -= vNonOrthFlux2;
			
			
			// // double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);;
			// // linAL.back() -= muF/dPN*area*tmpCoeff2;
			// // linAR.back() -= muF/dPN*area*tmpCoeff2;
			
			// // linA.at(face.owner)     += muF/dPN*area*tmpCoeff2;
			// // linA.at(face.neighbour) += muF/dPN*area*tmpCoeff2;
			
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA.at(face.owner) += (-linAL.back());

			// // convective term
			// linB0.at(face.owner) -= UF*RhoL*UnF*area;
			// linB1.at(face.owner) -= VF*RhoL*UnF*area;
			// linB2.at(face.owner) -= WF*RhoL*UnF*area;
			

			// // viscous term
			// double orgUL = mesh.cells[face.owner].var[controls.U];
			// double orgUR = UR;
			// double orgVL = mesh.cells[face.owner].var[controls.V];
			// double orgVR = VR;
			// double orgWL = mesh.cells[face.owner].var[controls.W];
			// double orgWR = WR;
			
			// linB0.at(face.owner) += muF*(orgUR-orgUL)/dPN_e*area;
			// linB1.at(face.owner) += muF*(orgVR-orgVL)/dPN_e*area;
			// linB2.at(face.owner) += muF*(orgWR-orgWL)/dPN_e*area;
		
			// // viscous term, properties
			// double delPhiECF;
			
			// double gradUf[3];
			// gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradUx_recv[proc_num]);
			// gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradUy_recv[proc_num]);
			// gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradUz_recv[proc_num]);
			// // delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
			
			// // gradUf[0] += nvec[0]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// // gradUf[1] += nvec[1]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// // gradUf[2] += nvec[2]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			
			// double gradVf[3];
			// gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradVx_recv[proc_num]);
			// gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradVy_recv[proc_num]);
			// gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradVz_recv[proc_num]);
			// // delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
			
			// // gradVf[0] += nvec[0]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// // gradVf[1] += nvec[1]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// // gradVf[2] += nvec[2]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			
			// double gradWf[3];
			// gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradWx_recv[proc_num]);
			// gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradWy_recv[proc_num]);
			// gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradWz_recv[proc_num]);
			// // delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
			
			// // gradWf[0] += nvec[0]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// // gradWf[1] += nvec[1]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// // gradWf[2] += nvec[2]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			
			// // viscous term, non-orthogonal
			// double vNonOrthFlux0 = 
				// nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			// double vNonOrthFlux1 = 
				// nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			// double vNonOrthFlux2 = 
				// nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
			// linB0[face.owner] += vNonOrthFlux0;
			// linB1[face.owner] += vNonOrthFlux1;
			// linB2[face.owner] += vNonOrthFlux2;
			
			
			// // double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);;
			// // linAL.back() -= muF/dPN*area*tmpCoeff2;
			
			// // linA.at(face.owner)     += muF/dPN*area*tmpCoeff2;
			
			
			
			
			// ++proc_num;
			
		// }
	// }
	
	
	// // // pressure term
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // linB0[i] -= cell.volume*gradP[i][0];
		// // linB1[i] -= cell.volume*gradP[i][1];
		// // linB2[i] -= cell.volume*gradP[i][2];
		
		// // // cout << cell.volume*gradP[i][1] << endl;
	// // }
	
	
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(
			// boundary.type[1] == "noSlip" ||
			// boundary.type[1] == "slip" ){
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec;
					// nvec.push_back(face.unitNormals[0]);
					// nvec.push_back(face.unitNormals[1]);
					// nvec.push_back(face.unitNormals[2]);
					
					// vector<double> distanceCells;
					// distanceCells.push_back(face.distCells[0]);
					// distanceCells.push_back(face.distCells[1]);
					// distanceCells.push_back(face.distCells[2]);
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
		
					// double UnL = UL*face.unitNormals[0]+
					             // VL*face.unitNormals[1]+
								 // WL*face.unitNormals[2];
					// double UnR = UR*face.unitNormals[0]+
					             // VR*face.unitNormals[1]+
								 // WR*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// // // convective term
					// // linA.at(face.owner) += RhoL*UnL*area;
					
					// linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					// linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					// linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					
				// }
			// }
			// else if(
			// boundary.type[1] == "fixedValue" ||
			// boundary.type[1] == "surfaceNormalFixedValue"
			// ){
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec;
					// nvec.push_back(face.unitNormals[0]);
					// nvec.push_back(face.unitNormals[1]);
					// nvec.push_back(face.unitNormals[2]);
					
					// vector<double> distanceCells;
					// distanceCells.push_back(face.distCells[0]);
					// distanceCells.push_back(face.distCells[1]);
					// distanceCells.push_back(face.distCells[2]);
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double orgUL = mesh.cells[face.owner].var[controls.U];
					// double orgVL = mesh.cells[face.owner].var[controls.V];
					// double orgWL = mesh.cells[face.owner].var[controls.W];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
					
					// double orgUnL = orgUL*face.unitNormals[0]+
					             // orgVL*face.unitNormals[1]+
								 // orgWL*face.unitNormals[2];
		
					// double UnL = UL*face.unitNormals[0]+
					             // VL*face.unitNormals[1]+
								 // WL*face.unitNormals[2];
					// double UnR = UR*face.unitNormals[0]+
					             // VR*face.unitNormals[1]+
								 // WR*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// // convective term
					// // linA.at(face.owner) += 0.5*RhoL*UnL*area;
					
					// // linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					// // linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					// // linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					// linB0[face.owner] -= UR*RhoR*UnF*area;
					// linB1[face.owner] -= VR*RhoR*UnF*area;
					// linB2[face.owner] -= WR*RhoR*UnF*area;
					// // linB0[face.owner] -= 0.5*(orgUL+UR)*RhoL*0.5*(orgUnL+UnR)*area;
					// // linB1[face.owner] -= 0.5*(orgVL+VR)*RhoL*0.5*(orgUnL+UnR)*area;
					// // linB2[face.owner] -= 0.5*(orgWL+WR)*RhoL*0.5*(orgUnL+UnR)*area;
					
				// }
			// }
			// else if(boundary.type[1] == "zeroGradient"){
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec;
					// nvec.push_back(face.unitNormals[0]);
					// nvec.push_back(face.unitNormals[1]);
					// nvec.push_back(face.unitNormals[2]);
					
					// vector<double> distanceCells;
					// distanceCells.push_back(face.distCells[0]);
					// distanceCells.push_back(face.distCells[1]);
					// distanceCells.push_back(face.distCells[2]);
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
		
					// double UnL = UL*face.unitNormals[0]+
					             // VL*face.unitNormals[1]+
								 // WL*face.unitNormals[2];
					// double UnR = UR*face.unitNormals[0]+
					             // VR*face.unitNormals[1]+
								 // WR*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// double orgUL = mesh.cells[face.owner].var[controls.U];
					// double orgVL = mesh.cells[face.owner].var[controls.V];
					// double orgWL = mesh.cells[face.owner].var[controls.W];
					
					// // convective term
					// linA[face.owner] += RhoL*UnF*area;
						
					// linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					// linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					// linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					
				// }
			// }
			// else if(boundary.type[1] == "inletOutlet"){
				
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					
					// double area = face.area;
					// vector<double> nvec;
					// nvec.push_back(face.unitNormals[0]);
					// nvec.push_back(face.unitNormals[1]);
					// nvec.push_back(face.unitNormals[2]);
					
					// vector<double> distanceCells;
					// distanceCells.push_back(face.distCells[0]);
					// distanceCells.push_back(face.distCells[1]);
					// distanceCells.push_back(face.distCells[2]);
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double orgUL = mesh.cells[face.owner].var[controls.U];
					// double orgVL = mesh.cells[face.owner].var[controls.V];
					// double orgWL = mesh.cells[face.owner].var[controls.W];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
					
					// double orgUnL = orgUL*face.unitNormals[0]+
					             // orgVL*face.unitNormals[1]+
								 // orgWL*face.unitNormals[2];
		
					// double UnL = UL*face.unitNormals[0]+
					             // VL*face.unitNormals[1]+
								 // WL*face.unitNormals[2];
					// double UnR = UR*face.unitNormals[0]+
					             // VR*face.unitNormals[1]+
								 // WR*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// // linB0[face.owner] -= UR*RhoR*UnF*area;
					// // linB1[face.owner] -= VR*RhoR*UnF*area;
					// // linB2[face.owner] -= WR*RhoR*UnF*area;
					
					// //=========
					// if( UnF > 0.0 ){
						// // outflow
						// // zeroGradient
						// linA[face.owner] += RhoL*UnF*area;
							
						// linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
						// linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
						// linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					// }
					// else{
						// // inflow
						// // fixedValue
						
						// linB0[face.owner] -= UR*RhoR*UnF*area;
						// linB1[face.owner] -= VR*RhoR*UnF*area;
						// linB2[face.owner] -= WR*RhoR*UnF*area;
					// }
					
					
						
					
					
				// }
					
			// }
			// else{
				
				// if(rank==0) cerr << "| #Error : not defined velocity B.C. name" << endl;
				// MPI_Barrier(MPI_COMM_WORLD);
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				
			// }
			
			
			
			// // diffusion term
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				
				// double area = face.area;
				// vector<double> nvec;
				// nvec.push_back(face.unitNormals[0]);
				// nvec.push_back(face.unitNormals[1]);
				// nvec.push_back(face.unitNormals[2]);
				
				// vector<double> distanceCells;
				// distanceCells.push_back(face.distCells[0]);
				// distanceCells.push_back(face.distCells[1]);
				// distanceCells.push_back(face.distCells[2]);
				
				// double UL = face.varL[controls.fU];
				// double VL = face.varL[controls.fV];
				// double WL = face.varL[controls.fW];
				
				// double UR = face.varR[controls.fU];
				// double VR = face.varR[controls.fV];
				// double WR = face.varR[controls.fW];
	
				// // diffusion term
				// double muF = face.varL[controls.fmu];
				// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
				// double dPN = sqrt(pow(distanceCells[0],2.0) + 
								  // pow(distanceCells[1],2.0) + 
								  // pow(distanceCells[2],2.0));
				// // double Ef[3];
				// // Ef[0] = distanceCells[0]/dPN;
				// // Ef[1] = distanceCells[1]/dPN;
				// // Ef[2] = distanceCells[2]/dPN;
				 
				// // // over-relaxed approach
				// // double overRelaxCoeff = 1.0;
				// // // double overRelaxCoeff = Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2];
				// // // double overRelaxCoeff = 1.0 / (Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2]);
				
				// // // double dPNEff = dPN * overRelaxCoeff;
				
				// // // double Tf[3];
				// // // Tf[0] = nvec[0] - Ef[0];
				// // // Tf[1] = nvec[1] - Ef[1];
				// // // Tf[2] = nvec[2] - Ef[2];
					
				// double orgUL = mesh.cells[face.owner].var[controls.U];
				// double orgVL = mesh.cells[face.owner].var[controls.V];
				// double orgWL = mesh.cells[face.owner].var[controls.W];
				
				// linA[face.owner] += muF/dPN_e*area;
				
				

				// double corrUL = orgUL;
				// double corrVL = orgVL;
				// double corrWL = orgWL;
				// // double tmpVec[3];
				// // tmpVec[0] = distanceCells[0]*(1.0-
					// // (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
				// // tmpVec[1] = distanceCells[1]*(1.0-
					// // (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
				// // tmpVec[2] = distanceCells[2]*(1.0-
					// // (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
					
				// // corrUL += gradP[face.owner][0]*tmpVec[0];
				// // corrUL += gradP[face.owner][1]*tmpVec[1];
				// // corrUL += gradP[face.owner][2]*tmpVec[2];
				
				// // corrVL += gradP[face.owner][0]*tmpVec[0];
				// // corrVL += gradP[face.owner][1]*tmpVec[1];
				// // corrVL += gradP[face.owner][2]*tmpVec[2];
				
				// // corrWL += gradP[face.owner][0]*tmpVec[0];
				// // corrWL += gradP[face.owner][1]*tmpVec[1];
				// // corrWL += gradP[face.owner][2]*tmpVec[2];
				
				
				// linB0[face.owner] += muF*(UR-corrUL)/dPN_e*area;
				// linB1[face.owner] += muF*(VR-corrVL)/dPN_e*area;
				// linB2[face.owner] += muF*(WR-corrWL)/dPN_e*area;

				// // // viscous term, properties
				// // double delPhiECF;
				
				// // double gradUf[3];
				// // gradUf[0] = (gradU[face.owner][0]);
				// // gradUf[1] = (gradU[face.owner][1]);
				// // gradUf[2] = (gradU[face.owner][2]);
				// // // delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
				
				// // // gradUf[0] += nvec[0]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				// // // gradUf[1] += nvec[1]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				// // // gradUf[2] += nvec[2]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				
				// // double gradVf[3];
				// // gradVf[0] = (gradV[face.owner][0]);
				// // gradVf[1] = (gradV[face.owner][1]);
				// // gradVf[2] = (gradV[face.owner][2]);
				// // // delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
				
				// // // gradVf[0] += nvec[0]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				// // // gradVf[1] += nvec[1]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				// // // gradVf[2] += nvec[2]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				
				// // double gradWf[3];
				// // gradWf[0] = (gradW[face.owner][0]);
				// // gradWf[1] = (gradW[face.owner][1]);
				// // gradWf[2] = (gradW[face.owner][2]);
				// // // delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
				
				// // // gradWf[0] += nvec[0]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
				// // // gradWf[1] += nvec[1]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
				// // // gradWf[2] += nvec[2]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
						
						
				// // // viscous term, non-orthogonal
				// // linB0[face.owner] += muF*
					// // (gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
				// // linB1[face.owner] += muF*
					// // (gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
				// // linB2[face.owner] += muF*
					// // (gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
				
				// // double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);
			
				// // linA[face.owner] += muF/dPN*area*tmpCoeff2;
				
				
				
			// }
			
			
		// }
		
	// }
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar0(mesh.cells.size(),0.0);
	// vector<double> resiVar1(mesh.cells.size(),0.0);
	// vector<double> resiVar2(mesh.cells.size(),0.0);
	
	// // cout << "AAAAAAA" << endl;
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
	
		// // solveHYPRE(mesh, resiVar0, linA, linAL, linAR, linB0,
			// // controls.solverU, controls.toleranceU, 
			// // controls.relTolU, controls.preconditionerU,
			// // controls.maxIterU);
		// // solveHYPRE(mesh, resiVar1, linA, linAL, linAR, linB1,
			// // controls.solverU, controls.toleranceU, 
			// // controls.relTolU, controls.preconditionerU,
			// // controls.maxIterU);
		// // solveHYPRE(mesh, resiVar2, linA, linAL, linAR, linB2,
			// // controls.solverU, controls.toleranceU, 
			// // controls.relTolU, controls.preconditionerU,
			// // controls.maxIterU);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
	
		// // solveHYPRE(mesh, resiVar0, linA, linAL, linAR, linB0,
			// // controls.solverFinalU, controls.toleranceFinalU, 
			// // controls.relTolFinalU, controls.preconditionerFinalU,
			// // controls.maxIterFinalU);
		// // solveHYPRE(mesh, resiVar1, linA, linAL, linAR, linB1,
			// // controls.solverFinalU, controls.toleranceFinalU, 
			// // controls.relTolFinalU, controls.preconditionerFinalU,
			// // controls.maxIterFinalU);
		// // solveHYPRE(mesh, resiVar2, linA, linAL, linAR, linB2,
			// // controls.solverFinalU, controls.toleranceFinalU, 
			// // controls.relTolFinalU, controls.preconditionerFinalU,
			// // controls.maxIterFinalU);
	// }
	
	// // cout << "BBBBBBB" << endl;
		// // // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	// double test0 = 0.0;

	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// // resiVar0[i] = linB0[i]/linA[i];
		// // resiVar1[i] = linB1[i]/linA[i];
		// // resiVar2[i] = linB2[i]/linA[i];
		
		// cell.var[controls.U] += controls.momVelURF * resiVar0[i];
		// cell.var[controls.V] += controls.momVelURF * resiVar1[i];
		// cell.var[controls.W] += controls.momVelURF * resiVar2[i];
		
		
		// test0 += controls.momVelURF * resiVar0[i];
		// test0 += controls.momVelURF * resiVar1[i];
		// test0 += controls.momVelURF * resiVar2[i];
		
		// // cout << resiVar0[i] << endl;
		
		
		// // if( cell.var[controls.P] <= controls.minP ) 
			// // cell.var[controls.P] = controls.minP;
		// // if( cell.var[controls.P] >= controls.maxP ) 
			// // cell.var[controls.P] = controls.maxP;
		
		// if( cell.var[controls.U] <= controls.minU ) 
			// cell.var[controls.U] = controls.minU;
		// if( cell.var[controls.U] >= controls.maxU ) 
			// cell.var[controls.U] = controls.maxU;
		
		// if( cell.var[controls.V] <= controls.minV ) 
			// cell.var[controls.V] = controls.minV;
		// if( cell.var[controls.V] >= controls.maxV ) 
			// cell.var[controls.V] = controls.maxV;
		
		// if( cell.var[controls.W] <= controls.minW ) 
			// cell.var[controls.W] = controls.minW;
		// if( cell.var[controls.W] >= controls.maxW ) 
			// cell.var[controls.W] = controls.maxW;
		
		
		// // cell.var.at(controls.UDV[0]) += gradP[i][0];
		// // cell.var.at(controls.UDV[1]) += gradP[i][1];
		// // cell.var.at(controls.UDV[2]) += gradP[i][2];
		
		
		// residuals[1] += pow(resiVar0[i],2.0)*cell.volume;
		// residuals[2] += pow(resiVar1[i],2.0)*cell.volume;
		// residuals[3] += pow(resiVar2[i],2.0)*cell.volume;
		
		
		// linAD[i] = linA[i] / cell.volume;
		// linAOD[i] = 0.0;
		// for(auto& j : cell.faces){
			// linAOD[i] += linAR[j] / cell.volume;
		// }
		
		
		
		// gradP[i][0] /= cell.volume;
		// gradP[i][1] /= cell.volume;
		// gradP[i][2] /= cell.volume;
		
		
	// }
	
	
	// // cout << "A : " << test0 << endl;
	
	
	// //========================
	// // calc Un 
	// // gradient P
	// // vector<vector<double>> gradP;
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	// // vector<double> limGradP;
	// // math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // for(int j=0; j<3; ++j){
			// // gradP[i][j] *= limGradP[i]; 
		// // }
	// // }
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	

	// // vector<double> gradPx_recv;
	// // vector<double> gradPy_recv;
	// // vector<double> gradPz_recv;
	// vector<double> resiVar0_recv;
	// vector<double> resiVar1_recv;
	// vector<double> resiVar2_recv;
	// vector<double> linAD_recv;
	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradPx_send;
		// vector<double> gradPy_send;
		// vector<double> gradPz_send;
		// vector<double> resiVar0_send;
		// vector<double> resiVar1_send;
		// vector<double> resiVar2_send;
		// vector<double> linAD_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradPx_send.push_back(gradP[face.owner][0]);
				// gradPy_send.push_back(gradP[face.owner][1]);
				// gradPz_send.push_back(gradP[face.owner][2]);
				// resiVar0_send.push_back(resiVar0[face.owner]);
				// resiVar1_send.push_back(resiVar1[face.owner]);
				// resiVar2_send.push_back(resiVar2[face.owner]);
				// linAD_recv.push_back(linAD[face.owner]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// gradPx_send, gradPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPy_send, gradPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradPz_send, gradPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		
		// mpi.setProcsFaceDatas(
					// resiVar0_send, resiVar0_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// resiVar1_send, resiVar1_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// resiVar2_send, resiVar2_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
					
		// mpi.setProcsFaceDatas(
					// linAD_send, linAD_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
					
		// // gradPx_send.clear();
		// // gradPy_send.clear();
		// // gradPz_send.clear();
	// }
	
	
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];

		// linAFL[i] = linAL[i];
		// linAFR[i] = linAR[i];
		
		
		
		// double area = face.area;
		// vector<double> nvec(3,0.0);
		// nvec[0] = face.unitNormals[0];
		// nvec[1] = face.unitNormals[1];
		// nvec[2] = face.unitNormals[2];
		
		// vector<double> distanceCells(3,0.0);
		// distanceCells[0] = face.distCells[0];
		// distanceCells[1] = face.distCells[1];
		// distanceCells[2] = face.distCells[2];
		
		// double wCL = face.wC;
		// double wCR = 1.0-face.wC;
		
		// double RhoL = face.varL[controls.fRho];
		
		// double RhoR = face.varR[controls.fRho];
		
		// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		// double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  // pow(distanceCells[1],2.0) + 
						  // pow(distanceCells[2],2.0));
				
		// double nonOrtholimiter = 1.0;
		// double cosAlpha = dPN_e/dPN;
		// if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			// nonOrtholimiter = 0.5;
		// }
		// else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			// nonOrtholimiter = 0.333;
		// }
		// else if( cosAlpha < 0.342 ){
			// nonOrtholimiter = 0.0;
		// }
				
		// double Ef[3];
		// Ef[0] = distanceCells[0]/dPN;
		// Ef[1] = distanceCells[1]/dPN;
		// Ef[2] = distanceCells[2]/dPN;
				
		// double Tf[3];
		// Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		// Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		// Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		// double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		

		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				
			// face.var[controls.Un] += controls.momVelURF * (
				// (wCL*resiVar0[face.owner]+wCR*resiVar0[face.neighbour])*nvec[0] +
				// (wCL*resiVar1[face.owner]+wCR*resiVar1[face.neighbour])*nvec[1] +
				// (wCL*resiVar2[face.owner]+wCR*resiVar2[face.neighbour])*nvec[2] );
				
				
			// // // pressure correction
			
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = mesh.cells[face.neighbour].var[controls.P];
				
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			// // // face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// // // wCL*1.0/linAD[face.owner]*( 
					// // // gradP[face.owner][0]*nvec[0] +
					// // // gradP[face.owner][1]*nvec[1] +
					// // // gradP[face.owner][2]*nvec[2] ) +
				// // // wCR*1.0/linAD[face.neighbour]*( 
					// // // gradP[face.neighbour][0]*nvec[0] +
					// // // gradP[face.neighbour][1]*nvec[1] +
					// // // gradP[face.neighbour][2]*nvec[2] ) );
			// // face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// // wCL*controls.timeStep/RhoL*( 
					// // gradP[face.owner][0]*nvec[0] +
					// // gradP[face.owner][1]*nvec[1] +
					// // gradP[face.owner][2]*nvec[2] ) +
				// // wCR*controls.timeStep/RhoR*( 
					// // gradP[face.neighbour][0]*nvec[0] +
					// // gradP[face.neighbour][1]*nvec[1] +
					// // gradP[face.neighbour][2]*nvec[2] ) );
				
			// // face.var[controls.Un] -= coeffDiss * controls.momVelURF * tmp1*(orgPR-orgPL)/dPN_e;
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			// // face.var[controls.Un] -= coeffDiss * controls.momVelURF 
				// // * nonOrtholimiter * tmp1*
				// // (gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				
			// face.var[controls.Un] += controls.momVelURF * (
				// (wCL*resiVar0[face.owner]+wCR*resiVar0_recv[proc_num])*nvec[0] +
				// (wCL*resiVar1[face.owner]+wCR*resiVar1_recv[proc_num])*nvec[1] +
				// (wCL*resiVar2[face.owner]+wCR*resiVar2_recv[proc_num])*nvec[2] );
				
				
			// // // pressure correction
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = face.varR[controls.fP];
				
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			// // // face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// // // wCL*1.0/linAD[face.owner]*( 
					// // // gradP[face.owner][0]*nvec[0] +
					// // // gradP[face.owner][1]*nvec[1] +
					// // // gradP[face.owner][2]*nvec[2] ) +
				// // // wCR*1.0/linAD_recv[proc_num]*( 
					// // // gradPx_recv[proc_num]*nvec[0] +
					// // // gradPy_recv[proc_num]*nvec[1] +
					// // // gradPz_recv[proc_num]*nvec[2] ) );
			// // face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// // wCL*controls.timeStep/RhoL*( 
					// // gradP[face.owner][0]*nvec[0] +
					// // gradP[face.owner][1]*nvec[1] +
					// // gradP[face.owner][2]*nvec[2] ) +
				// // wCR*controls.timeStep/RhoR*( 
					// // gradPx_recv[proc_num]*nvec[0] +
					// // gradPy_recv[proc_num]*nvec[1] +
					// // gradPz_recv[proc_num]*nvec[2] ) );
				
				
			// // face.var[controls.Un] -= coeffDiss * controls.momVelURF * tmp1*(orgPR-orgPL)/dPN_e;
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			// // face.var[controls.Un] -= coeffDiss * controls.momVelURF 
				// // * nonOrtholimiter * tmp1*
				// // (gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
				
			// ++proc_num;
			
		// }
		

	// }
	
	
	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB0.clear();
	// linB1.clear();
	// linB2.clear();
	// resiVar0.clear();
	// resiVar1.clear();
	// resiVar2.clear();
	
	// return 0;
	
	
// }



