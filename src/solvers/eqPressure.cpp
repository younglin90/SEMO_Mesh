
#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include <ctime>

double SEMO_Solvers_Builder::calcPressureEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	int iterNonOrthogonality){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 1;
    int A_n = B_n * B_n;
	
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	bool boolLimiters = true;
	// bool boolLimiters = false;
	int gradIterMax = 1;
	int gradBoundaryIterMax = 0;
	int gradIterMax2 = 2;
	
	
	// linear solver
	
	vector<double> A_p_vals(mesh.cells.size(), 0.0);
	
	
	SEMO_Utility_Math math;
	
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
		
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumU, mesh.cellsProcVar[controls.maximumU]);
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumV, mesh.cellsProcVar[controls.maximumV]);
		mpi.sendRecvTemporaryCellData(mesh, controls.maximumW, mesh.cellsProcVar[controls.maximumW]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumU, mesh.cellsProcVar[controls.minimumU]);
		mpi.sendRecvTemporaryCellData(mesh, controls.minimumV, mesh.cellsProcVar[controls.minimumV]);
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



	vector<int> A_rows(mesh.cells.size()*A_n, 0);
	vector<int> A_cols(mesh.cells.size()*A_n, 0);
	vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
	vector<double> B(mesh.cells.size()*B_n, 0.0);
	

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];

		int id = 0;
		int i_glo = 0;
		int j_glo = 0;
		int i_loc = i;
		int str_glo = mesh.startCellGlobal*B_n;
		int step_loc = mesh.cells.size();
		int step_glo = mesh.ncellsTotal;
		
		double volume = cell.volume;
		double timeStep = controls.timeStep;
		
		
		
		id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
		A_rows[id] = i_glo; A_cols[id] = j_glo;
		
		
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
		double conv_flux = 0.0;
		double diff_flux = 0.0;
		
		
		//------------------------
		// continuity
		id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
		A_vals[id] += ( + alpha * tmp1 / dPN * area );
		A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		A_vals.push_back(( - alpha * tmp1 / dPN * area ));

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( - alpha * tmp1 / dPN * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( + alpha * tmp1 / dPN * area ));
		}
		
 
		// ----------------------------
		B[step_loc_L*0 + i_loc_L] -= ( UnF * area );
		
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			B[step_loc_R*0 + i_loc_R] += ( UnF * area );
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
				
				// int ijStartL_local = B_n*(face.owner) - 1;
				// vector<int> id(B_n,0);
				// id[0] = A_n*(face.owner) - 1;
				// for(int j=1; j<B_n; ++j){
					// id[j] = id[j-1] + 4;
				// }

				
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
			
				
				// pressure coefficient
				double coeffP = 0.0;
				double coeffP_diff = 0.0;
				if( boundary.type[controls.P] == "fixedValue" ){
					coeffP = 0.0;
					coeffP_diff = 1.0;
				}
				else if( boundary.type[controls.P] == "zeroGradient" ){
					coeffP = 1.0;
					coeffP_diff = 0.0;
				}
				else if( boundary.type[controls.P] == "switch" ){
					double machNum = 
						sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
							 pow(mesh.cells[face.owner].var[controls.V],2.0)+
							 pow(mesh.cells[face.owner].var[controls.W],2.0))/
							  mesh.cells[face.owner].var[controls.C];
					coeffP = (machNum > 1.0) ? 0.0 : 1.0;
					coeffP_diff = (machNum > 1.0) ? 1.0 : 0.0;
				}
				else {
					cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				// face->cell interpolation
				double dPN = 0.5 * face.magPN;
				double alpha = face.alphaF;
				
				double orgPL = mesh.cells[face.owner].var[controls.P];
				double orgPR = face.varR[controls.fP];
				
				if( boundary.type[controls.P] == "fixedValue" ){
					UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
					for(int ii=0; ii<3; ++ii){
						UnF += alpha * tmp1 * 
							mesh.cellsGradientVar[controls.P][face.owner][ii]*face.unitNomalsPN[ii];
					}
				}
				
				int str_glo = mesh.startCellGlobal*B_n;
				int step_loc = mesh.cells.size();
				int step_glo = mesh.ncellsTotal;
				int i_loc_L = face.owner;
				int id=-1;
				
				
				// continuity
				id = step_loc*(B_n*0+0) + i_loc_L;
				A_vals[id] += coeffP_diff * ( alpha * tmp1 / dPN * area );
				
				// convective term
				B[step_loc*0 + i_loc_L] -= ( UnF * area );
				
			}
			
		}
		
	}
	
	
	vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	solveAMGCL("pressure", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);

	
	
	vector<double> resiVar_recv;
	if(size>1){
		mpi.sendRecvTemporaryData(mesh, resiVar, resiVar_recv);
	}
	
	vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
	for(int iter=0; iter<gradIterMax2; ++iter){
		
		vector<vector<double>> gradResiP_recv;
		if(size>1) mpi.sendRecvTemporaryVectorData(mesh, gradResiP, gradResiP_recv);
		
		vector<vector<double>> gradResiP_tmp(mesh.cells.size(),vector<double>(3,0.0));
		
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
			
			double wCL = face.wC;
			double wCR = 1.0-wCL;
			// double wCL = 0.5;
			// double wCR = 0.5;
			
			double varF = 0.0;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				varF = wCL * resiVar[face.owner] + wCR * resiVar[face.neighbour];
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				varF = wCL * resiVar[face.owner] + wCR * resiVar_recv[ip];
			}
			
			// skewness correction
			if(boolSkewnessCorrection){
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					for(int ii=0; ii<3; ++ii){
						varF += wCL * gradResiP[face.owner][ii]*face.vecSkewness[ii];
						varF += wCR * gradResiP[face.neighbour][ii]*face.vecSkewness[ii];
					}
				}
				else{
					for(int ii=0; ii<3; ++ii){
						varF += wCL * gradResiP[face.owner][ii]*face.vecSkewness[ii];
						varF += wCR * gradResiP_recv[ip][ii]*face.vecSkewness[ii];
					}
				}
			}
			
			
			for(int j=0; j<3; ++j){
				gradResiP_tmp[face.owner][j] += varF*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
			}
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				for(int j=0; j<3; ++j){
					gradResiP_tmp[face.neighbour][j] -=  varF*face.unitNormals[j]*face.area / mesh.cells[face.neighbour].volume;
				}
			}
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
			
		}
		
		// boundary face's nodes
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					auto& cell = mesh.cells[face.owner];
					
					double varF = resiVar[face.owner];
					if(boundary.type[0] == "fixedValue"){
						varF = 0.0;//0.5*varF;
					}
					
					for(int j=0; j<3; ++j){
						gradResiP_tmp[face.owner][j] += varF*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
					}
				}
				
			}
		}
		
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<3; ++j){
				gradResiP[i][j] = gradResiP_tmp[i][j];
			}
		}
		
		
		
		

		// vector<vector<double>> gradResiP_tmp(mesh.cells.size(),vector<double>(3,0.0));
		
		// // math.calcGaussGreen(mesh, -1, resiVar, gradResiP_tmp);
		// int dummy0;
		// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
			// -1, dummy0, resiVar, gradResiP_tmp);
		
		// // boundary face's nodes
		// for(auto& boundary : mesh.boundary){
			// if(boundary.neighbProcNo == -1){
				// if(boundary.type[0] != "fixedValue") continue;
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					// double varF = -0.5*resiVar[face.owner];
					// for(int j=0; j<3; ++j){
						// gradResiP_tmp[face.owner][j] += varF*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
					// }
				// }
			// }
		// }
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// for(int j=0; j<3; ++j){
				// gradResiP[i][j] = gradResiP_tmp[i][j];
			// }
		// }
		
		
		
	}
	
	
	
	double maxP = -1.e9;
	vector<double> maxP_xyz(3,0.0);
	double maxUn = -1.e9;
	vector<double> maxUn_xyz(3,0.0);
	
	double normDelPrim = 0.0;
	double normPrim = 0.0;
	
	// update
	double normP = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		cell.var[controls.P] += controls.prePreURF * resiVar[i];
		
		double tmp1 = controls.timeStep/cell.var[controls.Rho];
		// tmp1 = 1.0/(A_p_vals[i] + cell.var[controls.Rho]/controls.timeStep);
		
		double resiU = controls.preVelURF * tmp1*gradResiP[i][0];
		double resiV = controls.preVelURF * tmp1*gradResiP[i][1];
		double resiW = controls.preVelURF * tmp1*gradResiP[i][2];
		
		cell.var[controls.U] -= resiU;
		cell.var[controls.V] -= resiV;
		cell.var[controls.W] -= resiW;
		

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
		
		normDelPrim += pow(resiVar[i],2.0);
		normDelPrim += pow(resiU,2.0);
		normDelPrim += pow(resiV,2.0);
		normDelPrim += pow(resiW,2.0);
		
		normPrim += pow(cell.var[controls.P],2.0);
		normPrim += pow(cell.var[controls.U],2.0);
		normPrim += pow(cell.var[controls.V],2.0);
		normPrim += pow(cell.var[controls.W],2.0);
	}
	


	
	
	
	vector<double> maxP_glo(size,0.0);
	vector<double> maxP_x_glo(size,0.0);
	vector<double> maxP_y_glo(size,0.0);
	vector<double> maxP_z_glo(size,0.0);
	
	vector<double> maxUn_glo(size,0.0);
	vector<double> maxUn_x_glo(size,0.0);
	vector<double> maxUn_y_glo(size,0.0);
	vector<double> maxUn_z_glo(size,0.0);
	
	MPI_Allgather(&maxP, 1, MPI_DOUBLE, maxP_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxP_xyz[0], 1, MPI_DOUBLE, maxP_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxP_xyz[1], 1, MPI_DOUBLE, maxP_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxP_xyz[2], 1, MPI_DOUBLE, maxP_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	MPI_Allgather(&maxUn, 1, MPI_DOUBLE, maxUn_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[0], 1, MPI_DOUBLE, maxUn_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[1], 1, MPI_DOUBLE, maxUn_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather(&maxUn_xyz[2], 1, MPI_DOUBLE, maxUn_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	

	
	if(rank==0){
		maxP = maxP_glo[0];
		maxP_xyz[0] = maxP_x_glo[0]; maxP_xyz[1] = maxP_y_glo[0]; maxP_xyz[2] = maxP_z_glo[0];
		for(int i=1; i<size; ++i){
			if(maxP<maxP_glo[i]){
				maxP = maxP_glo[i];
				maxP_xyz[0] = maxP_x_glo[i]; maxP_xyz[1] = maxP_y_glo[i]; maxP_xyz[2] = maxP_z_glo[i];
			}
		}
		maxUn = maxUn_glo[0];
		maxUn_xyz[0] = maxUn_x_glo[0]; maxUn_xyz[1] = maxUn_y_glo[0]; maxUn_xyz[2] = maxUn_z_glo[0];
		for(int i=1; i<size; ++i){
			if(maxUn<maxUn_glo[i]){
				maxUn = maxUn_glo[i];
				maxUn_xyz[0] = maxUn_x_glo[i]; maxUn_xyz[1] = maxUn_y_glo[i]; maxUn_xyz[2] = maxUn_z_glo[i];
			}
		}
		
		cout << "-----------------" << endl;
		cout << " maximum P = " << maxP << ", xyz = " << maxP_xyz[0] << ", " << maxP_xyz[1] << ", " << maxP_xyz[2] << endl;
		cout << " maximum Un = " << maxUn << ", xyz = " << maxUn_xyz[0] << ", " << maxUn_xyz[1] << ", " << maxUn_xyz[2] << endl;
		cout << "-----------------" << endl;
	}
	

	
	
	return sqrt(normDelPrim)/sqrt(normPrim);
	
	
}













double SEMO_Solvers_Builder::calcPressureEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<double>& linAD,
	vector<double>& linAOD,
	vector<double>& linAFL,
	vector<double>& linAFR,
	vector<double>& residuals){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	// int proc_num=0;
	
	// vector<clock_t> startTime;
	
	// SEMO_Utility_Math math;

	// // gradient P
	// vector<vector<double>> gradP;
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
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


	// // diagonal terms
	// vector<double> linA(mesh.cells.size(),0.0);
	// vector<double> linB(mesh.cells.size(),0.0);

	// // diagonal, off-diagonal terms
	// vector<double> linAL(mesh.faces.size(),0.0);
	// vector<double> linAR(mesh.faces.size(),0.0);
	
	// proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
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
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double diff_flux = tmp1/dPN_e * area;
		
		// // own, ngb
		// linAR[i] = (-diff_flux);
		// // ngb, own
		// linAL[i] = (-diff_flux);
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // own, ngb
			// linA[face.owner] += (+diff_flux);
			// // ngb, own
			// linA[face.neighbour] -= (-diff_flux);
			
			// // convective term
			// linB[face.owner] -= UnF*area;
			// linB[face.neighbour] += UnF*area; 
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // own, ngb
			// linA[face.owner] += (+diff_flux);
			
			// // convective term
			// linB[face.owner] -= UnF*area;
			
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
								  
				// dPN = dPN_e;
							
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
			
				
				
				// // pressure coefficient
				// double coeffP = 0.0;
				// double coeffP_diff = 0.0;
				// if( boundary.type[controls.P] == "fixedValue" ){
					// coeffP = 0.0;
					// coeffP_diff = 1.0;
				// }
				// else if( boundary.type[controls.P] == "zeroGradient" ){
					// coeffP = 1.0;
					// coeffP_diff = 0.0;
					
				// }
				// else if( boundary.type[controls.P] == "switch" ){
					// double machNum = 
						// sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
							 // pow(mesh.cells[face.owner].var[controls.V],2.0)+
							 // pow(mesh.cells[face.owner].var[controls.W],2.0))/
							  // mesh.cells[face.owner].var[controls.C];
					// coeffP = (machNum > 1.0) ? 0.0 : 1.0;
					// coeffP_diff = (machNum > 1.0) ? 1.0 : 0.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				// // continuity
				// linA[face.owner] += coeffP_diff * ( tmp1 / dPN * area );
				
				// // convective term
				// linB[face.owner] -= ( UnF * area );
			// }
			
		// }
		
	// }
	
	
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar(mesh.cells.size(),0.0);
	
	// string solver;
	// double relTol;
	// double tolerance;
	// string preconditioner;
	// int maxIter;
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
		// solver = controls.solverP;
		// relTol = controls.toleranceP;
		// tolerance = controls.relTolP;
		// preconditioner = controls.preconditionerP;
		// maxIter = controls.maxIterP;
	
	// }
	// else{
		// solver = controls.solverFinalP;
		// relTol = controls.toleranceFinalP;
		// tolerance = controls.relTolFinalP;
		// preconditioner = controls.preconditionerFinalP;
		// maxIter = controls.maxIterFinalP;
	
	// }
	
	// // solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
		// // solver, tolerance, relTol, preconditioner, maxIter);
		
	// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
		// solver, tolerance, relTol, preconditioner, maxIter);
		
	
	
	

	
	// vector<vector<double>> gradResiP;
	// // math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
	// // math.calcLeastSquare(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
	
	
	
	// double maxP = -1.e9;
	// vector<double> maxP_xyz(3,0.0);
	// double maxUn = -1.e9;
	// vector<double> maxUn_xyz(3,0.0);
	
	// double normDelPrim = 0.0;
	// double normPrim = 0.0;
	
	// // update
	// double normP = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// cell.var[controls.P] += controls.prePreURF*resiVar[i];
		
		// double tmp1 = controls.timeStep/cell.var[controls.Rho];
		
		// double resiU = controls.preVelURF*tmp1*gradResiP[i][0];
		// double resiV = controls.preVelURF*tmp1*gradResiP[i][1];
		// double resiW = controls.preVelURF*tmp1*gradResiP[i][2];
		
		// cell.var[controls.U] -= resiU;
		// cell.var[controls.V] -= resiV;
		// cell.var[controls.W] -= resiW;
		
		// if( cell.var[controls.P] <= controls.minP ) 
			// cell.var[controls.P] = controls.minP;
		// if( cell.var[controls.P] >= controls.maxP ) 
			// cell.var[controls.P] = controls.maxP;
		
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
		
		// residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		// residuals[1] += pow(resiU,2.0)*cell.volume;
		// residuals[2] += pow(resiV,2.0)*cell.volume;
		// residuals[3] += pow(resiW,2.0)*cell.volume;
		
		// normDelPrim += pow(resiVar[i],2.0);
		// normDelPrim += pow(resiU,2.0);
		// normDelPrim += pow(resiV,2.0);
		// normDelPrim += pow(resiW,2.0);
		
		// normPrim += pow(cell.var[controls.P],2.0);
		// normPrim += pow(cell.var[controls.U],2.0);
		// normPrim += pow(cell.var[controls.V],2.0);
		// normPrim += pow(cell.var[controls.W],2.0);
	// }
	


	
	
	
	// vector<double> maxP_glo(size,0.0);
	// vector<double> maxP_x_glo(size,0.0);
	// vector<double> maxP_y_glo(size,0.0);
	// vector<double> maxP_z_glo(size,0.0);
	
	// vector<double> maxUn_glo(size,0.0);
	// vector<double> maxUn_x_glo(size,0.0);
	// vector<double> maxUn_y_glo(size,0.0);
	// vector<double> maxUn_z_glo(size,0.0);
	
	// MPI_Allgather(&maxP, 1, MPI_DOUBLE, maxP_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[0], 1, MPI_DOUBLE, maxP_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[1], 1, MPI_DOUBLE, maxP_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[2], 1, MPI_DOUBLE, maxP_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// MPI_Allgather(&maxUn, 1, MPI_DOUBLE, maxUn_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[0], 1, MPI_DOUBLE, maxUn_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[1], 1, MPI_DOUBLE, maxUn_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[2], 1, MPI_DOUBLE, maxUn_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	

	
	// if(rank==0){
		// maxP = maxP_glo[0];
		// maxP_xyz[0] = maxP_x_glo[0]; maxP_xyz[1] = maxP_y_glo[0]; maxP_xyz[2] = maxP_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(maxP<maxP_glo[i]){
				// maxP = maxP_glo[i];
				// maxP_xyz[0] = maxP_x_glo[i]; maxP_xyz[1] = maxP_y_glo[i]; maxP_xyz[2] = maxP_z_glo[i];
			// }
		// }
		// maxUn = maxUn_glo[0];
		// maxUn_xyz[0] = maxUn_x_glo[0]; maxUn_xyz[1] = maxUn_y_glo[0]; maxUn_xyz[2] = maxUn_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(maxUn<maxUn_glo[i]){
				// maxUn = maxUn_glo[i];
				// maxUn_xyz[0] = maxUn_x_glo[i]; maxUn_xyz[1] = maxUn_y_glo[i]; maxUn_xyz[2] = maxUn_z_glo[i];
			// }
		// }
		
		// cout << "-----------------" << endl;
		// cout << " maximum P = " << maxP << ", xyz = " << maxP_xyz[0] << ", " << maxP_xyz[1] << ", " << maxP_xyz[2] << endl;
		// cout << " maximum Un = " << maxUn << ", xyz = " << maxUn_xyz[0] << ", " << maxUn_xyz[1] << ", " << maxUn_xyz[2] << endl;
		// cout << "-----------------" << endl;
	// }
	

	
	
	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	// return sqrt(normP);
	
	
}














































// #include "build.h"
// #include <cmath>
// #include <array>
// #include <numeric>
// #include <ctime>




// double SEMO_Solvers_Builder::calcPressureEq(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& linAFL,
	// vector<double>& linAFR,
	// vector<double>& residuals){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi; 
	
	// vector<clock_t> startTime;
	
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
	
	
	
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	
	
	
	

	
	// // //=========================================
	// // vector<double> maxXPF(mesh.cells.size(),0.0);
	// // vector<double> maxS(mesh.cells.size(),0.0);	// internal faces
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double xPF = sqrt(
				// // pow(face.x - mesh.cells[face.owner].x,2.0)+
				// // pow(face.y - mesh.cells[face.owner].y,2.0)+
				// // pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// // maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// // double xNF = sqrt(
				// // pow(face.x - mesh.cells[face.neighbour].x,2.0)+
				// // pow(face.y - mesh.cells[face.neighbour].y,2.0)+
				// // pow(face.z - mesh.cells[face.neighbour].z,2.0));
				
			// // maxXPF[face.neighbour] = max(maxXPF[face.neighbour],xNF);
			
			// // maxS[face.owner] = max(maxS[face.owner],face.area);
			// // maxS[face.neighbour] = max(maxS[face.neighbour],face.area);
			
			
		// // }
		// // else{
			// // double xPF = sqrt(
				// // pow(face.x - mesh.cells[face.owner].x,2.0)+
				// // pow(face.y - mesh.cells[face.owner].y,2.0)+
				// // pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// // maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// // maxS[face.owner] = max(maxS[face.owner],face.area);
			
		// // }
	// // }
	// // //=========================================
	
	
	// // //=========================================
	// // vector<vector<double>> gradPLSQ;
	// // vector<vector<double>> gradPGG;
	// // math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradPLSQ);
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradPGG);
	// // vector<vector<double>> gradP;
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		// // double beta0 = min(1.0, 2.0/AR);
		
		// // vector<double> tmpVec;
		// // tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][0]);
		// // tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][1]);
		// // tmpVec.push_back(beta0*gradPLSQ[i][0]+(1.0-beta0)*gradPGG[i][2]);
		
		// // gradP.push_back(tmpVec);
		
	// // }
	// // //=========================================
	
	
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
	




	// // vector<vector<double>> gradU;
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // vector<vector<double>> gradV;
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // vector<vector<double>> gradW;
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);


	

	// // // calc gauss-green gradient
	// // vector<vector<double>> gradGrav(mesh.cells.size(),vector<double>(3,0.0));
	
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // // double RhoF = 2.0/(1.0/face.varL[controls.fRho]+1.0/face.varR[controls.fRho]);
			// // double RhoF = 0.5*(face.varL[controls.fRho]+face.varR[controls.fRho]);
			// // gradGrav[face.owner][1] += 
				// // 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.owner].volume;
			// // gradGrav[face.neighbour][1] -= 
				// // 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.neighbour].volume;
		// // }
		// // else{
			// // // double RhoF = 2.0/(1.0/face.varL[controls.fRho]+1.0/face.varR[controls.fRho]);
			// // double RhoF = 0.5*(face.varL[controls.fRho]+face.varR[controls.fRho]);
			// // gradGrav[face.owner][1] += 
				// // 0.5*RhoF*(-9.8)*face.unitNormals[1]*face.distCells[1]*face.area/mesh.cells[face.owner].volume;
		// // }
	// // }




	// vector<double> cellVolume_recv;
	// if(size>1){
		// vector<double> cellVolume_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// cellVolume_send.push_back(mesh.cells[face.owner].volume);
			// }
		// }
		// mpi.setProcsFaceDatas(
					// cellVolume_send, cellVolume_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// cellVolume_send.clear();
		
	// }





	// // linAD
	// vector<double> linAD_send;
	// vector<double> linAD_recv;
	// if(size>1){
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// linAD_send.push_back(linAD[face.owner]);
			// }
		// }
		
		// mpi.setProcsFaceDatas(
					// linAD_send, linAD_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// linAD_send.clear();
		
	// }







	// vector<double> linA;
	// vector<double> linB;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		
		// linA.push_back(0.0);
		// linB.push_back(0.0);
			
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
			// // nonOrtholimiter = 0.333;
		// }
		// else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			// nonOrtholimiter = 0.333;
			// // nonOrtholimiter = 0.1;
			// // nonOrtholimiter = 0.4;
		// }
		// else if( cosAlpha < 0.342 ){
			// nonOrtholimiter = 0.0;
			// // nonOrtholimiter = 0.3;
		// }
		// // nonOrtholimiter = 0.333;
		
		// // cout << cosAlpha << " " << nonOrtholimiter << endl;
				
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
		
		// // double coeffP = area/dPNEff * dPN/dPNEff;
		
		
		
		
		
		
		
		
		
		
		
		
		
		// double UnF = wCL*UnL+wCR*UnR;
		
		
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
		
		// // double RhoF = RhoL*weightL + RhoR*weightR;
		// double RhoF = 1.0;//RhoL*weightL + RhoR*weightR;
		
		// // linAL.push_back(-tmp1/dPN * area*overRelaxCoeff);
		// // linAR.push_back(-tmp1/dPN * area*overRelaxCoeff);
		
		// linAL.push_back(-tmp1/dPN_e * area);
		// linAR.push_back(-tmp1/dPN_e * area);
		
		// // linAL.back() += (-tmp1/dPN * area*overRelaxCoeff * (nvec[0]*Tf[0]+nvec[1]*Tf[1]+nvec[2]*Tf[2]));
		// // linAR.back() += (-tmp1/dPN * area*overRelaxCoeff * (nvec[0]*Tf[0]+nvec[1]*Tf[1]+nvec[2]*Tf[2]));
		
		// // linAL.push_back(-0.5*tmp1*coeffP);
		// // linAR.push_back(-0.5*tmp1*coeffP);
			
		// // convective term
		// // linB[face.owner] -= RhoF*UnF*area;
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA[face.owner]     -= linAL[i];
			// linA[face.neighbour] -= linAR[i];
			
			// // convective term
			// linB[face.owner] -= RhoF*UnF*area;
			// linB[face.neighbour] += RhoF*UnF*area; 
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA[face.owner] -= linAL[i];
			
			// // convective term
			// linB[face.owner] -= RhoF*UnF*area;
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
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
					
					// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
					// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
					
					// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
					// double dPN = sqrt(pow(distanceCells[0],2.0) + 
									  // pow(distanceCells[1],2.0) + 
									  // pow(distanceCells[2],2.0));
								
					// double UnF = 0.5*UnL+0.5*UnR;
				
					// double RhoL = face.varL[controls.fRho];
					// double tmp1 = controls.timeStep/RhoL;
					// // double tmp1 = 1.0/linAD[face.owner];
							
					// linA[face.owner] += tmp1*area/dPN_e;
					
					// double RhoF = 1.0;
					
					// // // convective term
					// // double orgPL = mesh.cells[face.owner].var[controls.P];
					// // double PR = face.varR[controls.fP];
					
					// // // UnF = 
						// // // ( 
						// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
						// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
						// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] );
					// // UnF = 
						// // ( 
						// // (UL + controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
						// // (VL + controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
						// // (WL + controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] );
						
					// // double corrPL = orgPL;
					// // UnF = UnF - tmp1*(PR-corrPL)/dPN_e;
					
					// linB[face.owner] -= RhoF*UnF*area;
				// }
			// }
			// else{
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
					
					// double UL = face.varL[controls.fU];
					// double VL = face.varL[controls.fV];
					// double WL = face.varL[controls.fW];
					
					// double UR = face.varR[controls.fU];
					// double VR = face.varR[controls.fV];
					// double WR = face.varR[controls.fW];
					
					// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
					// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
					
					// double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
					
					// double UnF = 0.5*UnL+0.5*UnR;
				
					// double RhoL = face.varL[controls.fRho];
					// // double tmp1 = controls.timeStep/RhoL;
					// double RhoF = 1.0;
					
					// // // convective term
					// // double orgPL = mesh.cells[face.owner].var[controls.P];
					// // double PR = orgPL 
						// // + gradP[face.owner][0]*distanceCells[0]
						// // + gradP[face.owner][1]*distanceCells[1]
						// // + gradP[face.owner][2]*distanceCells[2];
					
					// // // UnF = 
						// // // ( 
						// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
						// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
						// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] );
					// // UnF = 
						// // ( 
						// // (UL + controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
						// // (VL + controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
						// // (WL + controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] );
						
					// // double corrPL = orgPL;
					// // UnF = UnF - tmp1*(PR-corrPL)/dPN_e;
					
					// // convective term
					// linB[face.owner] -= RhoF*UnF*area;
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
			
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverP, controls.toleranceP, 
			// // controls.relTolP, controls.preconditionerP,
			// // controls.maxIterP);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalP, controls.toleranceFinalP, 
			// controls.relTolFinalP, controls.preconditionerFinalP,
			// controls.maxIterFinalP);
			
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverFinalP, controls.toleranceFinalP, 
			// // controls.relTolFinalP, controls.preconditionerFinalP,
			// // controls.maxIterFinalP);
	// }
	

	

	
	// vector<vector<double>> gradResiP;
	// // math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
	// // math.calcLeastSquare(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
	
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradResiP);
	
	// // //==========================================
	// // vector<vector<double>> gradResiPLSQ;
	// // vector<vector<double>> gradResiPGG;
	// // math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiPLSQ);
	// // math.calcGaussGreen(mesh, 0, resiVar, gradResiPGG);
	// // vector<vector<double>> gradResiP;
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		// // double beta0 = min(1.0, 2.0/AR);
		
		// // vector<double> tmpVec;
		// // tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][0]);
		// // tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][1]);
		// // tmpVec.push_back(beta0*gradResiPLSQ[i][0]+(1.0-beta0)*gradResiPGG[i][2]);
		
		// // gradResiP.push_back(tmpVec);
		
	// // }
	// // //==========================================
	
	
	// // update
	// double test0 = 0.0;
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// cell.var[controls.P] += controls.prePreURF*resiVar[i];
		
		// // double resiU = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][0];
		// // double resiV = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][1];
		// // double resiW = controls.preVelURF/(linAD[i]+linAOD[i])*gradResiP[i][2];
		// // double resiU = controls.preVelURF/linAD[i]*gradResiP[i][0];
		// // double resiV = controls.preVelURF/linAD[i]*gradResiP[i][1];
		// // double resiW = controls.preVelURF/linAD[i]*gradResiP[i][2];
		// double resiU = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0];
		// double resiV = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1];
		// double resiW = controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2];
		
		// cell.var[controls.U] -= resiU;
		// cell.var[controls.V] -= resiV;
		// cell.var[controls.W] -= resiW;
		
		

		// if( cell.var[controls.P] <= controls.minP ) 
			// cell.var[controls.P] = controls.minP;
		// if( cell.var[controls.P] >= controls.maxP ) 
			// cell.var[controls.P] = controls.maxP;
		
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
		
		// test0 += controls.prePreURF*resiVar[i];
		// // test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][0];
		// // test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][1];
		// // test0 += controls.preVelURF*controls.timeStep/cell.var[controls.Rho]*gradResiP[i][2];
		
		
		// // cell.var[controls.UDV[0]] = resiVar[i];
		// // cell.var[controls.UDV[1]] = gradResiP[i][1];
		
		
		// // cell.var.at(controls.UDV[0]) += gradP[i][0];
		// // cell.var.at(controls.UDV[1]) += gradP[i][1];
		// // cell.var.at(controls.UDV[2]) += gradP[i][2];
		
		// // cell.var[controls.U] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][0];
		// // cell.var[controls.V] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][1];
		// // cell.var[controls.W] -= controls.timeStep/cell.var[controls.Rho]*gradP[i][2];
		
		
				
		// residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		// residuals[1] += pow(resiU,2.0)*cell.volume;
		// residuals[2] += pow(resiV,2.0)*cell.volume;
		// residuals[3] += pow(resiW,2.0)*cell.volume;
		
	// }
	
	
	
	// // cout << "A : " << test0 << endl;
	
	


	// vector<double> gradResiPx_recv;
	// vector<double> gradResiPy_recv;
	// vector<double> gradResiPz_recv;
	// vector<double> resiVar_recv;
	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradResiPx_send;
		// vector<double> gradResiPy_send;
		// vector<double> gradResiPz_send;
		// vector<double> resiVar_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradResiPx_send.push_back(gradResiP[face.owner][0]);
				// gradResiPy_send.push_back(gradResiP[face.owner][1]);
				// gradResiPz_send.push_back(gradResiP[face.owner][2]);
				// resiVar_send.push_back(resiVar[face.owner]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// gradResiPx_send, gradResiPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradResiPy_send, gradResiPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradResiPz_send, gradResiPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		
		// mpi.setProcsFaceDatas(
					// resiVar_send, resiVar_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
					
		// // gradPx_send.clear();
		// // gradPy_send.clear();
		// // gradPz_send.clear();
	// }
	
	
	
	// proc_num = 0;
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
				
		// double Tf[3];
		// Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		// Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		// Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		// double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;

		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // tmp1 = (wCL/linAD[face.owner] + wCR/linAD[face.neighbour]);
			
			// face.var[controls.Un] -= 
				// controls.preVelURF * tmp1*(resiVar[face.neighbour]-resiVar[face.owner])/dPN_e;
				
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradResiP[face.owner][0]+wCR*gradResiP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradResiP[face.owner][1]+wCR*gradResiP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradResiP[face.owner][2]+wCR*gradResiP[face.neighbour][2]);
			
			// // face.var[controls.Un] -= controls.preVelURF * 
				// // nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // tmp1 = (wCL/linAD[face.owner] + wCR/linAD_recv[proc_num]);
			
			// face.var[controls.Un] -= 
				// controls.preVelURF * tmp1*(resiVar_recv[proc_num]-resiVar[face.owner])/dPN_e;
				
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradResiP[face.owner][0]+wCR*gradResiPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradResiP[face.owner][1]+wCR*gradResiPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradResiP[face.owner][2]+wCR*gradResiPz_recv[proc_num]);
			
			// // face.var[controls.Un] -= controls.preVelURF * 
				// // nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			// ++proc_num;
			
		// }
		
	// }
	
	
	
	
	
	
	
	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	// return 0;
	
	
// }








