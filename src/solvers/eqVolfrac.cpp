#include "build.h"
#include <cmath>
#include <array>
#include <numeric>




double SEMO_Solvers_Builder::calcVolfracEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 1;
    int A_n = B_n * B_n;
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	int gradIterMax = 1;
	
	SEMO_Utility_Math math;
	
	
	setCellVarMinMax(mesh, controls.U, controls.maximumU, controls.minimumU);
	setCellVarMinMax(mesh, controls.V, controls.maximumV, controls.minimumV);
	setCellVarMinMax(mesh, controls.W, controls.maximumW, controls.minimumW);
	
	
	
	// vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradU(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradV(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradW(mesh.cells.size(),vector<double>(3,0.0));
	
	// mesh.cellsGradientVar[controls.P].clear();
	// mesh.cellsGradientVar[controls.U].clear();
	// mesh.cellsGradientVar[controls.V].clear();
	// mesh.cellsGradientVar[controls.W].clear();
	// mesh.cellsGradientVar[controls.T].clear();
	// mesh.cellsGradientVar[controls.VF[0]].clear();
	
	mesh.cellsGradientVar[controls.P].resize(mesh.cells.size(),vector<double>(9,0.0));
	mesh.cellsGradientVar[controls.U].resize(mesh.cells.size(),vector<double>(9,0.0));
	mesh.cellsGradientVar[controls.V].resize(mesh.cells.size(),vector<double>(9,0.0));
	mesh.cellsGradientVar[controls.W].resize(mesh.cells.size(),vector<double>(9,0.0));
	mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(9,0.0));
	
	for(int iter=0; iter<gradIterMax; ++iter){
		
		// math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		// math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		// math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		
		// math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		// math.calcGaussGreen(mesh, controls.U, controls.fU, mesh.cellsGradientVar[controls.U]);
		// math.calcGaussGreen(mesh, controls.V, controls.fV, mesh.cellsGradientVar[controls.V]);
		// math.calcGaussGreen(mesh, controls.W, controls.fW, mesh.cellsGradientVar[controls.W]);
		// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], mesh.cellsGradientVar[controls.VF[0]]);
			
		vector<double> dummyVec;
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.P, controls.fP, dummyVec, mesh.cellsGradientVar[controls.P]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.U, controls.fU, dummyVec, mesh.cellsGradientVar[controls.U]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.V, controls.fV, dummyVec, mesh.cellsGradientVar[controls.V]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.W, controls.fW, dummyVec, mesh.cellsGradientVar[controls.W]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.VF[0], controls.fVF[0], dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
	}
	
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

        int ijStart_local = B_n*(i) - 1;
        int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
        int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		int id = 0;
		int i_glo = 0;
		int j_glo = 0;
		int i_loc = i;
		int str_glo = mesh.startCellGlobal*B_n;
		int step_loc = mesh.cells.size();
		int step_glo = mesh.ncellsTotal;
		
		double volume = cell.volume;
		double timeStep = controls.timeStep;
		
		double VF = cell.var[controls.VF[0]];
		double oldVF = cell.var[controls.oldVF[0]];
		double old2VF = cell.var[controls.Qm[3]];
		

		double beta_phi = 1.0;
		double phi_SOU = ( oldVF - old2VF ) /
					(abs(VF-old2VF)+1.e-200)*( VF>old2VF ? 1.0 : -1.0 );
		if(0.0 < phi_SOU && phi_SOU > 1.0){
			beta_phi = 0.0;
		}
		else{
			beta_phi = 1.0;
		}
		double coeff_old1 = controls.oldTimeStep/(controls.timeStep+controls.oldTimeStep);
		double coeff_old2 = controls.old2TimeStep/(controls.oldTimeStep+controls.old2TimeStep);
		double cons_np12 = VF + beta_phi * coeff_old1 * (VF-oldVF);
		double cons_nm12 = oldVF + beta_phi * coeff_old2 * (oldVF-old2VF);
		
		// volume fraction
		id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
		A_rows[id] = i_glo; A_cols[id] = j_glo;
		
		A_vals[id] = volume/timeStep;
		// A_vals[id] = 1.5*volume/timeStep;
		// A_vals[id] = (1.0+beta_phi*coeff_old1)*volume/timeStep;
		
		B[step_loc*0 + i_loc] = -(VF - oldVF)*volume/timeStep;
		// B[step_loc*0 + i_loc] = -(1.5*VF - 2.0*oldVF + 0.5*old2VF)*volume/timeStep;
		// B[step_loc*0 + i_loc] = -(cons_np12 - cons_nm12)*volume/timeStep;
		
		
		
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
	
		double VFL_HO = face.varL[controls.fVF[0]];
		double VFR_HO = face.varR[controls.fVF[0]];
		
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
		
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		double UnF = wCL*UnL+wCR*UnR;
		
		
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
				
				// UL_skew += gradUL[ii]*face.vecPF[ii];
				// UR_skew += gradUR[ii]*face.vecNF[ii];
				// VL_skew += gradVL[ii]*face.vecPF[ii];
				// VR_skew += gradVR[ii]*face.vecNF[ii];
				// WL_skew += gradWL[ii]*face.vecPF[ii];
				// WR_skew += gradWR[ii]*face.vecNF[ii];
			}
			UL_skew = max(minimumUL,min(maximumUL,UL_skew));
			UR_skew = max(minimumUR,min(maximumUR,UR_skew));
			VL_skew = max(minimumVL,min(maximumVL,VL_skew));
			VR_skew = max(minimumVR,min(maximumVR,VR_skew));
			WL_skew = max(minimumWL,min(maximumWL,WL_skew));
			WR_skew = max(minimumWR,min(maximumWR,WR_skew));
			
			UnL = UL_skew*nvec[0] + VL_skew*nvec[1] + WL_skew*nvec[2];
			UnR = UR_skew*nvec[0] + VR_skew*nvec[1] + WR_skew*nvec[2];
			UnF = wCL*UnL+wCR*UnR;
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
		// calcInterpolVelSurfTens(UnF, srcSFT_L, srcSFT_R,
			// RhoL, RhoR, wCL, wCR, 
			// controls.timeStep,
			// species[0].sigma, kappaL, kappaR,
			// VFL, VFR,
			// alpha, dPN,
			// nvec,
			// face.unitNomalsPN);
		// calcInterpolVelUnsteady(UnF, 
			// wCL, wCR, 
			// face.var[controls.Un],
			// UL_old, UR_old,
			// VL_old, VR_old,
			// WL_old, WR_old,
			// nvec);
		
		
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = 1.0 - weightL;
		
		double VFF_HO = weightL*VFL_HO + weightR*VFR_HO;
		
		// if(boolSkewnessCorrection){
			// for(int ii=0; ii<3; ++ii){
				// VFF_HO += weightL * gradVFL[ii]*face.vecSkewness[ii];
				// VFF_HO += weightR * gradVFR[ii]*face.vecSkewness[ii];
			// }
			// VFF_HO = max(0.0,min(1.0,VFF_HO));
		// }
		
		
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
        // VF'
		id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
		A_vals[id] += ( weightL * UnF * area );
		A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		A_vals.push_back(( weightR * UnF * area ));
		
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( weightR * UnF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( weightL * UnF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
        // ----------------------------
        B[step_loc_L*0 + i_loc_L] -= ( VFF_HO * UnF * area );

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			B[step_loc_R*0 + i_loc_R] += ( VFF_HO * UnF * area );
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
				
				double area = face.area;
		
				double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
							 face.varL[controls.fV]*face.unitNormals[1] + 
							 face.varL[controls.fW]*face.unitNormals[2];
				double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
							 face.varR[controls.fV]*face.unitNormals[1] + 
							 face.varR[controls.fW]*face.unitNormals[2];
							
				double UnF = 0.5*UnL+0.5*UnR;
				
				double RhoL = face.varL[controls.fRho];
				double tmp1 = 1.0/RhoL*controls.timeStep;
				
				// for(int ii=0; ii<3; ++ii){
					// UnF += tmp1 * gradP[face.owner][ii]*face.unitNormals[ii];
				// }
				
				double VFF = 0.5*face.varL[controls.fVF[0]] + 0.5*face.varR[controls.fVF[0]];
			
				// volume-fraction coefficient
				double coeffVF = 0.0;
				if( boundary.type[controls.VF[0]] == "fixedValue" ){
					coeffVF = 0.0;
				}
				else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
					coeffVF = 1.0;
					
				}
				else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
					double ownNorVel =  
						mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					coeffVF = (ownNorVel > 1.0) ? 1.0 : 0.0;
				}
				else {
					cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				int str_glo = mesh.startCellGlobal*B_n;
				int step_loc = mesh.cells.size();
				int step_glo = mesh.ncellsTotal;
				int i_loc_L = face.owner;
				int id=-1;
				
				// continuity
				id = step_loc*(B_n*0+0) + i_loc_L;
				A_vals[id] += coeffVF * ( UnF * area );
				
				// convective term
				B[step_loc*0 + i_loc_L] -= ( VFF * UnF * area );
			}
			
			
		}
		
	}
	
	
	vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	solveAMGCL("volume_fraction", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);
	
	
	double normDelVF = 0.0;
	double normVF = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		// resiVar[i] = B[i]/A_vals[i];
		
		cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
		normDelVF += pow(resiVar[i],2.0);
		
		normVF += pow(cell.var[controls.VF[0]],2.0);
		
	}
	
	return sqrt(normDelVF);// /(sqrt(normVF+1.d-50));
	
}








double SEMO_Solvers_Builder::calcVolfracEq(
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
	
	
	// gradient P
	vector<vector<double>> gradP;
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
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
	
	
	vector<double> gradPx_recv;
	vector<double> gradPy_recv;
	vector<double> gradPz_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> gradPx_send;
		vector<double> gradPy_send;
		vector<double> gradPz_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				gradPx_send.push_back(gradP[face.owner][0]);
				gradPy_send.push_back(gradP[face.owner][1]);
				gradPz_send.push_back(gradP[face.owner][2]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					gradPx_send, gradPx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradPy_send, gradPy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradPz_send, gradPz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradPx_send.clear();
		gradPy_send.clear();
		gradPz_send.clear();
	}
	





	vector<double> linA(mesh.cells.size(),0.0);
	vector<double> linAL(mesh.faces.size(),0.0);
	vector<double> linAR(mesh.faces.size(),0.0);
	vector<double> linB(mesh.cells.size(),0.0);

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		
		// time marching -> first order euler
		linA[i] = cell.volume/controls.timeStep;
		
		linB[i] = -(
			1.0*cell.var[controls.VF[0]] - 1.0*cell.var[controls.oldVF[0]]
			)*cell.volume/controls.timeStep;
			
		// time marching -> second order upwind euler
		// linA[i] = 1.5*cell.volume/controls.timeStep;
		
		// linB[i] = -(
			// 1.5*cell.var[controls.VF[0]] - 2.0*cell.var[controls.oldVF[0]] + 0.5*cell.var[controls.Qm[3]]
			// )*cell.volume/controls.timeStep;
		
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
		
		double Ef[3];
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
				
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		double PL = face.varL[controls.fP];
		double PR = face.varR[controls.fP];
		
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			UnF += (
				wCL*( 
				(controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				(controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				(controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				wCR*( 
				(controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				(controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				(controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = PR;
			
			UnF += (
				wCL*( 
				(controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				(controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				(controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				wCR*( 
				(controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				(controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				(controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			++proc_num;
		}
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		double YiL = face.varL[controls.fVF[0]];
		double YiR = face.varR[controls.fVF[0]];
		double YiF =YiL*weightL+YiR*weightR;
		
		double conv_flux_L = weightL*UnF*area;
		double conv_flux_R = weightR*UnF*area;
		
		// own, ngb
		linA[face.owner] += (+conv_flux_L);
		linAR[i] = (+conv_flux_R);
		// ngb, own
		linAL[i] = (-conv_flux_L);
		
		// convective term
		linB[face.owner] -= YiF*UnF*area;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// ngb, own
			linA[face.neighbour] += (-conv_flux_R);
			
			// convective term
			linB[face.neighbour] += YiF*UnF*area;
		}
		
	}
	
	
	// linear solver : PETSc library
	vector<double> resiVar(mesh.cells.size(),0.0);

	string solver;
	double relTol;
	double tolerance;
	string preconditioner;
	int maxIter;
	
	if(controls.iterPBs != controls.iterPBsMax-1){
		solver = controls.solverVF[0];
		relTol = controls.toleranceVF[0];
		tolerance = controls.relTolVF[0];
		preconditioner = controls.preconditionerVF[0];
		maxIter = controls.maxIterVF[0];
	
	}
	else{
		solver = controls.solverFinalU;
		relTol = controls.toleranceFinalVF[0];
		tolerance = controls.relTolFinalVF[0];
		preconditioner = controls.preconditionerFinalVF[0];
		maxIter = controls.maxIterFinalVF[0];
	
	}
	
	solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
		solver, tolerance, relTol, preconditioner, maxIter);
		
	// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
		// solver, tolerance, relTol, preconditioner, maxIter); 
		
		
	// update
	double normVF = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		
		cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
		normVF += pow(resiVar[i],2.0);
		
	}
	
	
	linA.clear();
	linAL.clear();
	linAR.clear();
	linB.clear();
	resiVar.clear();
	
	return sqrt(normVF);
	
	
}









double SEMO_Solvers_Builder::calcVolfracExplicitEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	SEMO_MPI_Builder mpi;
	
	SEMO_Utility_Math math;
	
	
	// gradient P
	vector<vector<double>> gradP;
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
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
	
	
	vector<double> gradPx_recv;
	vector<double> gradPy_recv;
	vector<double> gradPz_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> gradPx_send;
		vector<double> gradPy_send;
		vector<double> gradPz_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				gradPx_send.push_back(gradP[face.owner][0]);
				gradPy_send.push_back(gradP[face.owner][1]);
				gradPz_send.push_back(gradP[face.owner][2]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					gradPx_send, gradPx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradPy_send, gradPy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradPz_send, gradPz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradPx_send.clear();
		gradPy_send.clear();
		gradPz_send.clear();
	}
	





	vector<double> linA(mesh.cells.size(),0.0);
	vector<double> linAL(mesh.faces.size(),0.0);
	vector<double> linAR(mesh.faces.size(),0.0);
	vector<double> linB(mesh.cells.size(),0.0);

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		
		// time marching -> first order euler
		linA[i] = cell.volume/controls.timeStep;
		
		linB[i] = -(
			1.0*cell.var[controls.VF[0]] - 1.0*cell.var[controls.oldVF[0]]
			)*cell.volume/controls.timeStep;
			
		// time marching -> second order upwind euler
		// linA[i] = 1.5*cell.volume/controls.timeStep;
		
		// linB[i] = -(
			// 1.5*cell.var[controls.VF[0]] - 2.0*cell.var[controls.oldVF[0]] + 0.5*cell.var[controls.Qm[3]]
			// )*cell.volume/controls.timeStep;
		
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
		
		double Ef[3];
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
				
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		double PL = face.varL[controls.fP];
		double PR = face.varR[controls.fP];
		
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			UnF += (
				wCL*( 
				(controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				(controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				(controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				wCR*( 
				(controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				(controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				(controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = PR;
			
			UnF += (
				wCL*( 
				(controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				(controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				(controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				wCR*( 
				(controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				(controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				(controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			++proc_num;
		}
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		double YiL = face.varL[controls.fVF[0]];
		double YiR = face.varR[controls.fVF[0]];
		double YiF =YiL*weightL+YiR*weightR;
		
		double conv_flux_L = weightL*UnF*area;
		double conv_flux_R = weightR*UnF*area;
		
		// convective term
		linB[face.owner] -= YiF*UnF*area;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// convective term
			linB[face.neighbour] += YiF*UnF*area;
		}
		
	}
	
		
	// update
	double normVF = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.VF[0]] += linB[i]/linA[i];
		cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
		normVF += pow(linB[i]/linA[i],2.0);
		
	}
	
	
	linA.clear();
	linAL.clear();
	linAR.clear();
	linB.clear();
	
	return sqrt(normVF);
	
	
}









// #include "build.h"
// #include <cmath>
// #include <array>
// #include <numeric>

// void SEMO_Solvers_Builder::calcVolfracEq(
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
	
	// SEMO_Utility_Math math;
	
	
	// // gradient P
	// vector<vector<double>> gradP;
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// // math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
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
	

	
	// // vector<double> gradPx_recv;
	// // vector<double> gradPy_recv;
	// // vector<double> gradPz_recv;
	// // if(size>1){
		// // // processor faces
		// // // gradP , 
		// // vector<double> gradPx_send;
		// // vector<double> gradPy_send;
		// // vector<double> gradPz_send;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // gradPx_send.push_back(gradP[face.owner][0]);
				// // gradPy_send.push_back(gradP[face.owner][1]);
				// // gradPz_send.push_back(gradP[face.owner][2]);
			// // }
		// // }
		// // SEMO_MPI_Builder mpi;
		
		// // mpi.setProcsFaceDatasDouble(
					// // gradPx_send, gradPx_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // mpi.setProcsFaceDatasDouble(
					// // gradPy_send, gradPy_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // mpi.setProcsFaceDatasDouble(
					// // gradPz_send, gradPz_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // gradPx_send.clear();
		// // gradPy_send.clear();
		// // gradPz_send.clear();
	// // }
	





	// // vector<double> cellVolume_recv;
	// // if(size>1){
		// // vector<double> cellVolume_send;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // cellVolume_send.push_back(mesh.cells[face.owner].volume);
			// // }
		// // }
		// // mpi.setProcsFaceDatasDouble(
					// // cellVolume_send, cellVolume_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // cellVolume_send.clear();
		
	// // }


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






	// vector<double> linA(mesh.cells.size(),0.0);
	// vector<double> linAL(mesh.faces.size(),0.0);
	// vector<double> linAR(mesh.faces.size(),0.0);
	// vector<double> linB(mesh.cells.size(),0.0);

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		
		// // time marching -> first order euler
		// // linA[i] = cell.volume/controls.timeStep;
		
		// // linB[i] = -(
			// // 1.0*cell.var[controls.VF[0]] - 1.0*cell.var[controls.oldVF[0]]
			// // )*cell.volume/controls.timeStep;
			
		// // time marching -> second order upwind euler
		// linA[i] = 1.5*cell.volume/controls.timeStep;
		
		// linB[i] = -(
			// 1.5*cell.var[controls.VF[0]] - 2.0*cell.var[controls.oldVF[0]] + 0.5*cell.var[controls.Qm[3]]
			// )*cell.volume/controls.timeStep;
		
	// }
	
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
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double RhoL = face.varL[controls.fRho];
		// double RhoR = face.varR[controls.fRho];
	
		// // double tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
	
		// // double dPN = 
			// // sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
			                        // // std::begin(distanceCells), 0.0));
		// // double dPNEff = area/dPN/dPN*
			// // (std::inner_product(std::begin(nvec), std::end(nvec), 
			                    // // std::begin(distanceCells), 0.0));
		
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
		
		
		// double PL = face.varL[controls.fP];
		// double PR = face.varR[controls.fP];
		
		
		// double UnF = wCL*UnL+wCR*UnR;
		// // double UnF = 0.5*UnL+0.5*UnR;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// UnF = face.var[controls.Un];
		// }
		
		
		
		

		// // double UnF = wCL*UnL+wCR*UnR;

		// // double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			// // UnF += (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			// // UnF -= nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		// // }
		// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = PR;
			
			// // UnF += (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			// // UnF -= tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			// // UnF -=  nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			// // ++proc_num;
		// // }
		
		
		
		
		
		
		
		
		
		// // // double RhieChowCoeff = 0.8;
		// // double RhieChowCoeff = 1.0;
		
		// // double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		// // // double tmp1 = (0.5/RhoL + 0.5/RhoR)*controls.timeStep;
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			
			// // // UnF = 
				// // // wCL*( 
				// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // // wCR*( 
				// // // (UR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][0])*nvec[0] +
				// // // (VR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][1])*nvec[1] +
				// // // (WR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][2])*nvec[2] );
				
			// // UnF += RhieChowCoeff * (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // // UnF += RhieChowCoeff * (
				// // // 0.5*( 
				// // // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // // 0.5*( 
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPN;
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPN * overRelaxCoeff;
			// // UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			// // // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradP[face.neighbour][0]);
			// // // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradP[face.neighbour][1]);
			// // // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradP[face.neighbour][2]);
			
			// // // gradPf[0] = (newWCL*gradP[face.owner][0]+newWCR*gradP[face.neighbour][0]);
			// // // gradPf[1] = (newWCL*gradP[face.owner][1]+newWCR*gradP[face.neighbour][1]);
			// // // gradPf[2] = (newWCL*gradP[face.owner][2]+newWCR*gradP[face.neighbour][2]);
			
			// // // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// // UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		// // }
		// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // // double cellVolL = mesh.cells[face.owner].volume;
			// // // double cellVolR = cellVolume_recv[proc_num];
			
			// // // // tmp1 = (cellVolL+cellVolR)/
				// // // // (cellVolL*RhoL/controls.timeStep + cellVolR*RhoR/controls.timeStep);
				
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = PR;
			
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			
			// // // UnF = 
				// // // wCL*( 
				// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // // wCR*( 
				// // // (UR + 1.0/linAD_recv[proc_num]*gradPx_recv[proc_num])*nvec[0] +
				// // // (VR + 1.0/linAD_recv[proc_num]*gradPy_recv[proc_num])*nvec[1] +
				// // // (WR + 1.0/linAD_recv[proc_num]*gradPz_recv[proc_num])*nvec[2] );
			
			// // UnF += RhieChowCoeff * (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			// // // UnF += RhieChowCoeff * (
				// // // 0.5*( 
				// // // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // // 0.5*( 
				// // // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
				
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPNEff;
			// // UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			// // // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradPx_recv[proc_num]);
			// // // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradPy_recv[proc_num]);
			// // // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradPz_recv[proc_num]);
			// // // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// // UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			// // ++proc_num;
		// // }
		// // else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// // // tmp1 = 1.0/linAD[face.owner];
			// // tmp1 = controls.timeStep/RhoL;
			
		// // }
		
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double YiL = face.varL[controls.fVF[0]];
		// double YiR = face.varR[controls.fVF[0]];
		// // double YiF = 0.5*(YiL+YiR);
		// double YiF =YiL*weightL+YiR*weightR;
		
		
		
		
		
		// // linAL[i] = -0.5*weightL*UnF*area;
		// // linAR[i] = +0.5*weightR*UnF*area;
		// linAL[i] = -weightL*UnF*area;
		// linAR[i] = +weightR*UnF*area;
		
		// linB[face.owner] -= YiF*UnF*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			// linA[face.neighbour] += (-linAR[i]);
			
			// // convective term
			// linB[face.neighbour] += YiF*UnF*area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			
		// }
		
	// }
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar(mesh.cells.size(),0.0);


	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverVF[0], controls.toleranceVF[0], 
			// controls.relTolVF[0], controls.preconditionerVF[0],
			// controls.maxIterVF[0]);
	
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverVF[0], controls.toleranceVF[0], 
			// // controls.relTolVF[0], controls.preconditionerVF[0],
			// // controls.maxIterVF[0]);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			// controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			// controls.maxIterFinalVF[0]);
	
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			// // controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			// // controls.maxIterFinalVF[0]);
	// }
	
	
	
	// // update
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// // resiVar[i] = linB[i]/linA[i];
		
		// cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		
		
		// residuals[5] += pow(resiVar[i],2.0)*cell.volume;
		
		
		// cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
	// }
	


	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	
// }
























































// #include "build.h"
// #include <cmath>
// #include <array>
// #include <numeric>

// double SEMO_Solvers_Builder::calcVolfracEq(
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
	
	// SEMO_Utility_Math math;
	
	
	// // // gradient P
	// // vector<vector<double>> gradP;
	// // // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // // math.calcMGG(mesh, controls.P, controls.fP, 1, 1.e-8, gradP);
	// // // math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// // math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
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
	
	
	// // vector<double> gradPx_recv;
	// // vector<double> gradPy_recv;
	// // vector<double> gradPz_recv;
	// // if(size>1){
		// // // processor faces
		// // // gradP , 
		// // vector<double> gradPx_send;
		// // vector<double> gradPy_send;
		// // vector<double> gradPz_send;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // gradPx_send.push_back(gradP[face.owner][0]);
				// // gradPy_send.push_back(gradP[face.owner][1]);
				// // gradPz_send.push_back(gradP[face.owner][2]);
			// // }
		// // }
		// // SEMO_MPI_Builder mpi;
		
		// // mpi.setProcsFaceDatasDouble(
					// // gradPx_send, gradPx_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // mpi.setProcsFaceDatasDouble(
					// // gradPy_send, gradPy_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // mpi.setProcsFaceDatasDouble(
					// // gradPz_send, gradPz_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // gradPx_send.clear();
		// // gradPy_send.clear();
		// // gradPz_send.clear();
	// // }
	





	// // vector<double> cellVolume_recv;
	// // if(size>1){
		// // vector<double> cellVolume_send;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // cellVolume_send.push_back(mesh.cells[face.owner].volume);
			// // }
		// // }
		// // mpi.setProcsFaceDatasDouble(
					// // cellVolume_send, cellVolume_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // cellVolume_send.clear();
		
	// // }


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






	// vector<double> linA(mesh.cells.size(),0.0);
	// vector<double> linAL(mesh.faces.size(),0.0);
	// vector<double> linAR(mesh.faces.size(),0.0);
	// vector<double> linB(mesh.cells.size(),0.0);

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		
		// // time marching -> first order euler
		// // linA[i] = cell.volume/controls.timeStep;
		
		// // linB[i] = -(
			// // 1.0*cell.var[controls.VF[0]] - 1.0*cell.var[controls.oldVF[0]]
			// // )*cell.volume/controls.timeStep;
			
		// // time marching -> second order upwind euler
		// linA[i] = 1.5*cell.volume/controls.timeStep;
		
		// linB[i] = -(
			// 1.5*cell.var[controls.VF[0]] - 2.0*cell.var[controls.oldVF[0]] + 0.5*cell.var[controls.Qm[3]]
			// )*cell.volume/controls.timeStep;
		
	// }
	
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
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double RhoL = face.varL[controls.fRho];
		// double RhoR = face.varR[controls.fRho];
	
		// // double tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
	
		// // double dPN = 
			// // sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
			                        // // std::begin(distanceCells), 0.0));
		// // double dPNEff = area/dPN/dPN*
			// // (std::inner_product(std::begin(nvec), std::end(nvec), 
			                    // // std::begin(distanceCells), 0.0));
		
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
		
		
		// double PL = face.varL[controls.fP];
		// double PR = face.varR[controls.fP];
		
		
		// double UnF = wCL*UnL+wCR*UnR;
		// // double UnF = 0.5*UnL+0.5*UnR;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// UnF = face.var[controls.Un];
		// }
		
		// // // double RhieChowCoeff = 0.8;
		// // double RhieChowCoeff = 1.0;
		
		// // double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		// // // double tmp1 = (0.5/RhoL + 0.5/RhoR)*controls.timeStep;
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = mesh.cells[face.neighbour].var[controls.P];
			
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			
			// // // UnF = 
				// // // wCL*( 
				// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // // wCR*( 
				// // // (UR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][0])*nvec[0] +
				// // // (VR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][1])*nvec[1] +
				// // // (WR + 1.0/linAD[face.neighbour]*gradP[face.neighbour][2])*nvec[2] );
				
			// // UnF += RhieChowCoeff * (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // // UnF += RhieChowCoeff * (
				// // // 0.5*( 
				// // // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // // 0.5*( 
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][0])*nvec[0] +
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][1])*nvec[1] +
				// // // (controls.timeStep/RhoR*gradP[face.neighbour][2])*nvec[2] ) );
				
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPN;
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPN * overRelaxCoeff;
			// // UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			// // // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradP[face.neighbour][0]);
			// // // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradP[face.neighbour][1]);
			// // // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradP[face.neighbour][2]);
			
			// // // gradPf[0] = (newWCL*gradP[face.owner][0]+newWCR*gradP[face.neighbour][0]);
			// // // gradPf[1] = (newWCL*gradP[face.owner][1]+newWCR*gradP[face.neighbour][1]);
			// // // gradPf[2] = (newWCL*gradP[face.owner][2]+newWCR*gradP[face.neighbour][2]);
			
			// // // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// // UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
		// // }
		// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // // double cellVolL = mesh.cells[face.owner].volume;
			// // // double cellVolR = cellVolume_recv[proc_num];
			
			// // // // tmp1 = (cellVolL+cellVolR)/
				// // // // (cellVolL*RhoL/controls.timeStep + cellVolR*RhoR/controls.timeStep);
				
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = PR;
			
			// // // tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			
			// // // UnF = 
				// // // wCL*( 
				// // // (UL + 1.0/linAD[face.owner]*gradP[face.owner][0])*nvec[0] +
				// // // (VL + 1.0/linAD[face.owner]*gradP[face.owner][1])*nvec[1] +
				// // // (WL + 1.0/linAD[face.owner]*gradP[face.owner][2])*nvec[2] ) +
				// // // wCR*( 
				// // // (UR + 1.0/linAD_recv[proc_num]*gradPx_recv[proc_num])*nvec[0] +
				// // // (VR + 1.0/linAD_recv[proc_num]*gradPy_recv[proc_num])*nvec[1] +
				// // // (WR + 1.0/linAD_recv[proc_num]*gradPz_recv[proc_num])*nvec[2] );
			
			// // UnF += RhieChowCoeff * (
				// // wCL*( 
				// // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // wCR*( 
				// // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
			
			// // // UnF += RhieChowCoeff * (
				// // // 0.5*( 
				// // // (controls.timeStep/RhoL*gradP[face.owner][0])*nvec[0] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][1])*nvec[1] +
				// // // (controls.timeStep/RhoL*gradP[face.owner][2])*nvec[2] ) +
				// // // 0.5*( 
				// // // (controls.timeStep/RhoR*gradPx_recv[proc_num])*nvec[0] +
				// // // (controls.timeStep/RhoR*gradPy_recv[proc_num])*nvec[1] +
				// // // (controls.timeStep/RhoR*gradPz_recv[proc_num])*nvec[2] ) );
				
			// // // UnF = UnF - tmp1*(orgPR-orgPL)/dPNEff;
			// // UnF -= RhieChowCoeff * tmp1*(orgPR-orgPL)/dPN_e;
			
			
			// // // non-orthogonal, over-relaxed approach
			// // double gradPf[3];
			// // gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			// // gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			// // gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			// // // gradPf[0] = (0.5*gradP[face.owner][0]+0.5*gradPx_recv[proc_num]);
			// // // gradPf[1] = (0.5*gradP[face.owner][1]+0.5*gradPy_recv[proc_num]);
			// // // gradPf[2] = (0.5*gradP[face.owner][2]+0.5*gradPz_recv[proc_num]);
			// // // double delPhiECF = gradPf[0]*Ef[0]+gradPf[1]*Ef[1]+gradPf[2]*Ef[2];
			
			// // // gradPf[0] += nvec[0]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[1] += nvec[1]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			// // // gradPf[2] += nvec[2]*overRelaxCoeff*((orgPR-orgPL)/dPN - delPhiECF);
			
			// // UnF -= RhieChowCoeff * nonOrtholimiter * tmp1*(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
			
			
			// // ++proc_num;
		// // }
		// // else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// // // tmp1 = 1.0/linAD[face.owner];
			// // tmp1 = controls.timeStep/RhoL;
			
		// // }
		
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double YiL = face.varL[controls.fVF[0]];
		// double YiR = face.varR[controls.fVF[0]];
		// // double YiF = 0.5*(YiL+YiR);
		// double YiF =YiL*weightL+YiR*weightR;
		
		
		
		
		
		// // linAL[i] = -0.5*weightL*UnF*area;
		// // linAR[i] = +0.5*weightR*UnF*area;
		// linAL[i] = -weightL*UnF*area;
		// linAR[i] = +weightR*UnF*area;
		
		// linB[face.owner] -= YiF*UnF*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			// linA[face.neighbour] += (-linAR[i]);
			
			// // convective term
			// linB[face.neighbour] += YiF*UnF*area;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			
		// }
		
	// }
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar(mesh.cells.size(),0.0);


	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverVF[0], controls.toleranceVF[0], 
			// controls.relTolVF[0], controls.preconditionerVF[0],
			// controls.maxIterVF[0]);
	
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverVF[0], controls.toleranceVF[0], 
			// // controls.relTolVF[0], controls.preconditionerVF[0],
			// // controls.maxIterVF[0]);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			// controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			// controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			// controls.maxIterFinalVF[0]);
	
		// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
			// // controls.solverFinalVF[0], controls.toleranceFinalVF[0], 
			// // controls.relTolFinalVF[0], controls.preconditionerFinalVF[0],
			// // controls.maxIterFinalVF[0]);
	// }
	
	
	
	// // update
	// double normVF = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// // resiVar[i] = linB[i]/linA[i];
		
		// cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		
		// normVF += pow(resiVar[i],2.0);
		
		
		// residuals[5] += pow(resiVar[i],2.0)*cell.volume;
		
		
		// cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
		
	// }
	


	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB.clear();
	// resiVar.clear();
	
	
	// return sqrt(normVF);
	
	
// }


















