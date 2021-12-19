#include "../build.h"
#include <cmath>
#include <array>
#include <numeric>
#include <ctime>


double SEMO_Solvers_Builder::calcCompMomentumEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	
	// // diagonal terms
    // int B_n = 1;
    // int A_n = B_n * B_n;
	
	
	// // vector<clock_t> startTime;
	
	// SEMO_Utility_Math math;
	
	// int proc_num=0;

	// // gradient P
	// vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	// // math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
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
	
	
	
	
	
	// vector<double> P_recv, U_recv, V_recv, W_recv;
	// if(size>1){
		// // processor faces
		// vector<double> P_send, U_send, V_send, W_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// P_send.push_back(mesh.cells[face.owner].var[controls.P]);
				// U_send.push_back(mesh.cells[face.owner].var[controls.U]);
				// V_send.push_back(mesh.cells[face.owner].var[controls.V]);
				// W_send.push_back(mesh.cells[face.owner].var[controls.W]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// P_send, P_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// U_send, U_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// V_send, V_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// W_send, W_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// P_send.clear();
		// U_send.clear();
		// V_send.clear();
		// W_send.clear();
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
	
	
	
	
	
	
	// // =======================================
	// // surface tension
	// vector<vector<double>> gradAi;
	// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// // math.calcMGG(mesh, controls.VF[0], controls.fVF[0], 1, 1.e-8, gradAi);
	// // math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradAi);
	

	// vector<double> gradAix_recv;
	// vector<double> gradAiy_recv;
	// vector<double> gradAiz_recv;
	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradAix_send;
		// vector<double> gradAiy_send;
		// vector<double> gradAiz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradAix_send.push_back(gradAi[face.owner][0]);
				// gradAiy_send.push_back(gradAi[face.owner][1]);
				// gradAiz_send.push_back(gradAi[face.owner][2]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// gradAix_send, gradAix_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradAiy_send, gradAiy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradAiz_send, gradAiz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradAix_send.clear();
		// gradAiy_send.clear();
		// gradAiz_send.clear();
	// }
	
	
	// vector<double> kappa;
	// this->calcCurvature(mesh, controls.VF[0], kappa);
	
	// vector<double> kappa_recv;
	// if(size>1){
		// // processor faces
		// vector<double> kappa_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// kappa_send.push_back(kappa[face.owner]);
			// }
		// }
		// mpi.setProcsFaceDatas(
					// kappa_send, kappa_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// kappa_send.clear();
	// }
	// // =======================================



	
	
	

	// // diagonal terms
    // vector<int> A_rows(mesh.cells.size()*A_n, 0);
    // vector<int> A_cols(mesh.cells.size()*A_n, 0);
    // vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
    // // vector<double> B(mesh.cells.size()*B_n, 0.0);
    // vector<vector<double>> B(3,vector<double>(mesh.cells.size(), 0.0));
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];

        // int ijStart_local = B_n*(i) - 1;
        // int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
        // int Astart = A_n*(i) - 1;
        // // int id = Astart;
		
		// double P = cell.var[controls.P];
		// double Rho = cell.var[controls.Rho];
		// double U = cell.var[controls.U];
		// double V = cell.var[controls.V];
		// double W = cell.var[controls.W];
		// double oldP = cell.var[controls.oldP];
		// double oldU = cell.var[controls.oldU];
		// double oldV = cell.var[controls.oldV];
		// double oldW = cell.var[controls.oldW];
		// double oldRho = cell.var[controls.oldRho];
		// double oldHt = cell.var[controls.oldHt];
		// vector<double> oldMF;
		// oldMF.push_back(cell.var[controls.oldMF[0]]);
		
		// double old2P = cell.var[controls.Qm[0]];
		// double old2U = cell.var[controls.Qm[1]];
		// double old2V = cell.var[controls.Qm[2]];
		// double old2W = cell.var[controls.Qm[3]];
		// double old2Rho = cell.var[controls.Qm[5]];
		// double old2Ht = cell.var[controls.Qm[6]];
		// vector<double> old2MF;
		// old2MF.push_back(cell.var[controls.Qm[7]]);
		
		// double volume = cell.volume;
		// double timeStep = controls.timeStep;
		
		// int id = 0;
		// int i_glo = 0;
		// int j_glo = 0;
		// int i_loc = i;
		// int str_glo = mesh.startCellGlobal*B_n;
		// int step_loc = mesh.cells.size();
		// int step_glo = mesh.ncellsTotal;
		
		// double cTime = 1.0;
		// // double cTime = 1.5;
		
		// double Ht = cell.var[controls.Ht];
		
		// double dHtDP = cell.var[controls.dHtDP];
		// double dHtDT = cell.var[controls.dHtDT];
		
		// double dRhoDP = cell.var[controls.dRhoDP];
		// double dRhoDT = cell.var[controls.dRhoDT];
		// vector<double> dRhoDMF;
		// vector<double> dHtDMF;
		// vector<double> MF;
		// for(int ns=0; ns<controls.nSp-1; ++ns){
			// MF.push_back(cell.var[controls.MF[ns]]);
			// dRhoDMF.push_back(cell.var[controls.dRhoDMF[ns]]);
			// dHtDMF.push_back(cell.var[controls.dHtDMF[ns]]);
		// }


        // id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        // A_rows[id] = i_glo; A_cols[id] = j_glo;
        // A_vals[id] = cTime * ( Rho )*volume/timeStep;
        
		
        // B[0][step_loc*0 + i_loc] = -(Rho*U - oldRho*oldU)*volume / timeStep;
        // B[1][step_loc*0 + i_loc] = -(Rho*V - oldRho*oldV)*volume / timeStep;
        // B[2][step_loc*0 + i_loc] = -(Rho*W - oldRho*oldW)*volume / timeStep;
		
		// // gravity force terms
        // B[0][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[0]*volume;
        // B[1][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[1]*volume;
        // B[2][step_loc*0 + i_loc] += Rho*controls.gravityAcceleration[2]*volume;


		// // surface tension force terms
        // B[0][step_loc*0 + i_loc] += volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
        // B[1][step_loc*0 + i_loc] += volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
        // B[2][step_loc*0 + i_loc] += volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
		
	// }
	
	
	
	
	// proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		// int ijStartL_local = B_n*(face.owner) - 1;
		// int ijStartR_local = B_n*(face.neighbour) - 1;
		
        // int ijStartL = B_n*mesh.startCellGlobal + ijStartL_local;
        // int ijStartR = B_n*mesh.startCellGlobal + ijStartR_local;
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// ijStartR = B_n*mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] + 
					   // B_n*(mesh.procNeighbCellNo[proc_num]) - 1;
		// }

		
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
		// // double wCL = face.wVC;
		// // double wCL = 0.5;
		// double wCR = 1.0-wCL;
		
		// double wCPL = wCL;
		// double wCPR = wCR;
		// // double wCPL = 0.5;
		// // double wCPR = 0.5;
		
		// double UL = face.varL[controls.fU];
		// double VL = face.varL[controls.fV];
		// double WL = face.varL[controls.fW];
		// double PL = face.varL[controls.fP];
		// double RhoL = face.varL[controls.fRho];
		// double CL = face.varL[controls.fC];
		// double HtL = face.varL[controls.fHt];
		// double dRhoDPL = face.varL[controls.fdRhoDP];
		// double dRhoDTL = face.varL[controls.fdRhoDT];
		// double dHtDPL = face.varL[controls.fdHtDP];
		// double dHtDTL = face.varL[controls.fdHtDT];
		// vector<double> MFL, dRhoDMFL, dHtDMFL;
		// MFL.push_back(face.varL[controls.fMF[0]]);
		// dRhoDMFL.push_back(face.varL[controls.fdRhoDMF[0]]);
		// dHtDMFL.push_back(face.varL[controls.fdHtDMF[0]]);
		// double muL = face.varL[controls.fmu];
		
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		// double PR = face.varR[controls.fP];
		// double RhoR = face.varR[controls.fRho];
		// double CR = face.varR[controls.fC];
		// double HtR = face.varR[controls.fHt];
		// double dRhoDPR = face.varR[controls.fdRhoDP];
		// double dRhoDTR = face.varR[controls.fdRhoDT];
		// double dHtDPR = face.varR[controls.fdHtDP];
		// double dHtDTR = face.varR[controls.fdHtDT];
		// vector<double> MFR, dRhoDMFR, dHtDMFR;
		// MFR.push_back(face.varR[controls.fMF[0]]);
		// dRhoDMFR.push_back(face.varR[controls.fdRhoDMF[0]]);
		// dHtDMFR.push_back(face.varR[controls.fdHtDMF[0]]);
		// double muR = face.varR[controls.fmu];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double UnF = wCL*UnL+wCR*UnR;
		// double muF = wCL*muL + wCR*muR;
		
		// double Rho_star = 1.0 / (wCPL/RhoL + wCPR/RhoR);
		// double tmp1 = controls.timeStep / Rho_star;
		// // // F.Denner et al., 2013
		// // double d_F = 0.0;
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
			// // d_F += wCR*( A_p_vals[face.neighbour] + RhoR/controls.timeStep );
		// // }
		// // else{
			// // d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
			// // d_F += wCR*( A_p_vals_recv[proc_num] + RhoR/controls.timeStep );
		// // }
		// // tmp1 = 1.0/d_F;
	
	
		// double dPN = face.magPN;
		// // vector<double> Ef(3,0.0);
		// // Ef[0] = face.vecPN[0]/dPN;
		// // Ef[1] = face.vecPN[1]/dPN;
		// // Ef[2] = face.vecPN[2]/dPN;
		// // double alpha = nvec[0]*Ef[0] + nvec[1]*Ef[1] + nvec[2]*Ef[2];
		// double alpha = face.alphaF;
		


		// double orgPL = mesh.cells[face.owner].var[controls.P];
		// double orgUL = mesh.cells[face.owner].var[controls.U];
		// double orgVL = mesh.cells[face.owner].var[controls.V];
		// double orgWL = mesh.cells[face.owner].var[controls.W];
		// double orgPR = 0.0;
		// double orgUR = 0.0;
		// double orgVR = 0.0;
		// double orgWR = 0.0;
		// vector<double> gradPL(3,0.0);
		// vector<double> gradUL(3,0.0);
		// vector<double> gradVL(3,0.0);
		// vector<double> gradWL(3,0.0);
		// vector<double> gradPR(3,0.0);
		// vector<double> gradUR(3,0.0);
		// vector<double> gradVR(3,0.0);
		// vector<double> gradWR(3,0.0);
		// for(int ii=0; ii<3; ++ii){
			// gradPL[ii] = gradP[face.owner][ii];
			// gradUL[ii] = gradU[face.owner][ii];
			// gradVL[ii] = gradV[face.owner][ii];
			// gradWL[ii] = gradW[face.owner][ii];
		// }
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// orgPR = mesh.cells[face.neighbour].var[controls.P];
			// orgUR = mesh.cells[face.neighbour].var[controls.U];
			// orgVR = mesh.cells[face.neighbour].var[controls.V];
			// orgWR = mesh.cells[face.neighbour].var[controls.W];
			// for(int ii=0; ii<3; ++ii){
				// gradPR[ii] = gradP[face.neighbour][ii];
				// gradUR[ii] = gradU[face.neighbour][ii];
				// gradVR[ii] = gradV[face.neighbour][ii];
				// gradWR[ii] = gradW[face.neighbour][ii];
			// }
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// orgPR = P_recv[proc_num];
			// orgUR = U_recv[proc_num];
			// orgVR = V_recv[proc_num];
			// orgWR = W_recv[proc_num];
			// gradPR[0] = gradPx_recv[proc_num]; gradPR[1] = gradPy_recv[proc_num]; gradPR[2] = gradPz_recv[proc_num];
			// gradUR[0] = gradUx_recv[proc_num]; gradUR[1] = gradUy_recv[proc_num]; gradUR[2] = gradUz_recv[proc_num];
			// gradVR[0] = gradVx_recv[proc_num]; gradVR[1] = gradVy_recv[proc_num]; gradVR[2] = gradVz_recv[proc_num];
			// gradWR[0] = gradWx_recv[proc_num]; gradWR[1] = gradWy_recv[proc_num]; gradWR[2] = gradWz_recv[proc_num];
		// }
		
		
		
		
		// // pressure correction (cell->face)
		// double UnF_RC = 0.0;
		// for(int ii=0; ii<3; ++ii){
			// UnF_RC += alpha * tmp1 * Rho_star * wCPL/RhoL*gradPL[ii]*face.unitNomalsPN[ii];
			// UnF_RC += alpha * tmp1 * Rho_star * wCPR/RhoR*gradPR[ii]*face.unitNomalsPN[ii];
		// }
		
		// UnF_RC -= alpha * tmp1*(orgPR-orgPL)/dPN;
		// // for(int ii=0; ii<3; ++ii){
			// // UnF_RC -= tmp1 * (gradPR[ii]*face.vecNdN[ii] - gradPL[ii]*face.vecPdP[ii])/dPN;
		// // }
		
		// // YYL riemann
		// double Cbar = 0.5*(CL + CR);
		// double ML = UnL/Cbar;
		// double MR = UnR/Cbar;
		// double MLP = M_func(ML,1.0,0.125);
		// double MRM = M_func(MR,-1.0,0.125);
		// double PLP = pre_func(ML,1.0,0.1875);
		// double PRM = pre_func(MR,-1.0,0.1875);
		// double KLR = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR));
		// double Mdash = min(1.0,KLR/Cbar);
		// double G_coeff = -max(min(ML,0.0),-1.0)*min(max(MR,0.0),1.0);
		// double Vn = (RhoL*abs(UnL)+RhoR*abs(UnR)) / (RhoL + RhoR);
		// double VnP = (1.0-G_coeff)*Vn + G_coeff*abs(UnL);
		// double VnM = (1.0-G_coeff)*Vn + G_coeff*abs(UnR);
		// double Mdot = 0.5*(RhoL*(UnL+VnP)+RhoR*(UnR-VnM)-pow((1.0-Mdash),2.0)/Cbar*(PR-PL));
		// double MdotL = 0.5*(Mdot+abs(Mdot));
		// double MdotR = 0.5*(Mdot-abs(Mdot));
		
		// double chi = pow((1.0-Mdash),2.0);
		// double rhohat = 0.5*(RhoL+RhoR);
		
		
		// // UnF = (weightL/RhoL + weightR/RhoR)*Mdot + UnF_RC;
		// // UnF = (weightL/RhoL + weightR/RhoR)*Mdot;
        // // UnF = (MLP+MRM)*Cbar + UnF_RC - 0.5*chi/rhohat/Cbar*(PR-PL);
        // UnF = (MLP+MRM)*Cbar + UnF_RC - 0.5*chi/rhohat/Cbar*(orgPR-orgPL);
        // // UnF = 0.5*(UnL+UnR) + UnF_RC - 0.5*chi/rhohat/Cbar*(orgPR-orgPL);
        // // UnF = 0.5*(UnL+UnR) + UnF_RC;
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = 1.0 - weightL;
		
		// double WUL = 0.5;
		// double WUR = 0.5;
        // if(abs(ML) > 1.0){ 
            // WUL = 0.5*(1.0 + ((ML > 0.0) ? 1.0 : -1.0) );
		// }
        // else{
            // WUL = 0.5*(ML + 1.0) + 0.125*2.0*2.0*(ML*ML-1.0);
        // }
        // if(abs(MR) > 1.0){ 
            // WUR = 0.5*(1.0 - ((MR > 0.0) ? 1.0 : -1.0) );
		// }
        // else{
            // WUR = -0.5*(MR - 1.0) - 0.125*2.0*2.0*(MR*MR-1.0);
        // }
		
		// double WPL = 0.5 - 0.5*Mdash*PLP*PRM*0.5/Cbar*(UnR-UnL) + 0.5*Mdash*(PLP+PRM-1.0) + 0.5*(PLP-PRM);//PLP;
		// double WPR = 0.5 - 0.5*Mdash*PLP*PRM*0.5/Cbar*(UnR-UnL) + 0.5*Mdash*(PLP+PRM-1.0) - 0.5*(PLP-PRM);//PRM;
		
		// double diffDP = alpha * tmp1/dPN + 0.5*chi/rhohat/Cbar; 
		// // double diffDP = alpha * tmp1/dPN; 
		
		// double UF   = weightL*UL   + weightR*UR;
		// double VF   = weightL*VL   + weightR*VR;
		// double WF   = weightL*WL   + weightR*WR;
		// double RhoF = weightL*RhoL + weightR*RhoR;
		// double HtF  = weightL*HtL  + weightR*HtR;
		// vector<double> MFF;
		// MFF.push_back(weightL*MFL[0] + weightR*MFR[0]);
		
		// double PF = WPL*PL + WPR*PR;
		
		
		// // double diff2nd_U = 0.0;
		// // double diff2nd_V = 0.0;
		// // double diff2nd_W = 0.0;
		// // // for(int ii=0; ii<3; ++ii){
			// // // diff2nd_U += muF*(gradUR[ii]*face.vecNdN[ii] - gradUL[ii]*face.vecPdP[ii])/dPN*area;
			// // // diff2nd_V += muF*(gradVR[ii]*face.vecNdN[ii] - gradVL[ii]*face.vecPdP[ii])/dPN*area;
			// // // diff2nd_W += muF*(gradWR[ii]*face.vecNdN[ii] - gradWL[ii]*face.vecPdP[ii])/dPN*area;
		// // // }
		
		
		
		// int str_glo_L = mesh.startCellGlobal*B_n;
		// int str_glo_R = mesh.startCellGlobal*B_n;
		// int step_loc_L = mesh.cells.size();
		// int step_loc_R = mesh.cells.size();
		// int i_loc_L = face.owner;
		// int i_loc_R = face.neighbour;
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// str_glo_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]]*B_n;
			// step_loc_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]+1] - 
			             // mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
			// i_loc_R = mesh.procNeighbCellNo[proc_num];
		// }
		
		// int id = -1;
		// int i_glo = -1;
		// int j_glo = -1;
		// double conv_flux = 0.0;
		// double diff_flux = 0.0;


		// // viscous terms
		// vector<double> dUdxyzF(3,0.0);
		// vector<double> dVdxyzF(3,0.0);
		// vector<double> dWdxyzF(3,0.0);
		
		// double dUdnF = 0.0;
		// double dVdnF = 0.0;
		// double dWdnF = 0.0;

		// for(int ii=0; ii<3; ++ii){
			// dUdxyzF[ii] += wCL * gradUL[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dUdxyzF[ii] += wCR * gradUR[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			
			// dVdxyzF[ii] += wCL * gradVL[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dVdxyzF[ii] += wCR * gradVR[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			
			// dWdxyzF[ii] += wCL * gradWL[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dWdxyzF[ii] += wCR * gradWR[ii];// * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			
			// dUdnF += wCL * gradUL[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dUdnF += wCR * gradUR[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			
			// dVdnF += wCL * gradVL[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dVdnF += wCR * gradVR[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			
			// dWdnF += wCL * gradWL[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
			// dWdnF += wCR * gradWR[ii] * ( nvec[ii] - alpha * face.unitNomalsPN[ii] );
		// }
		// dUdnF += alpha * (orgUR-orgUL)/dPN;
		// dVdnF += alpha * (orgVR-orgVL)/dPN;
		// dWdnF += alpha * (orgWR-orgWL)/dPN;
		
		
		// // Ref : Blazek's book, pp.20-21
		// double txx = dUdxyzF[0] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
		// double txy = dVdxyzF[0];
		// double txz = dWdxyzF[0];
		// double tyx = dUdxyzF[1];
		// double tyy = dVdxyzF[1] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
		// double tyz = dWdxyzF[1];
		// double tzx = dUdxyzF[2];
		// double tzy = dVdxyzF[2];
		// double tzz = dWdxyzF[2] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
		// // txx += dUdxyzF[0];
		// // txy += dUdxyzF[1];
		// // txz += dUdxyzF[2];
		// // tyx += dVdxyzF[0];
		// // tyy += dVdxyzF[1];
		// // tyz += dVdxyzF[2];
		// // tzx += dWdxyzF[0];
		// // tzy += dWdxyzF[1];
		// // tzz += dWdxyzF[2];
		
		
		// double diff_viscous = muF * alpha / dPN;
		
		// vector<double> viscous_flux(3,0.0);
		// viscous_flux[0] = muF * (txx*nvec[0] + txy*nvec[1] + txz*nvec[2] + dUdnF) * area;
		// // flux[1] -= nx * 2.0/3.0 * rhoF * tkei;
		// viscous_flux[1] = muF * (tyx*nvec[0] + tyy*nvec[1] + tyz*nvec[2] + dVdnF) * area;
		// // flux[2] -= ny * 2.0/3.0 * rhoF * tkei;
		// viscous_flux[2] = muF * (tzx*nvec[0] + tzy*nvec[1] + tzz*nvec[2] + dWdnF) * area;
		// // flux[3] -= nz * 2.0/3.0 * rhoF * tkei;
		// // viscous_flux[4] = muF * (txx*UF + txy*VF + txz*WF + dUdxyzF[0]*UF + dUdxyzF[1]*VF + dUdxyzF[2]*WF) * nvec[0] * area +
                          // // muF * (tyx*UF + tyy*VF + tyz*WF + dVdxyzF[0]*UF + dVdxyzF[1]*VF + dVdxyzF[2]*WF) * nvec[1] * area +
                          // // muF * (tzx*UF + tzy*VF + tzz*WF + dWdxyzF[0]*UF + dWdxyzF[1]*VF + dWdxyzF[2]*WF) * nvec[2] * area;
		

        // //------------------------
		// id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
		// A_vals[id] += ( weightL * RhoL * UnF * area + diff_viscous * area );
		// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		// A_vals.push_back(( weightR * RhoR * UnF * area - diff_viscous * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ) {
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			// A_vals[id] -= ( weightR * RhoR * UnF * area - diff_viscous * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( weightL * RhoL * UnF * area + diff_viscous * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ) {
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
		

		// vector<double> convective_flux(3,0.0);
		// convective_flux[0] = RhoF * UnF * UF * area + PF * nvec[0] * area;
		// convective_flux[1] = RhoF * UnF * VF * area + PF * nvec[1] * area;
		// convective_flux[2] = RhoF * UnF * WF * area + PF * nvec[2] * area;


        // // ----------------------------
		// for(int ii=0; ii<3; ++ii){
			// B[ii][step_loc_L*0 + i_loc_L] -= ( convective_flux[ii] - viscous_flux[ii] );
		// }
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// for(int ii=0; ii<3; ++ii){
				// B[ii][step_loc_R*0 + i_loc_R] += ( convective_flux[ii] - viscous_flux[ii] );
			// }
		// }
		
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
		
	// }
	
	
	
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				
				// // int ijStartL_local = B_n*(face.owner) - 1;
				// // vector<int> id(B_n,0);
				// // id[0] = A_n*(face.owner) - 1;
				// // for(int j=1; j<B_n; ++j){
					// // id[j] = id[j-1] + 4;
				// // }

				
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
				
				// double UnF = 0.5*UnL+0.5*UnR;
				
				// double RhoF = 0.5*face.varL[controls.fRho] + 0.5*face.varR[controls.fRho];
				
				// double UF = 0.5*face.varL[controls.fU] + 0.5*face.varR[controls.fU];
				// double VF = 0.5*face.varL[controls.fV] + 0.5*face.varR[controls.fV];
				// double WF = 0.5*face.varL[controls.fW] + 0.5*face.varR[controls.fW];
				// double PF = 0.5*face.varL[controls.fP] + 0.5*face.varR[controls.fP];
			
				
				// vector<double> MFF, dRhoDMFF, dHtDMFF;
				// double dRhoDPF = 0.5*face.varL[controls.fdRhoDP] + 0.5*face.varR[controls.fdRhoDP];
				// double dRhoDTF = 0.5*face.varL[controls.fdRhoDT] + 0.5*face.varR[controls.fdRhoDT];
				// double dHtDPF = 0.5*face.varL[controls.fdHtDP] + 0.5*face.varR[controls.fdHtDP];
				// double dHtDTF = 0.5*face.varL[controls.fdHtDT] + 0.5*face.varR[controls.fdHtDT];
				// double HtF = 0.5*face.varL[controls.fHt] + 0.5*face.varR[controls.fHt];
				// MFF.push_back(0.5*face.varL[controls.fMF[0]] + 0.5*face.varR[controls.fMF[0]]);
				// dRhoDMFF.push_back(0.5*face.varL[controls.fdRhoDMF[0]] + 0.5*face.varR[controls.fdRhoDMF[0]]);
				// dHtDMFF.push_back(0.5*face.varL[controls.fdHtDMF[0]] + 0.5*face.varR[controls.fdHtDMF[0]]);
				
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
				
				// // velocity coefficient
				// double coeffUVW = 0.0;
				// double coeffUVW_diff = 0.0;
				// if( boundary.type[controls.U] == "fixedValue" ){
					// coeffUVW = 0.0;
					// coeffUVW_diff = 1.0;
				// }
				// else if( boundary.type[controls.U] == "zeroGradient" ){
					// coeffUVW = 1.0;
					// coeffUVW_diff = 0.0;
				// }
				// else if( boundary.type[controls.U] == "slip" ){
					// coeffUVW = 0.0;
					// coeffUVW_diff = 1.0;
				// }
				// else if( boundary.type[controls.U] == "noSlip" ){
					// coeffUVW = 0.0;
					// coeffUVW_diff = 1.0;
				// }
				// else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
					// coeffUVW = 0.0;
					// coeffUVW_diff = 1.0;
				// }
				// else if( boundary.type[controls.U] == "inletOutlet" ){
					// double ownNorVel =  
						// mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// coeffUVW = (ownNorVel > 1.0) ? 1.0 : 0.0;
					// coeffUVW_diff = (ownNorVel > 1.0) ? 0.0 : 1.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.U << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				// // temperature coefficient
				// double coeffT = 0.0;
				// double coeffT_diff = 0.0;
				// if( boundary.type[controls.T] == "fixedValue" ){
					// coeffT = 0.0;
					// coeffT_diff = 1.0;
				// }
				// else if( boundary.type[controls.T] == "zeroGradient" ){
					// coeffT = 1.0;
					// // coeffT = 0.5;
					// coeffT_diff = 0.0;
					
				// }
				// else if( boundary.type[controls.T] == "switch" ){
					// double machNum = 
						// sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
							 // pow(mesh.cells[face.owner].var[controls.V],2.0)+
							 // pow(mesh.cells[face.owner].var[controls.W],2.0))/
							  // mesh.cells[face.owner].var[controls.C];
					// coeffT = (machNum > 1.0) ? 0.0 : 1.0;
					// coeffT_diff = (machNum > 1.0) ? 1.0 : 0.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.T << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				// // mass-fraction coefficient
				// double coeffMF = 0.0;
				// double coeffMF_diff = 0.0;
				// if( boundary.type[controls.VF[0]] == "fixedValue" ){
					// coeffMF = 0.0;
					// coeffMF_diff = 1.0;
				// }
				// else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
					// coeffMF = 1.0;
					// coeffMF_diff = 0.0;
				// }
				// else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
					// double ownNorVel =  
						// mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// coeffMF = (ownNorVel > 1.0) ? 1.0 : 0.0;
					// coeffMF_diff = (ownNorVel > 1.0) ? 0.0 : 1.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				
				

				// int str_glo = mesh.startCellGlobal*B_n;
				// int step_loc = mesh.cells.size();
				// int step_glo = mesh.ncellsTotal;
				// int i_loc_L = face.owner;
				// int id=-1;
				
				
				// double dPN = face.magPN;
				// double RhoL = face.varL[controls.fRho];
				// double tmp1 = 1.0/RhoL*controls.timeStep;
				
				// // properties, viscous term
				// double muF = 0.5*face.varL[controls.fmu] + 0.5*face.varR[controls.fmu];
				
				// double orgUL = mesh.cells[face.owner].var[controls.U];
				// double orgUR = face.varR[controls.fU];
				// double orgVL = mesh.cells[face.owner].var[controls.V];
				// double orgVR = face.varR[controls.fV];
				// double orgWL = mesh.cells[face.owner].var[controls.W];
				// double orgWR = face.varR[controls.fW];


				// // viscous terms
				// vector<double> dUdxyzF(3,0.0);
				// vector<double> dVdxyzF(3,0.0);
				// vector<double> dWdxyzF(3,0.0);
				
				// double dUdnF = 0.0;
				// double dVdnF = 0.0;
				// double dWdnF = 0.0;

				// for(int ii=0; ii<3; ++ii){
					// dUdxyzF[ii] += gradU[face.owner][ii];
					// dVdxyzF[ii] += gradV[face.owner][ii];
					// dWdxyzF[ii] += gradW[face.owner][ii];
				// }
				// dUdnF += (orgUR-orgUL)/dPN;
				// dVdnF += (orgVR-orgVL)/dPN;
				// dWdnF += (orgWR-orgWL)/dPN;

				// // Ref : Blazek's book, pp.20-21
				// double txx = dUdxyzF[0] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
				// double txy = dVdxyzF[0];
				// double txz = dWdxyzF[0];
				// double tyx = dUdxyzF[1];
				// double tyy = dVdxyzF[1] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
				// double tyz = dWdxyzF[1];
				// double tzx = dUdxyzF[2];
				// double tzy = dVdxyzF[2];
				// double tzz = dWdxyzF[2] - 2.0/3.0 * (dUdxyzF[0] + dVdxyzF[1] + dWdxyzF[2]);
				// // txx += dUdxyzF[0];
				// // txy += dUdxyzF[1];
				// // txz += dUdxyzF[2];
				// // tyx += dVdxyzF[0];
				// // tyy += dVdxyzF[1];
				// // tyz += dVdxyzF[2];
				// // tzx += dWdxyzF[0];
				// // tzy += dWdxyzF[1];
				// // tzz += dWdxyzF[2];
						
				
				// vector<double> viscous_flux(3,0.0);
				// viscous_flux[0] = muF * (txx*nvec[0] + txy*nvec[1] + txz*nvec[2] + dUdnF) * area;
				// // flux[1] -= nx * 2.0/3.0 * rhoF * tkei;
				// viscous_flux[1] = muF * (tyx*nvec[0] + tyy*nvec[1] + tyz*nvec[2] + dVdnF) * area;
				// // flux[2] -= ny * 2.0/3.0 * rhoF * tkei;
				// viscous_flux[2] = muF * (tzx*nvec[0] + tzy*nvec[1] + tzz*nvec[2] + dWdnF) * area;
				// // flux[3] -= nz * 2.0/3.0 * rhoF * tkei;

				
				// id = step_loc*(B_n*0+0) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoF * UnF * area );
				// A_vals[id] += coeffUVW_diff *( muF / dPN * area );
				

				// vector<double> convective_flux(3,0.0);
				// convective_flux[0] = RhoF * UnF * UF * area + PF * nvec[0] * area;
				// convective_flux[1] = RhoF * UnF * VF * area + PF * nvec[1] * area;
				// convective_flux[2] = RhoF * UnF * WF * area + PF * nvec[2] * area;


				// // ----------------------------
				// for(int ii=0; ii<3; ++ii){
					// B[ii][step_loc*0 + i_loc_L] -= ( convective_flux[ii] - viscous_flux[ii] );
				// }
				
				
			// }
		// }
	// }
	
	
	// // linear solver : PETSc library
	// vector<vector<double>> resiVar(3,vector<double>(B_n*mesh.cells.size(),0.0));
	
	// solveAMGCL("momentum", mesh, B_n, A_rows, A_cols, A_vals, B[0], B[1], B[2], resiVar[0], resiVar[1], resiVar[2]);
	
	
	// // update
	// double relaxP = 1.0; //controls.prePreURF;
	// double relaxUVW = 1.0; //controls.momVelURF;
	// double relaxT = 1.0; //controls.momVelURF;
	// double relaxMF = 1.0; //controls.momVelURF;
	// // double relaxP = 0.3;
	// // double relaxUVW = 0.3;
	// // double relaxT = 0.3;
	// // double relaxMF = 0.3;
	
	
	// double maxP = -1.e9;
	// double minP = 1.e9;
	// vector<double> maxP_xyz(3,0.0);
	// vector<double> minP_xyz(3,0.0);
	// double maxUn = -1.e9;
	// vector<double> maxUn_xyz(3,0.0);
	// double maxT = -1.e9;
	// double minT = 1.e9;
	// vector<double> maxT_xyz(3,0.0);
	// vector<double> minT_xyz(3,0.0);
	
	// double normDelPrim = 0.0;
	// double normPrim = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
        // int ijStart = B_n*(i) - 1;
        // int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		// int step_loc = mesh.cells.size();
		// int i_loc = i;
		
		// cell.var[controls.U] +=     relaxUVW * resiVar[0][step_loc*0 + i_loc];
		// cell.var[controls.V] +=     relaxUVW * resiVar[1][step_loc*0 + i_loc];
		// cell.var[controls.W] +=     relaxUVW * resiVar[2][step_loc*0 + i_loc];
		
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
		
		// if( cell.var[controls.T] <= controls.minT ) 
			// cell.var[controls.T] = controls.minT;
		// if( cell.var[controls.T] >= controls.maxT ) 
			// cell.var[controls.T] = controls.maxT; 
		 

		// if(maxP<cell.var[controls.P]){
			// maxP = cell.var[controls.P];
			// maxP_xyz[0] = cell.x; maxP_xyz[1] = cell.y; maxP_xyz[2] = cell.z;
		// }
		// if(minP>cell.var[controls.P]){
			// minP = cell.var[controls.P];
			// minP_xyz[0] = cell.x; minP_xyz[1] = cell.y; minP_xyz[2] = cell.z;
		// }
		// double Un = sqrt(
		// cell.var[controls.U]*cell.var[controls.U]+
		// cell.var[controls.V]*cell.var[controls.V]+
		// cell.var[controls.W]*cell.var[controls.W]);
		// if(maxUn<Un){
			// maxUn = Un;
			// maxUn_xyz[0] = cell.x; maxUn_xyz[1] = cell.y; maxUn_xyz[2] = cell.z;
		// }
		// if(maxT<cell.var[controls.T]){
			// maxT = cell.var[controls.T];
			// maxT_xyz[0] = cell.x; maxT_xyz[1] = cell.y; maxT_xyz[2] = cell.z;
		// }
		// if(minT>cell.var[controls.T]){
			// minT = cell.var[controls.T];
			// minT_xyz[0] = cell.x; minT_xyz[1] = cell.y; minT_xyz[2] = cell.z;
		// }
		
		
		// normDelPrim += pow(resiVar[0][step_loc*0 + i_loc],2.0);
		// normDelPrim += pow(resiVar[1][step_loc*0 + i_loc],2.0);
		// normDelPrim += pow(resiVar[2][step_loc*0 + i_loc],2.0);
		
		// normPrim += pow(cell.var[controls.P],2.0);
		// normPrim += pow(cell.var[controls.U],2.0);
		// normPrim += pow(cell.var[controls.V],2.0);
		// normPrim += pow(cell.var[controls.W],2.0);
		// normPrim += pow(cell.var[controls.T],2.0);
		// normPrim += pow(cell.var[controls.MF[0]],2.0);
	// }
	

	
	// vector<double> maxP_glo(size,0.0);
	// vector<double> maxP_x_glo(size,0.0);
	// vector<double> maxP_y_glo(size,0.0);
	// vector<double> maxP_z_glo(size,0.0);
	
	// vector<double> minP_glo(size,0.0);
	// vector<double> minP_x_glo(size,0.0);
	// vector<double> minP_y_glo(size,0.0);
	// vector<double> minP_z_glo(size,0.0);
	
	// vector<double> maxUn_glo(size,0.0);
	// vector<double> maxUn_x_glo(size,0.0);
	// vector<double> maxUn_y_glo(size,0.0);
	// vector<double> maxUn_z_glo(size,0.0);
	
	// vector<double> maxT_glo(size,0.0);
	// vector<double> maxT_x_glo(size,0.0);
	// vector<double> maxT_y_glo(size,0.0);
	// vector<double> maxT_z_glo(size,0.0);
	
	// vector<double> minT_glo(size,0.0);
	// vector<double> minT_x_glo(size,0.0);
	// vector<double> minT_y_glo(size,0.0);
	// vector<double> minT_z_glo(size,0.0);
	
	// MPI_Allgather(&maxP, 1, MPI_DOUBLE, maxP_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[0], 1, MPI_DOUBLE, maxP_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[1], 1, MPI_DOUBLE, maxP_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxP_xyz[2], 1, MPI_DOUBLE, maxP_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// MPI_Allgather(&minP, 1, MPI_DOUBLE, minP_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minP_xyz[0], 1, MPI_DOUBLE, minP_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minP_xyz[1], 1, MPI_DOUBLE, minP_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minP_xyz[2], 1, MPI_DOUBLE, minP_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// MPI_Allgather(&maxUn, 1, MPI_DOUBLE, maxUn_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[0], 1, MPI_DOUBLE, maxUn_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[1], 1, MPI_DOUBLE, maxUn_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxUn_xyz[2], 1, MPI_DOUBLE, maxUn_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// MPI_Allgather(&maxT, 1, MPI_DOUBLE, maxT_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxT_xyz[0], 1, MPI_DOUBLE, maxT_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxT_xyz[1], 1, MPI_DOUBLE, maxT_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&maxT_xyz[2], 1, MPI_DOUBLE, maxT_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	
	// MPI_Allgather(&minT, 1, MPI_DOUBLE, minT_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minT_xyz[0], 1, MPI_DOUBLE, minT_x_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minT_xyz[1], 1, MPI_DOUBLE, minT_y_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	// MPI_Allgather(&minT_xyz[2], 1, MPI_DOUBLE, minT_z_glo.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
	

	
	// if(rank==0){
		// maxP = maxP_glo[0];
		// maxP_xyz[0] = maxP_x_glo[0]; maxP_xyz[1] = maxP_y_glo[0]; maxP_xyz[2] = maxP_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(maxP<maxP_glo[i]){
				// maxP = maxP_glo[i];
				// maxP_xyz[0] = maxP_x_glo[i]; maxP_xyz[1] = maxP_y_glo[i]; maxP_xyz[2] = maxP_z_glo[i];
			// }
		// }
		// minP = minP_glo[0];
		// minP_xyz[0] = minP_x_glo[0]; minP_xyz[1] = minP_y_glo[0]; minP_xyz[2] = minP_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(minP>minP_glo[i]){
				// minP = minP_glo[i];
				// minP_xyz[0] = minP_x_glo[i]; minP_xyz[1] = minP_y_glo[i]; minP_xyz[2] = minP_z_glo[i];
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
		// maxT = maxT_glo[0];
		// maxT_xyz[0] = maxT_x_glo[0]; maxT_xyz[1] = maxT_y_glo[0]; maxT_xyz[2] = maxT_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(maxT<maxT_glo[i]){
				// maxT = maxT_glo[i];
				// maxT_xyz[0] = maxT_x_glo[i]; maxT_xyz[1] = maxT_y_glo[i]; maxT_xyz[2] = maxT_z_glo[i];
			// }
		// }
		// minT = minT_glo[0];
		// minT_xyz[0] = minT_x_glo[0]; minT_xyz[1] = minT_y_glo[0]; minT_xyz[2] = minT_z_glo[0];
		// for(int i=1; i<size; ++i){
			// if(minT>minT_glo[i]){
				// minT = minT_glo[i];
				// minT_xyz[0] = minT_x_glo[i]; minT_xyz[1] = minT_y_glo[i]; minT_xyz[2] = minT_z_glo[i];
			// }
		// }
		
		
		// // double resi_print = sqrt(normDelPrim)/sqrt(normPrim)/(double)B_n;
		// double resi_print = sqrt(normDelPrim)/(double)B_n;
		
		// cout << endl;
		// cout << " └ maximum P  = " << maxP << ", locations = (" << maxP_xyz[0] << " " << maxP_xyz[1] << " " << maxP_xyz[2] << ")" << endl;
		// cout << " └ minimum P  = " << minP << ", locations = (" << minP_xyz[0] << " " << minP_xyz[1] << " " << minP_xyz[2] << ")" << endl;
		// cout << " └ maximum Un = " << maxUn << ", locations = (" << maxUn_xyz[0] << " " << maxUn_xyz[1] << " " << maxUn_xyz[2] << ")" << endl;
		// cout << " └ maximum T  = " << maxT << ", locations = (" << maxT_xyz[0] << " " << maxT_xyz[1] << " " << maxT_xyz[2] << ")" << endl;
		// cout << " └ minimum T  = " << minT << ", locations = (" << minT_xyz[0] << " " << minT_xyz[1] << " " << minT_xyz[2] << ")" << endl;
		// cout << " └ total residuals = " << resi_print << endl;
		// cout << endl;
		
	// }
	

	// return sqrt(normDelPrim)/sqrt(normPrim)/(double)B_n;
	
	
}

