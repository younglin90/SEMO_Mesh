#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include <ctime>

//**********************************************
//
// avoid ill-conditioned matrix
//
//**********************************************
double SEMO_Solvers_Builder::calcCoupledEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 4;
    int A_n = B_n * B_n;
	
	
	// vector<clock_t> startTime;
	
	SEMO_Utility_Math math;
	
	int proc_num=0;

	// gradient P
	vector<vector<double>> gradP;
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
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
	
	
	// F. Denner, 2013
	vector<double> A_p_vals(mesh.cells.size(), 0.0);
	proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		int ijStartL_local = B_n*(face.owner) - 1;
		int ijStartR_local = B_n*(face.neighbour) - 1;
		
        int ijStartL = B_n*mesh.startCellGlobal + ijStartL_local;
        int ijStartR = B_n*mesh.startCellGlobal + ijStartR_local;
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			ijStartR = B_n*mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] + 
					   B_n*(mesh.procNeighbCellNo[proc_num]) - 1;
		}

		
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
		// double wCL = face.wVC;
		// double wCL = 0.5;
		double wCR = 1.0-wCL;
		
		double wCPL = 0.5;
		double wCPR = 0.5;
		
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
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
		double tmp1 = controls.timeStep / Rho_star;
			
		double dPN = face.magPN;
				
		double orgPL = mesh.cells[face.owner].var[controls.P];
		double orgPR = 0.0;
		vector<double> gradPL(3,0.0);
		gradPL[0] = gradP[face.owner][0];
		gradPL[1] = gradP[face.owner][1];
		gradPL[2] = gradP[face.owner][2];
		
		vector<double> gradPR(3,0.0);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			orgPR = mesh.cells[face.neighbour].var[controls.P];
			gradPR[0] = gradP[face.neighbour][0];
			gradPR[1] = gradP[face.neighbour][1];
			gradPR[2] = gradP[face.neighbour][2];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			orgPR = PR;
			gradPR[0] = gradPx_recv[proc_num];
			gradPR[1] = gradPy_recv[proc_num];
			gradPR[2] = gradPz_recv[proc_num];
		}
		
		
		for(int ii=0; ii<3; ++ii){
			UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*face.vecPN[ii]/dPN;
			UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*face.vecPN[ii]/dPN;
			// UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			// UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
		}
			
		UnF -= tmp1*(orgPR-orgPL)/dPN;
		// for(int ii=0; ii<3; ++ii){
			// UnF -= tmp1 * (gradPR[ii]*face.vecNdN[ii] - gradPL[ii]*face.vecPdP[ii])/dPN;
		// }
		
		
			 
			 
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = 1.0 - weightL;
		
		double RhoF = 
			weightL*face.varL[controls.fRho] +
			weightR*face.varR[controls.fRho];
		
		A_p_vals[face.owner] += weightL * RhoF*UnF*face.area / mesh.cells[face.owner].volume;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			A_p_vals[face.neighbour] -=  weightR * RhoF*UnF*face.area / mesh.cells[face.neighbour].volume;
		}
		
		
	}
	
	
	vector<double> A_p_vals_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> A_p_vals_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				A_p_vals_send.push_back(A_p_vals[face.owner]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					A_p_vals_send, A_p_vals_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		A_p_vals_send.clear();
	}



	// diagonal terms
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
		
		double Rho = cell.var[controls.Rho];
		double U = cell.var[controls.U];
		double V = cell.var[controls.V];
		double W = cell.var[controls.W];
		double oldRho = cell.var[controls.oldRho];
		double oldU = cell.var[controls.oldU];
		double oldV = cell.var[controls.oldV];
		double oldW = cell.var[controls.oldW];
		double volume = cell.volume;
		double timeStep = controls.timeStep;
		
		int id = 0;
		int i_glo = 0;
		int j_glo = 0;
		int i_loc = i;
		int str_glo = mesh.startCellGlobal*B_n;
		int step_loc = mesh.cells.size();
		int step_glo = mesh.ncellsTotal;

        // continuity
        id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        
        id = step_loc*(B_n*0+1) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*1 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        
        id = step_loc*(B_n*0+2) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*2 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        
        id = step_loc*(B_n*0+3) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*3 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        // x-momentum
        id = step_loc*(B_n*1+0) + i_loc; i_glo = str_glo + step_loc*1 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        id = step_loc*(B_n*1+1) + i_loc; i_glo = str_glo + step_loc*1 + i_loc; j_glo = str_glo + step_loc*1 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        A_vals[id] = Rho*volume/timeStep;
		
        id = step_loc*(B_n*1+2) + i_loc; i_glo = str_glo + step_loc*1 + i_loc; j_glo = str_glo + step_loc*2 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        
        id = step_loc*(B_n*1+3) + i_loc; i_glo = str_glo + step_loc*1 + i_loc; j_glo = str_glo + step_loc*3 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        // y-momentum
        id = step_loc*(B_n*2+0) + i_loc; i_glo = str_glo + step_loc*2 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        id = step_loc*(B_n*2+1) + i_loc; i_glo = str_glo + step_loc*2 + i_loc; j_glo = str_glo + step_loc*1 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
		
        id = step_loc*(B_n*2+2) + i_loc; i_glo = str_glo + step_loc*2 + i_loc; j_glo = str_glo + step_loc*2 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        A_vals[id] = Rho*volume/timeStep;
        
        id = step_loc*(B_n*2+3) + i_loc; i_glo = str_glo + step_loc*2 + i_loc; j_glo = str_glo + step_loc*3 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        // z-momentum
        id = step_loc*(B_n*3+0) + i_loc; i_glo = str_glo + step_loc*3 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;

        id = step_loc*(B_n*3+1) + i_loc; i_glo = str_glo + step_loc*3 + i_loc; j_glo = str_glo + step_loc*1 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
		
        id = step_loc*(B_n*3+2) + i_loc; i_glo = str_glo + step_loc*3 + i_loc; j_glo = str_glo + step_loc*2 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        
        id = step_loc*(B_n*3+3) + i_loc; i_glo = str_glo + step_loc*3 + i_loc; j_glo = str_glo + step_loc*3 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        A_vals[id] = Rho*volume/timeStep;


        
        B[step_loc*0 + i_loc] = 0.0;
        B[step_loc*1 + i_loc] = -Rho*(U - oldU)*volume / timeStep;
        B[step_loc*2 + i_loc] = -Rho*(V - oldV)*volume / timeStep;
        B[step_loc*3 + i_loc] = -Rho*(W - oldW)*volume / timeStep;
		
		
		// // gravity force terms
        // B[step_loc*1 + i_loc] += Rho*controls.gravityAcceleration[0]*volume;
        // B[step_loc*2 + i_loc] += Rho*controls.gravityAcceleration[1]*volume;
        // B[step_loc*3 + i_loc] += Rho*controls.gravityAcceleration[2]*volume;


		// // surface tension force terms
        // B[ijStart + 2] += volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
        // B[ijStart + 3] += volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
        // B[ijStart + 4] += volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
		
	}
	
	
	
	
	proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		int ijStartL_local = B_n*(face.owner) - 1;
		int ijStartR_local = B_n*(face.neighbour) - 1;
		
        int ijStartL = B_n*mesh.startCellGlobal + ijStartL_local;
        int ijStartR = B_n*mesh.startCellGlobal + ijStartR_local;
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			ijStartR = B_n*mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] + 
					   B_n*(mesh.procNeighbCellNo[proc_num]) - 1;
		}

		
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
		// double wCL = face.wVC;
		// double wCL = 0.5;
		double wCR = 1.0-wCL;
		
		double wCPL = wCL;
		double wCPR = wCR;
		
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
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
		double tmp1 = controls.timeStep / Rho_star;
		
		// // F.Denner et al., 2013
		// double d_F = 0.0;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// tmp1 = wCL*(A_p_vals[face.owner] + RhoL/controls.timeStep);
			// tmp1 += wCR*(A_p_vals[face.neighbour] + RhoR/controls.timeStep);
		// }
		// else{
			// tmp1 = wCL*(A_p_vals[face.owner] + RhoL/controls.timeStep);
			// tmp1 += wCR*(A_p_vals_recv[proc_num] + RhoR/controls.timeStep);
		// }
		// tmp1 = 1.0/tmp1;
		
		
		double dPN = face.magPN;
		// dPN = sqrt(pow(distanceCells[0],2.0) + 
					  // pow(distanceCells[1],2.0) + 
					  // pow(distanceCells[2],2.0));
				
		double orgPL = mesh.cells[face.owner].var[controls.P];
		double orgPR = 0.0;
		vector<double> gradPL(3,0.0);
		gradPL[0] = gradP[face.owner][0];
		gradPL[1] = gradP[face.owner][1];
		gradPL[2] = gradP[face.owner][2];
		
		vector<double> gradPR(3,0.0);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			orgPR = mesh.cells[face.neighbour].var[controls.P];
			gradPR[0] = gradP[face.neighbour][0];
			gradPR[1] = gradP[face.neighbour][1];
			gradPR[2] = gradP[face.neighbour][2];
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			orgPR = PR;
			gradPR[0] = gradPx_recv[proc_num];
			gradPR[1] = gradPy_recv[proc_num];
			gradPR[2] = gradPz_recv[proc_num];
		}
		
		
		for(int ii=0; ii<3; ++ii){
			UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*face.vecPN[ii]/dPN;
			UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*face.vecPN[ii]/dPN;
			// UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			// UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
		}
			
		UnF -= tmp1*(orgPR-orgPL)/dPN;
		// for(int ii=0; ii<3; ++ii){
			// UnF -= tmp1 * (gradPR[ii]*face.vecNdN[ii] - gradPL[ii]*face.vecPdP[ii])/dPN;
		// }
		
		
 
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = 1.0 - weightL;
		
		double UF = weightL*UL + weightR*UR;
		double VF = weightL*VL + weightR*VR;
		double WF = weightL*WL + weightR*WR;
		double RhoF = weightL*RhoL + weightR*RhoR;
		double PF = wCPL*PL + wCPR*PR;
		
		
		
		int str_glo_L = mesh.startCellGlobal*B_n;
		int str_glo_R = mesh.startCellGlobal*B_n;
		int step_loc_L = mesh.cells.size();
		int step_loc_R = mesh.cells.size();
		int i_loc_L = face.owner;
		int i_loc_R = face.neighbour;
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			str_glo_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]]*B_n;
			step_loc_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]+1] - 
			             mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
			i_loc_R = mesh.procNeighbCellNo[proc_num];
		}
		
		int id = -1;
		int i_glo = -1;
		int j_glo = -1;
		double conv_flux = 0.0;
		double diff_flux = 0.0;

        //------------------------
        // continuity
        // p'
		id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
		A_vals[id] += ( + tmp1 / dPN * area );
		A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		A_vals.push_back(( - tmp1 / dPN * area ));
		
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( - tmp1 / dPN * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( + tmp1 / dPN * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
			
        
        // // u'
        // id = step_loc_L*(B_n*0+1) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*1 + i_loc_R;
        // A_vals[id] += ( wCL * nvec[0] * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( wCR * nvec[0] * area ));
		
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*0+1) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*1 + i_loc_L;
			// A_vals[id] -= ( wCR * nvec[0] * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( wCL * nvec[0] * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
        
        // // v'
        // id = step_loc_L*(B_n*0+2) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*2 + i_loc_R;
        // A_vals[id] += ( wCL * nvec[1] * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( wCR * nvec[1] * area ));
		
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*0+2) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*2 + i_loc_L;
			// A_vals[id] -= ( wCR * nvec[1] * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( wCL * nvec[1] * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
        
        // // w'
        // id = step_loc_L*(B_n*0+3) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*3 + i_loc_R;
        // A_vals[id] += ( wCL * nvec[2] * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( wCR * nvec[2] * area ));
		
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*0+3) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*3 + i_loc_L;
			// A_vals[id] -= ( wCR * nvec[2] * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( wCL * nvec[2] * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
		
        

        //------------------------
        // x-momentum

        
        // p'
        id = step_loc_L*(B_n*1+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*1 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
        A_vals[id] += ( wCPL * nvec[0] * area + RhoL * tmp1 / dPN * UF * area );
        A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        A_vals.push_back(( wCPR * nvec[0] * area - RhoL * tmp1 / dPN * UF * area ));
	
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*1+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*1 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( wCPR * nvec[0] * area - RhoR * tmp1 / dPN * UF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( wCPL * nvec[0] * area + RhoR * tmp1 / dPN * UF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
        // u'
		id = step_loc_L*(B_n*1+1) + i_loc_L; i_glo = str_glo_L + step_loc_L*1 + i_loc_L; j_glo = str_glo_R + step_loc_R*1 + i_loc_R;
		// A_vals[id] += ( RhoL * wCL * nvec[0] * UF * area + weightL * RhoL * UnF * area );
		A_vals[id] += ( weightL * RhoL * UnF * area );
		A_rows.push_back(i_glo); A_cols.push_back(j_glo);
		// A_vals.push_back(( RhoL * wCR * nvec[0] * UF * area + weightR * RhoL * UnF * area ));
		A_vals.push_back(( weightR * RhoL * UnF * area ));
			
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*1+1) + i_loc_R; i_glo = str_glo_R + step_loc_R*1 + i_loc_R; j_glo = str_glo_L + step_loc_L*1 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[0] * UF * area + weightR * RhoR * UnF * area );
			A_vals[id] -= ( weightR * RhoR * UnF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[0] * UF * area + weightL * RhoR * UnF * area ));
			A_vals.push_back(-( weightL * RhoR * UnF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
		
        // // v'
        // id = step_loc_L*(B_n*1+2) + i_loc_L; i_glo = str_glo_L + step_loc_L*1 + i_loc_L; j_glo = str_glo_R + step_loc_R*2 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[1] * UF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * UF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*1+2) + i_loc_R; i_glo = str_glo_R + step_loc_R*1 + i_loc_R; j_glo = str_glo_L + step_loc_L*2 + i_loc_L;
			// A_vals[id] -= (( RhoR * wCR * nvec[1] * UF * area ));
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[1] * UF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
        
		
        // // w'
        // id = step_loc_L*(B_n*1+3) + i_loc_L; i_glo = str_glo_L + step_loc_L*1 + i_loc_L; j_glo = str_glo_R + step_loc_R*3 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[2] * UF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * UF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*1+3) + i_loc_R; i_glo = str_glo_R + step_loc_R*1 + i_loc_R; j_glo = str_glo_L + step_loc_L*3 + i_loc_L;
			// A_vals[id] -= (( RhoR * wCR * nvec[2] * UF * area ));
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[2] * UF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		

		
		
        //------------------------
        // y-momentum
        
		
        // p'
        id = step_loc_L*(B_n*2+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*2 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
        A_vals[id] += ( wCPL * nvec[1] * area + RhoL * tmp1 / dPN * VF * area );
        A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        A_vals.push_back(( wCPR * nvec[1] * area - RhoL * tmp1 / dPN * VF * area ));
		
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*2+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*2 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( wCPR * nvec[1] * area - RhoR * tmp1 / dPN * VF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( wCPL * nvec[1] * area + RhoR * tmp1 / dPN * VF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
		
		
        // // u'
        // id = step_loc_L*(B_n*2+1) + i_loc_L; i_glo = str_glo_L + step_loc_L*2 + i_loc_L; j_glo = str_glo_R + step_loc_R*1 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[0] * VF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[0] * VF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*2+1) + i_loc_R; i_glo = str_glo_R + step_loc_R*2 + i_loc_R; j_glo = str_glo_L + step_loc_L*1 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[0] * VF * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[0] * VF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
		
         
        // v'
        
        id = step_loc_L*(B_n*2+2) + i_loc_L; i_glo = str_glo_L + step_loc_L*2 + i_loc_L; j_glo = str_glo_R + step_loc_R*2 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[1] * VF * area + weightL * RhoL * UnF * area );
        A_vals[id] += ( weightL * RhoL * UnF * area );
        A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * VF * area + weightR * RhoL * UnF * area ));
        A_vals.push_back(( weightR * RhoL * UnF * area ));
		
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*2+2) + i_loc_R; i_glo = str_glo_R + step_loc_R*2 + i_loc_R; j_glo = str_glo_L + step_loc_L*2 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[1] * VF * area + weightR * RhoR * UnF * area );
			A_vals[id] -= ( weightR * RhoR * UnF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[1] * VF * area + weightL * RhoR * UnF * area ));
			A_vals.push_back(-( weightL * RhoR * UnF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
		
        
		
        // // w'
        // id = step_loc_L*(B_n*2+3) + i_loc_L; i_glo = str_glo_L + step_loc_L*2 + i_loc_L; j_glo = str_glo_R + step_loc_R*3 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[2] * VF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * VF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*2+3) + i_loc_R; i_glo = str_glo_R + step_loc_R*2 + i_loc_R; j_glo = str_glo_L + step_loc_L*3 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[2] * VF * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[2] * VF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
		
		
		
		
        // ------------------------
        // z-momentum
        
		
        // p'
        id = step_loc_L*(B_n*3+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*3 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
        A_vals[id] += ( wCPL * nvec[2] * area + RhoL * tmp1 / dPN * WF * area );
        A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        A_vals.push_back(( wCPR * nvec[2] * area - RhoL * tmp1 / dPN * WF * area ));
		
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*3+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*3 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
			A_vals[id] -= ( wCPR * nvec[2] * area - RhoR * tmp1 / dPN * WF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(-( wCPL * nvec[2] * area + RhoR * tmp1 / dPN * WF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
        
		
        // // u'
        // id = step_loc_L*(B_n*3+1) + i_loc_L; i_glo = str_glo_L + step_loc_L*3 + i_loc_L; j_glo = str_glo_R + step_loc_R*1 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[0] * WF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[0] * WF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*3+1) + i_loc_R; i_glo = str_glo_R + step_loc_R*3 + i_loc_R; j_glo = str_glo_L + step_loc_L*1 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[0] * WF * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[0] * WF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		

        
		
        // // v'
        // id = step_loc_L*(B_n*3+2) + i_loc_L; i_glo = str_glo_L + step_loc_L*3 + i_loc_L; j_glo = str_glo_R + step_loc_R*2 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[1] * WF * area );
        // A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * WF * area ));
			
		// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		// }

		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			// id = step_loc_R*(B_n*3+2) + i_loc_R; i_glo = str_glo_R + step_loc_R*3 + i_loc_R; j_glo = str_glo_L + step_loc_L*2 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[1] * WF * area );
			// A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[1] * WF * area ));
			
			// if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				// A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			// }
		// }
		
        
        // w'
        
        id = step_loc_L*(B_n*3+3) + i_loc_L; i_glo = str_glo_L + step_loc_L*3 + i_loc_L; j_glo = str_glo_R + step_loc_R*3 + i_loc_R;
        // A_vals[id] += ( RhoL * wCL * nvec[2] * WF * area + weightL * RhoL * UnF * area );
        A_vals[id] += ( weightL * RhoL * UnF * area );
        A_rows.push_back(i_glo); A_cols.push_back(j_glo);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * WF * area + weightR * RhoL * UnF * area ));
        A_vals.push_back(( weightR * RhoL * UnF * area ));
			
		if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
			A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
		}

		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			id = step_loc_R*(B_n*3+3) + i_loc_R; i_glo = str_glo_R + step_loc_R*3 + i_loc_R; j_glo = str_glo_L + step_loc_L*3 + i_loc_L;
			// A_vals[id] -= ( RhoR * wCR * nvec[2] * WF * area + weightR * RhoR * UnF * area );
			A_vals[id] -= ( weightR * RhoR * UnF * area );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(-( RhoR * wCL * nvec[2] * WF * area + weightL * RhoR * UnF * area ));
			A_vals.push_back(-( weightL * RhoR * UnF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}
		}
		
		
		
		
		
        // ----------------------------
        // B
        B[step_loc_L*0 + i_loc_L] -= ( UnF * area );
        B[step_loc_L*1 + i_loc_L] -= ( RhoL * UF * UnF * area + PF * nvec[0] * area );
        B[step_loc_L*2 + i_loc_L] -= ( RhoL * VF * UnF * area + PF * nvec[1] * area );
        B[step_loc_L*3 + i_loc_L] -= ( RhoL * WF * UnF * area + PF * nvec[2] * area );

		
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			B[step_loc_R*0 + i_loc_R] += ( UnF * area );
			B[step_loc_R*1 + i_loc_R] += ( RhoR * UF * UnF * area + PF * nvec[0] * area );
			B[step_loc_R*2 + i_loc_R] += ( RhoR * VF * UnF * area + PF * nvec[1] * area );
			B[step_loc_R*3 + i_loc_R] += ( RhoR * WF * UnF * area + PF * nvec[2] * area );
		}
		
		
		
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
		
	}
	
	
	
	
	
	
	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				
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
				
				double dPN = face.magPN;
					
				// for(int ii=0; ii<3; ++ii){
					// UnF += tmp1 * gradP[face.owner][ii]*nvec[ii];
				// }
				
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
				
				// velocity coefficient
				double coeffUVW = 0.0;
				if( boundary.type[controls.U] == "fixedValue" ){
					coeffUVW = 0.0;
					// coeffUVW = 0.5;
				}
				else if( boundary.type[controls.U] == "zeroGradient" ){
					coeffUVW = 1.0;
					// coeffUVW = 0.5;
				}
				else if( boundary.type[controls.U] == "slip" ){
					coeffUVW = 0.0;
				}
				else if( boundary.type[controls.U] == "noSlip" ){
					coeffUVW = 0.0;
				}
				else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
					coeffUVW = 0.0;
					// coeffUVW = 0.5;
				}
				else if( boundary.type[controls.U] == "inletOutlet" ){
					double ownNorVel =  
						mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					coeffUVW = (ownNorVel > 1.0) ? 1.0 : 0.0;
				}
				else {
					cerr << "| #Error : not defined B.C., var = " << controls.U << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
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
				
				
				
				
				
				
				// if(abs(UnF)>0.0){
					// UnF -= tmp1 * (face.varR[controls.fP]-face.varL[controls.fP])/dPN;
					// UnF += tmp1 * gradP[face.owner][0]*nvec[0];
					// UnF += tmp1 * gradP[face.owner][1]*nvec[1];
					// UnF += tmp1 * gradP[face.owner][2]*nvec[2];
				// }
				
				
				
				
				

				int str_glo = mesh.startCellGlobal*B_n;
				int step_loc = mesh.cells.size();
				int step_glo = mesh.ncellsTotal;
				int i_loc_L = face.owner;
				int id=-1;
				
				
				// coeffP = 50.0 * coeffP;
				// tmp1 = 1000.0 * tmp1;
				
				
				// continuity
				id = step_loc*(B_n*0+0) + i_loc_L;
				A_vals[id] += coeffP_diff * ( tmp1 / dPN * area );
				
				// id = step_loc*(B_n*0+1) + i_loc_L;
				// A_vals[id] += coeffUVW *( nvec[0] * area );
				
				// id = step_loc*(B_n*0+2) + i_loc_L;
				// A_vals[id] += coeffUVW *( nvec[1] * area );
				
				// id = step_loc*(B_n*0+3) + i_loc_L;
				// A_vals[id] += coeffUVW *( nvec[2] * area );
				
				
				
				// x-momentum
				id = step_loc*(B_n*1+0) + i_loc_L;
				A_vals[id] += ( coeffP * nvec[0] * area + coeffP_diff * RhoL * tmp1 / dPN * UF * area );
				// A_vals[id] += ( coeffP * nvec[0] * area );
				
				id = step_loc*(B_n*1+1) + i_loc_L;
				// A_vals[id] += coeffUVW * ( RhoL * nvec[0] * UF * area + RhoL * UnF * area );
				A_vals[id] += coeffUVW * ( RhoL * UnF * area );
				
				
				// id = step_loc*(B_n*1+2) + i_loc_L;
				// A_vals[id] += coeffUVW * ( RhoL * nvec[1] * UF * area );
				
				
				
				// id = step_loc*(B_n*1+3) + i_loc_L;
				// A_vals[id] += coeffUVW * ( RhoL * nvec[2] * UF * area );
				
				
				
				// y-momentum
				id = step_loc*(B_n*2+0) + i_loc_L;
				A_vals[id] += ( coeffP * nvec[1] * area + coeffP_diff * RhoL * tmp1 / dPN * VF * area );
				// A_vals[id] += ( coeffP * nvec[1] * area );
				
				
				// id = step_loc*(B_n*2+1) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[0] * VF * area );
				
				
				
				id = step_loc*(B_n*2+2) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[1] * VF * area + RhoL * UnF * area );
				A_vals[id] += coeffUVW *( RhoL * UnF * area );
				
				
				
				// id = step_loc*(B_n*2+3) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[2] * VF * area );
				
				

				
				// z-momentum
				id = step_loc*(B_n*3+0) + i_loc_L;
				A_vals[id] += ( coeffP * nvec[2] * area + coeffP_diff * RhoL * tmp1 / dPN * WF * area );
				// A_vals[id] += ( coeffP * nvec[2] * area );
				
				
				// id = step_loc*(B_n*3+1) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[0] * WF * area );
				
				
				
				// id = step_loc*(B_n*3+2) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[1] * WF * area );
				
				
				
				id = step_loc*(B_n*3+3) + i_loc_L;
				// A_vals[id] += coeffUVW *( RhoL * nvec[2] * WF * area + RhoL * UnF * area );
				A_vals[id] += coeffUVW *( RhoL * UnF * area );
				
				
				
				// convective term
				B[step_loc*0 + i_loc_L] -= ( UnF * area );
				B[step_loc*1 + i_loc_L] -= ( RhoL * UF * UnF * area + PF * face.unitNormals[0] * area );
				B[step_loc*2 + i_loc_L] -= ( RhoL * VF * UnF * area + PF * face.unitNormals[1] * area );
				B[step_loc*3 + i_loc_L] -= ( RhoL * WF * UnF * area + PF * face.unitNormals[2] * area );
			}
			
			
			// if(boundary.type[0] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double area = face.area;
					
					// double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 // face.varL[controls.fV]*face.unitNormals[1] + 
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 // face.varR[controls.fV]*face.unitNormals[1] + 
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							       // face.unitNormals[1]*face.distCells[1] + 
								   // face.unitNormals[2]*face.distCells[2];
								
					// double UnF = 0.5*UnL+0.5*UnR;
				
					// double RhoL = face.varL[controls.fRho];
					// double tmp1 = controls.timeStep/RhoL;
					
					// linA[face.owner] += tmp1*area/dPN_e;
					
					// linB[face.owner] -= UnF*area;
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double area = face.area;
					// double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 // face.varL[controls.fV]*face.unitNormals[1] + 
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 // face.varR[controls.fV]*face.unitNormals[1] + 
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							       // face.unitNormals[1]*face.distCells[1] + 
								   // face.unitNormals[2]*face.distCells[2];
								
					// double UnF = 0.5*UnL+0.5*UnR;
				
					// double RhoL = face.varL[controls.fRho];
					
					// // convective term
					// linB[face.owner] -= UnF*area;
				// }
			// }
			
		}
		
	}
	
	
	
	
	
	// linear solver : PETSc library
	vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	
	string solver = "gmres";
	double relTol = 0.0;
	double tolerance = 0.0;
	string preconditioner = "mg";
	int maxIter = 50;
	
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
	
	
	
	
	// solvePETSc_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
		// B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
		// "lgmres", 1.e-8, 1.e-6, preconditioner, controls.maxIterP);
		
	// solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
		// solver, tolerance, relTol, preconditioner, maxIter);
		
	solveAMGCL("coupled", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);
	
	
	
	// update
	// double relaxP = 0.2;
	// double relaxUVW = 0.2;
	double relaxP = 0.7; //controls.prePreURF;
	double relaxUVW = 0.7; //controls.momVelURF;
	
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
		
		// double save_relaxP = relaxP;
		// double save_relaxUVW = relaxUVW;
		// if(abs(resiVar[step_loc*0 + i_loc]) > abs(cell.var[controls.P])*1.e-1 ){
			// relaxP = 1.e-5;
		// }
		// if(
		// abs(resiVar[step_loc*1 + i_loc]) > abs(cell.var[controls.U])*1.e-1 ||
		// abs(resiVar[step_loc*2 + i_loc]) > abs(cell.var[controls.V])*1.e-1 ||
		// abs(resiVar[step_loc*3 + i_loc]) > abs(cell.var[controls.W])*1.e-1
		// ){
			// relaxUVW = 1.e-5;
		// }
		
		cell.var[controls.P] += relaxP * resiVar[step_loc*0 + i_loc];
		cell.var[controls.U] += relaxUVW * resiVar[step_loc*1 + i_loc];
		cell.var[controls.V] += relaxUVW * resiVar[step_loc*2 + i_loc];
		cell.var[controls.W] += relaxUVW * resiVar[step_loc*3 + i_loc];
		
		
		// relaxP = save_relaxP;
		// relaxUVW = save_relaxUVW ;
		
		
		
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
		
		
		normDelPrim += pow(resiVar[step_loc*0 + i_loc],2.0);
		normDelPrim += pow(resiVar[step_loc*1 + i_loc],2.0);
		normDelPrim += pow(resiVar[step_loc*2 + i_loc],2.0);
		normDelPrim += pow(resiVar[step_loc*3 + i_loc],2.0);
		
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
	

	// for(int i=0; i<size; ++i){
		// if( maxUn_x_glo[i] < 0.000465 && maxUn_x_glo[i] > 0.000462 &&
		// maxUn_y_glo[i] < -0.0680 && maxUn_y_glo[i]> -0.0683  &&
		// maxUn_z_glo[i] < 0.00179 && maxUn_z_glo[i] > 0.00176
		// ){
	
			// cout << maxUn_xyz[0] << " " << maxUn_xyz[1] << " " << maxUn_xyz[2] << " " << endl;
		// }
	// }
	
	
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
		
		cout.precision(3);
		cout << "|  └ Coupled | ";
		cout << scientific << "P max  = " << maxP << ", xyz = " << maxP_xyz[0] << ", " << maxP_xyz[1] << ", " << maxP_xyz[2] << " | ";
		cout << endl;
		cout << "             | ";
		cout << scientific << "Un max = " << maxUn << ", xyz = " << maxUn_xyz[0] << ", " << maxUn_xyz[1] << ", " << maxUn_xyz[2] << " | ";
		cout.unsetf(ios::scientific);
		cout << endl;
		// cout << "-----------------" << endl;
		// cout << " maximum P = " << maxP << ", xyz = " << maxP_xyz[0] << ", " << maxP_xyz[1] << ", " << maxP_xyz[2] << endl;
		// cout << " maximum Un = " << maxUn << ", xyz = " << maxUn_xyz[0] << ", " << maxUn_xyz[1] << ", " << maxUn_xyz[2] << endl;
		// cout << "-----------------" << endl;
	}
	

	// for(auto& cell : mesh.cells){
		// if( cell.x < 0.000465 && cell.x > 0.000462 &&
		// cell.y < -0.0680 && cell.y > -0.0683  &&
		// cell.z < 0.00179 && cell.z > 0.00176
		// ){
			
			// cout << cell.x << " " << cell.y << " " << cell.z << " " << endl;
		// }
		
	// }
	
	
	return sqrt(normDelPrim)/sqrt(normPrim);
	
	
}







// //**********************************************
// //
// // ill - conditioned
// //
// //**********************************************
// double SEMO_Solvers_Builder::calcCoupledEq(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<SEMO_Species>& species
	// ){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	
	// // diagonal terms
    // int B_n = 4;
    // int A_n = B_n * B_n;
	
	
	// // vector<clock_t> startTime;
	
	// SEMO_Utility_Math math;
	
	// int proc_num=0;

	// // gradient P
	// vector<vector<double>> gradP;
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
    // vector<int> A_rows(mesh.cells.size()*A_n, 0);
    // vector<int> A_cols(mesh.cells.size()*A_n, 0);
    // vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
    // vector<double> B(mesh.cells.size()*B_n, 0.0);
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];

        // int ijStart_local = B_n*(i) - 1;
        // int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
        // int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		// double Rho = cell.var[controls.Rho];
		// double U = cell.var[controls.U];
		// double V = cell.var[controls.V];
		// double W = cell.var[controls.W];
		// double oldU = cell.var[controls.oldU];
		// double oldV = cell.var[controls.oldV];
		// double oldW = cell.var[controls.oldW];
		// double volume = cell.volume;
		// double timeStep = controls.timeStep;

        // // continuity
        // id += 1;
        // A_rows[id] = ijStart_global + 1; A_cols[id] = ijStart_global + 1;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 1; A_cols[id] = ijStart_global + 2;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 1; A_cols[id] = ijStart_global + 3;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 1; A_cols[id] = ijStart_global + 4;

        // // x-momentum
        // id += 1;
        // A_rows[id] = ijStart_global + 2; A_cols[id] = ijStart_global + 1;

        // id += 1;
        // A_rows[id] = ijStart_global + 2; A_cols[id] = ijStart_global + 2;
        // A_vals[id] = Rho*volume/timeStep;
		
        // id += 1;
        // A_rows[id] = ijStart_global + 2; A_cols[id] = ijStart_global + 3;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 2; A_cols[id] = ijStart_global + 4;

        // // y-momentum

        // id += 1;
        // A_rows[id] = ijStart_global + 3; A_cols[id] = ijStart_global + 1;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 3; A_cols[id] = ijStart_global + 2;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 3; A_cols[id] = ijStart_global + 3;
        // A_vals[id] = Rho*volume/timeStep;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 3; A_cols[id] = ijStart_global + 4;

        // // z-momentum

        // id += 1;
        // A_rows[id] = ijStart_global + 4; A_cols[id] = ijStart_global + 1;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 4; A_cols[id] = ijStart_global + 2;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 4; A_cols[id] = ijStart_global + 3;
        
        // id += 1;
        // A_rows[id] = ijStart_global + 4; A_cols[id] = ijStart_global + 4;
        // A_vals[id] = Rho*volume/timeStep;


        
        // B[ijStart_local + 1] = 0.0;
        // B[ijStart_local + 2] = -Rho*(U - oldU)*volume / timeStep;
        // B[ijStart_local + 3] = -Rho*(V - oldV)*volume / timeStep;
        // B[ijStart_local + 4] = -Rho*(W - oldW)*volume / timeStep;
		
		
		// // gravity force terms
        // B[ijStart_local + 2] += Rho*controls.gravityAcceleration[0]*volume;
        // B[ijStart_local + 3] += Rho*controls.gravityAcceleration[1]*volume;
        // B[ijStart_local + 4] += Rho*controls.gravityAcceleration[2]*volume;


		// // // surface tension force terms
        // // B[ijStart + 2] += volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
        // // B[ijStart + 3] += volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
        // // B[ijStart + 4] += volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
		
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
				
		// // double nonOrtholimiter = 0.0;
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
		
		// vector<double> Ef(3,0.0);
		// Ef[0] = distanceCells[0]/dPN;
		// Ef[1] = distanceCells[1]/dPN;
		// Ef[2] = distanceCells[2]/dPN;
		// // // over-relaxed
		// // double over_Ef = Ef[0]*nvec[0] + Ef[1]*nvec[1] + Ef[2]*nvec[2];
		// // Ef[0] /= over_Ef; Ef[1] /= over_Ef; Ef[2] /= over_Ef;
		// // double magEf = sqrt(Ef[0]*Ef[0]+Ef[1]*Ef[1]+Ef[2]*Ef[2]);
		// // dPN /= magEf;
		// // dPN = dPN_e;
			
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
		// // tmp1 = d_F_hat;
		
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
		
		// double UF = weightL*UL + weightR*UR;
		// double VF = weightL*VL + weightR*VR;
		// double WF = weightL*WL + weightR*WR;
		// double RhoF = weightL*RhoL + weightR*RhoR;
		// double PF = wCL*PL + wCR*PR;
		
		
        // int iL = A_n*(face.owner) - 1;
        // int iR = A_n*(face.neighbour) - 1;

		
		
		
		

        // //------------------------
        // // continuity
        // // p'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( tmp1 / dPN * area );
        // A_rows.push_back(ijStartL + 1); A_cols.push_back(ijStartR + 1);
        // A_vals.push_back(( - tmp1 / dPN * area ));
        
        // if(iR>-1) A_vals[iR] -= ( - tmp1 / dPN * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 1); if(iR>-1) A_cols.push_back(ijStartL + 1);
        // if(iR>-1) A_vals.push_back(-( tmp1 / dPN * area ));
        
        // // u'
        // iL += 1; iR += 1;
        
        // A_vals[iL] += ( wCL * nvec[0] * area );
        // A_rows.push_back(ijStartL + 1); A_cols.push_back(ijStartR + 2);
        // A_vals.push_back(( wCR * nvec[0] * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[0] * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 1); if(iR>-1) A_cols.push_back(ijStartL + 2);
        // if(iR>-1) A_vals.push_back(-( wCL * nvec[0] * area ));
        
        // // v'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( wCL* nvec[1] * area );
        // A_rows.push_back(ijStartL + 1); A_cols.push_back(ijStartR + 3);
        // A_vals.push_back(( wCR * nvec[1] * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[1] * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 1); if(iR>-1) A_cols.push_back(ijStartL + 3);
        // if(iR>-1) A_vals.push_back(-( wCL * nvec[1] * area ));
        
        // // w'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( wCL * nvec[2] * area );
        // A_rows.push_back(ijStartL + 1); A_cols.push_back(ijStartR + 4);
        // A_vals.push_back(( wCR * nvec[2] * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[2] * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 1); if(iR>-1) A_cols.push_back(ijStartL + 4);
        // if(iR>-1) A_vals.push_back(-( wCL * nvec[2] * area ));
        

        // //------------------------
        // // x-momentum

        // // p'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( wCL * nvec[0] * area + RhoL * tmp1 / dPN * UF * area );
        // A_rows.push_back(ijStartL + 2); A_cols.push_back(ijStartR + 1);
        // A_vals.push_back(( wCR * nvec[0] * area - RhoL * tmp1 / dPN * UF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[0] * area - RhoR * tmp1 / dPN * UF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 2); if(iR>-1) A_cols.push_back(ijStartL + 1);
        // if(iR>-1) A_vals.push_back(-( wCL * nvec[0] * area + RhoR * tmp1 / dPN * UF * area ));
        
        // // u'
        // iL += 1; iR += 1;
        
        // A_vals[iL] += ( RhoL * wCL * nvec[0] * UF * area + weightL * RhoL * UnF * area );
        // A_rows.push_back(ijStartL + 2); A_cols.push_back(ijStartR + 2);
        // A_vals.push_back(( RhoL * wCR * nvec[0] * UF * area + weightR * RhoL * UnF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[0] * UF * area + weightR * RhoR * UnF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 2); if(iR>-1) A_cols.push_back(ijStartL + 2);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[0] * UF * area + weightL * RhoR * UnF * area ));
        
        // // v'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[1] * UF * area );
        // A_rows.push_back(ijStartL + 2); A_cols.push_back(ijStartR + 3);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * UF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[1] * UF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 2); if(iR>-1) A_cols.push_back(ijStartL + 3);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[1] * UF * area ));
        
        // // w'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[2] * UF * area );
        // A_rows.push_back(ijStartL + 2); A_cols.push_back(ijStartR + 4);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * UF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[2] * UF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 2); if(iR>-1) A_cols.push_back(ijStartL + 4);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[2] * UF * area ));
        
		
		
        // //------------------------
        // // y-momentum
        
        // // p'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( wCL * nvec[1] * area + RhoL * tmp1 / dPN * VF * area );
        // A_rows.push_back(ijStartL + 3); A_cols.push_back(ijStartR + 1);
        // A_vals.push_back(( wCR * nvec[1] * area - RhoL * tmp1 / dPN * VF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[1] * area - RhoR * tmp1 / dPN * VF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 3); if(iR>-1) A_cols.push_back(ijStartL + 1);
        // if(iR>-1) A_vals.push_back(-( wCL * nvec[1] * area + RhoR * tmp1 / dPN * VF * area ));
        
        // // u'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[0] * VF * area );
        // A_rows.push_back(ijStartL + 3); A_cols.push_back(ijStartR + 2);
        // A_vals.push_back(( RhoL * wCR * nvec[0] * VF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[0] * VF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 3); if(iR>-1) A_cols.push_back(ijStartL + 2);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[0] * VF * area ));
        
        // // v'
        // iL += 1; iR += 1;
        
        // A_vals[iL] += ( RhoL * wCL * nvec[1] * VF * area + weightL * RhoL * UnF * area );
        // A_rows.push_back(ijStartL + 3); A_cols.push_back(ijStartR + 3);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * VF * area + weightR * RhoL * UnF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[1] * VF * area + weightR * RhoR * UnF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 3); if(iR>-1) A_cols.push_back(ijStartL + 3);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[1] * VF * area + weightL * RhoR * UnF * area ));
        
        // // w'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[2] * VF * area );
        // A_rows.push_back(ijStartL + 3); A_cols.push_back(ijStartR + 4);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * VF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[2] * VF * area );
        // if(iR>-1) A_rows.push_back(ijStartR + 3); if(iR>-1) A_cols.push_back(ijStartL + 4);
        // if(iR>-1) A_vals.push_back(-( RhoR * wCL * nvec[2] * VF * area ));
		
		
		
        // //------------------------
        // // z-momentum
        
        // // p'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( wCL * nvec[2] * area + RhoL * tmp1 / dPN * WF * area );
        // A_rows.push_back(ijStartL + 4); A_cols.push_back(ijStartR + 1);
        // A_vals.push_back(( wCR * nvec[2] * area - RhoL * tmp1 / dPN * WF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( wCR * nvec[2] * area - RhoR * tmp1 / dPN * WF * area );
        // if(iR>-1)A_rows.push_back(ijStartR + 4); if(iR>-1) A_cols.push_back(ijStartL + 1);
        // if(iR>-1)A_vals.push_back(-( wCL * nvec[2] * area + RhoR * tmp1 / dPN * WF * area ));
        
        // // u'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[0] * WF * area );
        // A_rows.push_back(ijStartL + 4); A_cols.push_back(ijStartR + 2);
        // A_vals.push_back(( RhoL * wCR * nvec[0] * WF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[0] * WF * area );
        // if(iR>-1)A_rows.push_back(ijStartR + 4); if(iR>-1) A_cols.push_back(ijStartL + 2);
        // if(iR>-1)A_vals.push_back(-( RhoR * wCL * nvec[0] * WF * area ));
        
        // // v'
        // iL += 1; iR += 1;
        
        // A_vals[iL] += ( RhoL * wCL * nvec[1] * WF * area );
        // A_rows.push_back(ijStartL + 4); A_cols.push_back(ijStartR + 3);
        // A_vals.push_back(( RhoL * wCR * nvec[1] * WF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[1] * WF * area );
        // if(iR>-1)A_rows.push_back(ijStartR + 4); if(iR>-1) A_cols.push_back(ijStartL + 3);
        // if(iR>-1)A_vals.push_back(-( RhoR * wCL * nvec[1] * WF * area ));
        
        // // w'
        // iL += 1; iR += 1;

        // A_vals[iL] += ( RhoL * wCL * nvec[2] * WF * area + weightL * RhoL * UnF * area );
        // A_rows.push_back(ijStartL + 4); A_cols.push_back(ijStartR + 4);
        // A_vals.push_back(( RhoL * wCR * nvec[2] * WF * area + weightR * RhoL * UnF * area ));
        
        // if(iR>-1) A_vals[iR] -= ( RhoR * wCR * nvec[2] * WF * area + weightR * RhoR * UnF * area );
        // if(iR>-1)A_rows.push_back(ijStartR + 4); if(iR>-1) A_cols.push_back(ijStartL + 4);
        // if(iR>-1)A_vals.push_back(-( RhoR * wCL * nvec[2] * WF * area + weightL * RhoR * UnF * area ));
		
		
		
        // // ----------------------------
        // // B
        // B[ijStartL_local + 1] -= ( UnF * area );
        // B[ijStartL_local + 2] -= ( RhoL * UF * UnF * area + PF * nvec[0] * area );
        // B[ijStartL_local + 3] -= ( RhoL * VF * UnF * area + PF * nvec[1] * area );
        // B[ijStartL_local + 4] -= ( RhoL * WF * UnF * area + PF * nvec[2] * area );

		
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// B[ijStartR_local + 1] += ( UnF * area );
			// B[ijStartR_local + 2] += ( RhoR * UF * UnF * area + PF * nvec[0] * area );
			// B[ijStartR_local + 3] += ( RhoR * VF * UnF * area + PF * nvec[1] * area );
			// B[ijStartR_local + 4] += ( RhoR * WF * UnF * area + PF * nvec[2] * area );
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
				
				// int ijStartL_local = B_n*(face.owner) - 1;
				// vector<int> id(B_n,0);
				// id[0] = A_n*(face.owner) - 1;
				// for(int j=1; j<B_n; ++j){
					// id[j] = id[j-1] + 4;
				// }

				
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
				// if( boundary.type[controls.U] == "fixedValue" ){
					// coeffUVW = 0.0;
				// }
				// else if( boundary.type[controls.U] == "zeroGradient" ){
					// coeffUVW = 1.0;
				// }
				// else if( boundary.type[controls.U] == "slip" ){
					// coeffUVW = 0.0;
				// }
				// else if( boundary.type[controls.U] == "noSlip" ){
					// coeffUVW = 0.0;
				// }
				// else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
					// coeffUVW = 0.0;
				// }
				// else if( boundary.type[controls.U] == "inletOutlet" ){
					// double ownNorVel =  
						// mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// coeffUVW = (ownNorVel > 1.0) ? 1.0 : 0.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.U << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				// // volume-fraction coefficient
				// double coeffVF = 0.0;
				// if( boundary.type[controls.VF[0]] == "fixedValue" ){
					// coeffVF = 0.0;
				// }
				// else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
					// coeffVF = 1.0;
					
				// }
				// else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
					// double ownNorVel =  
						// mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// coeffVF = (ownNorVel > 1.0) ? 1.0 : 0.0;
				// }
				// else {
					// cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				
				// // continuity
				// A_vals[id[0]+1] += coeffP_diff * ( tmp1 / dPN * area );
				// A_vals[id[0]+2] += coeffUVW *( nvec[0] * area );
				// A_vals[id[0]+3] += coeffUVW *( nvec[1] * area );
				// A_vals[id[0]+4] += coeffUVW *( nvec[2] * area );

				
				// // x-momentum
				// A_vals[id[1]+1] += ( coeffP * nvec[0] * area + coeffP_diff * RhoL * tmp1 / dPN * UF * area );
				// A_vals[id[1]+2] += coeffUVW * ( RhoL * nvec[0] * UF * area + RhoL * UnF * area );
				// A_vals[id[1]+3] += coeffUVW * ( RhoL * nvec[1] * UF * area );
				// A_vals[id[1]+4] += coeffUVW * ( RhoL * nvec[2] * UF * area );
				
				
				// // y-momentum
				// A_vals[id[2]+1] += ( coeffP * nvec[1] * area + coeffP_diff * RhoL * tmp1 / dPN * VF * area );
				// A_vals[id[2]+2] += coeffUVW *( RhoL * nvec[0] * VF * area );
				// A_vals[id[2]+3] += coeffUVW *( RhoL * nvec[1] * VF * area + RhoL * UnF * area );
				// A_vals[id[2]+4] += coeffUVW *( RhoL * nvec[2] * VF * area );
				

				
				// // z-momentum
				// A_vals[id[3]+1] += ( coeffP * nvec[2] * area + coeffP_diff * RhoL * tmp1 / dPN * WF * area );
				// A_vals[id[3]+2] += coeffUVW *( RhoL * nvec[0] * WF * area );
				// A_vals[id[3]+3] += coeffUVW *( RhoL * nvec[1] * WF * area );
				// A_vals[id[3]+4] += coeffUVW *( RhoL * nvec[2] * WF * area + RhoL * UnF * area );

				
				// // convective term
				// B[ijStartL_local + 1] -= ( UnF * area );
				// B[ijStartL_local + 2] -= ( RhoL * UF * UnF * area + PF * face.unitNormals[0] * area );
				// B[ijStartL_local + 3] -= ( RhoL * VF * UnF * area + PF * face.unitNormals[1] * area );
				// B[ijStartL_local + 4] -= ( RhoL * WF * UnF * area + PF * face.unitNormals[2] * area );
			// }
			
			
			// // if(boundary.type[0] == "fixedValue"){
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double area = face.area;
					
					// // double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 // // face.varL[controls.fV]*face.unitNormals[1] + 
								 // // face.varL[controls.fW]*face.unitNormals[2];
					// // double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 // // face.varR[controls.fV]*face.unitNormals[1] + 
								 // // face.varR[controls.fW]*face.unitNormals[2];
					
					// // double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							       // // face.unitNormals[1]*face.distCells[1] + 
								   // // face.unitNormals[2]*face.distCells[2];
								
					// // double UnF = 0.5*UnL+0.5*UnR;
				
					// // double RhoL = face.varL[controls.fRho];
					// // double tmp1 = controls.timeStep/RhoL;
					
					// // linA[face.owner] += tmp1*area/dPN_e;
					
					// // linB[face.owner] -= UnF*area;
				// // }
			// // }
			// // else{
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double area = face.area;
					// // double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 // // face.varL[controls.fV]*face.unitNormals[1] + 
								 // // face.varL[controls.fW]*face.unitNormals[2];
					// // double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 // // face.varR[controls.fV]*face.unitNormals[1] + 
								 // // face.varR[controls.fW]*face.unitNormals[2];
					
					// // double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							       // // face.unitNormals[1]*face.distCells[1] + 
								   // // face.unitNormals[2]*face.distCells[2];
								
					// // double UnF = 0.5*UnL+0.5*UnR;
				
					// // double RhoL = face.varL[controls.fRho];
					
					// // // convective term
					// // linB[face.owner] -= UnF*area;
				// // }
			// // }
			
		// }
		
	// }
	
	
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	
	// string solver = "gmres";
	// double relTol = 0.0;
	// double tolerance = 0.0;
	// string preconditioner = "mg";
	// int maxIter = 50;
	
	// // if(controls.iterPBs != controls.iterPBsMax-1){
		// // solver = controls.solverP;
		// // relTol = controls.toleranceP;
		// // tolerance = controls.relTolP;
		// // preconditioner = controls.preconditionerP;
		// // maxIter = controls.maxIterP;
	
	// // }
	// // else{
		// // solver = controls.solverFinalP;
		// // relTol = controls.toleranceFinalP;
		// // tolerance = controls.relTolFinalP;
		// // preconditioner = controls.preconditionerFinalP;
		// // maxIter = controls.maxIterFinalP;
	
	// // }
	
	
	
	
	// // solvePETSc_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
		// // B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
		// // "lgmres", 1.e-8, 1.e-6, preconditioner, controls.maxIterP);
		
	// // solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
		// // solver, tolerance, relTol, preconditioner, maxIter);
		
	// solveAMGCL(mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);
	
	
	
	// // update
	// // double relaxP = 0.2;
	// // double relaxUVW = 0.2;
	// double relaxP = controls.prePreURF;
	// double relaxUVW = controls.momVelURF;
	
	// double normPrim = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
        // int ijStart = B_n*(i) - 1;
        // int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		// cell.var[controls.P] += relaxP * resiVar[ijStart+1];
		// cell.var[controls.U] += relaxUVW * resiVar[ijStart+2];
		// cell.var[controls.V] += relaxUVW * resiVar[ijStart+3];
		// cell.var[controls.W] += relaxUVW * resiVar[ijStart+4];
		
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
		
		
		// normPrim += pow(resiVar[ijStart+1],2.0);
		// normPrim += pow(resiVar[ijStart+2],2.0);
		// normPrim += pow(resiVar[ijStart+3],2.0);
		// normPrim += pow(resiVar[ijStart+4],2.0);
		
	// }
	
	
	// return sqrt(normPrim);
	
	
// }






// double SEMO_Solvers_Builder::calcCoupledEq(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<SEMO_Species>& species
	// ){
	

    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	
	// // diagonal terms
    // int B_n = 4;
    // int A_n = B_n * B_n;
	
	
	// // vector<clock_t> startTime;
	
	// SEMO_Utility_Math math;
	
	// int proc_num=0;

	// // gradient P
	// vector<vector<double>> gradP;
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
	// int nrows = mesh.cells.size();
	// int nnz = mesh.non_zeros;
    // vector<double> A_vals(nnz*A_n, 0.0);
    // vector<double> B_vals(nrows*B_n, 0.0);

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];

        // int ijStart_local = B_n*(i) - 1;
        // int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
        // int Astart = A_n*(i) - 1;
        // int id = Astart;
		
		// double Rho = cell.var[controls.Rho];
		// double U = cell.var[controls.U];
		// double V = cell.var[controls.V];
		// double W = cell.var[controls.W];
		// double oldU = cell.var[controls.oldU];
		// double oldV = cell.var[controls.oldV];
		// double oldW = cell.var[controls.oldW];
		// double volume = cell.volume;
		// double timeStep = controls.timeStep;
		
		// int str_row = mesh.CRS_ptr[i]*B_n*B_n;
		
		// int str_col_p = (mesh.CRS_ptr[i+1] - mesh.CRS_ptr[i])*0;
		// int str_col_u = (mesh.CRS_ptr[i+1] - mesh.CRS_ptr[i])*1;
		// int str_col_v = (mesh.CRS_ptr[i+1] - mesh.CRS_ptr[i])*2;
		// int str_col_w = (mesh.CRS_ptr[i+1] - mesh.CRS_ptr[i])*3;
		// int n_col = (mesh.CRS_ptr[i+1]-mesh.CRS_ptr[i])*B_n;
		
		// int i_col = mesh.CRS_col_ptr_dig[i] - mesh.CRS_ptr[i];
		
		// vector<int> dig_p_eq(B_n,0);
		// vector<int> dig_u_eq(B_n,0);
		// vector<int> dig_v_eq(B_n,0);
		// vector<int> dig_w_eq(B_n,0);
		// for(int j=0; j<B_n; ++j){
			// dig_p_eq[j] = str_row + n_col*j + str_col_p + i_col;
			// dig_u_eq[j] = str_row + n_col*j + str_col_u + i_col;
			// dig_v_eq[j] = str_row + n_col*j + str_col_v + i_col;
			// dig_w_eq[j] = str_row + n_col*j + str_col_w + i_col;
		// }
		
        // // continuity

        // // x-momentum
		// A_vals[dig_u_eq[1]] = Rho*volume/timeStep;

        // // y-momentum
		// A_vals[dig_v_eq[2]] = Rho*volume/timeStep;

        // // z-momentum
		// A_vals[dig_w_eq[3]] = Rho*volume/timeStep;


        // B_vals[i*B_n + 0] = 0.0;
        // B_vals[i*B_n + 1] = -Rho*(U - oldU)*volume / timeStep;
        // B_vals[i*B_n + 2] = -Rho*(V - oldV)*volume / timeStep;
        // B_vals[i*B_n + 3] = -Rho*(W - oldW)*volume / timeStep;
		
		// // gravity force terms
        // B_vals[i*B_n + 1] += Rho*controls.gravityAcceleration[0]*volume;
        // B_vals[i*B_n + 2] += Rho*controls.gravityAcceleration[1]*volume;
        // B_vals[i*B_n + 3] += Rho*controls.gravityAcceleration[2]*volume;


		// // // surface tension force terms
        // // B_vals[nrows*1 + i] += volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
        // // B_vals[nrows*2 + i] += volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
        // // B_vals[nrows*3 + i] += volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
		
	// }
	
	
	
	
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
		// // double wCL = face.wVC;
		// // double wCL = 0.5;
		// double wCR = 1.0-wCL;
		
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
				
		// // double nonOrtholimiter = 0.0;
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
		
		// double UF = weightL*UL + weightR*UR;
		// double VF = weightL*VL + weightR*VR;
		// double WF = weightL*WL + weightR*WR;
		// double RhoF = weightL*RhoL + weightR*RhoR;
		// double PF = wCL*PL + wCR*PR;
		
		
		
		
		// int row_own = -1;
		// int row_ngb = -1;
		// int n_col_p_own = -1;
		// int n_col_p_ngb = -1;
		// int col_i_own_own = -1;
		// int col_i_ngb_ngb = -1;
		// int col_i_own_ngb = -1;
		// int col_i_ngb_own = -1;
		// int str_a_ii = -1;
		// int str_a_jj = -1;
		// int str_a_ij = -1;
		// int str_a_ji = -1;
		// int irow = -1;
		// int icol = -1;
		// int str_row_own = -1;
		// int str_row_ngb = -1;
		
		// str_row_own = mesh.CRS_ptr[face.owner]*B_n*B_n;
		// n_col_p_own = (mesh.CRS_ptr[face.owner+1] - mesh.CRS_ptr[face.owner]);
		// row_own = n_col_p_own*B_n;
		// col_i_own_own = (mesh.CRS_col_ptr_dig[face.owner] - mesh.CRS_ptr[face.owner]);
		// col_i_own_ngb = (mesh.CRS_col_ptr_LR[i] - mesh.CRS_ptr[face.owner]);
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// str_row_ngb = mesh.CRS_ptr[face.neighbour]*B_n*B_n;
			// row_ngb = mesh.CRS_ptr[face.neighbour]*B_n;
			// n_col_p_ngb = (mesh.CRS_ptr[face.neighbour+1] - mesh.CRS_ptr[face.neighbour]);
			// col_i_ngb_ngb = (mesh.CRS_col_ptr_dig[face.neighbour] - mesh.CRS_ptr[face.neighbour]);
			// col_i_ngb_own = (mesh.CRS_col_ptr_RL[i] - mesh.CRS_ptr[face.neighbour]);
		// }
		
		// vector<vector<int>> a_ii(B_n,vector<int>(B_n,-1));
		// vector<vector<int>> a_jj(B_n,vector<int>(B_n,-1));
		// vector<vector<int>> a_ij(B_n,vector<int>(B_n,-1));
		// vector<vector<int>> a_ji(B_n,vector<int>(B_n,-1));
		// for(int ii=0; ii<B_n; ++ii){
			// str_a_ii = str_row_own + row_own*ii + col_i_own_own;
			// str_a_jj = str_row_ngb + row_ngb*ii + col_i_ngb_ngb;
			// str_a_ij = str_row_own + row_own*ii + col_i_own_ngb;
			// str_a_ji = str_row_ngb + row_ngb*ii + col_i_ngb_own;
			// for(int jj=0; jj<B_n; ++jj){
				// a_ii[ii][jj] = str_a_ii + n_col_p_own*jj;
				// a_jj[ii][jj] = str_a_jj + n_col_p_ngb*jj;
				// a_ij[ii][jj] = str_a_ij + n_col_p_own*jj;
				// a_ji[ii][jj] = str_a_ji + n_col_p_ngb*jj;
			// }
		// }

        // //------------------------
        // // continuity
		// irow = 0;
		
        // // p'
		// icol = 0;
		
		// A_vals[a_ii[irow][icol]] += ( tmp1 / dPN * area );
		// A_vals[a_ij[irow][icol]] = ( - tmp1 / dPN * area );

		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( - tmp1 / dPN * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( tmp1 / dPN * area );
		
        // // u'
		// icol = 1;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[0] * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[0] * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[0] * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[0] * area );
		
		
        // // v'
		// icol = 2;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[1] * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[1] * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[1] * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[1] * area );
		
        
        // // w'
		// icol = 3;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[2] * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[2] * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[2] * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[2] * area );

        

        // //------------------------
        // // x-momentum
		// irow = 1;

        // // p'
		// icol = 0;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[0] * area + RhoL * tmp1 / dPN * UF * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[0] * area - RhoL * tmp1 / dPN * UF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[0] * area - RhoR * tmp1 / dPN * UF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[0] * area + RhoR * tmp1 / dPN * UF * area );
		
		
        // // u'
		// icol = 1;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[0] * UF * area + weightL * RhoL * UnF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[0] * UF * area + weightR * RhoL * UnF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[0] * UF * area + weightR * RhoR * UnF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[0] * UF * area + weightL * RhoR * UnF * area );
		
		
        // // v'
		// icol = 2;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[1] * UF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[1] * UF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[1] * UF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[1] * UF * area );
		
		
        // // w'
		// icol = 3;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[2] * UF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[2] * UF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[2] * UF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[2] * UF * area );
		
		
        // //------------------------
        // // y-momentum
		// irow = 2;
        
        // // p'
		// icol = 0;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[1] * area + RhoL * tmp1 / dPN * VF * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[1] * area - RhoL * tmp1 / dPN * VF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[1] * area - RhoR * tmp1 / dPN * VF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[1] * area + RhoR * tmp1 / dPN * VF * area );
		
		
        // // u'
		// icol = 1;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[0] * VF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[0] * VF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[0] * VF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[0] * VF * area );
		
		
        // // v'
		// icol = 2;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[1] * VF * area + weightL * RhoL * UnF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[1] * VF * area + weightR * RhoL * UnF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[1] * VF * area + weightR * RhoR * UnF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[1] * VF * area + weightL * RhoR * UnF * area );
		
		
        // // w'
		// icol = 3;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[2] * VF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[2] * VF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[2] * VF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[2] * VF * area );
		
		
		
        // //------------------------
        // // z-momentum
		// irow = 3;
        
        // // p'
		// icol = 0;
		
		// A_vals[a_ii[irow][icol]] += ( wCL * nvec[2] * area + RhoL * tmp1 / dPN * WF * area );
		// A_vals[a_ij[irow][icol]] = ( wCR * nvec[2] * area - RhoL * tmp1 / dPN * WF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( wCR * nvec[2] * area - RhoR * tmp1 / dPN * WF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( wCL * nvec[2] * area + RhoR * tmp1 / dPN * WF * area );
		
		
        // // u'
		// icol = 1;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[0] * WF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[0] * WF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[0] * WF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[0] * WF * area );
		
		
        // // v'
		// icol = 2;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[1] * WF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[1] * WF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[1] * WF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[1] * WF * area );
		
		
        // // w'
		// icol = 3;
		
		// A_vals[a_ii[irow][icol]] += ( RhoL * wCL * nvec[2] * WF * area + weightL * RhoL * UnF * area );
		// A_vals[a_ij[irow][icol]] = ( RhoL * wCR * nvec[2] * WF * area + weightR * RhoL * UnF * area );
		
		// if(a_jj[irow][icol]>=0) A_vals[a_jj[irow][icol]] -= ( RhoR * wCR * nvec[2] * WF * area + weightR * RhoR * UnF * area );
		// if(a_ji[irow][icol]>=0) A_vals[a_ji[irow][icol]] = -( RhoR * wCL * nvec[2] * WF * area + weightL * RhoR * UnF * area );
		
		
        // // ----------------------------
        // // B
        // B_vals[face.owner*B_n + 0] -= ( UnF * area );
        // B_vals[face.owner*B_n + 1] -= ( RhoL * UF * UnF * area + PF * nvec[0] * area );
        // B_vals[face.owner*B_n + 2] -= ( RhoL * VF * UnF * area + PF * nvec[1] * area );
        // B_vals[face.owner*B_n + 3] -= ( RhoL * WF * UnF * area + PF * nvec[2] * area );

		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// B_vals[face.neighbour*B_n + 0] += ( UnF * area );
			// B_vals[face.neighbour*B_n + 1] += ( RhoR * UF * UnF * area + PF * nvec[0] * area );
			// B_vals[face.neighbour*B_n + 2] += ( RhoR * VF * UnF * area + PF * nvec[1] * area );
			// B_vals[face.neighbour*B_n + 3] += ( RhoR * WF * UnF * area + PF * nvec[2] * area );
		// }
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
		
	// }
	
	
	
	
	
	
	// // // boundary
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			
			// // for(int i=str; i<end; ++i){
				// // auto& face = mesh.faces[i];
				
				// // double area = face.area;
				// // vector<double> nvec(3,0.0);
				// // nvec[0] = face.unitNormals[0];
				// // nvec[1] = face.unitNormals[1];
				// // nvec[2] = face.unitNormals[2];
		
				// // double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
							 // // face.varL[controls.fV]*face.unitNormals[1] + 
							 // // face.varL[controls.fW]*face.unitNormals[2];
				// // double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
							 // // face.varR[controls.fV]*face.unitNormals[1] + 
							 // // face.varR[controls.fW]*face.unitNormals[2];
				
				// // double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							   // // face.unitNormals[1]*face.distCells[1] + 
							   // // face.unitNormals[2]*face.distCells[2];
				
				// // double dPN = sqrt(pow(face.distCells[0],2.0) + 
								  // // pow(face.distCells[1],2.0) + 
								  // // pow(face.distCells[2],2.0));
							
				// // double UnF = 0.5*UnL+0.5*UnR;
				
				// // double RhoL = face.varL[controls.fRho];
				// // double tmp1 = 1.0/RhoL*controls.timeStep;
				
				// // // for(int ii=0; ii<3; ++ii){
					// // // UnF += tmp1 * gradP[face.owner][ii]*nvec[ii];
				// // // }
				
				// // double UF = 0.5*face.varL[controls.fU] + 0.5*face.varR[controls.fU];
				// // double VF = 0.5*face.varL[controls.fV] + 0.5*face.varR[controls.fV];
				// // double WF = 0.5*face.varL[controls.fW] + 0.5*face.varR[controls.fW];
				// // double PF = 0.5*face.varL[controls.fP] + 0.5*face.varR[controls.fP];
			
				
				
				// // // pressure coefficient
				// // double coeffP = 0.0;
				// // double coeffP_diff = 0.0;
				// // if( boundary.type[controls.P] == "fixedValue" ){
					// // coeffP = 0.0;
					// // coeffP_diff = 1.0;
				// // }
				// // else if( boundary.type[controls.P] == "zeroGradient" ){
					// // coeffP = 1.0;
					// // coeffP_diff = 0.0;
					
				// // }
				// // else if( boundary.type[controls.P] == "switch" ){
					// // double machNum = 
						// // sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
							 // // pow(mesh.cells[face.owner].var[controls.V],2.0)+
							 // // pow(mesh.cells[face.owner].var[controls.W],2.0))/
							  // // mesh.cells[face.owner].var[controls.C];
					// // coeffP = (machNum > 1.0) ? 0.0 : 1.0;
					// // coeffP_diff = (machNum > 1.0) ? 1.0 : 0.0;
				// // }
				// // else {
					// // cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
					// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// // }
				
				// // // velocity coefficient
				// // double coeffUVW = 0.0;
				// // if( boundary.type[controls.U] == "fixedValue" ){
					// // coeffUVW = 0.0;
				// // }
				// // else if( boundary.type[controls.U] == "zeroGradient" ){
					// // coeffUVW = 1.0;
				// // }
				// // else if( boundary.type[controls.U] == "slip" ){
					// // coeffUVW = 0.0;
				// // }
				// // else if( boundary.type[controls.U] == "noSlip" ){
					// // coeffUVW = 0.0;
				// // }
				// // else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
					// // coeffUVW = 0.0;
				// // }
				// // else if( boundary.type[controls.U] == "inletOutlet" ){
					// // double ownNorVel =  
						// // mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// // mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// // mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// // coeffUVW = (ownNorVel > 1.0) ? 1.0 : 0.0;
				// // }
				// // else {
					// // cerr << "| #Error : not defined B.C., var = " << controls.U << endl;
					// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// // }
				
				// // // volume-fraction coefficient
				// // double coeffVF = 0.0;
				// // if( boundary.type[controls.VF[0]] == "fixedValue" ){
					// // coeffVF = 0.0;
				// // }
				// // else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
					// // coeffVF = 1.0;
					
				// // }
				// // else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
					// // double ownNorVel =  
						// // mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
						// // mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
						// // mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					// // coeffVF = (ownNorVel > 1.0) ? 1.0 : 0.0;
				// // }
				// // else {
					// // cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
					// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// // }
				
				
				
				
				// // int str_row_own = mesh.CRS_ptr[face.owner]*B_n*B_n;
				// // int row_own = mesh.CRS_ptr[face.owner]*B_n;
				// // int n_col_p_own = (mesh.CRS_ptr[face.owner+1] - mesh.CRS_ptr[face.owner]);
				// // int col_i_own_own = (mesh.CRS_col_ptr_dig[face.owner] - mesh.CRS_ptr[face.owner]);
				
				// // vector<vector<int>> a_ii(B_n,vector<int>(B_n,0));
				// // for(int ii=0; ii<B_n; ++ii){
					// // int str_a_ii = str_row_own + row_own*ii + col_i_own_own;
					// // for(int jj=0; jj<B_n; ++jj){
						// // a_ii[ii][jj] = str_a_ii + n_col_p_own*jj;
					// // }
				// // }
				
				
				// // // continuity
				// // A_vals[a_ii[0][0]] += coeffP_diff * ( tmp1 / dPN * area );
				// // A_vals[a_ii[0][1]] += coeffUVW *( nvec[0] * area );
				// // A_vals[a_ii[0][2]] += coeffUVW *( nvec[1] * area );
				// // A_vals[a_ii[0][3]] += coeffUVW *( nvec[2] * area );

				
				// // // x-momentum
				// // A_vals[a_ii[1][0]] += ( coeffP * nvec[0] * area + coeffP_diff * RhoL * tmp1 / dPN * UF * area );
				// // A_vals[a_ii[1][1]] += coeffUVW * ( RhoL * nvec[0] * UF * area + RhoL * UnF * area );
				// // A_vals[a_ii[1][2]] += coeffUVW * ( RhoL * nvec[1] * UF * area );
				// // A_vals[a_ii[1][3]] += coeffUVW * ( RhoL * nvec[2] * UF * area );
				
				
				// // // y-momentum
				// // A_vals[a_ii[2][0]] += ( coeffP * nvec[1] * area + coeffP_diff * RhoL * tmp1 / dPN * VF * area );
				// // A_vals[a_ii[2][1]] += coeffUVW *( RhoL * nvec[0] * VF * area );
				// // A_vals[a_ii[2][2]] += coeffUVW *( RhoL * nvec[1] * VF * area + RhoL * UnF * area );
				// // A_vals[a_ii[2][3]] += coeffUVW *( RhoL * nvec[2] * VF * area );
				

				
				// // // z-momentum
				// // A_vals[a_ii[3][0]] += ( coeffP * nvec[2] * area + coeffP_diff * RhoL * tmp1 / dPN * WF * area );
				// // A_vals[a_ii[3][1]] += coeffUVW *( RhoL * nvec[0] * WF * area );
				// // A_vals[a_ii[3][2]] += coeffUVW *( RhoL * nvec[1] * WF * area );
				// // A_vals[a_ii[3][3]] += coeffUVW *( RhoL * nvec[2] * WF * area + RhoL * UnF * area );

				
				// // // convective term
				// // B_vals[face.owner*B_n + 0] -= ( UnF * area );
				// // B_vals[face.owner*B_n + 1] -= ( RhoL * UF * UnF * area + PF * face.unitNormals[0] * area );
				// // B_vals[face.owner*B_n + 2] -= ( RhoL * VF * UnF * area + PF * face.unitNormals[1] * area );
				// // B_vals[face.owner*B_n + 3] -= ( RhoL * WF * UnF * area + PF * face.unitNormals[2] * area );
			// // }
			
			
			
		// // }
		
	// // }
	
	
	
	
	
	// vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	
	
	// solveAMGCL_Flows(mesh, B_n, A_vals, B_vals, resiVar);
	
	
	
	// // update
	// // double relaxP = 0.2;
	// // double relaxUVW = 0.2;
	// double relaxP = controls.prePreURF;
	// double relaxUVW = controls.momVelURF;
	
	// double normPrim = 0.0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// cell.var[controls.P] += relaxP * resiVar[nrows*0 + i];
		// cell.var[controls.U] += relaxUVW * resiVar[nrows*1 + i];
		// cell.var[controls.V] += relaxUVW * resiVar[nrows*2 + i];
		// cell.var[controls.W] += relaxUVW * resiVar[nrows*3 + i];
		
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
		
		
		// normPrim += pow(resiVar[nrows*0 + i],2.0);
		// normPrim += pow(resiVar[nrows*1 + i],2.0);
		// normPrim += pow(resiVar[nrows*2 + i],2.0);
		// normPrim += pow(resiVar[nrows*3 + i],2.0);
		
	// }
	
	
	// return sqrt(normPrim);
	
	
// }



