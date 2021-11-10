#include "build.h"
#include <cmath>
#include <array>
#include <numeric>




double SEMO_Solvers_Builder::calcVolfracEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 1;
    int A_n = B_n * B_n;
	

	
	SEMO_Utility_Math math;
	
	int proc_num=0;
	
	// gradient P
	vector<vector<double>> gradP;
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
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
	


	// gradient U, V, W
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	// math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// math.calcMGG(mesh, controls.U, controls.fU, 10, 1.e-8, gradU);
	// math.calcMGG(mesh, controls.V, controls.fV, 10, 1.e-8, gradV);
	// math.calcMGG(mesh, controls.W, controls.fW, 10, 1.e-8, gradW);
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	vector<double> gradUx_recv, gradUy_recv, gradUz_recv;
	vector<double> gradVx_recv, gradVy_recv, gradVz_recv;
	vector<double> gradWx_recv, gradWy_recv, gradWz_recv;
	if(size>1){
		// processor faces
		vector<double> gradUx_send, gradUy_send, gradUz_send;
		vector<double> gradVx_send, gradVy_send, gradVz_send;
		vector<double> gradWx_send, gradWy_send, gradWz_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				gradUx_send.push_back(gradU[face.owner][0]);
				gradUy_send.push_back(gradU[face.owner][1]);
				gradUz_send.push_back(gradU[face.owner][2]);
				
				gradVx_send.push_back(gradV[face.owner][0]);
				gradVy_send.push_back(gradV[face.owner][1]);
				gradVz_send.push_back(gradV[face.owner][2]);
				
				gradWx_send.push_back(gradW[face.owner][0]);
				gradWy_send.push_back(gradW[face.owner][1]);
				gradWz_send.push_back(gradW[face.owner][2]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					gradUx_send, gradUx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradUy_send, gradUy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradUz_send, gradUz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradUx_send.clear();
		gradUy_send.clear();
		gradUz_send.clear();
					
		mpi.setProcsFaceDatas(
					gradVx_send, gradVx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradVy_send, gradVy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradVz_send, gradVz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradVx_send.clear();
		gradVy_send.clear();
		gradVz_send.clear();
					
		mpi.setProcsFaceDatas(
					gradWx_send, gradWx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradWy_send, gradWy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradWz_send, gradWz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradWx_send.clear();
		gradWy_send.clear();
		gradWz_send.clear();
	}
	
	
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
		
		
        // volume fraction
        id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
        A_rows[id] = i_glo; A_cols[id] = j_glo;
        A_vals[id] = volume/timeStep;
		
		
        B[step_loc*0 + i_loc] = -(VF - oldVF)*volume/timeStep;
		
		
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
		// double wCPL = 0.5;
		// double wCPR = 0.5;
		
		double UL = face.varL[controls.fU];
		double VL = face.varL[controls.fV];
		double WL = face.varL[controls.fW];
		double PL = face.varL[controls.fP];
		double RhoL = face.varL[controls.fRho];
		double VFL = face.varL[controls.fVF[0]];
		
		double UR = face.varR[controls.fU];
		double VR = face.varR[controls.fV];
		double WR = face.varR[controls.fW];
		double PR = face.varR[controls.fP];
		double RhoR = face.varR[controls.fRho];
		double VFR = face.varR[controls.fVF[0]];
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double Rho_star = 1.0 / (wCPL/RhoL + wCPR/RhoR);
		double tmp1 = controls.timeStep / Rho_star;
		// // F.Denner et al., 2013
		// double d_F = 0.0;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
			// d_F += wCR*( A_p_vals[face.neighbour] + RhoR/controls.timeStep );
		// }
		// else{
			// d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
			// d_F += wCR*( A_p_vals_recv[proc_num] + RhoR/controls.timeStep );
		// }
		// tmp1 = 1.0/d_F;
		
		
		double dPN = face.magPN;
		// vector<double> Ef(3,0.0);
		// Ef[0] = face.vecPN[0]/dPN;
		// Ef[1] = face.vecPN[1]/dPN;
		// Ef[2] = face.vecPN[2]/dPN;
		// double alpha = nvec[0]*Ef[0] + nvec[1]*Ef[1] + nvec[2]*Ef[2];
		double alpha = face.alphaF;
		
		
		
		
		double orgPL = mesh.cells[face.owner].var[controls.P];
		double orgPR = 0.0;
		vector<double> gradPL(3,0.0);
		gradPL[0] = gradP[face.owner][0];
		gradPL[1] = gradP[face.owner][1];
		gradPL[2] = gradP[face.owner][2];
		vector<double> gradPR(3,0.0);
		// properties
		double orgUL = mesh.cells[face.owner].var[controls.U];
		double orgVL = mesh.cells[face.owner].var[controls.V];
		double orgWL = mesh.cells[face.owner].var[controls.W];
		double orgUR = 0.0;
		double orgVR = 0.0;
		double orgWR = 0.0;
		vector<double> gradUL(3,0.0);
		vector<double> gradVL(3,0.0);
		vector<double> gradWL(3,0.0);
		vector<double> gradUR(3,0.0);
		vector<double> gradVR(3,0.0);
		vector<double> gradWR(3,0.0);
		
		for(int ii=0; ii<3; ++ii){
			gradUL[ii] = gradU[face.owner][ii];
			gradVL[ii] = gradV[face.owner][ii];
			gradWL[ii] = gradW[face.owner][ii];
		}
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			orgPR = mesh.cells[face.neighbour].var[controls.P];
			gradPR[0] = gradP[face.neighbour][0];
			gradPR[1] = gradP[face.neighbour][1];
			gradPR[2] = gradP[face.neighbour][2];
			
			orgUR = mesh.cells[face.neighbour].var[controls.U];
			orgVR = mesh.cells[face.neighbour].var[controls.V];
			orgWR = mesh.cells[face.neighbour].var[controls.W];
			for(int ii=0; ii<3; ++ii){
				gradUR[ii] = gradU[face.neighbour][ii];
				gradVR[ii] = gradV[face.neighbour][ii];
				gradWR[ii] = gradW[face.neighbour][ii];
			}
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			orgPR = PR;
			gradPR[0] = gradPx_recv[proc_num];
			gradPR[1] = gradPy_recv[proc_num];
			gradPR[2] = gradPz_recv[proc_num];
			
			orgUR = UR;
			orgVR = VR;
			orgWR = WR;
			gradUR[0] = gradUx_recv[proc_num]; gradUR[1] = gradUy_recv[proc_num]; gradUR[2] = gradUz_recv[proc_num];
			gradVR[0] = gradVx_recv[proc_num]; gradVR[1] = gradVy_recv[proc_num]; gradVR[2] = gradVz_recv[proc_num];
			gradWR[0] = gradWx_recv[proc_num]; gradWR[1] = gradWy_recv[proc_num]; gradWR[2] = gradWz_recv[proc_num];
			
		}
		

		// // skewness correction
		// for(int ii=0; ii<3; ++ii){
			// UnF += wCL * gradUL[ii]*face.vecSkewness[ii]*nvec[0];
			// UnF += wCR * gradUR[ii]*face.vecSkewness[ii]*nvec[0];
			
			// UnF += wCL * gradVL[ii]*face.vecSkewness[ii]*nvec[1];
			// UnF += wCR * gradVR[ii]*face.vecSkewness[ii]*nvec[1];
			
			// UnF += wCL * gradWL[ii]*face.vecSkewness[ii]*nvec[2];
			// UnF += wCR * gradWR[ii]*face.vecSkewness[ii]*nvec[2];
		// }
		
		// pressure correction (cell->face)
		for(int ii=0; ii<3; ++ii){
			UnF += alpha * tmp1 * Rho_star * wCPL/RhoL*gradPL[ii]*face.unitNomalsPN[ii];
			UnF += alpha * tmp1 * Rho_star * wCPR/RhoR*gradPR[ii]*face.unitNomalsPN[ii];
		}
			
		UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
		// for(int ii=0; ii<3; ++ii){
			// UnF -= tmp1 * (gradPR[ii]*face.vecNdN[ii] - gradPL[ii]*face.vecPdP[ii])/dPN;
		// }
		
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = 1.0 - weightL;
		
		double VFF = weightL*VFL + weightR*VFR;
		
		
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
        B[step_loc_L*0 + i_loc_L] -= ( VFF * UnF * area );

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			B[step_loc_R*0 + i_loc_R] += ( VFF * UnF * area );
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
	

	// solvePETSc_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
		// B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
		// "bicgstab", 1.e-200, 1.e-12, "mg", controls.maxIterVF[0]);
		
	// solveHypre_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
		// B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
		// "bicgstab", 1.e-6, 1.e-6, "mg", 50);
		
	
	

	

	// update
	// double relaxVF = 0.99;
	
	double normDelVF = 0.0;
	double normVF = 0.0;
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.VF[0]] += controls.vofVofURF * resiVar[i];
		cell.var[controls.VF[0]] = max(1.e-12,min(1.0-1.e-12,cell.var[controls.VF[0]]));
		
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


















