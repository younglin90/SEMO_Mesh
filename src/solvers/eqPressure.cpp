

#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include <ctime>




double SEMO_Solvers_Builder::calcPressureEq(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int iterNonOrthogonality){
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	SEMO_MPI_Builder mpi;
	
	
	// diagonal terms
    int B_n = 1;
    int A_n = B_n * B_n;
	

	
	SEMO_Utility_Math math;
	
	int proc_num=0;
	
	// gradient P
	vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	vector<double> gradPx_recv;
	vector<double> gradPy_recv;
	vector<double> gradPz_recv;
	
	// linear solver
	vector<double> resiVar(B_n*mesh.cells.size(),0.0);
	vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
	vector<double> gradResiPx_recv;
	vector<double> gradResiPy_recv;
	vector<double> gradResiPz_recv;
	
	
	vector<double> A_p_vals(mesh.cells.size(), 0.0);
	
	
	// for(int iterNonOrtho=0; iterNonOrtho<iterNonOrthogonality; ++iterNonOrtho){
	for(int iterNonOrtho=0; iterNonOrtho<1; ++iterNonOrtho){ 
		
		
		// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
		math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
		
		if(iterNonOrtho!=0) math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
		// math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
		
		
		
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
		
		if(size>1){
			// processor faces
			// gradP , 
			vector<double> gradResiPx_send;
			vector<double> gradResiPy_send;
			vector<double> gradResiPz_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					gradResiPx_send.push_back(gradResiP[face.owner][0]);
					gradResiPy_send.push_back(gradResiP[face.owner][1]);
					gradResiPz_send.push_back(gradResiP[face.owner][2]);
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
			gradResiPx_send.clear();
			gradResiPy_send.clear();
			gradResiPz_send.clear();
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
	
	

		// //=============================================
		// // F. Denner, 2013
		// vector<double> A_p_vals(mesh.cells.size(), 0.0);
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
			
			// double UnF = wCL*UnL+wCR*UnR;
			
			// // double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
			// // double tmp1 = controls.timeStep / Rho_star;
			
			// // double dPN = face.magPN;
			
			// // double orgPL = mesh.cells[face.owner].var[controls.P];
			// // double orgPR = 0.0;
			// // vector<double> gradPL(3,0.0);
			// // gradPL[0] = gradP[face.owner][0];
			// // gradPL[1] = gradP[face.owner][1];
			// // gradPL[2] = gradP[face.owner][2];
			// // vector<double> gradPR(3,0.0);
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// // orgPR = mesh.cells[face.neighbour].var[controls.P];
				// // gradPR[0] = gradP[face.neighbour][0];
				// // gradPR[1] = gradP[face.neighbour][1];
				// // gradPR[2] = gradP[face.neighbour][2];
				
				
			// // }
			// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// // orgPR = PR;
				// // gradPR[0] = gradPx_recv[proc_num];
				// // gradPR[1] = gradPy_recv[proc_num];
				// // gradPR[2] = gradPz_recv[proc_num];
				
			// // }
			
			// // for(int ii=0; ii<3; ++ii){
				// // UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*face.vecPN[ii]/dPN;
				// // UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*face.vecPN[ii]/dPN;
			// // }
				
			// // UnF -= tmp1*(orgPR-orgPL)/dPN;
			// // for(int ii=0; ii<3; ++ii){
				// // UnF -= tmp1 * (gradPR[ii]*face.vecNdN[ii] - gradPL[ii]*face.vecPdP[ii])/dPN;
			// // }
				
			// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
			// double weightR = 1.0 - weightL;
			
			// double RhoF = 
				// weightL*face.varL[controls.fRho] +
				// weightR*face.varR[controls.fRho];
			
			// A_p_vals[face.owner] += weightL * RhoF*UnF*face.area;// / mesh.cells[face.owner].volume;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// A_p_vals[face.neighbour] -=  weightR * RhoF*UnF*face.area;// / mesh.cells[face.neighbour].volume;
			// }
			
			
		// }
		
		

		// // boundary
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					// double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 // face.varL[controls.fV]*face.unitNormals[1] + 
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 // face.varR[controls.fV]*face.unitNormals[1] + 
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double UnF = 0.5*UnL+0.5*UnR;
					// double RhoL = face.varL[controls.fRho];
					// double RhoF = RhoL;
					
					// A_p_vals[face.owner] += RhoF*UnF*face.area;// / mesh.cells[face.owner].volume;
				
				// }
				
				
			// }
			
		// }
		
		
		
		
		
		// vector<double> A_p_vals_recv;
		// if(size>1){
			// // processor faces
			// // gradP , 
			// vector<double> A_p_vals_send;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// A_p_vals_send.push_back(A_p_vals[face.owner]);
				// }
			// }
			// // SEMO_MPI_Builder mpi;
			
			// mpi.setProcsFaceDatas(
						// A_p_vals_send, A_p_vals_recv,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// A_p_vals_send.clear();
		// }
		// //=============================================



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
		
		
		// if(rank==0) cout << "AAAAAAAAA" << endl;

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
				// // d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
				// // d_F += wCR*( A_p_vals[face.neighbour] + RhoR/controls.timeStep );
				// d_F += wCL*( mesh.cells[face.owner].volume / (A_p_vals[face.owner]+1.e-50) );
				// d_F += wCR*( mesh.cells[face.neighbour].volume / (A_p_vals[face.neighbour]+1.e-50) );
			// }
			// else{
				// // d_F += wCL*( A_p_vals[face.owner] + RhoL/controls.timeStep );
				// // d_F += wCR*( A_p_vals_recv[proc_num] + RhoR/controls.timeStep );
				// d_F += wCL*( mesh.cells[face.owner].volume / (A_p_vals[face.owner]+1.e-50) );
				// d_F += wCR*( mesh.cells[face.owner].volume / (A_p_vals_recv[proc_num]+1.e-50) );
			// }
			// tmp1 = d_F / (1.0 + d_F * Rho_star / controls.timeStep);
			
				
				
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
			vector<double> gradResiPL(3,0.0);
			gradPL[0] = gradP[face.owner][0];
			gradPL[1] = gradP[face.owner][1];
			gradPL[2] = gradP[face.owner][2];
			
			gradResiPL[0] = gradResiP[face.owner][0];
			gradResiPL[1] = gradResiP[face.owner][1];
			gradResiPL[2] = gradResiP[face.owner][2];
			
			vector<double> gradPR(3,0.0);
			vector<double> gradResiPR(3,0.0);
			
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
			
			// // non-orthogonal
			// double coeff_non_orth = 1.0;
			// B[step_loc_L*0 + i_loc_L] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[0] * face.vecNdN[0] );
			// B[step_loc_L*0 + i_loc_L] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[1] * face.vecNdN[1] );
			// B[step_loc_L*0 + i_loc_L] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[2] * face.vecNdN[2] );
			// B[step_loc_L*0 + i_loc_L] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[0] * face.vecPdP[0] );
			// B[step_loc_L*0 + i_loc_L] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[1] * face.vecPdP[1] );
			// B[step_loc_L*0 + i_loc_L] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[2] * face.vecPdP[2] );
			
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				B[step_loc_R*0 + i_loc_R] += ( UnF * area );
			
				// // non-orthogonal
				// B[step_loc_R*0 + i_loc_R] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[0] * face.vecNdN[0] );
				// B[step_loc_R*0 + i_loc_R] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[1] * face.vecNdN[1] );
				// B[step_loc_R*0 + i_loc_R] -= coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPR[2] * face.vecNdN[2] );
				// B[step_loc_R*0 + i_loc_R] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[0] * face.vecPdP[0] );
				// B[step_loc_R*0 + i_loc_R] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[1] * face.vecPdP[1] );
				// B[step_loc_R*0 + i_loc_R] += coeff_non_orth * ( alpha * tmp1 / dPN * area * gradResiPL[2] * face.vecPdP[2] );
			}
			
			
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
			
			
		}
		
		
		
		// if(rank==0) cout << "BBBBBBBBBB" << endl;
		
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
					
					// double dPN_e = face.unitNormals[0]*face.distCells[0] + 
								   // face.unitNormals[1]*face.distCells[1] + 
								   // face.unitNormals[2]*face.distCells[2];
					
					// double dPN = sqrt(pow(face.distCells[0],2.0) + 
									  // pow(face.distCells[1],2.0) + 
									  // pow(face.distCells[2],2.0));
								
					double UnF = 0.5*UnL+0.5*UnR;
					
					double RhoL = face.varL[controls.fRho];
					double tmp1 = 1.0/RhoL*controls.timeStep;
					
					// tmp1 = 1.0 / (A_p_vals[face.owner] + RhoL/controls.timeStep);
					
					// for(int ii=0; ii<3; ++ii){
						// UnF += tmp1 * gradP[face.owner][ii]*nvec[ii];
					// }
					
					
					// dPN = dPN_e;
					
					double dPN = face.magPN;
					
					
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
					
					

					int str_glo = mesh.startCellGlobal*B_n;
					int step_loc = mesh.cells.size();
					int step_glo = mesh.ncellsTotal;
					int i_loc_L = face.owner;
					int id=-1;
					
					
					// continuity
					id = step_loc*(B_n*0+0) + i_loc_L;
					A_vals[id] += coeffP_diff * ( tmp1 / dPN * area );
					
					// convective term
					B[step_loc*0 + i_loc_L] -= ( UnF * area );
					
					// non-orthogonal
					vector<double> gradResiPL(3,0.0);
					gradResiPL[0] = gradResiP[face.owner][0];
					gradResiPL[1] = gradResiP[face.owner][1];
					gradResiPL[2] = gradResiP[face.owner][2];
					
					B[step_loc*0 + i_loc_L] += ( tmp1 / dPN * area * gradResiPL[0] * face.vecNdN[0] );
					B[step_loc*0 + i_loc_L] += ( tmp1 / dPN * area * gradResiPL[1] * face.vecNdN[1] );
					B[step_loc*0 + i_loc_L] += ( tmp1 / dPN * area * gradResiPL[2] * face.vecNdN[2] );
					B[step_loc*0 + i_loc_L] -= ( tmp1 / dPN * area * gradResiPL[0] * face.vecPdP[0] );
					B[step_loc*0 + i_loc_L] -= ( tmp1 / dPN * area * gradResiPL[1] * face.vecPdP[1] );
					B[step_loc*0 + i_loc_L] -= ( tmp1 / dPN * area * gradResiPL[2] * face.vecPdP[2] );
				}
				
			}
			
		}
		
		

		
		
		solveAMGCL("pressure", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);

		// solveHypre_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
			// B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
			// "bicgstab", 1.e-6, 1.e-6, "mg", 50);
		
		// solvePETSc_Test_Coupled(mesh, resiVar, A_rows, A_cols, A_vals, B,
			// B_n*mesh.cells.size(), B_n*mesh.ncellsTotal,
			// "bicgstab", 1.e-6, 1.e-6, "mg", 50);
		
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// SEMO_Cell& cell = mesh.cells[i];
			
			// cell.var[controls.P] += controls.prePreURF * resiVar[i];
		// }
		
	
	
	}
	
	

	
	// math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
	
	for(auto& var : gradResiP){
		std::fill(var.begin(), var.end(), 0.0);
	}
	
	vector<double> resiVar_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> resiVar_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				resiVar_send.push_back(resiVar[face.owner]);
			}
		}
		mpi.setProcsFaceDatas(
					resiVar_send, resiVar_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		resiVar_send.clear();
	}
	
	proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		double wCPL = face.wC;
		// double wCPL = 0.5;
		double wCPR = 1.0-wCPL;
		
		double PF = 0.0;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			PF = wCPL * resiVar[face.owner] + wCPR * resiVar[face.neighbour];
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			PF = wCPL * resiVar[face.owner] + wCPR * resiVar_recv[proc_num];
		}
		
		for(int j=0; j<3; ++j){
			gradResiP[face.owner][j] += PF*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
		}
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(int j=0; j<3; ++j){
				gradResiP[face.neighbour][j] -=  PF*face.unitNormals[j]*face.area / mesh.cells[face.neighbour].volume;
			}
		}
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
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
					varF = 0.5*varF;
				}
				
				for(int j=0; j<3; ++j){
					gradResiP[face.owner][j] += varF*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
				}
			}
			
		}
	}
	

	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradResiPx_send;
		// vector<double> gradResiPy_send;
		// vector<double> gradResiPz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradResiPx_send.push_back(gradResiP[face.owner][0]);
				// gradResiPy_send.push_back(gradResiP[face.owner][1]);
				// gradResiPz_send.push_back(gradResiP[face.owner][2]);
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
		// gradResiPx_send.clear();
		// gradResiPy_send.clear();
		// gradResiPz_send.clear();
	// }
	

	

	// // pressure skewness correction
	// proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		// double wCPL = face.wC;
		// // double wCPL = 0.5;
		// double wCPR = 1.0-wCPL;
		
		// vector<double> gradResiPL(3,0.0);
		// for(int ii=0; ii<3; ++ii){
			// gradResiPL[ii] = gradResiP[face.owner][ii];
		// }
		
		// vector<double> gradResiPR(3,0.0);
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// for(int ii=0; ii<3; ++ii){
				// gradResiPR[ii] = gradResiP[face.neighbour][ii];
			// }
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// gradResiPR[0] = gradResiPx_recv[proc_num];
			// gradResiPR[1] = gradResiPy_recv[proc_num];
			// gradResiPR[2] = gradResiPz_recv[proc_num];
		// }
		
		// double PF_Skewness = 0.0;
		// for(int ii=0; ii<3; ++ii){
			// PF_Skewness += wCPL * gradResiPL[ii]*face.vecSkewness[ii];
			// PF_Skewness += wCPR * gradResiPR[ii]*face.vecSkewness[ii];
		// }
		
		// for(int j=0; j<3; ++j){
			// gradResiP[face.owner][j] += PF_Skewness*face.unitNormals[j]*face.area / mesh.cells[face.owner].volume;
		// }
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// for(int j=0; j<3; ++j){
				// gradResiP[face.neighbour][j] -=  PF_Skewness*face.unitNormals[j]*face.area / mesh.cells[face.neighbour].volume;
			// }
		// }
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
	// }
	

		   
	
	
	
	
	
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
	

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	int proc_num=0;
	
	vector<clock_t> startTime;
	
	SEMO_Utility_Math math;

	// gradient P
	vector<vector<double>> gradP;
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
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
		
		double wCL = face.wC;
		// double wCL = face.wVC;
		// double wCL = 0.5;
		double wCR = 1.0-wCL;
		
		double UnF = 
		    (wCL*face.varL[controls.fU]*face.unitNormals[0] +
			 wCL*face.varL[controls.fV]*face.unitNormals[1] +
			 wCL*face.varL[controls.fW]*face.unitNormals[2] +
			 wCR*face.varR[controls.fU]*face.unitNormals[0] +
			 wCR*face.varR[controls.fV]*face.unitNormals[1] +
			 wCR*face.varR[controls.fW]*face.unitNormals[2]);
			 
		// if( UnF < std::numeric_limits<double>::min() ){
			// UnF = std::numeric_limits<double>::min();
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


	// diagonal terms
	vector<double> linA(mesh.cells.size(),0.0);
	vector<double> linB(mesh.cells.size(),0.0);

	// diagonal, off-diagonal terms
	vector<double> linAL(mesh.faces.size(),0.0);
	vector<double> linAR(mesh.faces.size(),0.0);
	
	proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
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
		
		// double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  // pow(distanceCells[1],2.0) + 
						  // pow(distanceCells[2],2.0));
		double dPN = sqrt(distanceCells[0]*distanceCells[0] + 
						  distanceCells[1]*distanceCells[1] + 
						  distanceCells[2]*distanceCells[2]);
				
		double nonOrtholimiter = 1.0;
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
			  
		dPN = dPN_e;
				
		vector<double> Ef(3,0.0);
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
		// // over-relaxed
		// double over_Ef = Ef[0]*nvec[0] + Ef[1]*nvec[1] + Ef[2]*nvec[2];
		// Ef[0] /= over_Ef; Ef[1] /= over_Ef; Ef[2] /= over_Ef;
		// double magEf = sqrt(Ef[0]*Ef[0]+Ef[1]*Ef[1]+Ef[2]*Ef[2]);
		// dPN /= magEf;
		// // dPN = dPN_e;
			
		vector<double> Tf(3,0.0);
		Tf[0] = nvec[0] - Ef[0];
		Tf[1] = nvec[1] - Ef[1];
		Tf[2] = nvec[2] - Ef[2];
		
		double UnF = wCL*UnL+wCR*UnR;
		
		double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
		double tmp1 = controls.timeStep / Rho_star;
		// F.Denner et al., 2013
		double d_F = 0.0;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			d_F = wCL/(A_p_vals[face.owner]) + wCR/(A_p_vals[face.neighbour]);
		}
		else{
			d_F = 1.0/(A_p_vals[face.owner]);
		}
		double d_F_hat = d_F / (2.0 + d_F / tmp1);
		if(d_F>1.e9) d_F_hat = tmp1;
		tmp1 = d_F_hat;
		
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
			UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
			// UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*Ef[ii];
			// UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*Ef[ii];
		}
			
		UnF -= tmp1*(orgPR-orgPL)/dPN;
		
		// // non-orthogonal, over-relaxed approach
		// for(int ii=0; ii<3; ++ii){
			// // UnF -= nonOrtholimiter * tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*Tf[ii];
			// // UnF -= nonOrtholimiter * tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*Tf[ii];
			// UnF -= nonOrtholimiter * tmp1 * wCL*gradPL[ii]*Tf[ii];
			// UnF -= nonOrtholimiter * tmp1 * wCR*gradPR[ii]*Tf[ii];
		// }
		
		
		// // F.Denner et al., 2018
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// UnF += face.var[controls.Un];
			// UnF -= 0.5*mesh.cells[face.owner].var[controls.oldU]*nvec[0];
			// UnF -= 0.5*mesh.cells[face.owner].var[controls.oldV]*nvec[1];
			// UnF -= 0.5*mesh.cells[face.owner].var[controls.oldW]*nvec[2];
			// UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldU]*nvec[0];
			// UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldV]*nvec[1];
			// UnF -= 0.5*mesh.cells[face.neighbour].var[controls.oldW]*nvec[2];
			// if(controls.iterPBs+1==controls.iterPBsMax){
				// face.var[controls.Un] = UnF;
			// }
		// }
		
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		double diff_flux = tmp1/dPN_e * area;
		
		// own, ngb
		linAR[i] = (-diff_flux);
		// ngb, own
		linAL[i] = (-diff_flux);
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// own, ngb
			linA[face.owner] += (+diff_flux);
			// ngb, own
			linA[face.neighbour] -= (-diff_flux);
			
			// convective term
			linB[face.owner] -= UnF*area;
			linB[face.neighbour] += UnF*area; 
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// own, ngb
			linA[face.owner] += (+diff_flux);
			
			// convective term
			linB[face.owner] -= UnF*area;
			
			++proc_num;
		}
		
	}
	
	
	
	
	

	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				
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
				
				double dPN_e = face.unitNormals[0]*face.distCells[0] + 
							   face.unitNormals[1]*face.distCells[1] + 
							   face.unitNormals[2]*face.distCells[2];
				
				double dPN = sqrt(pow(face.distCells[0],2.0) + 
								  pow(face.distCells[1],2.0) + 
								  pow(face.distCells[2],2.0));
								  
				dPN = dPN_e;
							
				double UnF = 0.5*UnL+0.5*UnR;
				
				double RhoL = face.varL[controls.fRho];
				double tmp1 = 1.0/RhoL*controls.timeStep;
				
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
				
				// continuity
				linA[face.owner] += coeffP_diff * ( tmp1 / dPN * area );
				
				// convective term
				linB[face.owner] -= ( UnF * area );
			}
			
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
		solver = controls.solverP;
		relTol = controls.toleranceP;
		tolerance = controls.relTolP;
		preconditioner = controls.preconditionerP;
		maxIter = controls.maxIterP;
	
	}
	else{
		solver = controls.solverFinalP;
		relTol = controls.toleranceFinalP;
		tolerance = controls.relTolFinalP;
		preconditioner = controls.preconditionerFinalP;
		maxIter = controls.maxIterFinalP;
	
	}
	
	// solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
		// solver, tolerance, relTol, preconditioner, maxIter);
		
	solveHYPRE(mesh, resiVar, linA, linAL, linAR, linB,
		solver, tolerance, relTol, preconditioner, maxIter);
		
	
	
	

	
	vector<vector<double>> gradResiP;
	// math.calcGaussGreen(mesh, 0, resiVar, gradResiP);
	// math.calcLeastSquare(mesh, 0, resiVar, gradResiP);
	math.calcLeastSquare2nd(mesh, 0, resiVar, gradResiP);
	
	
	
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
		
		cell.var[controls.P] += controls.prePreURF*resiVar[i];
		
		double tmp1 = controls.timeStep/cell.var[controls.Rho];
		
		double resiU = controls.preVelURF*tmp1*gradResiP[i][0];
		double resiV = controls.preVelURF*tmp1*gradResiP[i][1];
		double resiW = controls.preVelURF*tmp1*gradResiP[i][2];
		
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
		
		residuals[0] += pow(resiVar[i],2.0)*cell.volume;
		residuals[1] += pow(resiU,2.0)*cell.volume;
		residuals[2] += pow(resiV,2.0)*cell.volume;
		residuals[3] += pow(resiW,2.0)*cell.volume;
		
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
	

	
	
	linA.clear();
	linAL.clear();
	linAR.clear();
	linB.clear();
	resiVar.clear();
	
	return sqrt(normP);
	
	
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








