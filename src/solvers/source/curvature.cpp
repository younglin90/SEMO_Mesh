#include "../build.h"
#include <cmath>
#include <array>
#include <numeric>


void SEMO_Solvers_Builder::calcCurvature(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double>& kappa){
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	SEMO_Utility_Math math;
	
	
	int iterVFSmoothingMax = 6;
	int iterVFSmoothing_Eikonal_Max = 6;
	int iterFirstSmoothingMax = 0;
	int iterSecondSmoothingMax = 0;
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	// bool boolNonOrthoCorrection = true;
	bool boolNonOrthoCorrection = false;
	
	int gradIterMax_LS = 1;
	int gradIterMax_GG = 0;
	
	
	
	vector<double> smoothAi(mesh.cells.size());
	for(int i=0; i<mesh.cells.size(); ++i){
		smoothAi[i] = mesh.cells[i].var[cn];
	}
	
	
	
	// //================================================
	// // smoothing Ai
	// for(int iter=0; iter<iterVFSmoothingMax; ++iter){

		// vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
		// // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// int dummy0;
		// for(int i=0; i<gradIterMax_LS; ++i){
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
				// cn, dummy0, smoothAi, gradAi);
		// }
		
		// vector<double> AiUp(mesh.cells.size(),0.0);
		// vector<double> AiDown(mesh.cells.size(),0.0);
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC; double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
				// double AiF = wCL*smoothAi[face.owner] +
							 // wCR*smoothAi[face.neighbour];
					
				// // skewness correction
				// if(boolSkewnessCorrection){
					// for(int ii=0; ii<3; ++ii){
						// AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
						// AiF += wCR * gradAi[face.neighbour][ii]*face.vecSkewness[ii];
					// }
					// AiF = max(0.0,min(1.0,AiF));
				// }
				
				// AiUp[face.owner] += AiF*face.area;
				// AiUp[face.neighbour] += AiF*face.area;
				
				// AiDown[face.owner] += face.area;
				// AiDown[face.neighbour] += face.area;
			
			// }
			
		// }

		// // boundary
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double AiF = smoothAi[face.owner];
					// for(int ii=0; ii<3; ++ii){
						// AiF += gradAi[face.owner][ii]*face.vecPF[ii];
					// }
					// AiF = max(0.0,min(1.0,AiF));
					
					// AiUp[face.owner] += AiF*face.area;
				
					// AiDown[face.owner] += face.area;
					
				// }
			// }
		// }
		
		// if(size>1){
			// // processor faces
			// vector<double> sendValues;
			// vector<vector<double>> sendGrad(3,vector<double>());
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// sendValues.push_back(smoothAi[face.owner]);
					// for(int ii=0; ii<3; ++ii){
						// sendGrad[ii].push_back(gradAi[face.owner][ii]);
					// }
				// }
			// }
			// vector<double> recvValues;
			// vector<vector<double>> recvGrad(3,vector<double>());
			// mpi.setProcsFaceDatas(
						// sendValues, recvValues,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int i=0; i<3; ++i){
				// mpi.setProcsFaceDatas(
							// sendGrad[i], recvGrad[i],
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
			// }
			// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double AiF = wCL*smoothAi[face.owner]+wCR*recvValues[ip];

					// // skewness correction
					// if(boolSkewnessCorrection){
						// for(int ii=0; ii<3; ++ii){
							// AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
							// AiF += wCR * recvGrad[ii][ip]*face.vecSkewness[ii];
						// }
						// AiF = max(0.0,min(1.0,AiF));
					// }
					
					// AiUp[face.owner] += AiF*face.area;
				
					// AiDown[face.owner] += face.area;
					
					// ++ip;
				// }
			// }
		// }
	
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// smoothAi[i] = AiUp[i]/AiDown[i];
		// }
		
	// }
	// //================================================
	
	
	
	
	//================================================
	// smoothing Ai -> level-set : solving Eikonal equation
	
	vector<bool> faceSmoothAi(mesh.faces.size(),false);
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		double wCL = face.wC;
		double wCR = 1.0-wCL;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			faceSmoothAi[i] = wCL*smoothAi[face.owner] + wCR*smoothAi[face.neighbour];
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			faceSmoothAi[i] = smoothAi[face.owner];
		}
	}
	// processor faces
	if(size>1){
		int proc_num=0;
		vector<double> var_send, var_recv;
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				var_send.push_back(smoothAi[face.owner]);
				++proc_num;
			}
		}
		mpi.setProcsFaceDatas(var_send, var_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			double wCL = face.wC;
			double wCR = 1.0-wCL;
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				faceSmoothAi[i] = wCL*smoothAi[face.owner] + wCR*var_recv[ip];
				++ip;
			}
		}
	}
	
	vector<bool> interfacePresentCell(mesh.cells.size(),false);
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if((smoothAi[face.owner]-0.5)*(faceSmoothAi[i]-0.5)<=0.0){
			interfacePresentCell[face.owner] = true;
		}
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			if((smoothAi[face.neighbour]-0.5)*(faceSmoothAi[i]-0.5)<=0.0){
				interfacePresentCell[face.neighbour] = true;
			}
		}
	}
	
	
	
	
	// vector<bool> interfacePresentCell(mesh.cells.size(),false);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// for(auto j : cell.stencil){
			// auto& cellSten = mesh.cells[j];
			// if((smoothAi[i]-0.5)*(smoothAi[j]-0.5)<0.0){
				// interfacePresentCell[i] = true;
				// break;
			// }
		// }
	// }
	
	// // processor faces
	// if(size>1){
		// int proc_num=0;
		// vector<double> var_send, var_recv;
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// var_send.push_back(smoothAi[face.owner]);
				// ++proc_num;
			// }
		// }
		// mpi.setProcsFaceDatas(var_send, var_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// vector<int> interfacePresentCell_send(proc_num,0);
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// for(auto j : face.stencil){
					// auto& cellSten = mesh.cells[j];
					// if((var_recv[ip]-0.5)*(smoothAi[j]-0.5)<0.0){
						// interfacePresentCell_send[ip] = 1;
						// break;
					// }
				// }
				// ++ip;
			// }
		// }
		// vector<int> interfacePresentCell_recv;
		// mpi.setProcsFaceDatas(interfacePresentCell_send, interfacePresentCell_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// if(interfacePresentCell_recv[ip]==1){
					// interfacePresentCell[face.owner] = true;
				// }
				// ++ip;
			// }
		// }
	
	// }
	
	
	// vof 함수를 level-set 함수로 변경
	for(int i=0; i<mesh.cells.size(); ++i){
		smoothAi[i] = smoothAi[i]-0.5;
	}
	
	// interface cell 의 level-set 값 재계산
	{
		vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(9,0.0));
		{
			int dummy0;
			math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
				-1, dummy0, smoothAi, gradSmoothAi);
		}
		for(int i=0; i<mesh.cells.size(); ++i){
			if(interfacePresentCell[i]){
				double magGrad = 0.0;
				magGrad += pow(gradSmoothAi[i][0],2.0);
				magGrad += pow(gradSmoothAi[i][1],2.0);
				magGrad += pow(gradSmoothAi[i][2],2.0);
				// magGrad = sqrt(magGrad);
				if(magGrad>1.e-8){
					smoothAi[i] = smoothAi[i] / magGrad;
				}
			}
		}
	}
	
	// double tau_ls = 1.e15;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// tau_ls = min(tau_ls,pow(mesh.cells[i].volume,0.3));
	// }
	// double tau_ls_glob;
	// MPI_Allreduce(&tau_ls, &tau_ls_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	// tau_ls = 0.1*tau_ls_glob;
	
	// 로컬 시간스템 구하기
	vector<double> tau_ls_local(mesh.cells.size(),0.0);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		// tau_ls_local[i] = 0.1 * pow(mesh.cells[i].volume,0.3);
		double maxA=0.0;
		double minA=1000000.0;
		for(auto& j : cell.faces){
			maxA = max(maxA, mesh.faces[j].area);
			minA = min(minA, mesh.faces[j].area);
		}
		tau_ls_local[i] = cell.volume / maxA / 1.0;
		tau_ls_local[i] = min(tau_ls_local[i],pow(cell.volume,0.3));
	}
	
	double tau_ls = 1.e15;
	for(int i=0; i<mesh.cells.size(); ++i){
		tau_ls = min(tau_ls,pow(mesh.cells[i].volume,0.3));
	}
	double tau_ls_glob;
	MPI_Allreduce(&tau_ls, &tau_ls_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	tau_ls = tau_ls_glob;
	
	
	
	
	vector<double> velocity0(mesh.cells.size(),0.0);
	for(int i=0; i<mesh.cells.size(); ++i){
		velocity0[i] = smoothAi[i] / sqrt(pow(smoothAi[i],2.0) +
			pow(abs(smoothAi[i])*pow(mesh.cells[i].volume,0.333),2.0));
	}
	vector<double> velocity0_recv;
	if(size>1){
		// processor faces
		vector<double> velocity0_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				velocity0_send.push_back(velocity0[face.owner]);
			}
		}
		mpi.setProcsFaceDatas(
					velocity0_send, velocity0_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	}
	
	// std::copy(smoothAi.begin(),smoothAi.end(),ini_smoothAi.begin());
	

	// // 그레디언트 계산
	// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(3,0.0));
	// {
		// int dummy0;
		// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
			// -1, dummy0, smoothAi, gradSmoothAi);
	// }
			
			
	// // HLL 스킴으로 계산, A. Chiapolino et al., 2021
	// for(int iter=0; iter< 5; ++iter){
		
		// vector<double> smoothAi_old(mesh.cells.size(),0.0);
		// std::copy(smoothAi.begin(),smoothAi.end(),smoothAi_old.begin());
		
		
		// // for(int iterIter=0; iterIter< 1; ++iterIter)
		// {
			// for(int rk=0; rk<2; ++rk)
			// {
				// double rkCoeff = 1.0;
				// if(rk==0) rkCoeff = 0.5;
				// if(rk==1) rkCoeff = 1.0;
			
				
				// // 그레디언트 계산
				// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(3,0.0));
				// {
					// int dummy0;
					// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
						// -1, dummy0, smoothAi, gradSmoothAi);
				// }
				// vector<double> smoothAi_recv;
				// vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
				// if(size>1){
					// // processor faces
					// vector<double> smoothAi_send;
					// vector<vector<double>> gradSmoothAi_send(3,vector<double>());
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						
						// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
							// smoothAi_send.push_back(smoothAi[face.owner]);
							// gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
							// gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
							// gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
						// }
					// }
					// mpi.setProcsFaceDatas(
								// smoothAi_send, smoothAi_recv,
								// mesh.countsProcFaces, mesh.countsProcFaces, 
								// mesh.displsProcFaces, mesh.displsProcFaces);
					// for(int ii=0; ii<3; ++ii){
						// mpi.setProcsFaceDatas(
									// gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
									// mesh.countsProcFaces, mesh.countsProcFaces, 
									// mesh.displsProcFaces, mesh.displsProcFaces);
					// }
				// }
						
						
				// //==================================================
				// // gradient 업데이트
				// // vector<vector<double>> gradResiduals(mesh.cells.size(),vector<double>(3,0.0));
				// vector<double> phiResiduals(mesh.cells.size(),0.0);
				// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
					// auto& face = mesh.faces[i];
							
					// vector<double> nvec(3,0.0);
					// std::copy(face.unitNormals.begin(),face.unitNormals.end(),nvec.begin());
					
					// double wCL = face.wC; double wCR = 1.0-wCL;
					
					// double velocity0_L = velocity0[face.owner];
					// double velocity0_R = 0.0;
					
					// double phiL = smoothAi[face.owner];
					// double phiR = 0.0;
					
					// vector<double> gradPhiL(3,0.0);
					// for(int ii=0; ii<3; ++ii){
						// gradPhiL[ii] = gradSmoothAi[face.owner][ii];
					// }
					// vector<double> gradPhiR(3,0.0);
					
					// if(face.getType() == SEMO_Types::INTERNAL_FACE){
						// velocity0_R = velocity0[face.neighbour];
						// phiR = smoothAi[face.neighbour];
						// for(int ii=0; ii<3; ++ii){
							// gradPhiR[ii] = gradSmoothAi[face.neighbour][ii];
						// }
					// }
					// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						// velocity0_R = velocity0_recv[ip];
						// phiR = smoothAi_recv[ip];
						// for(int ii=0; ii<3; ++ii){
							// gradPhiR[ii] = gradSmoothAi_recv[ii][ip];
						// }
					// }
					// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
						// velocity0_R = velocity0[face.owner];
						// phiR = smoothAi[face.owner];
						// for(int ii=0; ii<3; ++ii){
							// gradPhiR[ii] = gradSmoothAi[face.owner][ii];
						// }
					// }
					
					// double magUnL = sqrt(pow(gradPhiL[0],2.0)+pow(gradPhiL[1],2.0)+pow(gradPhiL[2],2.0));
					// double UL, VL, WL; UL = VL = WL = 0.0;
					// if(magUnL>1.e-8){
						// UL = velocity0_L * gradPhiL[0]/magUnL;
						// VL = velocity0_L * gradPhiL[1]/magUnL;
						// WL = velocity0_L * gradPhiL[2]/magUnL;
					// }
					// double UnL = UL*nvec[0]+VL*nvec[1]+WL*nvec[2];
					
					// double magUnR = sqrt(pow(gradPhiR[0],2.0)+pow(gradPhiR[1],2.0)+pow(gradPhiR[2],2.0));
					// double UR, VR, WR; UR = VR = WR = 0.0;
					// if(magUnR>1.e-8){
						// UR = velocity0_R * gradPhiR[0]/magUnR;
						// VR = velocity0_R * gradPhiR[1]/magUnR;
						// WR = velocity0_R * gradPhiR[2]/magUnR;
					// }
					// double UnR = UR*nvec[0]+VR*nvec[1]+WR*nvec[2];
					
					// double UnF = wCL*UnL+wCR*UnR;
					
					// double consL = gradPhiL[0]*nvec[0]+gradPhiL[1]*nvec[1]+gradPhiL[2]*nvec[2];
					// double consR = gradPhiR[0]*nvec[0]+gradPhiR[1]*nvec[1]+gradPhiR[2]*nvec[2];
					
					// double fluxL = gradPhiL[0]*UL+gradPhiL[1]*VL+gradPhiL[2]*WL;
					// double fluxR = gradPhiR[0]*UR+gradPhiR[1]*VR+gradPhiR[2]*WR;
					
					// double SL = min(UnL,UnR);
					// double SR = max(UnL,UnR);
					// double U_HLL = 0.0;
					// if(abs(SR-SL)>1.e-6){
						// U_HLL = (fluxL-fluxR+SR*consR-SL*consL) / (SR-SL);
					// }
					
					// // double flux = 0.0;
					// // if(U_HLL>0.0){
						// // flux = fluxL;
					// // }
					// // else if(U_HLL<0.0){
						// // flux = fluxR;
					// // }
					// double nF = 0.0;
					// double phiF = 0.0;
					// if(U_HLL>0.0){
						// phiF = phiL;
						// if(magUnL>1.e-8){
							// nF += gradPhiL[0]/magUnL*nvec[0];
							// nF += gradPhiL[1]/magUnL*nvec[1];
							// nF += gradPhiL[2]/magUnL*nvec[2];
						// }
					// }
					// // else if(U_HLL<0.0){
					// else{
						// phiF = phiR;
						// if(magUnR>1.e-8){
							// nF += gradPhiR[0]/magUnR*nvec[0];
							// nF += gradPhiR[1]/magUnR*nvec[1];
							// nF += gradPhiR[2]/magUnR*nvec[2];
						// }
					// }
					
					// // phiResiduals[face.owner] += (phiF-smoothAi_old[face.owner])*nF*face.area / mesh.cells[face.owner].volume;
					// phiResiduals[face.owner] += phiF*nF*face.area / mesh.cells[face.owner].volume;
					
					// if(face.getType() == SEMO_Types::INTERNAL_FACE){
						// // phiResiduals[face.neighbour] -= (phiF-smoothAi_old[face.neighbour])*nF*face.area / mesh.cells[face.neighbour].volume;
						// phiResiduals[face.neighbour] -= phiF*nF*face.area / mesh.cells[face.neighbour].volume;
					// }
					
					// // for(int ii=0; ii<3; ++ii){
						// // gradResiduals[face.owner][ii] -= flux*nvec[ii]*face.area;
					// // }
					
					// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
						// // for(int ii=0; ii<3; ++ii){
							// // gradResiduals[face.neighbour][ii] += flux*nvec[ii]*face.area;
						// // }
					// // }
					
					// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
					
					
				// }
					
				// double residual_ls = 0.0;
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if(!interfacePresentCell[i]){
						// // double resi_tmp = tau_ls_local[i] * ini_smoothAi[i] * (phiResiduals[i] / mesh.cells[i].volume - 1.0) ;
						// // double resi_tmp = rkCoeff * tau_ls * (phiResiduals[i] / mesh.cells[i].volume + velocity0[i]) ;
						// // double resi_tmp = rkCoeff * tau_ls * velocity0[i] * phiResiduals[i] / mesh.cells[i].volume;
						// double resi_tmp = -(smoothAi[i]-smoothAi_old[i]) + 
							// rkCoeff * tau_ls * ( -velocity0[i] * phiResiduals[i] + velocity0[i] );
						// // double resi_tmp = rkCoeff * tau_ls * ( -velocity0[i] * phiResiduals[i] + velocity0[i] );
						// // smoothAi[i] = smoothAi[i] + resi_tmp;
						
						// smoothAi[i] += resi_tmp;
						// residual_ls += resi_tmp*resi_tmp;
					// }
				// }
				// double residual_ls_glob = 0.0;
				// MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				// if(rank==0) cout << "| level-set re-initialization residual = " << residual_ls_glob << endl;
				
			// }
		// }
		
	// }
	
	
	
	
	
	
	
	
	
	// // HLL 스킴으로 계산, A. Chiapolino et al., 2021
	// for(int iter=0; iter< 5; ++iter){
		
		// vector<double> smoothAi_old(mesh.cells.size(),0.0);
		// std::copy(smoothAi.begin(),smoothAi.end(),smoothAi_old.begin());
		
		
		// // for(int iterIter=0; iterIter< 1; ++iterIter)
		// {
			// // for(int rk=0; rk<2; ++rk)
			// {
				// double rkCoeff = 1.0;
				// // if(rk==0) rkCoeff = 0.5;
				// // if(rk==1) rkCoeff = 1.0;
			
				
				// // 그레디언트 계산
				// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(3,0.0));
				// {
					// int dummy0;
					// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
						// -1, dummy0, smoothAi, gradSmoothAi);
				// }
				// vector<double> smoothAi_recv;
				// vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
				// if(size>1){
					// // processor faces
					// vector<double> smoothAi_send;
					// vector<vector<double>> gradSmoothAi_send(3,vector<double>());
					// for(int i=0; i<mesh.faces.size(); ++i){
						// auto& face = mesh.faces[i];
						
						// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
							// smoothAi_send.push_back(smoothAi[face.owner]);
							// gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
							// gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
							// gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
						// }
					// }
					// mpi.setProcsFaceDatas(
								// smoothAi_send, smoothAi_recv,
								// mesh.countsProcFaces, mesh.countsProcFaces, 
								// mesh.displsProcFaces, mesh.displsProcFaces);
					// for(int ii=0; ii<3; ++ii){
						// mpi.setProcsFaceDatas(
									// gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
									// mesh.countsProcFaces, mesh.countsProcFaces, 
									// mesh.displsProcFaces, mesh.displsProcFaces);
					// }
				// }
						
								
				// vector<vector<double>> phiResiduals(mesh.cells.size(),vector<double>(3,0.0));
				// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
					// auto& face = mesh.faces[i];
							
					// vector<double> nvec(3,0.0);
					// std::copy(face.unitNormals.begin(),face.unitNormals.end(),nvec.begin());
					
					// double wCL = face.wC; double wCR = 1.0-wCL;
					
					// double dPN = face.magPN;
					// double alpha = face.alphaF;
				
					// double dPhidN = 0.0;
					// if(face.getType() == SEMO_Types::INTERNAL_FACE){
						// dPhidN = alpha*(smoothAi[face.neighbour]-smoothAi[face.owner])/dPN;
						// for(int ii=0; ii<3; ++ii){
							// dPhidN += wCL*gradSmoothAi[face.owner][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
							// dPhidN += wCR*gradSmoothAi[face.neighbour][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
						// }
					// }
					// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						// dPhidN = alpha*(smoothAi_recv[ip]-smoothAi[face.owner])/dPN;
						// for(int ii=0; ii<3; ++ii){
							// dPhidN += wCL*gradSmoothAi[face.owner][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
							// dPhidN += wCR*gradSmoothAi_recv[ii][ip]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
						// }
					// }
					// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
						// dPhidN = 0.0;
					// }
					
					// {
						// vector<double> ai2(3,0.0);
						// for(int ii=0; ii<3; ++ii){
							// if(
							// (smoothAi[face.owner]>0.0 && face.unitNomalsPN[ii]*nvec[ii]>0.0) ||
							// (smoothAi[face.owner]<0.0 && face.unitNomalsPN[ii]*nvec[ii]<0.0) ){
								// ai2[ii] = min(0.0,dPhidN*nvec[ii]);
							// }
							// else if(
							// (smoothAi[face.owner]>0.0 && face.unitNomalsPN[ii]*nvec[ii]<0.0) ||
							// (smoothAi[face.owner]<0.0 && face.unitNomalsPN[ii]*nvec[ii]>0.0) ){
								// ai2[ii] = max(0.0,dPhidN*nvec[ii]);
							// }
							// phiResiduals[face.owner][ii] = max(phiResiduals[face.owner][ii],ai2[ii]*ai2[ii]);
						// }
					// }
					
					// if(face.getType() == SEMO_Types::INTERNAL_FACE)
					// {
						// vector<double> ai2(3,0.0);
						// for(int ii=0; ii<3; ++ii){
							// if(
							// (smoothAi[face.neighbour]>0.0 && (-face.unitNomalsPN[ii])*nvec[ii]>0.0) ||
							// (smoothAi[face.neighbour]<0.0 && (-face.unitNomalsPN[ii])*nvec[ii]<0.0) ){
								// ai2[ii] = min(0.0,dPhidN*nvec[ii]);
							// }
							// else if(
							// (smoothAi[face.neighbour]>0.0 && (-face.unitNomalsPN[ii])*nvec[ii]<0.0) ||
							// (smoothAi[face.neighbour]<0.0 && (-face.unitNomalsPN[ii])*nvec[ii]>0.0) ){
								// ai2[ii] = max(0.0,dPhidN*nvec[ii]);
							// }
							// phiResiduals[face.neighbour][ii] = max(phiResiduals[face.neighbour][ii],ai2[ii]*ai2[ii]);
						// }
					// }
					
					// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
					
				// }
				
					
				// double residual_ls = 0.0;
				// for(int i=0; i<mesh.cells.size(); ++i){
					// if(!interfacePresentCell[i]){
						// double magDelPhi = sqrt(
							// phiResiduals[i][0]*phiResiduals[i][0]+
							// phiResiduals[i][1]*phiResiduals[i][1]+
							// phiResiduals[i][2]*phiResiduals[i][2]);
						// double resi_tmp = -(smoothAi[i]-smoothAi_old[i]) + 
							// rkCoeff * tau_ls * ( -velocity0[i] * magDelPhi + velocity0[i] );
						
						// smoothAi[i] += resi_tmp;
						// residual_ls += resi_tmp*resi_tmp;
					// }
				// }
				// double residual_ls_glob = 0.0;
				// MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				// if(rank==0) cout << "| level-set re-initialization iter = " << iter << ", residual = " << residual_ls_glob << endl;
				
			// }
		// }
		
	// }
	
	
	
	
	
	



	
	for(int iter=0; iter< 8; ++iter){

		// =======================================
		// 그레디언트 구하기
		vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(3,0.0));
		{
			int dummy0;
			math.calcLeastSquare(mesh, "face", "1st", "input", 
				-1, dummy0, smoothAi, gradSmoothAi);
		}
		vector<double> smoothAi_recv;
		vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
		if(size>1){
			// processor faces
			vector<double> smoothAi_send;
			vector<vector<double>> gradSmoothAi_send(3,vector<double>());
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					smoothAi_send.push_back(smoothAi[face.owner]);
					gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
					gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
					gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
				}
			}
			mpi.setProcsFaceDatas(
						smoothAi_send, smoothAi_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			for(int ii=0; ii<3; ++ii){
				mpi.setProcsFaceDatas(
							gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
			}
		}
		
		// =======================================
		// 셀 surface normal vector 구하기
		vector<vector<double>> sufNorVec(mesh.cells.size(),vector<double>(3,0.0));
		// internal cells
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			double magGrad = sqrt(
				pow(gradSmoothAi[i][0],2.0)+
				pow(gradSmoothAi[i][1],2.0)+
				pow(gradSmoothAi[i][2],2.0));
			if(magGrad != 0.0){
				for(int ii=0; ii<3; ++ii){
					sufNorVec[i][ii] = gradSmoothAi[i][ii]/magGrad;
				}
			}
		}
		vector<vector<double>> sufNorVec_recv(3,vector<double>());
		if(size>1){
			// processor faces
			vector<vector<double>> sufNorVec_send(3,vector<double>());
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					sufNorVec_send[0].push_back(sufNorVec[face.owner][0]);
					sufNorVec_send[1].push_back(sufNorVec[face.owner][1]);
					sufNorVec_send[2].push_back(sufNorVec[face.owner][2]);
				}
			}
			for(int ii=0; ii<3; ++ii){
				mpi.setProcsFaceDatas(
							sufNorVec_send[ii], sufNorVec_recv[ii],
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
			}
		}
		
		// =======================================
		// 리스트 스퀘어 weight 구하기
		vector<double> weightLower(mesh.cells.size(),0.0);
		// internal cells
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
			
			double distX = face.distCells[0];
			double distY = face.distCells[1];
			double distZ = face.distCells[2];
			
			
			double wk = 0.0;
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][0]*(-distX);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][1]*(-distY);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][2]*(-distZ);
			wk = max(0.0,wk);
			
			weightLower[face.owner] += wk;
				
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				wk = 0.0;
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][0]*(distX);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][1]*(distY);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][2]*(distZ);
				wk = max(0.0,wk);
				
				weightLower[face.neighbour] += wk;
			}

		}
		// processor faces
		vector<double> weightLower_recv;
		if(size>1){
			vector<double> weightLower_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					weightLower_send.push_back(weightLower[face.owner]);
				}
			}
			mpi.setProcsFaceDatas(
						weightLower_send, weightLower_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
		}
		
		
		
		// =======================================
		// 리스트 스퀘어 초기화
		vector<vector<double>> vsum(mesh.cells.size(),vector<double>(6,0.0));
		
		for(auto& face : mesh.faces){
			
			double distX = face.distCells[0];
			double distY = face.distCells[1];
			double distZ = face.distCells[2];
			
			// weighting function
			// double wk = 1.0 / sqrt(pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			double wk = 0.0;
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][0]*(-distX);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][1]*(-distY);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][2]*(-distZ);
			wk = max(0.0,wk);
			if(abs(weightLower[face.owner])>1.e-8){
				wk /= weightLower[face.owner];
			}
			else{
				wk = 0.0;
			}
			
			vector<double> vari(3,0.0);
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
			
			// 대칭행렬
			for(int ii=0, num=0; ii<3; ++ii){
				for(int jj=0; jj<3; ++jj){
					if(jj>=ii){
						vsum[face.owner][num++] += wk * vari[ii]*vari[jj];
					}
				}
			}
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				wk = 0.0;
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][0]*(distX);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][1]*(distY);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][2]*(distZ);
				wk = max(0.0,wk);
				if(abs(weightLower[face.neighbour])>1.e-8){
					wk /= weightLower[face.neighbour];
				}
				else{
					wk = 0.0;
				}
				
				vari[0] *= (-1.0); vari[1] *= (-1.0); vari[2] *= (-1.0);
				
				for(int ii=0, num=0; ii<3; ++ii){
					for(int jj=0; jj<3; ++jj){
						if(jj>=ii){
							vsum[face.neighbour][num++] += wk * vari[ii]*vari[jj];
						}
					}
				}
			}
			
		}
		
		
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			// 대칭행렬의 역행렬도 대칭행렬
			vector<vector<double>> U(3,vector<double>(3,0.0));
			U[0][0] = vsum[i][0]; 
			U[0][1] = vsum[i][1];
			U[0][2] = vsum[i][2];
			U[1][0] = vsum[i][1];
			U[1][1] = vsum[i][3];
			U[1][2] = vsum[i][4];
			U[2][0] = vsum[i][2];
			U[2][1] = vsum[i][4];
			U[2][2] = vsum[i][5];
		
			// SVD 수도 역행렬 구하기
			vector<double> D(U[0].size());
			vector<vector<double>> V(U[0].size(), vector<double>(U[0].size(), 0.0));
			
			math.svdcmp(U, D, V);
			
			vector<vector<double>> invD(D.size(), vector<double>(D.size(), 0.0));
			for (int ii = 0; ii < D.size(); ++ii) {
				if (abs(D[ii]) > 1.e-8 ) {
					invD[ii][ii] = 1.0/D[ii];
				}
			}
	
			vector<vector<double>> transU;
			math.transpose(U, transU);
			vector<vector<double>> out;
			math.matmul(invD, transU, out);
			vector<vector<double>> pseudoInvMatrix;
			math.matmul(V, out, pseudoInvMatrix);
			
			vsum[i][0] = pseudoInvMatrix[0][0];
			vsum[i][1] = pseudoInvMatrix[0][1];
			vsum[i][2] = pseudoInvMatrix[0][2];
			vsum[i][3] = pseudoInvMatrix[1][1];
			vsum[i][4] = pseudoInvMatrix[1][2];
			vsum[i][5] = pseudoInvMatrix[2][2];
			
		}
		
		
		// =======================================
		// 그레디언트 계산
		vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
		
		
		// internal cells
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
	

			double distX = face.distCells[0];
			double distY = face.distCells[1];
			double distZ = face.distCells[2];
			
			// double wk = 1.0 / sqrt(pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			double wk = 0.0;
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][0]*(-distX);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][1]*(-distY);
			wk += (velocity0[face.owner] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.owner][2]*(-distZ);
			wk = max(0.0,wk);
			if(abs(weightLower[face.owner])>1.e-8){
				wk /= weightLower[face.owner];
			}
			else{
				wk = 0.0;
			}
			
			double DVar = 0.0;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				DVar = smoothAi[face.neighbour] - smoothAi[face.owner];
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				DVar = smoothAi_recv[ip] - smoothAi[face.owner];
			}
		
			vector<double> vari(3,0.0);
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
	
			for(int ii=0; ii<3; ++ii){
				grad_tmp[face.owner][ii] += wk * vari[ii] * DVar;
			}
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				wk = 0.0;
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][0]*(distX);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][1]*(distY);
				wk += (velocity0[face.neighbour] > 0.0 ? 1.0 : -1.0) * sufNorVec[face.neighbour][2]*(distZ);
				wk = max(0.0,wk);
				if(abs(weightLower[face.neighbour])>1.e-8){
					wk /= weightLower[face.neighbour];
				}
				else{
					wk = 0.0;
				}
				
				vari[0] *= (-1.0); vari[1] *= (-1.0); vari[2] *= (-1.0);
				DVar *= (-1.0);
				for(int ii=0; ii<3; ++ii){
					grad_tmp[face.neighbour][ii] += wk * vari[ii] * DVar;
				}
			}
				
			if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
		}
		
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double tmp0 = 
				vsum[i][0] * grad_tmp[i][0] +
				vsum[i][1] * grad_tmp[i][1] +
				vsum[i][2] * grad_tmp[i][2];
				
			double tmp1 = 
				vsum[i][1] * grad_tmp[i][0] +
				vsum[i][3] * grad_tmp[i][1] +
				vsum[i][4] * grad_tmp[i][2];
				
			double tmp2 = 
				vsum[i][2] * grad_tmp[i][0] +
				vsum[i][4] * grad_tmp[i][1] +
				vsum[i][5] * grad_tmp[i][2];
				
			grad_tmp[i][0] = tmp0;
			grad_tmp[i][1] = tmp1;
			grad_tmp[i][2] = tmp2;
			
		}
		
		// =======================================
		// Hamilton-Jacobi 방정식 풀기
		double residual_ls = 0.0;
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];

			if(!interfacePresentCell[i]){
				double magGrad = sqrt(
					pow(grad_tmp[i][0],2.0)+
					pow(grad_tmp[i][1],2.0)+
					pow(grad_tmp[i][2],2.0));
				// double magGrad = sqrt(
					// pow(gradSmoothAi[i][0],2.0)+
					// pow(gradSmoothAi[i][1],2.0)+
					// pow(gradSmoothAi[i][2],2.0));
					
				double resi_tmp = 0.1 * tau_ls_local[i]*velocity0[i]*(1.0-magGrad);
				// double resi_tmp = 0.1 * tau_ls*velocity0[i]*(1.0-magGrad);
				smoothAi[i] += resi_tmp;
				residual_ls += resi_tmp*resi_tmp;
			}
		}
		double residual_ls_glob = 0.0;
		MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if(rank==0) cout << "| level-set re-initialization residual = " << residual_ls_glob << endl;

	}





	
	// // ============================================
	// for(int iter=0; iter< 500; ++iter){
		
		// // 올드 값 저장
		// vector<double> smoothAi_old(mesh.cells.size(),0.0);
		// std::copy(smoothAi.begin(),smoothAi.end(),smoothAi_old.begin());
		
		// // 그레디언트 구하기
		// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(9,0.0));
		// {
			// int dummy0;
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
				// -1, dummy0, smoothAi, gradSmoothAi);
		// }
		// vector<double> smoothAi_recv;
		// vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
		// if(size>1){
			// // processor faces
			// vector<double> smoothAi_send;
			// vector<vector<double>> gradSmoothAi_send(3,vector<double>());
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// smoothAi_send.push_back(smoothAi[face.owner]);
					// gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
					// gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
					// gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
				// }
			// }
			// mpi.setProcsFaceDatas(
						// smoothAi_send, smoothAi_recv,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int ii=0; ii<3; ++ii){
				// mpi.setProcsFaceDatas(
							// gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
			// }
		// }
				
		// // HLL 플럭스 계산
		// vector<double> residuals(mesh.cells.size(),0.0);
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
					
			// vector<double> nvec(3,0.0);
			// nvec[0] = face.unitNormals[0];
			// nvec[1] = face.unitNormals[1];
			// nvec[2] = face.unitNormals[2];
			
			// double wCL = face.wC; double wCR = 1.0-wCL;
			
			// double phiL = smoothAi[face.owner];
			// double phiR = 0.0;
			// double vel0L = velocity0[face.owner];
			// double vel0R = 0.0;
			// vector<double> gradPhiL(3,0.0);
			// vector<double> gradPhiR(3,0.0);
			// for(int ii=0; ii<3; ++ii){
				// gradPhiL[ii] = gradSmoothAi[face.owner][ii];
			// }
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// phiR = smoothAi[face.neighbour];
				// vel0R = velocity0[face.neighbour];
				// for(int ii=0; ii<3; ++ii){
					// gradPhiR[ii] = gradSmoothAi[face.neighbour][ii];
				// }
			// }
			// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// phiR = smoothAi_recv[ip];
				// vel0R = velocity0_recv[ip];
				// for(int ii=0; ii<3; ++ii){
					// gradPhiR[ii] = gradSmoothAi_recv[ii][ip];
				// }
			// }
			// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				// phiR = smoothAi[face.owner];
				// vel0R = velocity0[face.owner];
				// for(int ii=0; ii<3; ++ii){
					// gradPhiR[ii] = gradSmoothAi[face.owner][ii];
				// }
			// }
			
			// double magPhiL = sqrt(pow(gradPhiL[0],2.0)+pow(gradPhiL[1],2.0)+pow(gradPhiL[2],2.0));
			// double magPhiR = sqrt(pow(gradPhiR[0],2.0)+pow(gradPhiR[1],2.0)+pow(gradPhiR[2],2.0));
			
			// if(magPhiL>1.e-6 && magPhiR>1.e-6){
				// double surfNVecL = (gradPhiL[0]*nvec[0] + gradPhiL[1]*nvec[1] + gradPhiL[2]*nvec[2])/magPhiL;
				// double surfNVecR = (gradPhiR[0]*nvec[0] + gradPhiR[1]*nvec[1] + gradPhiR[2]*nvec[2])/magPhiR;
				// double UnL = vel0L*surfNVecL;
				// double UnR = vel0R*surfNVecR;
				
				// double consL = gradPhiL[0]*nvec[0] + gradPhiL[1]*nvec[1] + gradPhiL[2]*nvec[2];
				// double consR = gradPhiR[0]*nvec[0] + gradPhiR[1]*nvec[1] + gradPhiR[2]*nvec[2];
				// double fluxL = vel0L*(pow(gradPhiL[0],2.0)+pow(gradPhiL[1],2.0)+pow(gradPhiL[2],2.0))/magPhiL;
				// double fluxR = vel0R*(pow(gradPhiR[0],2.0)+pow(gradPhiR[1],2.0)+pow(gradPhiR[2],2.0))/magPhiR;
				
				// double SL = min(UnL,UnR);
				// double SR = max(UnL,UnR);
				// if(abs(SR-SL)>1.e-6) {
					// double UStar = (SR*UnR-SL*UnL)/(SR-SL)-(fluxR-fluxL)/(SR-SL);
					
					// double flux0=0.0;
					// double flux1=0.0;
					// if(UStar>0.0){
						// flux0 = phiL*surfNVecL*face.area;
						// flux1 = surfNVecL*face.area;
					// }
					// else{
						// flux0 = phiR*surfNVecR*face.area;
						// flux1 = surfNVecR*face.area;
					// }
					
					// residuals[face.owner] -= vel0L*(flux0-smoothAi_old[face.owner]*flux1);
					// if(face.getType() == SEMO_Types::INTERNAL_FACE){
						// residuals[face.neighbour] += vel0R*(flux0-smoothAi_old[face.neighbour]*flux1);
					// }
				// }
			// }
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
			
		// }
		
		
		// double residual_ls = 0.0;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if(!interfacePresentCell[i]){
				// double resi_tmp = tau_ls_local[i] * (residuals[i]/mesh.cells[i].volume + velocity0[i]);
				// smoothAi[i] += 0.1*resi_tmp;
				// residual_ls += resi_tmp*resi_tmp;
			// }
		// }
		// double residual_ls_glob = 0.0;
		// MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// if(rank==0) cout << "| level-set re-initialization residual = " << residual_ls_glob << endl;
		
	// }
	
	
	
	
	
	// // ============================================
	// // Godunov flux 계산
	// for(int iter=0; iter< 500; ++iter){
		// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(9,0.0));
		// {
			// int dummy0;
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
				// -1, dummy0, smoothAi, gradSmoothAi);
		// }
		// vector<double> smoothAi_recv;
		// vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
		// if(size>1){
			// // processor faces
			// vector<double> smoothAi_send;
			// vector<vector<double>> gradSmoothAi_send(3,vector<double>());
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// smoothAi_send.push_back(smoothAi[face.owner]);
					// gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
					// gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
					// gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
				// }
			// }
			// mpi.setProcsFaceDatas(
						// smoothAi_send, smoothAi_recv,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int ii=0; ii<3; ++ii){
				// mpi.setProcsFaceDatas(
							// gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
			// }
		// }
				
		// vector<vector<double>> grad2Ai(mesh.cells.size(),vector<double>(3,0.0));
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
					
			// vector<double> nvec(3,0.0);
			// nvec[0] = face.unitNormals[0];
			// nvec[1] = face.unitNormals[1];
			// nvec[2] = face.unitNormals[2];
			
			// double wCL = face.wC;
			// double wCR = 1.0-wCL;
			
			// double dPN = face.magPN;
			// double alpha = face.alphaF;
		
			// double smoothAiL = 0.0;
			// double smoothAiR = 0.0;
			// double dPhidN = 0.0;
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// smoothAiL = smoothAi[face.owner];
				// smoothAiR = smoothAi[face.neighbour];
				// dPhidN = alpha*(smoothAiR-smoothAiL)/dPN;
				// for(int ii=0; ii<3; ++ii){
					// smoothAiL += gradSmoothAi[face.owner][ii]*(face.vecPF[ii]);
					// smoothAiR += gradSmoothAi[face.neighbour][ii]*(face.vecNF[ii]);
					// dPhidN += wCL*gradSmoothAi[face.owner][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
					// dPhidN += wCR*gradSmoothAi[face.neighbour][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
				// }
			// }
			// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// smoothAiL = smoothAi[face.owner];
				// smoothAiR = smoothAi_recv[ip];
				// dPhidN = alpha*(smoothAiR-smoothAiL)/dPN;
				// for(int ii=0; ii<3; ++ii){
					// smoothAiL += gradSmoothAi[face.owner][ii]*(face.vecPF[ii]);
					// smoothAiR += gradSmoothAi_recv[ii][ip]*(face.vecNF[ii]);
					// dPhidN += wCL*gradSmoothAi[face.owner][ii]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
					// dPhidN += wCR*gradSmoothAi_recv[ii][ip]*(nvec[ii]-alpha*face.unitNomalsPN[ii]);
				// }
			// }
			 
			
			// {
				// vector<double> ai2(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// // if(
					// // (smoothAiL>0.0 && face.unitNomalsPN[ii]>0.0) ||
					// // (smoothAiL<0.0 && face.unitNomalsPN[ii]<0.0) ){
						// // ai2[ii] = min(0.0,dPhidN*nvec[ii]);
					// // }
					// // else if(
					// // (smoothAiL>0.0 && face.unitNomalsPN[ii]<0.0) ||
					// // (smoothAiL<0.0 && face.unitNomalsPN[ii]>0.0) ){
						// // ai2[ii] = max(0.0,dPhidN*nvec[ii]);
					// // }
					// if(
					// (smoothAiL>0.0 && face.unitNomalsPN[ii]>0.0) ||
					// (smoothAiL<0.0 && face.unitNomalsPN[ii]<0.0) ){
						// ai2[ii] = min(0.0,dPhidN*nvec[ii]);
					// }
					// else if(
					// (smoothAiL>0.0 && face.unitNomalsPN[ii]<0.0) ||
					// (smoothAiL<0.0 && face.unitNomalsPN[ii]>0.0) ){
						// ai2[ii] = max(0.0,dPhidN*nvec[ii]);
					// }
					// grad2Ai[face.owner][ii] = max(grad2Ai[face.owner][ii],ai2[ii]*ai2[ii]);
				// }
			// }
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE)
			// {
				// vector<double> ai2(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// if(
					// (smoothAiR>0.0 && (-face.unitNomalsPN[ii])>0.0) ||
					// (smoothAiR<0.0 && (-face.unitNomalsPN[ii])<0.0) ){
						// ai2[ii] = min(0.0,dPhidN*(nvec[ii]));
					// }
					// else if(
					// (smoothAiR>0.0 && (-face.unitNomalsPN[ii])<0.0) ||
					// (smoothAiR<0.0 && (-face.unitNomalsPN[ii])>0.0) ){
						// ai2[ii] = max(0.0,dPhidN*(nvec[ii]));
					// }
					// grad2Ai[face.neighbour][ii] = max(grad2Ai[face.neighbour][ii],ai2[ii]*ai2[ii]);
				// }
			// }
			
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
			
		// }
		
		
		// double residual_ls = 0.0;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if(!interfacePresentCell[i]){
				// double magGrad = sqrt(grad2Ai[i][0]+grad2Ai[i][1]+grad2Ai[i][2]);
				// double resi_tmp = tau_ls_local[i] * velocity0[i] * (1.0-magGrad);
				// smoothAi[i] += 0.1 * resi_tmp;
				// residual_ls += resi_tmp*resi_tmp;
			// }
		// }
		// double residual_ls_glob = 0.0;
		// MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// if(rank==0) cout << "| level-set re-initialization residual = " << residual_ls_glob << endl;
		
	// }
	
	
	
	
	
	
	

	
	// // ============================================
	// // Lax-Friedrichs flux 계산
	// for(int iter=0; iter< 500; ++iter){
		// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(9,0.0));
		// {
			// int dummy0;
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "input", 
				// -1, dummy0, smoothAi, gradSmoothAi);
		// }
		// vector<double> smoothAi_recv;
		// vector<vector<double>> gradSmoothAi_recv(3,vector<double>());
		// if(size>1){
			// // processor faces
			// vector<double> smoothAi_send;
			// vector<vector<double>> gradSmoothAi_send(3,vector<double>());
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// smoothAi_send.push_back(smoothAi[face.owner]);
					// gradSmoothAi_send[0].push_back(gradSmoothAi[face.owner][0]);
					// gradSmoothAi_send[1].push_back(gradSmoothAi[face.owner][1]);
					// gradSmoothAi_send[2].push_back(gradSmoothAi[face.owner][2]);
				// }
			// }
			// mpi.setProcsFaceDatas(
						// smoothAi_send, smoothAi_recv,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int ii=0; ii<3; ++ii){
				// mpi.setProcsFaceDatas(
							// gradSmoothAi_send[ii], gradSmoothAi_recv[ii],
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
			// }
		// }
				
		// vector<vector<double>> grad2Ai(mesh.cells.size(),vector<double>(3,0.0));
		// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
					
			// vector<double> nvec(3,0.0);
			// nvec[0] = face.unitNormals[0];
			// nvec[1] = face.unitNormals[1];
			// nvec[2] = face.unitNormals[2];
			
			// double wCL = face.wC;
			// double wCR = 1.0-wCL;
			
			// double dPN = face.magPN;
			// double alpha = face.alphaF;
		
			// double smoothAiL = 0.0;
			// double smoothAiR = 0.0;
			// double dPhidN = 0.0;
			// double HemiltonL = 0.0;
			// double HemiltonR = 0.0;
			// vector<double> gradL(3,0.0);
			// vector<double> gradR(3,0.0);
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// smoothAiL = smoothAi[face.owner];
				// smoothAiR = smoothAi[face.neighbour];
				// for(int ii=0; ii<3; ++ii){
					// gradL[ii] = gradSmoothAi[face.owner][ii];
					// gradR[ii] = gradSmoothAi[face.neighbour][ii];
					// dPhidN += pow(0.5*(gradSmoothAi[face.owner][ii]+gradSmoothAi[face.neighbour][ii]),2.0);
					// HemiltonL += pow(gradSmoothAi[face.owner][ii],2.0);
					// HemiltonR += pow(gradSmoothAi[face.neighbour][ii],2.0);
				// }
				// dPhidN = sqrt(dPhidN);
				// HemiltonL = sqrt(HemiltonL);
				// HemiltonR = sqrt(HemiltonR);
			// }
			// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// smoothAiL = smoothAi[face.owner];
				// smoothAiR = smoothAi_recv[ip];
				// for(int ii=0; ii<3; ++ii){
					// gradL[ii] = gradSmoothAi[face.owner][ii];
					// gradR[ii] = gradSmoothAi_recv[ii][ip];
					// dPhidN += pow(0.5*(gradSmoothAi[face.owner][ii]+gradSmoothAi_recv[ii][ip]),2.0);
					// HemiltonL += pow(gradSmoothAi[face.owner][ii],2.0);
					// HemiltonR += pow(gradSmoothAi_recv[ii][ip],2.0);
				// }
				// dPhidN = sqrt(dPhidN);
				// HemiltonL = sqrt(HemiltonL);
				// HemiltonR = sqrt(HemiltonR);
			// }
			// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				// smoothAiL = smoothAi[face.owner];
				// smoothAiR = smoothAiL;
				// for(int ii=0; ii<3; ++ii){
					// gradL[ii] = gradSmoothAi[face.owner][ii];
					// gradR[ii] = gradSmoothAi[face.owner][ii];
					// dPhidN += pow(0.5*(gradSmoothAi[face.owner][ii]+gradSmoothAi[face.owner][ii]),2.0);
					// HemiltonL += pow(gradSmoothAi[face.owner][ii],2.0);
					// HemiltonR += pow(gradSmoothAi[face.owner][ii],2.0);
				// }
				// dPhidN = sqrt(dPhidN);
				// HemiltonL = sqrt(HemiltonL);
				// HemiltonR = sqrt(HemiltonR);
			// }
			
			// double maxHemiltonX = max(abs(HemiltonL*nvec[0]),abs(HemiltonR*nvec[0]));
			// double maxHemiltonY = max(abs(HemiltonL*nvec[1]),abs(HemiltonR*nvec[1]));
			// double maxHemiltonZ = max(abs(HemiltonL*nvec[2]),abs(HemiltonR*nvec[2]));
			
			// dPhidN -= maxHemiltonX*0.5*(gradR[0]-gradL[0]);
			// dPhidN -= maxHemiltonY*0.5*(gradR[1]-gradL[1]);
			// dPhidN -= maxHemiltonZ*0.5*(gradR[2]-gradL[2]);
			
			// grad2Ai[face.owner][0] += dPhidN * face.area / mesh.cells[face.owner].volume;
					
			// if(face.getType() == SEMO_Types::INTERNAL_FACE)
			// {
				// grad2Ai[face.neighbour][0] -= dPhidN * face.area / mesh.cells[face.neighbour].volume;
			// }
			
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
			
		// }
		
		
		// double residual_ls = 0.0;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if(!interfacePresentCell[i]){
				// double resi_tmp = tau_ls_local[i] * velocity0[i] * (1.0-grad2Ai[i][0]);
				// smoothAi[i] += 0.1 * resi_tmp;
				// residual_ls += resi_tmp*resi_tmp;
			// }
		// }
		// double residual_ls_glob = 0.0;
		// MPI_Allreduce(&residual_ls, &residual_ls_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// if(rank==0) cout << "| level-set re-initialization residual = " << residual_ls_glob << endl;
		
	// }
	
	
	
	
	
	
	
	
	
	
	
	
	// for(int i=0; i< 10 ; ++i){
		// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(9,0.0));
		// {
			// // math.calcGaussGreen(mesh, cn, smoothAi, gradSmoothAi);
			// int dummy0;
			// math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
				// -1, dummy0, smoothAi, gradSmoothAi);
		// }
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if(!interfacePresentCell[i]){
				// double magGrad = 0.0;
				// magGrad += pow(gradSmoothAi[i][0],2.0);
				// magGrad += pow(gradSmoothAi[i][1],2.0);
				// magGrad += pow(gradSmoothAi[i][2],2.0);
				// magGrad = sqrt(magGrad);
				// smoothAi[i] += tau_ls_local[i] * ( ini_smoothAi[i] > 0.0 ? 1.0 : -1.0 ) * (1.0-magGrad);
			// }
		// }
	// }
	//================================================
	
	
	// double delH = 0.0;
	// int interCellNum = 0;
	// double tau_ls = 1.e15;
	// vector<double> ini_smoothAi(mesh.cells.size(),0.0);
	// std::copy(smoothAi.begin(),smoothAi.end(),ini_smoothAi.begin());
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(interfacePresentCell[i]){
			// delH += pow(mesh.cells[i].volume,0.3);
			// ++interCellNum;
		// }
		// tau_ls = min(tau_ls,pow(mesh.cells[i].volume,0.3));
	// }
	// delH /= (double)interCellNum;
	
	// double eta_ls = 0.1*delH;
	// for(int i=0; i<8; ++i){
		// vector<vector<double>> gradSmoothAi(mesh.cells.size(),vector<double>(3,0.0));
		// {
			// int dummy0;
			// math.calcLeastSquare(mesh, "face", "1st", "input", 
				// -1, dummy0, smoothAi, gradSmoothAi);
		// }
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if(!interfacePresentCell[i]){
				// double magGrad = 0.0;
				// magGrad += pow(gradSmoothAi[i][0],2.0);
				// magGrad += pow(gradSmoothAi[i][1],2.0);
				// magGrad += pow(gradSmoothAi[i][2],2.0);
				// magGrad = sqrt(magGrad);
				// double coeff_ls = smoothAi[i]*tau_ls/sqrt(smoothAi[i]*smoothAi[i]+eta_ls*eta_ls);
				// smoothAi[i] += coeff_ls*(1.0-magGrad);
			// }
		// }
	// }
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// vector<double> Ai_recv;
	// if(size>1){
		// // processor faces
		// vector<double> sendValues;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// sendValues.push_back(mesh.cells[face.owner].var[cn]);
			// }
		// }
		// mpi.setProcsFaceDatas(
					// sendValues, Ai_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
	// }
	
	
	
	
	vector<double> smoothAi_recv;
	if(size>1){
		// processor faces
		vector<double> sendValues;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues.push_back(smoothAi[face.owner]);
			}
		}
		mpi.setProcsFaceDatas(
					sendValues, smoothAi_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	}
	
	
	//================================================
	// calc gauss-green gradient
	// vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
	vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(9,0.0));
	// for(int i=0; i<gradIterMax_LS; ++i)
	{
		int dummy0;
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
			-1, dummy0, smoothAi, gradAi);
	}
	// for(int i=0; i<gradIterMax_GG; ++i){
		// math.calcGaussGreen(mesh, -1, smoothAi, gradAi);
	// }
	
	
	
	// // exact curvature formula
	// vector<double> kappa_test(mesh.cells.size(),0.0);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// vector<vector<double>> grad_phi(1,vector<double>(3,0.0));
		
		// grad_phi[0][0] = gradAi[i][0];
		// grad_phi[0][1] = gradAi[i][1];
		// grad_phi[0][2] = gradAi[i][2];
		
		// vector<vector<double>> trans_grad_phi(3,vector<double>(1,0.0));
		
		// trans_grad_phi[0][0] = gradAi[i][0];
		// trans_grad_phi[1][0] = gradAi[i][1];
		// trans_grad_phi[2][0] = gradAi[i][2];
		
		// // 3 4 5
		  // // 6 7
		    // // 8
		// vector<vector<double>> adjoint_hessian_phi(3,vector<double>(3,0.0));
		
		// adjoint_hessian_phi[0][0] = gradAi[i][6]*gradAi[i][8] - gradAi[i][7]*gradAi[i][7];
		// adjoint_hessian_phi[0][1] = gradAi[i][7]*gradAi[i][5] - gradAi[i][4]*gradAi[i][7];
		// adjoint_hessian_phi[0][2] = gradAi[i][4]*gradAi[i][7] - gradAi[i][6]*gradAi[i][5];
		
		// adjoint_hessian_phi[1][0] = adjoint_hessian_phi[0][1];
		// adjoint_hessian_phi[1][1] = gradAi[i][3]*gradAi[i][8] - gradAi[i][5]*gradAi[i][5];
		// adjoint_hessian_phi[1][2] = gradAi[i][4]*gradAi[i][5] - gradAi[i][3]*gradAi[i][7];
		
		// adjoint_hessian_phi[2][0] = adjoint_hessian_phi[0][2];
		// adjoint_hessian_phi[2][1] = adjoint_hessian_phi[1][2];
		// adjoint_hessian_phi[2][2] = gradAi[i][3]*gradAi[i][6] - gradAi[i][4]*gradAi[i][4];
		
		// double tmp0 = 0.0;
		// tmp0 += adjoint_hessian_phi[0][0]*trans_grad_phi[0][0];
		// tmp0 += adjoint_hessian_phi[0][1]*trans_grad_phi[1][0];
		// tmp0 += adjoint_hessian_phi[0][2]*trans_grad_phi[2][0];
		
		// double tmp1 = 0.0;
		// tmp1 += adjoint_hessian_phi[1][0]*trans_grad_phi[0][0];
		// tmp1 += adjoint_hessian_phi[1][1]*trans_grad_phi[1][0];
		// tmp1 += adjoint_hessian_phi[1][2]*trans_grad_phi[2][0];
		
		// double tmp2 = 0.0;
		// tmp2 += adjoint_hessian_phi[2][0]*trans_grad_phi[0][0];
		// tmp2 += adjoint_hessian_phi[2][1]*trans_grad_phi[1][0];
		// tmp2 += adjoint_hessian_phi[2][2]*trans_grad_phi[2][0];
		
		// double upper = grad_phi[0][0]*tmp0 + grad_phi[0][1]*tmp1 + grad_phi[0][2]*tmp2;
		
		// // // cout << "AA" << endl;
		// // // vector<vector<double>> trans_grad_phi;
		// // // math.transpose(grad_phi, trans_grad_phi);
		// // // cout << "BB" << endl;
		// // vector<vector<double>> out;
		// // math.matmul(adjoint_hessian_phi, trans_grad_phi, out);
		// // // cout << "CC" << endl;
		// // vector<vector<double>> out2;
		// // math.matmul(grad_phi, out, out2);
		// // // cout << out2[0][0] << endl;
		
		// double mag = 0.0;
		// mag += pow(grad_phi[0][0],2.0);
		// mag += pow(grad_phi[0][1],2.0);
		// mag += pow(grad_phi[0][2],2.0);
		
		// // cout << "DD" << endl;
		// if(mag > 1.e-6){
			// kappa_test[i] = upper / pow(mag,2.0);
		// }
		// else{
			// kappa_test[i] = 0.0;
		// }
		
		
		
		
	// }
			
	
	
	
	// vector<vector<double>> gradAi_LS(mesh.cells.size(),vector<double>(3,0.0));
	// math.calcLeastSquare2nd(mesh, -1, smoothAi, gradAi_LS);
	// vector<vector<double>> gradAi_LS_recv(3,vector<double>());
	// if(size>1){
		// // processor faces
		// vector<vector<double>> gradAi_LS_send(3,vector<double>());
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradAi_LS_send[0].push_back(gradAi_LS[face.owner][0]);
				// gradAi_LS_send[1].push_back(gradAi_LS[face.owner][1]);
				// gradAi_LS_send[2].push_back(gradAi_LS[face.owner][2]);
			// }
		// }
		// for(int i=0; i<3; ++i){
			// mpi.setProcsFaceDatas(
						// gradAi_LS_send[i], gradAi_LS_recv[i],
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
		// }
	// }
	// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double wCL = face.wC; double wCR = 1.0 - wCL;
		// double dPN = face.magPN;
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// for(int ii=0; ii<3; ++ii){
				// gradAi[face.owner][ii] += (smoothAi[face.neighbour]-smoothAi[face.owner])/dPN * 
					// face.vecPF[ii] * face.area / mesh.cells[face.owner].volume;
				// gradAi[face.neighbour][ii] += (smoothAi[face.owner]-smoothAi[face.neighbour])/dPN * 
					// face.vecNF[ii] * face.area / mesh.cells[face.neighbour].volume;
			// }
			// // non-orthogonal
			// if(boolNonOrthoCorrection){
				// vector<double> vecGrad(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// vecGrad[ii] += wCL*gradAi_LS[face.owner][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);
					// vecGrad[ii] += wCR*gradAi_LS[face.neighbour][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);	
				// }
				// for(int ii=0; ii<3; ++ii){
					// gradAi[face.owner][ii] += vecGrad[ii] * face.vecPF[ii] * face.area / mesh.cells[face.owner].volume;
					// gradAi[face.neighbour][ii] -= vecGrad[ii] * face.vecNF[ii] * face.area / mesh.cells[face.neighbour].volume;
				// }
			// }
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// for(int ii=0; ii<3; ++ii){
				// gradAi[face.owner][ii] += (smoothAi_recv[ip]-smoothAi[face.owner])/dPN * 
					// face.vecPF[ii] * face.area / mesh.cells[face.owner].volume;
			// }
			// // non-orthogonal
			// if(boolNonOrthoCorrection){
				// vector<double> vecGrad(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// vecGrad[ii] += wCL*gradAi_LS[face.owner][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);
					// vecGrad[ii] += wCR*gradAi_LS_recv[ii][ip]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);	
				// }
				// for(int ii=0; ii<3; ++ii){
					// gradAi[face.owner][ii] += vecGrad[ii] * face.vecPF[ii] * face.area / mesh.cells[face.owner].volume;
				// }
			// }
			
			// ++ip;
		// }
	// }



	//================================================
	// surfNormalVec = dAidx/|dAidx|
	vector<vector<double>> surfNormalVec(3,vector<double>(mesh.cells.size()));
	
	for(int i=0; i<mesh.cells.size(); ++i){
		double magGrad = 0.0;
		for(int ii=0; ii<3; ++ii){
			magGrad += gradAi[i][ii]*gradAi[i][ii];
		}
		magGrad = sqrt(magGrad);
			
		// if( magGrad > std::numeric_limits<double>::min() ){
		if( magGrad > 1.e-8 ){
			for(int ii=0; ii<3; ++ii){
				surfNormalVec[ii][i] = gradAi[i][ii]/magGrad;
			}
		}
		else{
			for(int ii=0; ii<3; ++ii) surfNormalVec[ii][i] = 0.0;
		}
		
	}
	vector<vector<double>> surfNormalVec_recv(3,vector<double>());
	if(size>1){
		// processor faces
		vector<vector<double>> sendValues(3,vector<double>());
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues[0].push_back(surfNormalVec[0][face.owner]);
				sendValues[1].push_back(surfNormalVec[1][face.owner]);
				sendValues[2].push_back(surfNormalVec[2][face.owner]);
			}
		}
		for(int i=0; i<3; ++i){
			mpi.setProcsFaceDatas(
						sendValues[i], surfNormalVec_recv[i],
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
		}
	}
	
	
	
	
	
	//================================================
	// calc surface normal vector gradients
	
	
	vector<vector<vector<double>>> gradSurfNor_LS;
	gradSurfNor_LS.resize(3);
	for(int i=0; i<3; ++i){
		gradSurfNor_LS[i].resize(mesh.cells.size(),vector<double>(9,0.0));
	}
	// for(int i=0; i<gradIterMax_LS; ++i)
	{
		int dummy0;
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
			-1, dummy0, surfNormalVec[0], gradSurfNor_LS[0]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
			-1, dummy0, surfNormalVec[1], gradSurfNor_LS[1]);
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "input", 
			-1, dummy0, surfNormalVec[2], gradSurfNor_LS[2]);
	}
	// for(int i=0; i<gradIterMax_GG; ++i){
		// math.calcGaussGreen(mesh, -1, surfNormalVec[0], gradSurfNor_LS[0]);
		// math.calcGaussGreen(mesh, -1, surfNormalVec[1], gradSurfNor_LS[1]);
		// math.calcGaussGreen(mesh, -1, surfNormalVec[2], gradSurfNor_LS[2]);
	// }
	kappa.clear();
	kappa.resize(mesh.cells.size(),0.0);
	for(int i=0; i<mesh.cells.size(); ++i){
		kappa[i] = -(gradSurfNor_LS[0][i][0] + gradSurfNor_LS[1][i][1] + gradSurfNor_LS[2][i][2]);
	}
	
	// vector<vector<double>> gradSurfNor_LS(mesh.cells.size(),vector<double>(3,0.0));
	// math.calcGaussGreen(mesh, -1, surfNormalVec[0], surfNormalVec[1], surfNormalVec[2], gradSurfNor_LS);
	// kappa.clear();
	// kappa.resize(mesh.cells.size(),0.0);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = -(gradSurfNor_LS[i][0] + gradSurfNor_LS[i][1] + gradSurfNor_LS[i][2]);
	// }
	
	
	// vector<double> gradSurfNorX(mesh.cells.size(),0.0);
	// vector<double> gradSurfNorY(mesh.cells.size(),0.0);
	// vector<double> gradSurfNorZ(mesh.cells.size(),0.0);
	// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double wCL = face.wC; double wCR = 1.0 - wCL;
		// double dPN = face.magPN;
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// vector<double> surfVarL(3,0.0);
			// vector<double> surfVarR(3,0.0);
			// surfVarL[0] = surfNormalVec[0][face.owner];
			// surfVarL[1] = surfNormalVec[1][face.owner];
			// surfVarL[2] = surfNormalVec[2][face.owner];
			// surfVarR[0] = surfNormalVec[0][face.neighbour];
			// surfVarR[1] = surfNormalVec[1][face.neighbour];
			// surfVarR[2] = surfNormalVec[2][face.neighbour];
			
			// // skewness
			// if(boolSkewnessCorrection){
				// for(int ii=0; ii<3; ++ii){
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][0]*face.skewness[0];
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][1]*face.skewness[1];
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][2]*face.skewness[2];
					
					// surfVarR[ii] += gradSurfNor_LS[ii][face.neighbour][0]*face.skewness[0];
					// surfVarR[ii] += gradSurfNor_LS[ii][face.neighbour][1]*face.skewness[1];
					// surfVarR[ii] += gradSurfNor_LS[ii][face.neighbour][2]*face.skewness[2];
				// }
			// }
			
			// gradSurfNorX[face.owner] += (surfVarR[0]-surfVarL[0])/dPN * face.vecPF[0] * face.area / mesh.cells[face.owner].volume;
			// gradSurfNorY[face.owner] += (surfVarR[1]-surfVarL[1])/dPN * face.vecPF[1] * face.area / mesh.cells[face.owner].volume;
			// gradSurfNorZ[face.owner] += (surfVarR[2]-surfVarL[2])/dPN * face.vecPF[2] * face.area / mesh.cells[face.owner].volume;
			
			// gradSurfNorX[face.neighbour] += (surfVarL[0]-surfVarR[0])/dPN * face.vecNF[0] * face.area / mesh.cells[face.neighbour].volume;
			// gradSurfNorY[face.neighbour] += (surfVarL[1]-surfVarR[1])/dPN * face.vecNF[1] * face.area / mesh.cells[face.neighbour].volume;
			// gradSurfNorZ[face.neighbour] += (surfVarL[2]-surfVarR[2])/dPN * face.vecNF[2] * face.area / mesh.cells[face.neighbour].volume;
				
			// // non-orthogonal
			// if(boolNonOrthoCorrection){
				// vector<double> vecGradSurNor(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// vecGradSurNor[ii] += wCL*gradSurfNor_LS[ii][face.owner][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);
					// vecGradSurNor[ii] += wCR*gradSurfNor_LS[ii][face.neighbour][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);	
				// }
				// gradSurfNorX[face.owner] += vecGradSurNor[0] * face.vecPF[0] * face.area / mesh.cells[face.owner].volume;
				// gradSurfNorY[face.owner] += vecGradSurNor[1] * face.vecPF[1] * face.area / mesh.cells[face.owner].volume;
				// gradSurfNorZ[face.owner] += vecGradSurNor[2] * face.vecPF[2] * face.area / mesh.cells[face.owner].volume;
				
				// gradSurfNorX[face.neighbour] -= vecGradSurNor[0] * face.vecNF[0] * face.area / mesh.cells[face.neighbour].volume;
				// gradSurfNorY[face.neighbour] -= vecGradSurNor[1] * face.vecNF[1] * face.area / mesh.cells[face.neighbour].volume;
				// gradSurfNorZ[face.neighbour] -= vecGradSurNor[2] * face.vecNF[2] * face.area / mesh.cells[face.neighbour].volume;
			// }
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// vector<double> surfVarL(3,0.0);
			// vector<double> surfVarR(3,0.0);
			// surfVarL[0] = surfNormalVec[0][face.owner];
			// surfVarL[1] = surfNormalVec[1][face.owner];
			// surfVarL[2] = surfNormalVec[2][face.owner];
			// surfVarR[0] = surfNormalVec_recv[0][ip];
			// surfVarR[1] = surfNormalVec_recv[1][ip];
			// surfVarR[2] = surfNormalVec_recv[2][ip];
			
			// // skewness
			// if(boolSkewnessCorrection){
				// for(int ii=0; ii<3; ++ii){
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][0]*face.skewness[0];
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][1]*face.skewness[1];
					// surfVarL[ii] += gradSurfNor_LS[ii][face.owner][2]*face.skewness[2];
					
					// surfVarR[ii] += gradSurfNor_LS_recv[ii][face.neighbour][0]*face.skewness[0];
					// surfVarR[ii] += gradSurfNor_LS_recv[ii][face.neighbour][1]*face.skewness[1];
					// surfVarR[ii] += gradSurfNor_LS_recv[ii][face.neighbour][2]*face.skewness[2];
				// }
			// }
			
			// gradSurfNorX[face.owner] += (surfVarR[0]-surfVarL[0])/dPN * 
				// face.vecPF[0] * face.area / mesh.cells[face.owner].volume;
			// gradSurfNorY[face.owner] += (surfVarR[1]-surfVarL[1])/dPN * 
				// face.vecPF[1] * face.area / mesh.cells[face.owner].volume;
			// gradSurfNorZ[face.owner] += (surfVarR[2]-surfVarL[2])/dPN * 
				// face.vecPF[2] * face.area / mesh.cells[face.owner].volume;
			
			// // non-orthogonal
			// if(boolNonOrthoCorrection){
				// vector<double> vecGradSurNor(3,0.0);
				// for(int ii=0; ii<3; ++ii){
					// vecGradSurNor[ii] += wCL*gradSurfNor_LS[ii][face.owner][ii]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);
					// vecGradSurNor[ii] += wCR*gradSurfNor_LS_recv[ii][face.neighbour]*(face.unitNormals[ii] - face.unitNomalsPN[ii]);	
				// }
				// gradSurfNorX[face.owner] += vecGradSurNor[0] * face.vecPF[0] * face.area / mesh.cells[face.owner].volume;
				// gradSurfNorY[face.owner] += vecGradSurNor[1] * face.vecPF[1] * face.area / mesh.cells[face.owner].volume;
				// gradSurfNorZ[face.owner] += vecGradSurNor[2] * face.vecPF[2] * face.area / mesh.cells[face.owner].volume;
			// }
			
			
			// ++ip;
		// }
	// }
	
	
	// //================================================
	// // // calc Curvature (green-gause)
	// kappa.clear();
	// kappa.resize(mesh.cells.size(),0.0);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = -(gradSurfNor_LS[i][0] + gradSurfNor_LS[i][1] + gradSurfNor_LS[i][2]);
	// }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = -(gradSurfNorX[i] + gradSurfNorY[i] + gradSurfNorZ[i]);
	// }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = -(gradSurfNor_LS[0][i][0] + gradSurfNor_LS[1][i][1] + gradSurfNor_LS[2][i][2]);
	// }
	// gradAi.clear();
	// gradAi.resize(mesh.cells.size(),vector<double>(9,0.0));
	// math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// // kappa[i] = -100.0*(surfNormalVec[0][i] + surfNormalVec[1][i] + surfNormalVec[2][i]);
		// kappa[i] = -100.0*(gradAi[i][0] + gradAi[i][1] + gradAi[i][2]);
		// // kappa[i] = -1000.0;
	// }

	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = 1.0*(gradSurfNor_LS[0][i][0]*gradSurfNor_LS[0][i][0] + gradSurfNor_LS[1][i][1]*gradSurfNor_LS[1][i][1] + gradSurfNor_LS[2][i][2]*gradSurfNor_LS[2][i][2]);
	// }
	
	// //================================================
	// // calc Curvature (green-gause)
	// kappa.clear();
	// kappa.resize(mesh.cells.size(),0.0);
	// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double wCL = face.wC; double wCR = 1.0 - wCL;
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double var = 0.0;
			// var += (wCL*surfNormalVec[0][face.owner] + wCR*surfNormalVec[0][face.neighbour])*face.unitNormals[0];
			// var += (wCL*surfNormalVec[1][face.owner] + wCR*surfNormalVec[1][face.neighbour])*face.unitNormals[1];
			// var += (wCL*surfNormalVec[2][face.owner] + wCR*surfNormalVec[2][face.neighbour])*face.unitNormals[2];
			
			// kappa[face.owner] -= var * face.area / mesh.cells[face.owner].volume;
			// kappa[face.neighbour] += var * face.area / mesh.cells[face.neighbour].volume;
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// double var = 0.0;
			// var += (wCL*surfNormalVec[0][face.owner] + wCR*surfNormalVec_recv[0][ip])*face.unitNormals[0];
			// var += (wCL*surfNormalVec[1][face.owner] + wCR*surfNormalVec_recv[1][ip])*face.unitNormals[1];
			// var += (wCL*surfNormalVec[2][face.owner] + wCR*surfNormalVec_recv[2][ip])*face.unitNormals[2];
			
			// kappa[face.owner] -= var * face.area / mesh.cells[face.owner].volume;
			
			// ++ip;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// double var = 0.0;
			// var += (surfNormalVec[0][face.owner])*face.unitNormals[0];
			// var += (surfNormalVec[1][face.owner])*face.unitNormals[1];
			// var += (surfNormalVec[2][face.owner])*face.unitNormals[2];
			
			// kappa[face.owner] -= var * face.area / mesh.cells[face.owner].volume;
			
		// }
	// }
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// //================================================
	// // smoothing process of Curvature (denner)
	
	// // first smoothing
	// for(int iter=0; iter<iterFirstSmoothingMax; ++iter)
	// {
		// vector<double> smoothUp(mesh.cells.size(),0.0);
		// vector<double> smoothDown(mesh.cells.size(),0.0);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double weight = pow(1.0-2.0*abs(0.5-smoothAi[i]), 8.0);
			// smoothUp[i] = kappa[i]*weight;
			// smoothDown[i] = weight;
		// }
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC; double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// double weightL = pow(1.0-2.0*abs(0.5-smoothAi[face.owner]), 8.0);
				// double weightR = pow(1.0-2.0*abs(0.5-smoothAi[face.neighbour]), 8.0);
			
				// smoothUp[face.owner] += weightR*kappa[face.neighbour];
				// smoothDown[face.owner] += weightR;
				
				// smoothUp[face.neighbour] += weightL*kappa[face.owner];
				// smoothDown[face.neighbour] += weightL;
				
				// // double avgVar = wCL*smoothAi[face.owner] + wCR*smoothAi[face.neighbour];
				// // double weight = pow(1.0-2.0*avgVar, 8.0) * face.area;
				// // double avgKappa = wCL*kappa[face.owner] + wCR*kappa[face.neighbour];
			
				// // smoothUp[face.owner] += weight*avgKappa;
				// // smoothDown[face.owner] += weight;
				
				// // smoothUp[face.neighbour] += weight*avgKappa;
				// // smoothDown[face.neighbour] += weight;
			
			// }
			
		// }
		
		// if(size>1){
			// // processor faces
			// vector<double> sendValues;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// sendValues.push_back(kappa[face.owner]);
				// }
			// }
			// vector<double> recvValues;
			// mpi.setProcsFaceDatas(
						// sendValues, recvValues,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[ip]), 8.0);
					
					// smoothUp[face.owner] += weightR*recvValues[ip];
					// smoothDown[face.owner] += weightR;
					
					// // double avgVar = wCL*smoothAi[face.owner] + wCR*smoothAi_recv[ip];
					// // double weight = pow(1.0-2.0*avgVar, 8.0) * face.area;
					// // double avgKappa = wCL*kappa[face.owner] + wCR*recvValues[ip];
					
					// // smoothUp[face.owner] += weight*avgKappa;
					// // smoothDown[face.owner] += weight;
					
					// ++ip;
				// }
			// }
			
		// }
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
				// kappa[i] = smoothUp[i]/smoothDown[i];
			// }
			// // else{
				// // kappa[i] = 0.0;
			// // }
		// }
		
	// }
	
	
	
	
	// // second smoothing
	// for(int iter=0; iter<iterSecondSmoothingMax; ++iter)
	// {
		// vector<double> smoothUp(mesh.cells.size(),0.0);
		// vector<double> smoothDown(mesh.cells.size(),0.0);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double weight = pow(1.0-2.0*abs(0.5-smoothAi[i]), 8.0);
			// smoothUp[i] = kappa[i]*weight;
			// smoothDown[i] = weight;
		// }
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC; double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// double weightL = pow(1.0-2.0*abs(0.5-smoothAi[face.owner]), 8.0);
				// double weightR = pow(1.0-2.0*abs(0.5-smoothAi[face.neighbour]), 8.0);
				// double weightL_m = 0.0;
				// weightL_m += surfNormalVec[0][face.neighbour]*face.unitNomalsPN[0];
				// weightL_m += surfNormalVec[1][face.neighbour]*face.unitNomalsPN[1];
				// weightL_m += surfNormalVec[2][face.neighbour]*face.unitNomalsPN[2];
				// weightL_m = abs(weightL_m);
				// weightL_m = pow(weightL_m, 8.0);
				// double weightR_m = 0.0;
				// weightR_m += surfNormalVec[0][face.owner]*face.unitNomalsPN[0];
				// weightR_m += surfNormalVec[1][face.owner]*face.unitNomalsPN[1];
				// weightR_m += surfNormalVec[2][face.owner]*face.unitNomalsPN[2];
				// weightR_m = abs(weightR_m);
				// weightR_m = pow(weightR_m, 8.0);
			
				// smoothUp[face.owner] += weightR*weightR_m*kappa[face.neighbour];
				// smoothUp[face.neighbour] += weightL*weightL_m*kappa[face.owner];
				
				// smoothDown[face.owner] += weightR*weightR_m;
				// smoothDown[face.neighbour] += weightL*weightL_m;
			
			// }
			
		// }
		
		// if(size>1){
			// // processor faces
			// vector<double> sendValues;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// sendValues.push_back(kappa[face.owner]);
				// }
			// }
			// vector<double> recvValues;
			// mpi.setProcsFaceDatas(
						// sendValues, recvValues,
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
			// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[ip]), 8.0);
					// double weightR_m = 0.0;
					// // weightR_m += surfNormalVec_recv[0][ip]*face.unitNomalsPN[0];
					// // weightR_m += surfNormalVec_recv[1][ip]*face.unitNomalsPN[1];
					// // weightR_m += surfNormalVec_recv[2][ip]*face.unitNomalsPN[2];
					// weightR_m += surfNormalVec[0][face.owner]*face.unitNomalsPN[0];
					// weightR_m += surfNormalVec[1][face.owner]*face.unitNomalsPN[1];
					// weightR_m += surfNormalVec[2][face.owner]*face.unitNomalsPN[2];
					// weightR_m = abs(weightR_m);
					// weightR_m = pow(weightR_m, 8.0);
					
					// smoothUp[face.owner] += weightR*weightR_m*recvValues[ip];
					// smoothDown[face.owner] += weightR*weightR_m;
					
					// ++ip;
				// }
			// }
			
		// }
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
				// kappa[i] = smoothUp[i]/smoothDown[i];
			// }
			// // else{
				// // kappa[i] = 0.0;
			// // }
		// }
		
	// }
	
	
	
	
	
	
	
	
	
	


// for(int i=0; i<mesh.cells.size(); ++i){
	// auto& cell = mesh.cells[i];
	// cell.var[1] = (gradAi[i][0]);
	// cell.var[2] = (gradAi[i][1]);
	// cell.var[3] = (gradAi[i][2]);
	// // cell.var[1] = (surfNormalVec[0][i]);
	// // cell.var[2] = (surfNormalVec[1][i]);
	// // cell.var[3] = (surfNormalVec[2][i]);
	// // cell.var[1] = kappa[i]*(gradAi[i][0]);
	// // cell.var[2] = kappa[i]*(gradAi[i][1]);
	// // cell.var[3] = kappa[i]*(gradAi[i][2]);
	// cell.var[4] = (kappa[i]);
	// cell.var[0] = (smoothAi[i]);
	// // cell.var[5] = kappa_test[i];
	

	// // if(interfacePresentCell[i]){
		// // cell.var[5] = 1.0;
	// // }
	// // else{
		// // cell.var[5] = 0.0;
	// // }
// }

	
	
	
	
	
	

	
}

