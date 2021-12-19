// #include "build.h"
// #include <cmath>
// #include <array>
// #include <numeric>


// void SEMO_Solvers_Builder::calcCurvature(
	// SEMO_Mesh_Builder& mesh,
	// int cn,
	// vector<double>& kappa){
		
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	// SEMO_Utility_Math math;
	
	
	// int iterVFSmoothingMax = 4;
	// int iterFirstSmoothingMax = 1;
	// int iterSecondSmoothingMax = 1;
	
	
	
	
	// vector<double> smoothAi(mesh.cells.size());
	// for(int i=0; i<mesh.cells.size(); ++i){
		// smoothAi[i] = mesh.cells[i].var[cn];
	// }
	
	// //================================================
	// // smoothing Ai
	// for(int iter=0; iter<iterVFSmoothingMax; ++iter){

		// vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
		// // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// // for(int i=0; i<3; ++i){
			// math.calcLeastSquare2nd(mesh, cn, smoothAi, gradAi);
		// // }
		
		// vector<double> AiUp(mesh.cells.size(),0.0);
		// vector<double> AiDown(mesh.cells.size(),0.0);
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC; double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
				// double AiF = wCL*smoothAi[face.owner] +
							 // wCR*smoothAi[face.neighbour];
					
				// // skewness correction
				// for(int ii=0; ii<3; ++ii){
					// AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
					// AiF += wCR * gradAi[face.neighbour][ii]*face.vecSkewness[ii];
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
			// int num=0;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double AiF = wCL*smoothAi[face.owner]+wCR*recvValues[num];

					// // skewness correction
					// for(int ii=0; ii<3; ++ii){
						// AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
						// AiF += wCR * recvGrad[ii][num]*face.vecSkewness[ii];
					// }
					
					// AiUp[face.owner] += AiF*face.area;
				
					// AiDown[face.owner] += face.area;
					
					// ++num;
				// }
			// }
		// }
	
		
		// for(int i=0; i<mesh.cells.size(); ++i){
			// smoothAi[i] = AiUp[i]/AiDown[i];
		// }
		
	// }
	
	
	
	
	
	
	// vector<double> smoothAi_recv;
	// if(size>1){
		// // processor faces
		// vector<double> sendValues;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// sendValues.push_back(smoothAi[face.owner]);
			// }
		// }
		// mpi.setProcsFaceDatas(
					// sendValues, smoothAi_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
	// }
	
	
	// //================================================
	// // calc gauss-green gradient
	// vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
		// // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// math.calcLeastSquare2nd(mesh, cn, smoothAi, gradAi);
	// // for(int i=0; i<3; ++i){
		// // // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// // math.calcLeastSquare2nd(mesh, cn, smoothAi, gradAi);
	// // }

	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // double wCL = face.wC;
		// // // double wCL = 0.5;
		// // double wCR = 1.0 - wCL;
			
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double dPN = face.magPN;
			
			// // for(int ii=0; ii<3; ++ii){
				// // gradAi[face.owner][ii] += (smoothAi[face.neighbour]-smoothAi[face.owner])/dPN * 
					// // face.vecPF[ii] * face.area / mesh.cells[face.owner].volume;
				// // gradAi[face.neighbour][ii] += (smoothAi[face.owner]-smoothAi[face.neighbour])/dPN * 
					// // face.vecNF[ii] * face.area / mesh.cells[face.neighbour].volume;
			// // }
			
		// // }
	// // }

	// //================================================
	// // surfNormalVec = dAidx/|dAidx|
	// vector<vector<double>> surfNormalVec(3,vector<double>(mesh.cells.size()));
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// double magGrad = sqrt(
			// pow(gradAi[i][0],2.0)+pow(gradAi[i][1],2.0)+pow(gradAi[i][2],2.0));
			
		// if( magGrad > std::numeric_limits<double>::min() ){
			// for(int ii=0; ii<3; ++ii){
				// surfNormalVec[ii][i] = gradAi[i][ii]/magGrad;
			// }
		// }
		// else{
			// for(int ii=0; ii<3; ++ii) surfNormalVec[ii][i] = 0.0;
		// }
		
	// }
	// vector<vector<double>> surfNormalVec_recv(3,vector<double>());
	// if(size>1){
		// // processor faces
		// vector<vector<double>> sendValues(3,vector<double>());
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// sendValues[0].push_back(surfNormalVec[0][face.owner]);
				// sendValues[1].push_back(surfNormalVec[1][face.owner]);
				// sendValues[2].push_back(surfNormalVec[2][face.owner]);
			// }
		// }
		// for(int i=0; i<3; ++i){
			// mpi.setProcsFaceDatas(
						// sendValues[i], surfNormalVec_recv[i],
						// mesh.countsProcFaces, mesh.countsProcFaces, 
						// mesh.displsProcFaces, mesh.displsProcFaces);
		// }
	// }
	
	
	
	
	
	// //================================================
	// // calc surface normal vector gradients
	// vector<vector<double>> gradSurfNorX(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradSurfNorY(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradSurfNorZ(mesh.cells.size(),vector<double>(3,0.0));
	// // math.calcGaussGreen(mesh, cn, surfNormalVec[0], gradSurfNorX);
	// // math.calcGaussGreen(mesh, cn, surfNormalVec[1], gradSurfNorY);
	// // math.calcGaussGreen(mesh, cn, surfNormalVec[2], gradSurfNorZ);
	// math.calcLeastSquare2nd(mesh, cn, surfNormalVec[0], gradSurfNorX);
	// math.calcLeastSquare2nd(mesh, cn, surfNormalVec[1], gradSurfNorY);
	// math.calcLeastSquare2nd(mesh, cn, surfNormalVec[2], gradSurfNorZ);

	// // vector<double> gradSurfNorX(mesh.cells.size(),0.0);
	// // vector<double> gradSurfNorY(mesh.cells.size(),0.0);
	// // vector<double> gradSurfNorZ(mesh.cells.size(),0.0);
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // double wCL = face.wC;
		// // // double wCL = 0.5;
		// // double wCR = 1.0 - wCL;
			
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double dPN = face.magPN;
			
			// // gradSurfNorX[face.owner] += (surfNormalVec[0][face.neighbour]-surfNormalVec[0][face.owner])/dPN * 
				// // face.vecPF[0] * face.area / mesh.cells[face.owner].volume;
			// // gradSurfNorY[face.owner] += (surfNormalVec[1][face.neighbour]-surfNormalVec[1][face.owner])/dPN * 
				// // face.vecPF[1] * face.area / mesh.cells[face.owner].volume;
			// // gradSurfNorZ[face.owner] += (surfNormalVec[2][face.neighbour]-surfNormalVec[2][face.owner])/dPN * 
				// // face.vecPF[2] * face.area / mesh.cells[face.owner].volume;
				
			// // gradSurfNorX[face.neighbour] += (surfNormalVec[0][face.owner]-surfNormalVec[0][face.neighbour])/dPN * 
				// // face.vecNF[0] * face.area / mesh.cells[face.neighbour].volume;
			// // gradSurfNorY[face.neighbour] += (surfNormalVec[1][face.owner]-surfNormalVec[1][face.neighbour])/dPN * 
				// // face.vecNF[1] * face.area / mesh.cells[face.neighbour].volume;
			// // gradSurfNorZ[face.neighbour] += (surfNormalVec[2][face.owner]-surfNormalVec[2][face.neighbour])/dPN * 
				// // face.vecNF[2] * face.area / mesh.cells[face.neighbour].volume;
			
			
		// // }
	// // }
	
	
	// // //================================================
	// // // calc Curvature (green-gause)
	// kappa.clear();
	// kappa.resize(mesh.cells.size(),0.0);
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // kappa[i] = (gradSurfNorX[i] + gradSurfNorY[i] + gradSurfNorZ[i]);
	// // }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// kappa[i] = -(gradSurfNorX[i][0] + gradSurfNorY[i][1] + gradSurfNorZ[i][2]);
	// }

	// // kappa.clear();
	// // kappa.resize(mesh.cells.size(),0.0);
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // double wCL = face.wC;
		// // // double wCL = 0.5;
		// // double wCR = 1.0 - wCL;
			
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double varF = 0.0;
			// // for(int ii=0; ii<3; ++ii){
				// // varF += ( wCL*surfNormalVec[ii][face.owner]+wCR*surfNormalVec[ii][face.neighbour] ) * face.unitNormals[ii];
			// // }
			// // kappa[face.owner] += varF*face.area/mesh.cells[face.owner].volume;
			// // kappa[face.neighbour] -= varF*face.area/mesh.cells[face.neighbour].volume;
		// // }
	// // }
	// // // boundary
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // for(int i=str; i<end; ++i){
				// // auto& face = mesh.faces[i];
				
				// // double varF = 0.0;
				// // varF += surfNormalVec[0][face.owner]*face.unitNormals[0];
				// // varF += surfNormalVec[1][face.owner]*face.unitNormals[1];
				// // varF += surfNormalVec[2][face.owner]*face.unitNormals[2];
				
				// // kappa[face.owner] += varF*face.area/mesh.cells[face.owner].volume;
				
			// // }
		// // }
	// // }
	
	
	
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
				// smoothUp[face.neighbour] += weightL*kappa[face.owner];
				
				// smoothDown[face.owner] += weightR;
				// smoothDown[face.neighbour] += weightL;
			
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
			// int num=0;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[num]), 8.0);
					
					// smoothUp[face.owner] += weightR*recvValues[num];
				
					// smoothDown[face.owner] += weightR;
					
					// ++num;
				// }
			// }
			
		// }
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
				// kappa[i] = smoothUp[i]/smoothDown[i];
			// }
			// else{
				// kappa[i] = 0.0;
			// }
		// }
		
	// }
	
	
	
	
	// // second smoothing
	// for(int iter=0; iter<iterSecondSmoothingMax; ++iter)
	// {
		// vector<double> smoothUp(mesh.cells.size(),0.0);
		// vector<double> smoothDown(mesh.cells.size(),0.0);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double weight = pow(1.0-2.0*abs(0.5-smoothAi[i]), 8.0);
			// // double weight = 0.0;
			// // weight += surfNormalVec_recv[0][i]
			// // weight = pow(weight, 8.0);
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
			// int num=0;
			// for(int i=0; i<mesh.faces.size(); ++i){
				// auto& face = mesh.faces[i];
			
				// double wCL = face.wC; double wCR = 1.0 - wCL;
				
				// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[num]), 8.0);
					// double weightR_m = 0.0;
					// weightR_m += surfNormalVec_recv[0][num]*face.unitNomalsPN[0];
					// weightR_m += surfNormalVec_recv[1][num]*face.unitNomalsPN[1];
					// weightR_m += surfNormalVec_recv[2][num]*face.unitNomalsPN[2];
					// weightR_m = abs(weightR_m);
					// weightR_m = pow(weightR_m, 8.0);
					
					// smoothUp[face.owner] += weightR*weightR_m*recvValues[num];
				
					// smoothDown[face.owner] += weightR*weightR_m;
					
					// ++num;
				// }
			// }
			
		// }
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// if( abs(smoothDown[i]) > std::numeric_limits<double>::min() ){
				// kappa[i] = smoothUp[i]/smoothDown[i];
			// }
			// else{
				// kappa[i] = 0.0;
			// }
		// }
		
	// }
	
	
	

	
// }


// // void SEMO_Solvers_Builder::calcCurvature(
	// // SEMO_Mesh_Builder& mesh,
	// // int cn,
	// // vector<double>& kappa){
		
    // // int rank = MPI::COMM_WORLD.Get_rank(); 
    // // int size = MPI::COMM_WORLD.Get_size();
	
	// // SEMO_MPI_Builder mpi;
	
	// // SEMO_Utility_Math math;
	
	
	// // int iterVFSmoothingMax = 0;
	// // int iterFirstSmoothingMax = 0;
	// // int iterSecondSmoothingMax = 0;
	
	
	
	// // vector<double> AiUp(mesh.cells.size());
	// // vector<double> AiDown(mesh.cells.size());
	
	// // vector<double> smoothAi(mesh.cells.size());
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // smoothAi[i] = mesh.cells[i].var[cn];
	// // }
	
	// // //================================================
	// // // smoothing Ai
	// // for(int iter=0; iter<iterVFSmoothingMax; ++iter){

		// // vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
		// // // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// // for(int i=0; i<3; ++i){
			// // math.calcLeastSquare2nd(mesh, cn, smoothAi, gradAi);
		// // }
		
		// // std::fill(AiUp.begin(), AiUp.end(), 0.0);
		// // std::fill(AiDown.begin(), AiDown.end(), 0.0);
		// // // for(int i=0; i<mesh.cells.size(); ++i){
			// // // AiUp[i] = 0.0;
			// // // AiDown[i] = 0.0;
		// // // }
		
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // double wCL = face.wC;
			// // // double wCL = 0.5;
			// // double wCR = 1.0 - wCL;
			
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
				// // double AiF = 
					// // wCL*smoothAi[face.owner] +
					// // wCR*smoothAi[face.neighbour];
					
				// // // skewness correction
				// // for(int ii=0; ii<3; ++ii){
					// // AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
					// // AiF += wCR * gradAi[face.neighbour][ii]*face.vecSkewness[ii];
				// // }
				
				// // AiUp[face.owner] += AiF*face.area;
				// // AiUp[face.neighbour] += AiF*face.area;
				
				// // AiDown[face.owner] += face.area;
				// // AiDown[face.neighbour] += face.area;
			
			// // }
			
		// // }

		// // // boundary
		// // for(auto& boundary : mesh.boundary){
			
			// // if(boundary.neighbProcNo == -1){
				
				// // int str = boundary.startFace;
				// // int end = str + boundary.nFaces;
				
				// // for(int i=str; i<end; ++i){
					// // auto& face = mesh.faces[i];
					
					// // double AiF = smoothAi[face.owner];
					
					// // AiUp[face.owner] += AiF*face.area;
				
					// // AiDown[face.owner] += face.area;
					
				// // }
			// // }
		// // }
		
		// // if(size>1){
			// // // processor faces
			// // vector<double> sendValues;
			// // vector<vector<double>> sendGrad(3,vector<double>());
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // sendValues.push_back(smoothAi[face.owner]);
					// // for(int ii=0; ii<3; ++ii){
						// // sendGrad[ii].push_back(gradAi[face.owner][ii]);
					// // }
				// // }
			// // }
			// // vector<double> recvValues;
			// // vector<vector<double>> recvGrad(3,vector<double>());
			// // mpi.setProcsFaceDatas(
						// // sendValues, recvValues,
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
			// // for(int i=0; i<3; ++i){
				// // mpi.setProcsFaceDatas(
							// // sendGrad[i], recvGrad[i],
							// // mesh.countsProcFaces, mesh.countsProcFaces, 
							// // mesh.displsProcFaces, mesh.displsProcFaces);
			// // }
			// // int num=0;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
			
				// // double wCL = face.wC;
				// // // double wCL = 0.5;
				// // double wCR = 1.0 - wCL;
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// // double AiF = wCL*smoothAi[face.owner]+wCR*recvValues[num];

					// // // skewness correction
					// // for(int ii=0; ii<3; ++ii){
						// // AiF += wCL * gradAi[face.owner][ii]*face.vecSkewness[ii];
						// // AiF += wCR * recvGrad[ii][num]*face.vecSkewness[ii];
					// // }
					
					// // AiUp[face.owner] += AiF*face.area;
				
					// // AiDown[face.owner] += face.area;
					
					// // ++num;
				// // }
			// // }
		// // }
	
		
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // smoothAi[i] = AiUp[i]/AiDown[i];
		// // }
		
	// // }
	
	
	
	
	
	
	// // vector<double> smoothAi_recv;
	// // if(size>1){
		// // // processor faces
		// // vector<double> sendValues;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // sendValues.push_back(smoothAi[face.owner]);
			// // }
		// // }
		// // mpi.setProcsFaceDatas(
					// // sendValues, smoothAi_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
	// // }
	
	
	// // //================================================
	// // // calc gauss-green gradient
	// // vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
	// // for(int i=0; i<3; ++i){
		// // // math.calcGaussGreen(mesh, cn, smoothAi, gradAi);
		// // math.calcLeastSquare2nd(mesh, cn, smoothAi, gradAi);
	// // }
	

	// // //================================================
	// // // surfNormalVec = dAidx/|dAidx|
	// // vector<vector<double>> surfNormalVec(3,vector<double>(mesh.cells.size()));
	
	// // for(int i=0; i<mesh.cells.size(); ++i){
		
		// // double magGrad = sqrt(
			// // pow(gradAi[i][0],2.0)+pow(gradAi[i][1],2.0)+pow(gradAi[i][2],2.0));
			
		// // for(int ii=0; ii<3; ++ii){
			// // surfNormalVec[ii][i] = gradAi[i][ii]/(magGrad+1.e-200);
		// // }
		
	// // }
	// // vector<vector<double>> surfNormalVec_recv(3,vector<double>());
	// // if(size>1){
		// // // processor faces
		// // vector<vector<double>> sendValues(3,vector<double>());
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // sendValues[0].push_back(surfNormalVec[0][face.owner]);
				// // sendValues[1].push_back(surfNormalVec[1][face.owner]);
				// // sendValues[2].push_back(surfNormalVec[2][face.owner]);
			// // }
		// // }
		// // for(int i=0; i<3; ++i){
			// // mpi.setProcsFaceDatas(
						// // sendValues[i], surfNormalVec_recv[i],
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // }
	// // }
	
	
	
	
	
	
	
	
	// // //================================================
	// // // calc surface normal vector gradients
	// // vector<vector<double>> gradSurfNorX(mesh.cells.size(),vector<double>(3,0.0));
	// // vector<vector<double>> gradSurfNorY(mesh.cells.size(),vector<double>(3,0.0));
	// // vector<vector<double>> gradSurfNorZ(mesh.cells.size(),vector<double>(3,0.0));
	// // for(int i=0; i<3; ++i){
		// // math.calcLeastSquare2nd(mesh, -1, surfNormalVec[0], gradSurfNorX);
		// // math.calcLeastSquare2nd(mesh, -1, surfNormalVec[1], gradSurfNorY);
		// // math.calcLeastSquare2nd(mesh, -1, surfNormalVec[2], gradSurfNorZ);
	// // }
	
	
	// // //================================================
	// // // calc Curvature (green-gause)
	// // kappa.clear();
	// // kappa.resize(mesh.cells.size(),0.0);
	
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // double wCL = face.wC;
		// // // double wCL = 0.5;
		// // double wCR = 1.0 - wCL;
			
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// // double varF = 0.0;
			// // for(int ii=0; ii<3; ++ii){
				// // varF += ( wCL*surfNormalVec[ii][face.owner]+wCR*surfNormalVec[ii][face.neighbour] ) * face.unitNormals[ii];
			// // }
			
			// // vector<double> gradXL(3,0.0); vector<double> gradXR(3,0.0);
			// // vector<double> gradYL(3,0.0); vector<double> gradYR(3,0.0);
			// // vector<double> gradZL(3,0.0); vector<double> gradZR(3,0.0);
			
			// // for(int ii=0; ii<3; ++ii){
				// // gradXL[ii] = gradSurfNorX[face.owner][ii]; gradXR[ii] = gradSurfNorX[face.neighbour][ii];
				// // gradYL[ii] = gradSurfNorY[face.owner][ii]; gradYR[ii] = gradSurfNorY[face.neighbour][ii];
				// // gradZL[ii] = gradSurfNorZ[face.owner][ii]; gradZR[ii] = gradSurfNorZ[face.neighbour][ii];
			// // }
			
			// // // skewness correction
			// // for(int ii=0; ii<3; ++ii){
				// // varF += wCL * gradXL[ii]*face.vecSkewness[ii] * face.unitNormals[0];
				// // varF += wCR * gradXR[ii]*face.vecSkewness[ii] * face.unitNormals[0];
				
				// // varF += wCL * gradYL[ii]*face.vecSkewness[ii] * face.unitNormals[1];
				// // varF += wCR * gradYR[ii]*face.vecSkewness[ii] * face.unitNormals[1];
				
				// // varF += wCL * gradZL[ii]*face.vecSkewness[ii] * face.unitNormals[2];
				// // varF += wCR * gradZR[ii]*face.vecSkewness[ii] * face.unitNormals[2];
			// // }
			
			// // kappa[face.owner] += varF*face.area/mesh.cells[face.owner].volume;
			// // kappa[face.neighbour] -= varF*face.area/mesh.cells[face.neighbour].volume;
		// // }
	// // }
	// // // boundary
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // for(int i=str; i<end; ++i){
				// // auto& face = mesh.faces[i];
				
				// // double varF = 0.0;
				// // varF += surfNormalVec[0][face.owner]*face.unitNormals[0];
				// // varF += surfNormalVec[1][face.owner]*face.unitNormals[1];
				// // varF += surfNormalVec[2][face.owner]*face.unitNormals[2];
				// // for(int j=0; j<3; ++j){
					// // kappa[face.owner] += varF*face.area/mesh.cells[face.owner].volume;
				// // }
				
			// // }
		// // }
	// // }
	// // if(size>1){
		// // // processor faces
		// // vector<vector<double>> sendGradSkewValues(3,vector<double>());
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // double varX = 0.0;
				// // double varY = 0.0;
				// // double varZ = 0.0;
				// // for(int ii=0; ii<3; ++ii){
					// // varX += gradSurfNorX[face.owner][ii]*face.vecSkewness[ii];
					// // varY += gradSurfNorY[face.owner][ii]*face.vecSkewness[ii];
					// // varZ += gradSurfNorZ[face.owner][ii]*face.vecSkewness[ii];
				// // }
				// // sendGradSkewValues[0].push_back(varX);
				// // sendGradSkewValues[1].push_back(varY);
				// // sendGradSkewValues[2].push_back(varZ);
			// // }
		// // }
		// // vector<vector<double>> recvGradSkewValues(3,vector<double>());
		// // for(int i=0; i<3; ++i){
			// // mpi.setProcsFaceDatas(
						// // sendGradSkewValues[i], recvGradSkewValues[i],
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
		// // }
		// // int num=0;
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // double wCL = face.wC;
			// // // double wCL = 0.5;
			// // double wCR = 1.0 - wCL;
			
			// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
				// // double varF = 0.0;
				// // varF += ( wCL*surfNormalVec[0][face.owner]+wCR*surfNormalVec_recv[0][num] ) * face.unitNormals[0];
				// // varF += ( wCL*surfNormalVec[1][face.owner]+wCR*surfNormalVec_recv[1][num] ) * face.unitNormals[1];
				// // varF += ( wCL*surfNormalVec[2][face.owner]+wCR*surfNormalVec_recv[2][num] ) * face.unitNormals[2];
					
				// // // skewness correction
				// // {
					// // for(int ii=0; ii<3; ++ii){
						// // varF += wCL * gradSurfNorX[face.owner][ii]*face.vecSkewness[ii] * face.unitNormals[0];
						// // varF += wCL * gradSurfNorY[face.owner][ii]*face.vecSkewness[ii] * face.unitNormals[1];
						// // varF += wCL * gradSurfNorZ[face.owner][ii]*face.vecSkewness[ii] * face.unitNormals[2];
						
						// // varF += wCR * recvGradSkewValues[ii][num] * face.unitNormals[ii];
					// // }
				// // }
				
				// // kappa[face.owner] += varF*face.area/mesh.cells[face.owner].volume;
				
				// // ++num;
			// // }
		// // }
		
	// // }
	
	
	
	
	
	
	
	
	
	// // //================================================
	// // // smoothing process of Curvature (denner)
	
	// // // first smoothing
	// // for(int iter=0; iter<iterFirstSmoothingMax; ++iter)
	// // {
		// // vector<double> smoothUp(mesh.cells.size(),0.0);
		// // vector<double> smoothDown(mesh.cells.size(),0.0);
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // double weight = pow(1.0-2.0*abs(0.5-smoothAi[i]), 8.0);
			// // smoothUp[i] = kappa[i]*weight;
			// // smoothDown[i] = weight;
		// // }
		
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // double wCL = face.wC;
			// // // double wCL = 0.5;
			// // double wCR = 1.0 - wCL;
			
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// // double weightL = pow(1.0-2.0*abs(0.5-smoothAi[face.owner]), 8.0);
				// // double weightR = pow(1.0-2.0*abs(0.5-smoothAi[face.neighbour]), 8.0);
			
				// // smoothUp[face.owner] += weightR*kappa[face.neighbour];
				// // smoothUp[face.neighbour] += weightL*kappa[face.owner];
				
				// // smoothDown[face.owner] += weightR;
				// // smoothDown[face.neighbour] += weightL;
			
			// // }
			
		// // }
		
		// // if(size>1){
			// // // processor faces
			// // vector<double> sendValues;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // sendValues.push_back(kappa[face.owner]);
				// // }
			// // }
			// // vector<double> recvValues;
			// // mpi.setProcsFaceDatas(
						// // sendValues, recvValues,
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
			// // int num=0;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
			
				// // double wCL = face.wC;
				// // // double wCL = 0.5;
				// // double wCR = 1.0 - wCL;
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// // double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[num]), 8.0);
					
					// // smoothUp[face.owner] += weightR*recvValues[num];
				
					// // smoothDown[face.owner] += weightR;
					
					// // ++num;
				// // }
			// // }
			
		// // }
	
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // kappa[i] = smoothUp[i]/( smoothDown[i] + 1.e-200 );
		// // }
		
	// // }
	
	
	
	
	// // // second smoothing
	// // for(int iter=0; iter<iterSecondSmoothingMax; ++iter)
	// // {
		// // vector<double> smoothUp(mesh.cells.size(),0.0);
		// // vector<double> smoothDown(mesh.cells.size(),0.0);
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // double weight = pow(1.0-2.0*abs(0.5-smoothAi[i]), 8.0);
			// // // double weight = 0.0;
			// // // weight += surfNormalVec_recv[0][i]
			// // // weight = pow(weight, 8.0);
			// // smoothUp[i] = kappa[i]*weight;
			// // smoothDown[i] = weight;
		// // }
		
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // auto& face = mesh.faces[i];
			
			// // double wCL = face.wC;
			// // // double wCL = 0.5;
			// // double wCR = 1.0 - wCL;
			
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// // double weightL = pow(1.0-2.0*abs(0.5-smoothAi[face.owner]), 8.0);
				// // double weightR = pow(1.0-2.0*abs(0.5-smoothAi[face.neighbour]), 8.0);
				// // double weightL_m = 0.0;
				// // weightL_m += surfNormalVec[0][face.neighbour]*face.unitNomalsPN[0];
				// // weightL_m += surfNormalVec[1][face.neighbour]*face.unitNomalsPN[1];
				// // weightL_m += surfNormalVec[2][face.neighbour]*face.unitNomalsPN[2];
				// // weightL_m = abs(weightL_m);
				// // weightL_m = pow(weightL_m, 8.0);
				// // double weightR_m = 0.0;
				// // weightR_m += surfNormalVec[0][face.owner]*face.unitNomalsPN[0];
				// // weightR_m += surfNormalVec[1][face.owner]*face.unitNomalsPN[1];
				// // weightR_m += surfNormalVec[2][face.owner]*face.unitNomalsPN[2];
				// // weightR_m = abs(weightR_m);
				// // weightR_m = pow(weightR_m, 8.0);
			
				// // smoothUp[face.owner] += weightR*weightR_m*kappa[face.neighbour];
				// // smoothUp[face.neighbour] += weightL*weightL_m*kappa[face.owner];
				
				// // smoothDown[face.owner] += weightR*weightR_m;
				// // smoothDown[face.neighbour] += weightL*weightL_m;
			
			// // }
			
		// // }
		
		// // if(size>1){
			// // // processor faces
			// // vector<double> sendValues;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // sendValues.push_back(kappa[face.owner]);
				// // }
			// // }
			// // vector<double> recvValues;
			// // mpi.setProcsFaceDatas(
						// // sendValues, recvValues,
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);
			// // int num=0;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
			
				// // double wCL = face.wC;
				// // // double wCL = 0.5;
				// // double wCR = 1.0 - wCL;
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					// // double weightR = pow(1.0-2.0*abs(0.5-smoothAi_recv[num]), 8.0);
					// // double weightR_m = 0.0;
					// // weightR_m += surfNormalVec_recv[0][num]*face.unitNomalsPN[0];
					// // weightR_m += surfNormalVec_recv[1][num]*face.unitNomalsPN[1];
					// // weightR_m += surfNormalVec_recv[2][num]*face.unitNomalsPN[2];
					// // weightR_m = abs(weightR_m);
					// // weightR_m = pow(weightR_m, 8.0);
					
					// // smoothUp[face.owner] += weightR*weightR_m*recvValues[num];
				
					// // smoothDown[face.owner] += weightR*weightR_m;
					
					// // ++num;
				// // }
			// // }
			
		// // }
	
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // kappa[i] = smoothUp[i]/( smoothDown[i] + 1.e-200 );
		// // }
		
	// // }
	
	
	
	
	
	
// // }