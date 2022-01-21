#include "build.h"
#include <cmath>
#include <array>


//================ Segregated ===============

void SEMO_Solvers_Builder::setCompValuesLeftRightFaceForSegregated(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconIncomZeroOrder(mesh, controls, species);
	
	// calc recon. 1st order
	// calc gradient
	SEMO_Utility_Math math;
	
	
	vector<vector<double>> gradP;
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	// double underCoeff = 0.4;
	double underCoeff = 0.9;
	// double underCoeff = 1.0;
	// boundary faces treatment : pressure
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo != -1) continue;
		
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		
		if( boundary.type[controls.P] == "zeroGradient" ){
		// if( boundary.type[controls.U] == "noSlip" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				face.varL[controls.fP] += underCoeff*
				 ( (gradP[face.owner][0])*face.distCells[0]
				  +(gradP[face.owner][1])*face.distCells[1]
				  +(gradP[face.owner][2])*face.distCells[2] );
				face.varR[controls.fP] = face.varL[controls.fP];
			}
		}
	}
	
	
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
	// this->calcOtherDataIncomFromEOS(mesh, controls, species);
	
	
}


void SEMO_Solvers_Builder::setCompValuesLeftRightFaceWithReconPVForSegregated(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconIncomZeroOrder(mesh, controls, species);
	
	
	// calc recon. 1st order
	// calc gradient
	SEMO_Utility_Math math;
	
	
	
	
	vector<vector<double>> gradP;
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	
	// double underCoeff = 0.4;
	double underCoeff = 0.9;
	// double underCoeff = 1.0;
	// boundary faces treatment : pressure
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo != -1) continue;
		
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		
		if( boundary.type[controls.P] == "zeroGradient" ){
		// if( boundary.type[controls.U] == "noSlip" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				face.varL[controls.fP] += underCoeff*
				 ( (gradP[face.owner][0])*face.distCells[0]
				  +(gradP[face.owner][1])*face.distCells[1]
				  +(gradP[face.owner][2])*face.distCells[2] );
				face.varR[controls.fP] = face.varL[controls.fP];
			}
		}
	}
	
	
	
	

	//=========================
	vector<double> limGradP;
	vector<double> limGradU;
	vector<double> limGradV;
	vector<double> limGradW;
	
	
	if( controls.gradScheme == "Gauss linear" ){
		math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<3; ++j){
				gradP[i][j] *= limGradP[i]; 
			}
		}
		// internal faces
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double distFace[3];
				
				distFace[0] = face.x-mesh.cells[face.owner].x;
				distFace[1] = face.y-mesh.cells[face.owner].y;
				distFace[2] = face.z-mesh.cells[face.owner].z;
				
				face.varL[controls.fP] +=  
					 (gradP[face.owner][0]*distFace[0]+
					  gradP[face.owner][1]*distFace[1]+
					  gradP[face.owner][2]*distFace[2]);
				
				
				distFace[0] = face.x-mesh.cells[face.neighbour].x;
				distFace[1] = face.y-mesh.cells[face.neighbour].y;
				distFace[2] = face.z-mesh.cells[face.neighbour].z;
				
				face.varR[controls.fP] += 
					 (gradP[face.neighbour][0]*distFace[0]+
					  gradP[face.neighbour][1]*distFace[1]+
					  gradP[face.neighbour][2]*distFace[2]);
				
			}
		}
	}
	else if( controls.gradScheme == "Gauss vanLeer" ){
		this->calcVanLeer(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss MUSCL" ){
		this->calcMUSCL(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss minmod" ){
		this->calcMINMOD(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss QUICK" ){
		this->calcQUICK(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss upwind" ){
	}
	else{
		cerr << "| #Error : not defined gradSchemes at system/fvSchemes file" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	if( controls.divScheme == "Gauss linear" ){
		math.calcLimiterGradient(mesh, controls.U, controls.fU, gradU, limGradU);
		math.calcLimiterGradient(mesh, controls.V, controls.fV, gradV, limGradV);
		math.calcLimiterGradient(mesh, controls.W, controls.fW, gradW, limGradW);
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<3; ++j){
				gradU[i][j] *= limGradU[i];
				gradV[i][j] *= limGradV[i];
				gradW[i][j] *= limGradW[i];
			}
		}
		// internal faces
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double distFace[3];
				
				distFace[0] = face.x-mesh.cells[face.owner].x;
				distFace[1] = face.y-mesh.cells[face.owner].y;
				distFace[2] = face.z-mesh.cells[face.owner].z;
				
				face.varL[controls.fU] += 
					 (gradU[face.owner][0]*distFace[0]+
					  gradU[face.owner][1]*distFace[1]+
					  gradU[face.owner][2]*distFace[2]);
				face.varL[controls.fV] +=  
					 (gradV[face.owner][0]*distFace[0]+
					  gradV[face.owner][1]*distFace[1]+
					  gradV[face.owner][2]*distFace[2]);
				face.varL[controls.fW] += 
					 (gradW[face.owner][0]*distFace[0]+
					  gradW[face.owner][1]*distFace[1]+
					  gradW[face.owner][2]*distFace[2]);
				
				
				distFace[0] = face.x-mesh.cells[face.neighbour].x;
				distFace[1] = face.y-mesh.cells[face.neighbour].y;
				distFace[2] = face.z-mesh.cells[face.neighbour].z;
				
				face.varR[controls.fU] +=  
					 (gradU[face.neighbour][0]*distFace[0]+
					  gradU[face.neighbour][1]*distFace[1]+
					  gradU[face.neighbour][2]*distFace[2]);
				face.varR[controls.fV] +=  
					 (gradV[face.neighbour][0]*distFace[0]+
					  gradV[face.neighbour][1]*distFace[1]+
					  gradV[face.neighbour][2]*distFace[2]);
				face.varR[controls.fW] +=  
					 (gradW[face.neighbour][0]*distFace[0]+
					  gradW[face.neighbour][1]*distFace[1]+
					  gradW[face.neighbour][2]*distFace[2]);
				
			}
		}
	
	}
	else if( controls.divScheme == "Gauss vanLeer" ){
		this->calcVanLeer(mesh, controls, controls.U, controls.fU, gradU);
		this->calcVanLeer(mesh, controls, controls.V, controls.fV, gradV);
		this->calcVanLeer(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss MUSCL" ){
		this->calcMUSCL(mesh, controls, controls.U, controls.fU, gradU);
		this->calcMUSCL(mesh, controls, controls.V, controls.fV, gradV);
		this->calcMUSCL(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss minmod" ){
		this->calcMINMOD(mesh, controls, controls.U, controls.fU, gradU);
		this->calcMINMOD(mesh, controls, controls.V, controls.fV, gradV);
		this->calcMINMOD(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss QUICK" ){
		this->calcQUICK(mesh, controls, controls.U, controls.fU, gradU);
		this->calcQUICK(mesh, controls, controls.V, controls.fV, gradV);
		this->calcQUICK(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss upwind" ){
	}
	else{
		cerr << "| #Error : not defined divSchemes at system/fvSchemes file" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	
	
	
	
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
	// this->calcOtherDataIncomFromEOS(mesh, controls, species);
	
	
}


//================ Coupled ===============

void SEMO_Solvers_Builder::setCompValuesLeftRightFace(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// cout << "aaaaa1" << endl;
	// calc recon. zero order
	this->reconZeroOrder(mesh, controls, species);
	
	// cout << "aaaaa2" << endl;
	
	
	
	
	
	{
		// calc gradient
		SEMO_Utility_Math math;
		
		mesh.cellsGradientVar[controls.P].resize(mesh.cells.size(),vector<double>(3,0.0));
		mesh.cellsGradientVar[controls.U].resize(mesh.cells.size(),vector<double>(3,0.0));
		mesh.cellsGradientVar[controls.V].resize(mesh.cells.size(),vector<double>(3,0.0));
		mesh.cellsGradientVar[controls.W].resize(mesh.cells.size(),vector<double>(3,0.0));
		mesh.cellsGradientVar[controls.T].resize(mesh.cells.size(),vector<double>(3,0.0));
		// mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
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
						mesh.cellsGradientVar[controls.T][face.owner][ii] = 0.0;
						// mesh.cellsGradientVar[controls.VF[0]][face.owner][ii] = 0.0;
					}
				}
			}
		}
		{
			vector<double> dummyVec;
			math.calcGaussGreen(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
			math.calcGaussGreen(mesh, controls.U, controls.fU, mesh.cellsGradientVar[controls.U]);
			math.calcGaussGreen(mesh, controls.V, controls.fV, mesh.cellsGradientVar[controls.V]);
			math.calcGaussGreen(mesh, controls.W, controls.fW, mesh.cellsGradientVar[controls.W]);
			math.calcGaussGreen(mesh, controls.T, controls.fT, mesh.cellsGradientVar[controls.T]);
				// controls.P, controls.fP, dummyVec, mesh.cellsGradientVar[controls.P]);
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				// controls.U, controls.fU, dummyVec, mesh.cellsGradientVar[controls.U]);
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				// controls.V, controls.fV, dummyVec, mesh.cellsGradientVar[controls.V]);
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				// controls.W, controls.fW, dummyVec, mesh.cellsGradientVar[controls.W]);
			// math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				// controls.T, controls.fT, dummyVec, mesh.cellsGradientVar[controls.T]);
		}
		
		
		
		// this->calcMINMOD(mesh, controls, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		this->calcVanLeer(mesh, controls, controls.U, controls.fU, mesh.cellsGradientVar[controls.U]);
		this->calcVanLeer(mesh, controls, controls.V, controls.fV, mesh.cellsGradientVar[controls.V]);
		this->calcVanLeer(mesh, controls, controls.W, controls.fW, mesh.cellsGradientVar[controls.W]);
		this->calcMINMOD(mesh, controls, controls.T, controls.fT, mesh.cellsGradientVar[controls.T]);
		
		
		// vector<double> limGradP;
		// vector<double> limGradU;
		// vector<double> limGradV;
		// vector<double> limGradW;
		// vector<double> limGradT;
		// math.calcLimiterGradient(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P], limGradP);
		// math.calcLimiterGradient(mesh, controls.U, controls.fU, mesh.cellsGradientVar[controls.U], limGradU);
		// math.calcLimiterGradient(mesh, controls.V, controls.fV, mesh.cellsGradientVar[controls.V], limGradV);
		// math.calcLimiterGradient(mesh, controls.W, controls.fW, mesh.cellsGradientVar[controls.W], limGradW);
		// math.calcLimiterGradient(mesh, controls.T, controls.fT, mesh.cellsGradientVar[controls.T], limGradT);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// for(int j=0; j<3; ++j){
				// mesh.cellsGradientVar[controls.P][i][j] *= limGradP[i];
				// mesh.cellsGradientVar[controls.U][i][j] *= limGradU[i];
				// mesh.cellsGradientVar[controls.V][i][j] *= limGradV[i];
				// mesh.cellsGradientVar[controls.W][i][j] *= limGradW[i];
				// mesh.cellsGradientVar[controls.T][i][j] *= limGradT[i];
			// }
		// }
		// // internal faces
		// for(auto& face : mesh.faces){
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// double distFace[3];
				
				// distFace[0] = face.x-mesh.cells[face.owner].x;
				// distFace[1] = face.y-mesh.cells[face.owner].y;
				// distFace[2] = face.z-mesh.cells[face.owner].z;
				
				// // face.varL[controls.fP] += 
					 // // (mesh.cellsGradientVar[controls.P][face.owner][0]*distFace[0]+
					  // // mesh.cellsGradientVar[controls.P][face.owner][1]*distFace[1]+
					  // // mesh.cellsGradientVar[controls.P][face.owner][2]*distFace[2]);
				// face.varL[controls.fU] += 
					 // (mesh.cellsGradientVar[controls.U][face.owner][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.U][face.owner][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.U][face.owner][2]*distFace[2]);
				// face.varL[controls.fV] +=  
					 // (mesh.cellsGradientVar[controls.V][face.owner][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.V][face.owner][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.V][face.owner][2]*distFace[2]);
				// face.varL[controls.fW] += 
					 // (mesh.cellsGradientVar[controls.W][face.owner][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.W][face.owner][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.W][face.owner][2]*distFace[2]);
				// // face.varL[controls.fT] += 
					 // // (mesh.cellsGradientVar[controls.T][face.owner][0]*distFace[0]+
					  // // mesh.cellsGradientVar[controls.T][face.owner][1]*distFace[1]+
					  // // mesh.cellsGradientVar[controls.T][face.owner][2]*distFace[2]);
				
				
				// distFace[0] = face.x-mesh.cells[face.neighbour].x;
				// distFace[1] = face.y-mesh.cells[face.neighbour].y;
				// distFace[2] = face.z-mesh.cells[face.neighbour].z;
				
				// // face.varR[controls.fP] +=  
					 // // (mesh.cellsGradientVar[controls.P][face.neighbour][0]*distFace[0]+
					  // // mesh.cellsGradientVar[controls.P][face.neighbour][1]*distFace[1]+
					  // // mesh.cellsGradientVar[controls.P][face.neighbour][2]*distFace[2]);
				// face.varR[controls.fU] +=  
					 // (mesh.cellsGradientVar[controls.U][face.neighbour][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.U][face.neighbour][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.U][face.neighbour][2]*distFace[2]);
				// face.varR[controls.fV] +=  
					 // (mesh.cellsGradientVar[controls.V][face.neighbour][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.V][face.neighbour][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.V][face.neighbour][2]*distFace[2]);
				// face.varR[controls.fW] +=  
					 // (mesh.cellsGradientVar[controls.W][face.neighbour][0]*distFace[0]+
					  // mesh.cellsGradientVar[controls.W][face.neighbour][1]*distFace[1]+
					  // mesh.cellsGradientVar[controls.W][face.neighbour][2]*distFace[2]);
				// // face.varR[controls.fT] +=  
					 // // (mesh.cellsGradientVar[controls.T][face.neighbour][0]*distFace[0]+
					  // // mesh.cellsGradientVar[controls.T][face.neighbour][1]*distFace[1]+
					  // // mesh.cellsGradientVar[controls.T][face.neighbour][2]*distFace[2]);
				
			// }
		// }
	
	}
	
	
	
	
	
	
	
	// MSTACS
	{
		// calc gradient
		SEMO_Utility_Math math;
		
		vector<vector<double>> gradMF(mesh.cells.size(),vector<double>(3,0.0));
		
		{
			vector<double> dummyVec;
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.MF[0], controls.fMF[0], dummyVec, gradMF);
		}
		

		// processor faces
		if(size>1){
			vector<double> phi_send0, phi_recv0;
			vector<double> phi_send1, phi_recv1;
			vector<double> phi_send2, phi_recv2;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					phi_send0.push_back(gradMF[face.owner][0]);
					phi_send1.push_back(gradMF[face.owner][1]);
					phi_send2.push_back(gradMF[face.owner][2]);
				}
			}
			
			SEMO_MPI_Builder mpi;
			
			mpi.setProcsFaceDatas(
						phi_send0, phi_recv0,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						phi_send1, phi_recv1,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						phi_send2, phi_recv2,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
						
			int proc_num=0;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					vector<double> tmp;
					tmp.push_back(phi_recv0[proc_num]);
					tmp.push_back(phi_recv1[proc_num]);
					tmp.push_back(phi_recv2[proc_num]);
					gradMF.push_back(tmp);
					++proc_num;
				}
			}
		}
		
		
		// NVD, MSTACS
		this->calcMSTACS(mesh, controls, controls.MF[0], controls.fMF[0], gradMF, controls.fMF[0]);
	}
	
	
	
	
	// rho, C, Ht from EOS
	// this->calcOtherDataFromEOS(mesh, controls, species);
	this->calcOtherDataFromEOSMF(mesh, controls, species);
	// cout << "aaaaa3" << endl;
	
	
}




void SEMO_Solvers_Builder::setCompValuesLeftRightFaceWithVfMSTACS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconZeroOrder(mesh, controls, species);
	
	
	
	// calc gradient
	SEMO_Utility_Math math;
	
	vector<vector<double>> gradMF;
	
	math.calcGaussGreen(mesh, controls.MF[0], controls.fMF[0], gradMF);
	// math.calcLeastSquare2nd(mesh, controls.MF[0], controls.fMF[0], gradMF);
	
	

	// processor faces
	if(size>1){
		vector<double> phi_send0, phi_recv0;
		vector<double> phi_send1, phi_recv1;
		vector<double> phi_send2, phi_recv2;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send0.push_back(gradMF[face.owner][0]);
				phi_send1.push_back(gradMF[face.owner][1]);
				phi_send2.push_back(gradMF[face.owner][2]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send2, phi_recv2,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				vector<double> tmp;
				tmp.push_back(phi_recv0[proc_num]);
				tmp.push_back(phi_recv1[proc_num]);
				tmp.push_back(phi_recv2[proc_num]);
				gradMF.push_back(tmp);
				++proc_num;
			}
		}
	}
	
	
	
	// rho, C, Ht from EOS
	// this->calcOtherDataFromEOS(mesh, controls, species);
	this->calcOtherDataFromEOSMF(mesh, controls, species);
	
	
	
	
	
	// NVD, MSTACS
	for(auto& face : mesh.faces){
		face.varL[controls.fMF_HO[0]] = face.varL[controls.fMF[0]];
		face.varR[controls.fMF_HO[0]] = face.varR[controls.fMF[0]];
	}
	
	this->calcMSTACS(mesh, controls, controls.MF[0], controls.fMF[0], gradMF, controls.fMF_HO[0]);
	
	// for(auto& face : mesh.faces){
		// face.varL[controls.fMF_HO[0]] = max(1.e-8,min(1.0-1.e-8,face.varL[controls.fMF_HO[0]]));
		// face.varR[controls.fMF_HO[0]] = max(1.e-8,min(1.0-1.e-8,face.varR[controls.fMF_HO[0]]));
	// }
	 
	
	
	
	
	double dummy;
	for(auto& face : mesh.faces){
		
		vector<double> massFractions;
		vector<double> volumeFractions(controls.nSp,0.0);
		vector<double> dummyVec(controls.nSp,0.0);
		// vector<double> dHtDMF(controls.nSp,0.0);
		double MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			massFractions.push_back(face.varL[controls.fMF_HO[ns]]);
			MFnSp += face.varL[controls.fMF_HO[ns]];
		}
		massFractions.push_back(1.0 - MFnSp);
		
		
		this->getValuesFromEOSMF(
			species,
			face.varL[controls.fP], 
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW], 
			face.varL[controls.fT], massFractions,
			volumeFractions, face.varL[controls.fRho_HO], dummy, face.varL[controls.fHt_HO],
			dummy, dummy, 
			dummy, dummy,
			dummyVec, dummyVec );
			
			
		// for(int ns=0; ns<controls.nSp-1; ++ns){
			// face.varL[controls.fVF[ns]] = volumeFractions[ns];
			// face.varL[controls.fdRhoDMF[ns]] = dRhoDMF[ns];
			// face.varL[controls.fdHtDMF[ns]] = dHtDMF[ns];
		// }


		massFractions.clear();
		MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			massFractions.push_back(face.varR[controls.fMF_HO[ns]]);
			MFnSp += face.varR[controls.fMF_HO[ns]];
		}
		massFractions.push_back(1.0 - MFnSp);
		
		
		
		this->getValuesFromEOSMF(
			species,
			face.varR[controls.fP], 
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW], 
			face.varR[controls.fT], massFractions,
			volumeFractions, face.varR[controls.fRho_HO], dummy, face.varR[controls.fHt_HO],
			dummy, dummy, 
			dummy, dummy,
			dummyVec, dummyVec );
			
		// for(int ns=0; ns<controls.nSp-1; ++ns){
			// face.varR[controls.fVF[ns]] = volumeFractions[ns];
			// face.varR[controls.fdRhoDMF[ns]] = dRhoDMF[ns];
			// face.varR[controls.fdHtDMF[ns]] = dHtDMF[ns];
		// }
		
	}

	
}



//=========================================






void SEMO_Solvers_Builder::setCompValuesLeftRightFaceWithRecon(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconZeroOrder(mesh, controls, species);
	
	
	// calc recon. 1st order
	// calc gradient
	SEMO_Utility_Math math;
	
	
	
	
	vector<vector<double>> gradP;
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	
	// // double underCoeff = 0.4;
	// double underCoeff = 0.9;
	// // double underCoeff = 1.0;
	// // boundary faces treatment : pressure
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo != -1) continue;
		
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		
		// if( boundary.type[controls.P] == "zeroGradient" ){
		// // if( boundary.type[controls.U] == "noSlip" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				// face.varL[controls.fP] += underCoeff*
				 // ( (gradP[face.owner][0])*face.distCells[0]
				  // +(gradP[face.owner][1])*face.distCells[1]
				  // +(gradP[face.owner][2])*face.distCells[2] );
				// face.varR[controls.fP] = face.varL[controls.fP];
			// }
		// }
	// }
	
	
	
	

	//=========================
	vector<double> limGradP;
	vector<double> limGradU;
	vector<double> limGradV;
	vector<double> limGradW;
	
	
	if( controls.gradScheme == "Gauss linear" ){
		math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<3; ++j){
				gradP[i][j] *= limGradP[i]; 
			}
		}
		// internal faces
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double distFace[3];
				
				distFace[0] = face.x-mesh.cells[face.owner].x;
				distFace[1] = face.y-mesh.cells[face.owner].y;
				distFace[2] = face.z-mesh.cells[face.owner].z;
				
				face.varL[controls.fP] +=  
					 (gradP[face.owner][0]*distFace[0]+
					  gradP[face.owner][1]*distFace[1]+
					  gradP[face.owner][2]*distFace[2]);
				
				
				distFace[0] = face.x-mesh.cells[face.neighbour].x;
				distFace[1] = face.y-mesh.cells[face.neighbour].y;
				distFace[2] = face.z-mesh.cells[face.neighbour].z;
				
				face.varR[controls.fP] += 
					 (gradP[face.neighbour][0]*distFace[0]+
					  gradP[face.neighbour][1]*distFace[1]+
					  gradP[face.neighbour][2]*distFace[2]);
				
			}
		}
	}
	else if( controls.gradScheme == "Gauss vanLeer" ){
		this->calcVanLeer(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss MUSCL" ){
		this->calcMUSCL(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss minmod" ){
		this->calcMINMOD(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss QUICK" ){
		this->calcQUICK(mesh, controls, controls.P, controls.fP, gradP);
	}
	else if( controls.gradScheme == "Gauss upwind" ){
	}
	else{
		cerr << "| #Error : not defined gradSchemes at system/fvSchemes file" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	if( controls.divScheme == "Gauss linear" ){
		math.calcLimiterGradient(mesh, controls.U, controls.fU, gradU, limGradU);
		math.calcLimiterGradient(mesh, controls.V, controls.fV, gradV, limGradV);
		math.calcLimiterGradient(mesh, controls.W, controls.fW, gradW, limGradW);
		for(int i=0; i<mesh.cells.size(); ++i){
			for(int j=0; j<3; ++j){
				gradU[i][j] *= limGradU[i];
				gradV[i][j] *= limGradV[i];
				gradW[i][j] *= limGradW[i];
			}
		}
		// internal faces
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double distFace[3];
				
				distFace[0] = face.x-mesh.cells[face.owner].x;
				distFace[1] = face.y-mesh.cells[face.owner].y;
				distFace[2] = face.z-mesh.cells[face.owner].z;
				
				face.varL[controls.fU] += 
					 (gradU[face.owner][0]*distFace[0]+
					  gradU[face.owner][1]*distFace[1]+
					  gradU[face.owner][2]*distFace[2]);
				face.varL[controls.fV] +=  
					 (gradV[face.owner][0]*distFace[0]+
					  gradV[face.owner][1]*distFace[1]+
					  gradV[face.owner][2]*distFace[2]);
				face.varL[controls.fW] += 
					 (gradW[face.owner][0]*distFace[0]+
					  gradW[face.owner][1]*distFace[1]+
					  gradW[face.owner][2]*distFace[2]);
				
				
				distFace[0] = face.x-mesh.cells[face.neighbour].x;
				distFace[1] = face.y-mesh.cells[face.neighbour].y;
				distFace[2] = face.z-mesh.cells[face.neighbour].z;
				
				face.varR[controls.fU] +=  
					 (gradU[face.neighbour][0]*distFace[0]+
					  gradU[face.neighbour][1]*distFace[1]+
					  gradU[face.neighbour][2]*distFace[2]);
				face.varR[controls.fV] +=  
					 (gradV[face.neighbour][0]*distFace[0]+
					  gradV[face.neighbour][1]*distFace[1]+
					  gradV[face.neighbour][2]*distFace[2]);
				face.varR[controls.fW] +=  
					 (gradW[face.neighbour][0]*distFace[0]+
					  gradW[face.neighbour][1]*distFace[1]+
					  gradW[face.neighbour][2]*distFace[2]);
				
			}
		}
	
	}
	else if( controls.divScheme == "Gauss vanLeer" ){
		this->calcVanLeer(mesh, controls, controls.U, controls.fU, gradU);
		this->calcVanLeer(mesh, controls, controls.V, controls.fV, gradV);
		this->calcVanLeer(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss MUSCL" ){
		this->calcMUSCL(mesh, controls, controls.U, controls.fU, gradU);
		this->calcMUSCL(mesh, controls, controls.V, controls.fV, gradV);
		this->calcMUSCL(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss minmod" ){
		this->calcMINMOD(mesh, controls, controls.U, controls.fU, gradU);
		this->calcMINMOD(mesh, controls, controls.V, controls.fV, gradV);
		this->calcMINMOD(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss QUICK" ){
		this->calcQUICK(mesh, controls, controls.U, controls.fU, gradU);
		this->calcQUICK(mesh, controls, controls.V, controls.fV, gradV);
		this->calcQUICK(mesh, controls, controls.W, controls.fW, gradW);
	}
	else if( controls.divScheme == "Gauss upwind" ){
	}
	else{
		cerr << "| #Error : not defined divSchemes at system/fvSchemes file" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	
	vector<vector<double>> gradVF;
	
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	

	// processor faces
	if(size>1){
		vector<double> phi_send0, phi_recv0;
		vector<double> phi_send1, phi_recv1;
		vector<double> phi_send2, phi_recv2;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send0.push_back(gradVF[face.owner][0]);
				phi_send1.push_back(gradVF[face.owner][1]);
				phi_send2.push_back(gradVF[face.owner][2]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send2, phi_recv2,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				vector<double> tmp;
				tmp.push_back(phi_recv0[proc_num]);
				tmp.push_back(phi_recv1[proc_num]);
				tmp.push_back(phi_recv2[proc_num]);
				gradVF.push_back(tmp);
				++proc_num;
			}
		}
	}
	
	
	// NVD, MSTACS
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF, controls.fVF[0]);
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
}




void SEMO_Solvers_Builder::setCompValuesLeftRightFaceForSegregatedWithVfMSTACS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconIncomZeroOrder(mesh, controls, species);
	
	
	

	
	// calc gradient
	SEMO_Utility_Math math;
	
	vector<vector<double>> gradVF;
	
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	

	// processor faces
	if(size>1){
		vector<double> phi_send0, phi_recv0;
		vector<double> phi_send1, phi_recv1;
		vector<double> phi_send2, phi_recv2;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send0.push_back(gradVF[face.owner][0]);
				phi_send1.push_back(gradVF[face.owner][1]);
				phi_send2.push_back(gradVF[face.owner][2]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi_send2, phi_recv2,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				vector<double> tmp;
				tmp.push_back(phi_recv0[proc_num]);
				tmp.push_back(phi_recv1[proc_num]);
				tmp.push_back(phi_recv2[proc_num]);
				gradVF.push_back(tmp);
				++proc_num;
			}
		}
	}
	
	
	// NVD, MSTACS
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF, controls.fVF[0]);
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
}






void SEMO_Solvers_Builder::reconZeroOrder(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
			// (*mesh.boundary[0].setFunctionP)(1.0, 1.0, 1.0, 1.0);

	// internal faces
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
			face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
			face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
			face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
			face.varL[controls.fT] = mesh.cells[face.owner].var[controls.T];
			face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
			face.varL[controls.fMF[0]] = mesh.cells[face.owner].var[controls.MF[0]];
			face.varL[controls.fmu] = mesh.cells[face.owner].var[controls.mu];
			
			
			face.varR[controls.fP] = mesh.cells[face.neighbour].var[controls.P];
			face.varR[controls.fU] = mesh.cells[face.neighbour].var[controls.U];
			face.varR[controls.fV] = mesh.cells[face.neighbour].var[controls.V];
			face.varR[controls.fW] = mesh.cells[face.neighbour].var[controls.W];
			face.varR[controls.fT] = mesh.cells[face.neighbour].var[controls.T];
			face.varR[controls.fVF[0]] = mesh.cells[face.neighbour].var[controls.VF[0]];
			face.varR[controls.fMF[0]] = mesh.cells[face.neighbour].var[controls.MF[0]];
			face.varR[controls.fmu] = mesh.cells[face.neighbour].var[controls.mu];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
			face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
			face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
			face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
			face.varL[controls.fT] = mesh.cells[face.owner].var[controls.T];
			face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
			face.varL[controls.fMF[0]] = mesh.cells[face.owner].var[controls.MF[0]];
			face.varL[controls.fmu] = mesh.cells[face.owner].var[controls.mu];
		}
	}
	
	
	if(size>1){
		// processor faces's right values
		mpi.setCellDatasToFaceRight(mesh, 
					controls.P, controls.fP,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.U, controls.fU,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh,
					controls.V, controls.fV,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.W, controls.fW,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.T, controls.fT,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.VF[0], controls.fVF[0],
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.MF[0], controls.fMF[0],
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setCellDatasToFaceRight(mesh, 
					controls.mu, controls.fmu,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	}	

	// // boundary faces
	for(auto& boundary : mesh.boundary){
		// if(rank==0) cout << boundary.name << endl;
		
		if(boundary.neighbProcNo != -1) continue;
		
		
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		
		if( boundary.type[controls.P] == "fixedValue" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fP] = boundary.var[controls.P];
				face.varL[controls.fP] = 0.5*boundary.var[controls.P] + 
										 0.5*mesh.cells[face.owner].var[controls.P];
				face.varR[controls.fP] = face.varL[controls.fP];
			}
		}
		else if( boundary.type[controls.P] == "zeroGradient" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				face.varR[controls.fP] = face.varL[controls.fP];
			}
		}
		else if( boundary.type[controls.P] == "switch" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				double machNum = 
					sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
					     pow(mesh.cells[face.owner].var[controls.V],2.0)+
						 pow(mesh.cells[face.owner].var[controls.W],2.0))/
						  mesh.cells[face.owner].var[controls.C];
				if( machNum > 1.0 ){
					face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				}
				else{
					face.varL[controls.fP] = boundary.var[controls.P];
				}
				face.varR[controls.fP] = face.varL[controls.fP];
			}
		}
		// else if( boundary.type[controls.P] == "function" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// // boundary.setFunctionP(mesh, face, controls.fP, controls.P);
			// }
		// }
		else{
			cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// if(rank==0)cout << boundary.values[controls.VF[0]] << endl;
		
		if( boundary.type[controls.U] == "fixedValue" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fU] = boundary.var[controls.U];
				// face.varL[controls.fV] = boundary.var[controls.V];
				// face.varL[controls.fW] = boundary.var[controls.W];
				face.varL[controls.fU] = 0.5*boundary.var[controls.U] + 
										 0.5*mesh.cells[face.owner].var[controls.U];
				face.varL[controls.fV] = 0.5*boundary.var[controls.V] + 
										 0.5*mesh.cells[face.owner].var[controls.V];
				face.varL[controls.fW] = 0.5*boundary.var[controls.W] + 
										 0.5*mesh.cells[face.owner].var[controls.W];
				
				face.varR[controls.fU] = boundary.var[controls.U];
				face.varR[controls.fV] = boundary.var[controls.V];
				face.varR[controls.fW] = boundary.var[controls.W];
			}
		}
		else if( boundary.type[controls.U] == "zeroGradient" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
				face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
				face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
			}
		}
		else if( boundary.type[controls.U] == "slip" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				double U = mesh.cells[face.owner].var[controls.U];
				double V = mesh.cells[face.owner].var[controls.V];
				double W = mesh.cells[face.owner].var[controls.W];
				
				double norVel = U * face.unitNormals[0] + 
								V * face.unitNormals[1] + 
								W * face.unitNormals[2];
				
				double invU = U - norVel * face.unitNormals[0];
				double invV = V - norVel * face.unitNormals[1];
				double invW = W - norVel * face.unitNormals[2];
				
				face.varL[controls.fU] = invU;
				face.varL[controls.fV] = invV;
				face.varL[controls.fW] = invW;
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
			}
		}
		else if( boundary.type[controls.U] == "noSlip" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fU] = 0.0;
				face.varL[controls.fV] = 0.0;
				face.varL[controls.fW] = 0.0;
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
			}
		}
		else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
			double norVel = boundary.var[controls.U];
				// cout << boundary.var[controls.U] << endl;
			
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
			
				
				// face.varL[controls.fU] = norVel*face.unitNormals[0];
				// face.varL[controls.fV] = norVel*face.unitNormals[1];
				// face.varL[controls.fW] = norVel*face.unitNormals[2];
				face.varL[controls.fU] = 0.5*norVel*face.unitNormals[0]; + 
										 0.5*mesh.cells[face.owner].var[controls.U];
				face.varL[controls.fV] = 0.5*norVel*face.unitNormals[1]; + 
										 0.5*mesh.cells[face.owner].var[controls.V];
				face.varL[controls.fW] = 0.5*norVel*face.unitNormals[2]; + 
										 0.5*mesh.cells[face.owner].var[controls.W];
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
			}
			
		}
		else if( boundary.type[controls.U] == "inletOutlet" ){
			double norVel = boundary.var[controls.U];
				// cout << boundary.var[controls.U] << endl;
			
					
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				double ownNorVel =  
					mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
					mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
					mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					
				if( ownNorVel > 0.0 ){
					// outflow
					// zeroGradient
					face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
					face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
					face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
					
					face.varR[controls.fU] = face.varL[controls.fU];
					face.varR[controls.fV] = face.varL[controls.fV];
					face.varR[controls.fW] = face.varL[controls.fW];
					
				}
				else{
					// inflow
					// fixedValue
					
					// face.varL[controls.fU] = norVel*face.unitNormals[0];
					// face.varL[controls.fV] = norVel*face.unitNormals[1];
					// face.varL[controls.fW] = norVel*face.unitNormals[2];
					face.varL[controls.fU] = 0.5*norVel*face.unitNormals[0]; + 
											 0.5*mesh.cells[face.owner].var[controls.U];
					face.varL[controls.fV] = 0.5*norVel*face.unitNormals[1]; + 
											 0.5*mesh.cells[face.owner].var[controls.V];
					face.varL[controls.fW] = 0.5*norVel*face.unitNormals[2]; + 
											 0.5*mesh.cells[face.owner].var[controls.W];
					
					face.varR[controls.fU] = face.varL[controls.fU];
					face.varR[controls.fV] = face.varL[controls.fV];
					face.varR[controls.fW] = face.varL[controls.fW];
				}
			
				
			}
			
		}
		else if( boundary.type[controls.U] == "function" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				double phi;
				(*boundary.setFunctionVariables[1])(controls.time, face.x, face.y, face.z, phi);
				face.varL[controls.fU] = phi;
				face.varR[controls.fU] = phi;
				(*boundary.setFunctionVariables[2])(controls.time, face.x, face.y, face.z, phi);
				// cout << phi << endl;
				face.varL[controls.fV] = phi;
				face.varR[controls.fV] = phi;
				(*boundary.setFunctionVariables[3])(controls.time, face.x, face.y, face.z, phi);
				face.varL[controls.fW] = phi;
				face.varR[controls.fW] = phi;
			}
		}
		else{
			cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		
		if( boundary.type[controls.T] == "fixedValue" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fT] = boundary.var[controls.T];
				face.varR[controls.fT] = face.varL[controls.fT];
			}
		}
		else if( boundary.type[controls.T] == "zeroGradient" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fT] = mesh.cells[face.owner].var[controls.T];
				face.varR[controls.fT] = face.varL[controls.fT];
			}
		}
		else if( boundary.type[controls.T] == "function" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				// boundary.setFunctionT(mesh, face, controls.fT, controls.T);
			}
		}
		else{
			cerr << "| #Error : not defined B.C., var = " << controls.T << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		
		// if( boundary.type[controls.MF[0]] == "fixedValue" ){
		if( boundary.type[controls.VF[0]] == "fixedValue" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				face.varL[controls.fVF[0]] = boundary.var[controls.VF[0]];
				face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
				
				// face.varL[controls.fMF[0]] = boundary.var[controls.MF[0]];
				face.varL[controls.fMF[0]] = boundary.var[controls.VF[0]];
				face.varR[controls.fMF[0]] = face.varL[controls.fMF[0]];
			}
		}
		// else if( boundary.type[controls.MF[0]] == "zeroGradient" ){
		else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
				face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
				
				face.varL[controls.fMF[0]] = mesh.cells[face.owner].var[controls.MF[0]];
				face.varR[controls.fMF[0]] = face.varL[controls.fMF[0]];
			}
		}
		else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				double ownNorVel =  
					mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
					mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
					mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
					
				if( ownNorVel > 0.0 ){
					// outflow
					face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
					face.varL[controls.fMF[0]] = mesh.cells[face.owner].var[controls.MF[0]];
				}
				else{
					face.varL[controls.fVF[0]] = boundary.var[controls.VF[0]];
					face.varL[controls.fMF[0]] = boundary.var[controls.VF[0]];
				}
				
				face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
				face.varR[controls.fMF[0]] = face.varL[controls.fMF[0]];
			}
		}
		// else if( boundary.type[controls.VF[0]] == "function" ){
			
			// // CSV format
			// vector<vector<double>> readDatas;
			// {
				// string filename("./VF.csv");
				// string file_contents;
				// vector<vector<string>> csv_contents;
				// char delimiter = ',';
				// {
					// auto ss = ostringstream{};
					// ifstream input_file(filename);
					// if (!input_file.is_open()) {
						// cerr << "Could not open the file - '"
							 // << filename << "'" << endl;
						// exit(EXIT_FAILURE);
					// }
					// ss << input_file.rdbuf();
					// file_contents = ss.str();
				// }
				// istringstream sstream(file_contents);
				// vector<string> items;
				// string record;

				// int counter = 0;
				// while (std::getline(sstream, record)) {
					// istringstream line(record);
					// while (std::getline(line, record, delimiter)) {
						// items.push_back(record);
					// }

					// csv_contents.push_back(items);
					// items.clear();
					// counter += 1;
				// }
				// counter = 0;
				// for(auto& contents : csv_contents){
					// ++counter;
					// if(counter==1) continue;
					// vector<double> tmp_vec;
					// bool itemIsNumber = true;
					// for(auto& item : contents){
						// // for (char const &c : item) {
							// // if (std::isdigit(c) == 0) {
								// // itemIsNumber = false;
								// // break;
							// // }
						// // }
						// // if (std::isdigit(item[0]) != 0) {
						// // if (counttt != 0) {
						// // if (item.find_first_not_of("-+0123456789") == string::npos) {
							// tmp_vec.push_back(stod(item));
							// // cout << item << endl;
						// // }
					// }
					// readDatas.push_back(tmp_vec);
				// }
			// }
			
			
			
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				
				// double dist_save = 1.e9;
				// double VF_save;
				// for(auto& item : readDatas){
					// double VF_tmp = item[3];
					// double x_tmp = item[4];
					// double y_tmp = item[5];
					// double z_tmp = item[6];
					
					// double distance2 = (pow(face.x-x_tmp,2.0) + pow(face.y-y_tmp,2.0) + pow(face.z-z_tmp,2.0));
					// if(distance2 < dist_save){
						// dist_save = distance2;
						// VF_save = VF_tmp;
					// }
				// }
				
				// face.varL[controls.fVF[0]] = VF_save;
				// face.varL[controls.fMF[0]] = VF_save;
				
				// face.varR[controls.fVF[0]] = VF_save;
				// face.varR[controls.fMF[0]] = VF_save;
			// }
			
			
		// }
		else if( boundary.type[controls.VF[0]] == "function" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				double phi;
				(*boundary.setFunctionVariables[5])(controls.time, face.x, face.y, face.z, phi);
				face.varL[controls.fMF[0]] = phi;
				face.varR[controls.fMF[0]] = phi;
				face.varL[controls.fVF[0]] = phi;
				face.varR[controls.fVF[0]] = phi;
				
				// boundary.setFromFunction("./lib", mesh, face, controls.fVF[0], controls.VF[0]);
				// boundary.setFunctionMF(mesh, face, controls.fMF[0], controls.MF[0]);
			}
		}
		else{
			cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		// mu
		for(int i=str; i<end; ++i){
			SEMO_Face& face = mesh.faces[i];
			face.varL[controls.fmu] = mesh.cells[face.owner].var[controls.mu];
			face.varR[controls.fmu] = face.varL[controls.fmu];
		}
		
		
	}

}










void SEMO_Solvers_Builder::calcOtherDataFromEOS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){

	for(auto& face : mesh.faces){
		
		vector<double> massFractions;
		vector<double> volumeFractions;
		double MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			volumeFractions.push_back(face.varL[controls.fVF[ns]]);
			MFnSp += face.varL[controls.fVF[ns]];
		}
		volumeFractions.push_back(1.0 - MFnSp);
		
		this->getValuesFromEOSVF(
			species,
			face.varL[controls.fP], 
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW], 
			face.varL[controls.fT], volumeFractions, 
			face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
			massFractions );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varL[controls.fMF[ns]] = massFractions[ns];
		}


		massFractions.clear();
		volumeFractions.clear();
		MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			volumeFractions.push_back(face.varR[controls.fVF[ns]]);
			MFnSp += face.varR[controls.fVF[ns]];
		}
		volumeFractions.push_back(1.0 - MFnSp);
		
		this->getValuesFromEOSVF(
			species,
			face.varR[controls.fP], 
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW], 
			face.varR[controls.fT], volumeFractions, 
			face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			massFractions );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varR[controls.fMF[ns]] = massFractions[ns];
		}
		
	}


}






void SEMO_Solvers_Builder::calcOtherDataFromEOSMF(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){

	for(auto& face : mesh.faces){
		
		vector<double> massFractions;
		vector<double> volumeFractions(controls.nSp,0.0);
		vector<double> dRhoDMF(controls.nSp,0.0);
		vector<double> dHtDMF(controls.nSp,0.0);
		double MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			massFractions.push_back(face.varL[controls.fMF[ns]]);
			MFnSp += face.varL[controls.fMF[ns]];
		}
		massFractions.push_back(1.0 - MFnSp);
		
		
		this->getValuesFromEOSMF(
			species,
			face.varL[controls.fP], 
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW], 
			face.varL[controls.fT], massFractions,
			volumeFractions, face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
			face.varL[controls.fdRhoDP], face.varL[controls.fdHtDP], 
			face.varL[controls.fdRhoDT], face.varL[controls.fdHtDT],
			dRhoDMF, dHtDMF );
			
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varL[controls.fVF[ns]] = volumeFractions[ns];
			face.varL[controls.fdRhoDMF[ns]] = dRhoDMF[ns];
			face.varL[controls.fdHtDMF[ns]] = dHtDMF[ns];
		}


		massFractions.clear();
		MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			massFractions.push_back(face.varR[controls.fMF[ns]]);
			MFnSp += face.varR[controls.fMF[ns]];
		}
		massFractions.push_back(1.0 - MFnSp);
		
		
		
		this->getValuesFromEOSMF(
			species,
			face.varR[controls.fP], 
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW], 
			face.varR[controls.fT], massFractions,
			volumeFractions, face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			face.varR[controls.fdRhoDP], face.varR[controls.fdHtDP], 
			face.varR[controls.fdRhoDT], face.varR[controls.fdHtDT],
			dRhoDMF, dHtDMF );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varR[controls.fVF[ns]] = volumeFractions[ns];
			face.varR[controls.fdRhoDMF[ns]] = dRhoDMF[ns];
			face.varR[controls.fdHtDMF[ns]] = dHtDMF[ns];
		}
		
	}


}



