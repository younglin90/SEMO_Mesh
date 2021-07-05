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
	math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
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
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// // math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// // math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	
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



	// calc recon. zero order
	this->reconZeroOrder(mesh, controls, species);
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
}




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
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// // math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// // math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// // math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// // math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// // math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// // math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
	
	
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
		
		mpi.setProcsFaceDatasDouble(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
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
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF);
	
	
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
		
		mpi.setProcsFaceDatasDouble(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
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
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF);
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataFromEOS(mesh, controls, species);
	
	
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
		
		mpi.setProcsFaceDatasDouble(
					phi_send0, phi_recv0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					phi_send1, phi_recv1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
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
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF);
	
	
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
		vector<double> volumeFractions;
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
			face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
			volumeFractions );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varL[controls.fVF[ns]] = volumeFractions[ns];
		}


		massFractions.clear();
		volumeFractions.clear();
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
			face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			volumeFractions );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varR[controls.fVF[ns]] = volumeFractions[ns];
		}
		
	}


}



