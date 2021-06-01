#include "build.h"
#include <cmath>
#include <array>

#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::setIncomValuesLeftRightFace(
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
	
	
	
	
	// vector<vector<double>> gradP;
	// vector<vector<double>> gradU;
	// vector<vector<double>> gradV;
	// vector<vector<double>> gradW;
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	

	// NVD, QUICK
	// this->calcMINMOD(mesh, controls, controls.P, controls.fP, gradP);
	// this->calcMINMOD(mesh, controls, controls.U, controls.fU, gradU);
	// this->calcMINMOD(mesh, controls, controls.V, controls.fV, gradV);
	// this->calcMINMOD(mesh, controls, controls.W, controls.fW, gradW);
	
	
	
	// vector<double> limGradP;
	// vector<double> limGradU;
	// vector<double> limGradV;
	// vector<double> limGradW;
	// math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// math.calcLimiterGradient(mesh, controls.U, controls.fU, gradU, limGradU);
	// math.calcLimiterGradient(mesh, controls.V, controls.fV, gradV, limGradV);
	// math.calcLimiterGradient(mesh, controls.W, controls.fW, gradW, limGradW);
	
	// // internal faces
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double distFace[3];
			
			// distFace[0] = face.x-mesh.cells[face.owner].x;
			// distFace[1] = face.y-mesh.cells[face.owner].y;
			// distFace[2] = face.z-mesh.cells[face.owner].z;
			
			// face.varL[controls.fP] += limGradP[face.owner]
				// *(gradP[face.owner][0]*distFace[0]+
				  // gradP[face.owner][1]*distFace[1]+
				  // gradP[face.owner][2]*distFace[2]);
			// face.varL[controls.fU] += limGradU[face.owner]
				// *(gradU[face.owner][0]*distFace[0]+
				  // gradU[face.owner][1]*distFace[1]+
				  // gradU[face.owner][2]*distFace[2]);
			// face.varL[controls.fV] += limGradV[face.owner]
				// *(gradV[face.owner][0]*distFace[0]+
				  // gradV[face.owner][1]*distFace[1]+
				  // gradV[face.owner][2]*distFace[2]);
			// face.varL[controls.fW] += limGradW[face.owner]
				// *(gradW[face.owner][0]*distFace[0]+
				  // gradW[face.owner][1]*distFace[1]+
				  // gradW[face.owner][2]*distFace[2]);
			
			
			// distFace[0] = face.x-mesh.cells[face.neighbour].x;
			// distFace[1] = face.y-mesh.cells[face.neighbour].y;
			// distFace[2] = face.z-mesh.cells[face.neighbour].z;
			
			// face.varR[controls.fP] += limGradP[face.neighbour]
				// *(gradP[face.neighbour][0]*distFace[0]+
				  // gradP[face.neighbour][1]*distFace[1]+
				  // gradP[face.neighbour][2]*distFace[2]);
			// face.varR[controls.fU] += limGradU[face.neighbour]
				// *(gradU[face.neighbour][0]*distFace[0]+
				  // gradU[face.neighbour][1]*distFace[1]+
				  // gradU[face.neighbour][2]*distFace[2]);
			// face.varR[controls.fV] += limGradV[face.neighbour]
				// *(gradV[face.neighbour][0]*distFace[0]+
				  // gradV[face.neighbour][1]*distFace[1]+
				  // gradV[face.neighbour][2]*distFace[2]);
			// face.varR[controls.fW] += limGradW[face.neighbour]
				// *(gradW[face.neighbour][0]*distFace[0]+
				  // gradW[face.neighbour][1]*distFace[1]+
				  // gradW[face.neighbour][2]*distFace[2]);
			
		// }
	// }
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataIncomFromEOS(mesh, controls, species);
	
	
}







void SEMO_Solvers_Builder::setIncomValuesLeftRightFaceWithVelQUICK(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();



	// calc recon. zero order
	this->reconZeroOrder(mesh, controls, species);
	
	
	
	
	// calc gradient
	SEMO_Utility_Math math;
	
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	
	math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	
	// NVD, QUICK
	this->calcQUICK(mesh, controls, controls.U, controls.fU, gradU);
	this->calcQUICK(mesh, controls, controls.V, controls.fV, gradV);
	this->calcQUICK(mesh, controls, controls.W, controls.fW, gradW);
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataIncomFromEOS(mesh, controls, species);
	
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}




void SEMO_Solvers_Builder::setIncomValuesLeftRightFaceWithVfMSTACS(
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
	
	// NVD, MSTACS
	this->calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF);
	
	
	
	// rho, C, Ht from EOS
	this->calcOtherDataIncomFromEOS(mesh, controls, species);
	
	
}









void SEMO_Solvers_Builder::calcOtherDataIncomFromEOS(
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
		
		double constP = 101325.0;
		double constT = 300.0;
		
		this->getValuesFromEOSVF(
			species,
			constP, 
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW], 
			constT, volumeFractions, 
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
			constP, 
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW], 
			constT, volumeFractions, 
			face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			massFractions );
			
		for(int ns=0; ns<controls.nSp-1; ++ns){
			face.varR[controls.fMF[ns]] = massFractions[ns];
		}
		
	}


}

// void SEMO_Solvers_Builder::setIncomValuesLeftRightFace(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<SEMO_Species>& species){
	
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	

	// // internal faces
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
			// face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
			// face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
			// face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
			// face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
			
			
			// face.varR[controls.fP] = mesh.cells[face.neighbour].var[controls.P];
			// face.varR[controls.fU] = mesh.cells[face.neighbour].var[controls.U];
			// face.varR[controls.fV] = mesh.cells[face.neighbour].var[controls.V];
			// face.varR[controls.fW] = mesh.cells[face.neighbour].var[controls.W];
			// face.varR[controls.fVF[0]] = mesh.cells[face.neighbour].var[controls.VF[0]];
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
			// face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
			// face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
			// face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
			// face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
		// }
	// }
	

	// // processor faces's right values
	// mpi.setCellDatasToFaceRight(mesh, 
				// controls.P, controls.fP,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
				
	// mpi.setCellDatasToFaceRight(mesh, 
				// controls.U, controls.fU,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
				
	// mpi.setCellDatasToFaceRight(mesh,
				// controls.V, controls.fV,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
				
	// mpi.setCellDatasToFaceRight(mesh, 
				// controls.W, controls.fW,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
				
	// mpi.setCellDatasToFaceRight(mesh, 
				// controls.VF[0], controls.fVF[0],
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
				

	// // // boundary faces
	// for(auto& boundary : mesh.boundary){
		// // if(rank==0) cout << boundary.name << endl;
		
		// if(boundary.neighbProcNo != -1) continue;
		
		
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		
		// if( boundary.type[controls.P] == "fixedValue" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fP] = boundary.var[controls.P];
				// face.varR[controls.fP] = face.varL[controls.fP];
			// }
		// }
		// else if( boundary.type[controls.P] == "zeroGradient" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				// face.varR[controls.fP] = face.varL[controls.fP];
			// }
		// }
		// else if( boundary.type[controls.P] == "switch" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// double machNum = 
					// sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
					     // pow(mesh.cells[face.owner].var[controls.V],2.0)+
						 // pow(mesh.cells[face.owner].var[controls.W],2.0))/
						  // mesh.cells[face.owner].var[controls.C];
				// if( machNum > 1.0 ){
					// face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
				// }
				// else{
					// face.varL[controls.fP] = boundary.var[controls.P];
				// }
				// face.varR[controls.fP] = face.varL[controls.fP];
			// }
		// }
		// else{
			// cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		// // if(rank==0)cout << boundary.values[controls.VF[0]] << endl;
		
		// if( boundary.type[controls.U] == "fixedValue" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fU] = boundary.var[controls.U];
				// face.varL[controls.fV] = boundary.var[controls.V];
				// face.varL[controls.fW] = boundary.var[controls.W];
				
				// face.varR[controls.fU] = face.varL[controls.fU];
				// face.varR[controls.fV] = face.varL[controls.fV];
				// face.varR[controls.fW] = face.varL[controls.fW];
			// }
		// }
		// else if( boundary.type[controls.U] == "zeroGradient" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
				// face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
				// face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
				
				// face.varR[controls.fU] = face.varL[controls.fU];
				// face.varR[controls.fV] = face.varL[controls.fV];
				// face.varR[controls.fW] = face.varL[controls.fW];
			// }
		// }
		// else if( boundary.type[controls.U] == "slip" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// double U = mesh.cells[face.owner].var[controls.U];
				// double V = mesh.cells[face.owner].var[controls.V];
				// double W = mesh.cells[face.owner].var[controls.W];
				
				// double norVel = U * face.unitNormals[0] + 
								// V * face.unitNormals[1] + 
								// W * face.unitNormals[2];
				
				// double invU = U - norVel * face.unitNormals[0];
				// double invV = V - norVel * face.unitNormals[1];
				// double invW = W - norVel * face.unitNormals[2];
				
				// face.varL[controls.fU] = invU;
				// face.varL[controls.fV] = invV;
				// face.varL[controls.fW] = invW;
				
				// face.varR[controls.fU] = face.varL[controls.fU];
				// face.varR[controls.fV] = face.varL[controls.fV];
				// face.varR[controls.fW] = face.varL[controls.fW];
			// }
		// }
		// else if( boundary.type[controls.U] == "noSlip" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fU] = 0.0;
				// face.varL[controls.fV] = 0.0;
				// face.varL[controls.fW] = 0.0;
				
				// face.varR[controls.fU] = face.varL[controls.fU];
				// face.varR[controls.fV] = face.varL[controls.fV];
				// face.varR[controls.fW] = face.varL[controls.fW];
			// }
		// }
		// else if( boundary.type[controls.U] == "surfaceNormalFixedValue" ){
			// double norVel = boundary.var[controls.U];
				// // cout << boundary.var[controls.U] << endl;
			
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
			
				
				// face.varL[controls.fU] = norVel*face.unitNormals[0];
				// face.varL[controls.fV] = norVel*face.unitNormals[1];
				// face.varL[controls.fW] = norVel*face.unitNormals[2];
				
				// face.varR[controls.fU] = face.varL[controls.fU];
				// face.varR[controls.fV] = face.varL[controls.fV];
				// face.varR[controls.fW] = face.varL[controls.fW];
			// }
			
		// }
		// else{
			// cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		
		
		// if( boundary.type[controls.T] == "fixedValue" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fT] = boundary.var[controls.T];
				// face.varR[controls.fT] = face.varL[controls.fT];
			// }
		// }
		// else if( boundary.type[controls.T] == "zeroGradient" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fT] = mesh.cells[face.owner].var[controls.T];
				// face.varR[controls.fT] = face.varL[controls.fT];
			// }
		// }
		// else{
			// cerr << "| #Error : not defined B.C., var = " << controls.T << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		
		
		// if( boundary.type[controls.VF[0]] == "fixedValue" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fVF[0]] = boundary.var[controls.VF[0]];
				// face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
			// }
		// }
		// else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
			// for(int i=str; i<end; ++i){
				// SEMO_Face& face = mesh.faces[i];
				// face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
				// face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
			// }
		// }
		// else{
			// cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		
	// }

	
	
			

			
	// // rho, C, Ht from EOS
	// for(auto& face : mesh.faces){
		
		// // vector<double> massFractions;
		// // vector<double> volumeFractions;
		// // double VFnSp = 0.0;
		// // for(int ns=0; ns<controls.nSp-1; ++ns){
			// // volumeFractions.push_back(face.varL[controls.fVF[ns]]);
			// // VFnSp += face.varL[controls.fVF[ns]];
		// // }
		// // volumeFractions.push_back(1.0 - VFnSp);
		
		// // this->getValuesFromEOSVF(
			// // species,
			// // face.varL[controls.fP], 
			// // face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW], 
			// // face.varL[controls.fT], volumeFractions, 
			// // face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
			// // massFractions);



		// // massFractions.clear();
		// // volumeFractions.clear();
		// // VFnSp = 0.0;
		// // for(int ns=0; ns<controls.nSp-1; ++ns){
			// // volumeFractions.push_back(face.varR[controls.fVF[ns]]);
			// // VFnSp += face.varR[controls.fVF[ns]];
		// // }
		// // volumeFractions.push_back(1.0 - VFnSp);
		
		// // this->getValuesFromEOSVF(
			// // species,
			// // face.varR[controls.fP], 
			// // face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW], 
			// // face.varR[controls.fT], volumeFractions, 
			// // face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
			// // massFractions);

		// // for(int ns=0; ns<controls.nSp-1; ++ns){
			// // face.varR[controls.fMF[ns]] = massFractions[ns];
		// // }
		
		
		
		// face.varL[controls.fRho] = 
			// 1000.0*face.varL[controls.fVF[0]] + 
			// 1.0*(1.0-face.varL[controls.fVF[0]]);
		// face.varR[controls.fRho] = 
			// 1000.0*face.varR[controls.fVF[0]] + 
			// 1.0*(1.0-face.varR[controls.fVF[0]]);
			

			
	// }

	
// }