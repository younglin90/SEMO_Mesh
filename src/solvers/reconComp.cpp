#include "build.h"
#include <cmath>
#include <array>

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
	
	math.calcLeastSquare(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	// NVD
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
			face.varL[controls.fmu] = mesh.cells[face.owner].var[controls.mu];
			
			
			face.varR[controls.fP] = mesh.cells[face.neighbour].var[controls.P];
			face.varR[controls.fU] = mesh.cells[face.neighbour].var[controls.U];
			face.varR[controls.fV] = mesh.cells[face.neighbour].var[controls.V];
			face.varR[controls.fW] = mesh.cells[face.neighbour].var[controls.W];
			face.varR[controls.fT] = mesh.cells[face.neighbour].var[controls.T];
			face.varR[controls.fVF[0]] = mesh.cells[face.neighbour].var[controls.VF[0]];
			face.varR[controls.fmu] = mesh.cells[face.neighbour].var[controls.mu];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.varL[controls.fP] = mesh.cells[face.owner].var[controls.P];
			face.varL[controls.fU] = mesh.cells[face.owner].var[controls.U];
			face.varL[controls.fV] = mesh.cells[face.owner].var[controls.V];
			face.varL[controls.fW] = mesh.cells[face.owner].var[controls.W];
			face.varL[controls.fT] = mesh.cells[face.owner].var[controls.T];
			face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
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
				face.varL[controls.fP] = boundary.var[controls.P];
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
				face.varL[controls.fU] = boundary.var[controls.U];
				face.varL[controls.fV] = boundary.var[controls.V];
				face.varL[controls.fW] = boundary.var[controls.W];
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
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
			
				
				face.varL[controls.fU] = norVel*face.unitNormals[0];
				face.varL[controls.fV] = norVel*face.unitNormals[1];
				face.varL[controls.fW] = norVel*face.unitNormals[2];
				
				face.varR[controls.fU] = face.varL[controls.fU];
				face.varR[controls.fV] = face.varL[controls.fV];
				face.varR[controls.fW] = face.varL[controls.fW];
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
			}
		}
		// else if( boundary.type[controls.MF[0]] == "zeroGradient" ){
		else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				face.varL[controls.fVF[0]] = mesh.cells[face.owner].var[controls.VF[0]];
				face.varR[controls.fVF[0]] = face.varL[controls.fVF[0]];
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