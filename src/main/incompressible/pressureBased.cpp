#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>

#include "parmetis.h" 
#include "scotch.h" 

#include "../../mesh/build.h" 
#include "../../mesh/load.h" 
#include "../../mesh/geometric.h" 

#include "../../controls/build.h" 
#include "../../variables/build.h" 
#include "../../solvers/build.h" 


int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// cout << "AAAAAAAA" << endl;

	vector<SEMO_Species> species;
	
	SEMO_Controls_Builder controls;
		
	controls.readSpecies(species);
	
	// cout << "BB" << endl;
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	// cout << "CC" << endl;
	
	SEMO_Solvers_Builder solvers;
	
	SEMO_Mesh_Builder mesh;
	
	SEMO_MPI_Builder mpi;
	
	
	bool boolLoad = true;
	bool boolPartitioning = false;
	bool boolAMR = false;
	bool boolGeometric = true;
	
	if(boolLoad){
		
		SEMO_Mesh_Load load;
		double starttime = stod(controls.startFrom);
		string foldername;
		std::ostringstream streamObj;
		streamObj << starttime;
		foldername = "./save/" + streamObj.str() + "/";
		
		if(starttime == 0.0){
			foldername = "./save/0/";
		}
		
		load.vtu(foldername, mesh, controls, species);
		
		
		solvers.calcCellEOSVF(mesh, controls, species);
		
		for(auto& cell : mesh.cells){
			vector<double> Q(controls.nEq,0.0);
			Q[0] = cell.var[controls.Rho];
			Q[1] = cell.var[controls.Rho] * cell.var[controls.U];
			Q[2] = cell.var[controls.Rho] * cell.var[controls.V];
			Q[3] = cell.var[controls.Rho] * cell.var[controls.W];
			Q[4] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			for(int ns=0; ns<controls.nSp-1; ++ns){
				Q[5+ns] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
				
				// cout << cell.var[controls.MF[ns]] << endl;
			}
			for(int i=0; i<controls.nEq; ++i){
				cell.var[controls.Qm[i]] = Q[i];
				cell.var[controls.Qn[i]] = Q[i];
			}
		}
	
	}
	else {
		
		// load serial OpenFoam files
		mesh.loadFile("OpenFoam", "./grid/");
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
	
	// partitioning
	if(boolPartitioning){
		
		mesh.distributeOneToAll("EVENLY");
		
	}
	
	
	
	if(boolAMR){
		
		mesh.hexaOctAMR();
		
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	

	
	// geometric
	if(boolGeometric){
		
		SEMO_Utility_Math math;
		
		SEMO_Mesh_Geometric geometric;
		
		geometric.init(mesh);
		
		// math.initLeastSquare(mesh);
		math.initLeastSquare2nd(mesh);
	
	}
	
	

	bool initFlow = false;
	
	// flow initialization
	if(initFlow){
		
		solvers.setInitValues(mesh, controls);
		
		mesh.saveFile("vtu", "./save/0/", controls);
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
	


	bool calcFlow = true;
	
	// flow calculation
	if(calcFlow){
		
		SEMO_Mesh_Save save;
		
		
		

	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// auto& cell = mesh.cells[i];
		
		// if(cell.y < -0.157){
			// cell.var[controls.VF[0]] = 0.0;
		// }
		
	// }
		
				// string foldername2;
				// // if( controls.time < 1.e-6 ){
					// std::ostringstream streamObj2;
					// streamObj2 << controls.time;
					// foldername2 = "./save/" + streamObj2.str() + "/";
	// save.vtu(foldername2, mesh, controls, species);
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		
		while(
		controls.iterReal<controls.iterRealMax ||
		controls.time<stod(controls.stopAt) ){
			
	
			if(rank==0) {
				cout << "| real-time step = " << controls.iterReal 
				<< " | time = " << controls.time;
				// << " | timeStep = " << controls.timeStep << endl;
			}
		
			solvers.incompressiblePressureBased(mesh, controls, species);
				
			controls.time += controls.timeStep;
			
			++controls.iterReal;
			
			if(controls.iterReal % controls.saveInterval == 0){
				string foldername;
				// if( controls.time < 1.e-6 ){
					std::ostringstream streamObj;
					streamObj << controls.time;
					foldername = "./save/" + streamObj.str() + "/";
				// }
				// else{
					// foldername = "./save/" + to_string(controls.time) + "/";
				// }
				
				solvers.setCompValuesLeftRightFace(mesh, controls, species);
				
				save.vtu(foldername, mesh, controls, species);
			}
		}
		
	}
	
	
	if(rank==0) cout << "| End Program" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	return EXIT_SUCCESS;
}










// int main(int argc, char* argv[]) {


	// MPI::Init(); 
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// // cout << "AAAAAAAA" << endl;

	// vector<SEMO_Species> species;
	
	// SEMO_Controls_Builder controls;
		
	// controls.readSpecies(species);
	
	// // cout << "BB" << endl;
	
	// controls.readConfigures();
	
	// controls.setValues(species);
	
	// // cout << "CC" << endl;
	
	// SEMO_Solvers_Builder solvers;
	
	// SEMO_Mesh_Builder mesh;
	
	// SEMO_MPI_Builder mpi;
	
	
	// bool boolLoad = true;
	// bool boolPartitioning = false;
	// bool boolAMR = false;
	// bool boolGeometric = true;
	
	// if(boolLoad){
		
		// SEMO_Mesh_Load load;
		// double starttime = stod(controls.startFrom);
		// string foldername;
		// std::ostringstream streamObj;
		// streamObj << starttime;
		// foldername = "./save/" + streamObj.str() + "/";
		
		// if(starttime == 0.0){
			// foldername = "./save/0/";
		// }
		
		// load.vtu(foldername, mesh, controls, species);
		
		
		// solvers.calcCellEOSVF(mesh, controls, species);
		
		// for(auto& cell : mesh.cells){
			// vector<double> Q(controls.nEq,0.0);
			// Q[0] = cell.var[controls.Rho];
			// Q[1] = cell.var[controls.Rho] * cell.var[controls.U];
			// Q[2] = cell.var[controls.Rho] * cell.var[controls.V];
			// Q[3] = cell.var[controls.Rho] * cell.var[controls.W];
			// Q[4] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			// for(int ns=0; ns<controls.nSp-1; ++ns){
				// Q[5+ns] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
				
				// // cout << cell.var[controls.MF[ns]] << endl;
			// }
			// for(int i=0; i<controls.nEq; ++i){
				// cell.var[controls.Qm[i]] = Q[i];
				// cell.var[controls.Qn[i]] = Q[i];
			// }
		// }
	
	// }
	// else {
		
		// // load serial OpenFoam files
		// mesh.loadFile("OpenFoam", "./grid/");
	// }
	
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
	
	// // partitioning
	// if(boolPartitioning){
		
		// mesh.distributeOneToAll("EVENLY");
		
	// }
	
	
	
	// if(boolAMR){
		
		// mesh.hexaOctAMR();
		
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	

	
	// // geometric
	// if(boolGeometric){
		
		// SEMO_Utility_Math math;
		
		// SEMO_Mesh_Geometric geometric;
		
		// geometric.init(mesh);
		
		// math.initLeastSquare(mesh);
	
	// }
	
	

	// bool initFlow = false;
	
	// // flow initialization
	// if(initFlow){
		
		// solvers.setInitValues(mesh, controls);
		
		// mesh.saveFile("vtu", "./save/0/", controls);
		
		// MPI_Barrier(MPI_COMM_WORLD);
	// }
	
	
	


	// bool calcFlow = true;
	
	// // flow calculation
	// if(calcFlow){
		
		// SEMO_Mesh_Save save;
		
		// while(
		// controls.iterReal<controls.iterRealMax ||
		// controls.time<stod(controls.stopAt) ){
			
	
			// if(rank==0) {
				// cout << "| real-time step = " << controls.iterReal 
				// << " | time = " << controls.time << endl;
			// }
		
			// // if(controls.application=="incomPBased"){
				
				// solvers.incompressiblePressureBased(mesh, controls, species);
				
			// // }
			// // else if(controls.application=="comRhoBasedSingle"){
				
				// // solvers.compressibleDensityBasedSingleTime(mesh, controls, species);
				
			// // }
			// // else if(controls.application=="comRhoBasedDual"){
				
				// // solvers.compressibleDensityBasedDualTime(mesh, controls, species);
				
			// // }
			// // else if(controls.application=="incomPBased+comRhoBasedDual"){
				
				// // solvers.hybridBased(mesh, controls, species);
				
			// // }
			
			
			
			
			// for(auto& i : mesh.cells){
				
				// i.var[0] = i.x*exp(-i.x*i.x - i.y*i.y);
				
			// }
			
			

			
			// // // calc gradient
			// // vector<vector<double>> gradResiP;
			
			// // gradResiP.clear();
			// // gradResiP.resize(mesh.cells.size(),vector<double>(3,0.0));
			
			
			// // for(auto& face : mesh.faces){
				
				// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
					
					// // double DVar = 
						// // mesh.cells[face.neighbour].var[0] - mesh.cells[face.owner].var[0];
					
					// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
					// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
					// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
					
					// // gradResiP[face.neighbour][0] += face.distCells[0] * DVar;
					// // gradResiP[face.neighbour][1] += face.distCells[1] * DVar;
					// // gradResiP[face.neighbour][2] += face.distCells[2] * DVar;
					
				// // }
				
			// // }
			
			
			// // // boundary
			// // for(auto& boundary : mesh.boundary){
				
				// // if(boundary.neighbProcNo == -1){
					
					// // int str = boundary.startFace;
					// // int end = str + boundary.nFaces;
					
					// // if(boundary.type[0] == "fixedValue"){
						// // for(int i=str; i<end; ++i){
							// // auto& face = mesh.faces[i];
							
							// // double DVar = 
								// // 0.0 - mesh.cells[face.owner].var[0];
							
							// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
							// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
							// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
							
							
						// // }
					// // }
				// // }
			// // }
			
			
			// // // processor faces
			// // vector<double> sendValues;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // sendValues.push_back(mesh.cells[face.owner].var[0]);
				// // }
			// // }
			
			// // vector<double> recvValues;
			// // mpi.setProcsFaceDatasDouble(
						// // sendValues, recvValues,
						// // mesh.countsProcFaces, mesh.countsProcFaces, 
						// // mesh.displsProcFaces, mesh.displsProcFaces);

			// // int num=0;
			// // for(int i=0; i<mesh.faces.size(); ++i){
				// // auto& face = mesh.faces[i];
				
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // double DVar = 
						// // recvValues[num] - mesh.cells[face.owner].var[0];
					
					// // gradResiP[face.owner][0] += face.distCells[0] * DVar;
					// // gradResiP[face.owner][1] += face.distCells[1] * DVar;
					// // gradResiP[face.owner][2] += face.distCells[2] * DVar;
							
					// // ++num;
				// // }
			// // }
			
			
			
			// // for(int i=0; i<mesh.cells.size(); ++i){
				// // SEMO_Cell& cell = mesh.cells[i];
				
				// // double tmp0 = 
					// // cell.coeffLeastSquare[0] * gradResiP[i][0] +
					// // cell.coeffLeastSquare[1] * gradResiP[i][1] +
					// // cell.coeffLeastSquare[2] * gradResiP[i][2];
					
				// // double tmp1 = 
					// // cell.coeffLeastSquare[1] * gradResiP[i][0] +
					// // cell.coeffLeastSquare[3] * gradResiP[i][1] +
					// // cell.coeffLeastSquare[4] * gradResiP[i][2];
					
				// // double tmp2 = 
					// // cell.coeffLeastSquare[2] * gradResiP[i][0] +
					// // cell.coeffLeastSquare[4] * gradResiP[i][1] +
					// // cell.coeffLeastSquare[5] * gradResiP[i][2];
					
				// // cell.var[1] = tmp0;
				// // cell.var[2] = tmp1;
				// // cell.var[3] = tmp2;
				
			// // }
			
					
			
			
			

	
	// // velocities correction
	// vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double PF = 0.5*(mesh.cells[face.owner].var[0]+mesh.cells[face.neighbour].var[0]);
			
			// for(int j=0; j<3; ++j){
				// gradResiP[face.owner][j] += 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// gradResiP[face.neighbour][j] -= 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			// }
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
					
					// double PF = 0.5*(mesh.cells[face.owner].var[0]);
					// for(int j=0; j<3; ++j){
						// gradResiP[face.owner][j] += 
							// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// }
					
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					
					// double PF = mesh.cells[face.owner].var[0];
					// for(int j=0; j<3; ++j){
						// gradResiP[face.owner][j] += 
							// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// }
					
				// }
			// }
		// }
	// }
	
	// // processor faces
	// vector<double> sendValues;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// sendValues.push_back(mesh.cells[face.owner].var[0]);
		// }
	// }
	// vector<double> recvValues;
	// mpi.setProcsFaceDatasDouble(
				// sendValues, recvValues,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// int num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
			// double PF = 0.5*(mesh.cells[face.owner].var[0]+recvValues[num]);
			
			// for(int j=0; j<3; ++j){
				// gradResiP[face.owner][j] += 
					// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			// }
			// ++num;
		// }
	// }
			
			
			// for(int i=0; i<mesh.cells.size(); ++i){
				// SEMO_Cell& cell = mesh.cells[i];
				
				// cell.var[1] = gradResiP[i][0];
				// cell.var[2] = gradResiP[i][1];
				// cell.var[3] = gradResiP[i][2];
				
			// }
			
			
			
			
			// controls.time += controls.timeStep;
			
			// ++controls.iterReal;
			
			// if(controls.iterReal % controls.saveInterval == 0){
				// string foldername;
				// // if( controls.time < 1.e-6 ){
					// std::ostringstream streamObj;
					// streamObj << controls.time;
					// foldername = "./save/" + streamObj.str() + "/";
				// // }
				// // else{
					// // foldername = "./save/" + to_string(controls.time) + "/";
				// // }
				
				// solvers.setCompValuesLeftRightFace(mesh, controls, species);
				
				// save.vtu(foldername, mesh, controls, species);
			// }
		// }
		
	// }
	
	
	// if(rank==0) cout << "| End Program" << endl;
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	// return EXIT_SUCCESS;
// }














