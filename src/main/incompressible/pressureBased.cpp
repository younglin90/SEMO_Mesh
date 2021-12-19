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
#include "../../load/load.h" 
#include "../../mesh/geometric.h" 
#include "../../mesh/polyAMR.h"

#include "../../controls/build.h" 
#include "../../variables/build.h" 
#include "../../solvers/build.h" 


int main(int argc, char* argv[]) {


	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// return EXIT_SUCCESS;
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
	bool boolAMR = true;
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
			cell.var[controls.Qm[0]] = cell.var[controls.U];
			cell.var[controls.Qm[1]] = cell.var[controls.V];
			cell.var[controls.Qm[2]] = cell.var[controls.W];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qm[3+i]] = cell.var[controls.VF[i]];
			}
			
			cell.var[controls.Qn[0]] = cell.var[controls.U];
			cell.var[controls.Qn[1]] = cell.var[controls.V];
			cell.var[controls.Qn[2]] = cell.var[controls.W];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qn[3+i]] = cell.var[controls.VF[i]];
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
	
	
	
	// if(boolAMR){
		
		// // mesh.hexaOctAMR();
		
		// int tmpNum = 0;
		// for(auto& face : mesh.faces){
			// face.group = tmpNum;
			// ++tmpNum;
		// }
		
		// tmpNum = 0;
		// for(auto& cell : mesh.cells){
			// cell.group = tmpNum;
			// ++tmpNum;
		// }
		
		
		// // SEMO_Poly_AMR_Builder AMR;
		// // AMR.polyAMR(mesh, 0);
		// // AMR.polyAMR(mesh, 1);
		// // AMR.polyAMR(mesh, 1);
		
		// // MPI_Barrier(MPI_COMM_WORLD);
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	
	
	
	// geometric
	if(boolGeometric){
		
		SEMO_Utility_Math math;
		
		SEMO_Mesh_Geometric geometric;
		
		geometric.init(mesh);
		
		math.initLeastSquare(mesh);
		
		
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
		
		
		solvers.calcCellTransport(mesh, controls, species);
		
		
		
		SEMO_Mesh_Save save;
		
		
		while(
		controls.iterReal<controls.iterRealMax &&
		controls.time<stod(controls.stopAt) ){
			
	
			if(rank==0) {
				cout << "| real-time step = " << controls.iterReal 
				<< " | time = " << controls.time;
				// << " | timeStep = " << controls.timeStep << endl;
			}
		
			solvers.incompressiblePressureBased(mesh, controls, species);


			//==============================
			// AMR
			// if(controls.iterReal % 10 == 0 && controls.iterReal != 0){
			// if(controls.iterReal % 1 == 0){ 
			if
			( (controls.iterReal+1) % controls.intervalRefine == 0 ||
			  (controls.iterReal+1) % controls.intervalUnrefine == 0)
			{ 
			int numnum=0; 
			// while(1){
				// cout << numnum << endl;
				SEMO_Poly_AMR_Builder AMR;
				AMR.polyAMR(mesh, controls, species, 0);
				mesh.cellsGlobal();
				solvers.calcIncomCellEOSVF(mesh, controls, species);
				solvers.calcCellTransport(mesh, controls, species);
				// ++numnum;
			// }
			} 
			//==============================
			
				
			controls.time += controls.timeStep;
			
			++controls.iterReal;
			
			if(controls.saveControl == "timeStep"){
				if(controls.iterReal % (int)controls.saveInterval == 0){
					string foldername;
					std::ostringstream streamObj;
					streamObj << controls.time;
					foldername = "./save/" + streamObj.str() + "/";
					
					solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					save.vtu(foldername, mesh, controls, species);
					
					// save.particles(foldername, mesh, controls, species);
					
				}
			}
			else if(controls.saveControl == "runTime"){
				int jung = controls.time / controls.saveInterval;
				double namuji = controls.time - (double)jung * controls.saveInterval;
				if(
				namuji < controls.timeStep &&
				namuji >= 0.0
				){
					string foldername;
					std::ostringstream streamObj;
					streamObj << controls.time;
					foldername = "./save/" + streamObj.str() + "/";
					
					solvers.setCompValuesLeftRightFace(mesh, controls, species);
					
					save.vtu(foldername, mesh, controls, species);
					
					// save.particles(foldername, mesh, controls, species);
					
				}
			}
			
		}
		
		
		
	}
	
	
	if(rank==0) cout << "| End Program" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	return EXIT_SUCCESS;
}






