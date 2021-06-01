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
	

	vector<SEMO_Species> species;
	
	SEMO_Controls_Builder controls;
		
	controls.readSpecies(species);
	
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	
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
		
		SEMO_Mesh_Save save;
		
		while(
		controls.iterReal<controls.iterRealMax ||
		controls.time<stod(controls.stopAt) ){
			
	
			if(rank==0) {
				cout << "| real-time step = " << controls.iterReal 
				<< " | time = " << controls.time << endl;
			}
		
			// if(controls.application=="incomPBased"){
				
				// solvers.incompressiblePressureBased(mesh, controls, species);
				
			// }
			// else if(controls.application=="comRhoBasedSingle"){
				
				// solvers.compressibleDensityBasedSingleTime(mesh, controls, species);
				
			// }
			// else if(controls.application=="comRhoBasedDual"){
				
				// solvers.compressibleDensityBasedDualTime(mesh, controls, species);
				
			// }
			// else if(controls.application=="incomPBased+comRhoBasedDual"){
				
				solvers.hybridBased(mesh, controls, species);
				
			// }
			
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
				save.vtu(foldername, mesh, controls, species);
			}
		}
		
	}
	
	
	if(rank==0) cout << "| End Program" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	return EXIT_SUCCESS;
}














