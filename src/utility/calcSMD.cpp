#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>

#include <zlib.h>

using namespace std;

#include "../mesh/build.h" 
#include "../load/load.h" 
#include "../mesh/geometric.h" 

#include "../solvers/build.h" 
#include "../controls/build.h" 
#include "../variables/build.h" 
#include "../mesh/polyAMR.h"


void saveInitialField(string folder);
	
vector<SEMO_Species> species;
SEMO_Controls_Builder controls;
SEMO_MPI_Builder mpi;
SEMO_Mesh_Load load;
SEMO_Utility_Math math;
SEMO_Mesh_Geometric geometric;
SEMO_Mesh_Save save;
SEMO_Solvers_Builder solvers;


	
int main(int argc, char* argv[]) {
	

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	
	map<string,string> mapArgv;
	
	for(int i=1; i<argc; i+=2){
		string first = argv[i];
		string second;
		if(i+1==argc){
			second = "nan";
		}
		else{
			second = argv[i+1];
		}
		mapArgv.insert(make_pair(first, second));
		// cout <<  first<< " " << second << endl;
	}
	
	if( 
	mapArgv.find("-help") != mapArgv.end() ||
	mapArgv.find("-h") != mapArgv.end()
	){
		if(rank==0){
			cout << endl;
			cout << endl;
			cout << "┌─────── Potential Solver's helper ─────────────────────────────── " << endl;
			cout << "| -s \"real num.\" : start time" << endl;
			cout << "| -e \"real num.\" : end time" << endl;
			cout << "| -t \"real num.\" : save time step" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		MPI_Finalize();
		// exit(EXIT_FAILURE);
		exit(EXIT_SUCCESS);
		return EXIT_SUCCESS;
	}
	
	
	
	double srtTime = stod(mapArgv["-s"]);
	double endTime = stod(mapArgv["-e"]);
	double stepTime = stod(mapArgv["-t"]);
	// double cflInput = stod(mapArgv["-c"]);
	// int isave = stoi(mapArgv["-n"]);
	// inpTimeStep
	
	
	// 초기화
	controls.readSpecies(species);
	controls.readConfigures();
	controls.setValues(species);
	
	
	// 메쉬 및 데이터 읽기
	
	vector<double> strBox(3,0.0);
	vector<double> endBox(3,0.0);
	
	strBox[0] = -0.0012;
	strBox[1] = -0.1575;
	strBox[2] = -0.0012;
	
	endBox[0] = 0.0012;
	endBox[1] = -0.1571;
	endBox[2] = 0.0012;
	
	
	
	
	
	double D3 = 0.0;
	double D2 = 0.0;
	double time=srtTime;
	while(time<=endTime)
	{
		
		SEMO_Mesh_Builder mesh;
		{
			
			SEMO_Mesh_Load load;
			string foldername;
			std::ostringstream streamObj;
			streamObj << time;
			foldername = "./save/" + streamObj.str() + "/";
			if(time == 0.0) foldername = "./save/0/";
			
			string openFileName = foldername + "/plot.0.vtu";
			ifstream inputFile;
			inputFile.open(openFileName);
			if(inputFile.fail()){
				time += stepTime;
				continue;
			}
			
			load.vtu(foldername, mesh, controls, species);
			
			solvers.calcCellEOSVF(mesh, controls, species);
			
		}
		
		// geometric 계산
		{
			geometric.init(mesh);
			math.initLeastSquare(mesh);
		}

		
		// SMD 계산 
		{
			
			mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
			
			vector<double> dummyVec;
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.VF[0], controls.fVF[0], dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
			
			// if(size>1){
				// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.VF[0]], mesh.cellsProcGradientVar[controls.VF[0]]);
			// }
			
			
			for(int i=0, ip=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				
				if(
				cell.x >= strBox[0] && cell.x <= endBox[0] &&
				cell.y >= strBox[1] && cell.y <= endBox[1] &&
				cell.z >= strBox[2] && cell.z <= endBox[2])
				{
				
					
					vector<double> gradVF(3,0.0);
					for(int ii=0; ii<3; ++ii){
						gradVF[ii] = mesh.cellsGradientVar[controls.VF[0]][i][ii];
					}
					
					double magGrad=0.0;
					for(int ii=0; ii<3; ++ii){
						magGrad += gradVF[ii]*gradVF[ii];
					}
					magGrad = sqrt(magGrad);
					D2 += magGrad*cell.volume;
					
					D3 += cell.var[controls.VF[0]]*cell.volume;
				}
				
			}
			
		
		}
	
		// if(rank==0) cout << 
		// " | time = " << time << endl;
		
		double D3Reduced;
		double D2Reduced;
		MPI_Allreduce(&D3, &D3Reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&D2, &D2Reduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// print
		if(rank==0)
		{
			cout << " | time = " << time << endl;
			cout << "| Result : D3 = " << D3Reduced << ", D2 = " << D2Reduced << ", D32(SMD) = " << D3Reduced/D2Reduced << endl;
			
		}
			
		time += stepTime;
	}

		
	
	// {
		// string foldername;
		// std::ostringstream streamObj;
		// streamObj << controls.time;
		// foldername = "./save/" + streamObj.str() + "_SMD/";
		
		// save.vtu(foldername, mesh, controls, species);
	// }



	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	exit(EXIT_SUCCESS);
	return EXIT_SUCCESS;
}














