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
#include <sys/stat.h>
#include <cstring>
#include <vector>

using namespace std;
#include "../mesh/build.h" 
#include "../load/load.h" 
#include "../mesh/geometric.h" 

#include "../controls/build.h" 

// void mkdirs(char *dir_path);
// void saveExtractedData(
	// string folder,
	// string fileName,
	// vector<double> X,
	// vector<double> Y,
	// vector<double> var); 
	

	
int main(int argc, char* argv[]) {
	
	// MPI::Init(); 
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	

	// vector<SEMO_Species> species;
	// SEMO_Controls_Builder controls;
		
	// controls.readSpecies(species);
	
	// controls.readConfigures();
	
	// controls.setValues(species);
	
	
	// vector<double> targetTime;
	
	// double timeInitial = 0.00047;
	// for(int i=0; i<40; ++i){
		// double inpTime = timeInitial + 0.000002*(double)i;
		// targetTime.push_back(inpTime);
	// }
	

	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
	// }
	
	
	// for(int iTime=0; iTime<targetTime.size(); ++iTime){
		
		// double starttime = targetTime[iTime];
		
		// bool boolLoad = true;

		// SEMO_Mesh_Builder mesh;
		
		// if(boolLoad){
			
			// SEMO_Mesh_Load load;
			// // double starttime = stod(controls.startFrom);
			// string foldername;
			// std::ostringstream streamObj;
			// streamObj << starttime;
			// foldername = "./save/" + streamObj.str() + "/";
			
			// // if(starttime == 0.0){
				// // foldername = "./save/0/";
			// // }
			
			// load.vtu(foldername, mesh, controls, species);
			
			// SEMO_Mesh_Geometric geometric;
			
			// geometric.init(mesh);
		
		// }
		
		
		// vector<vector<double>> extDataRankX(3,vector<double>());
		// vector<vector<double>> extDataRankY(3,vector<double>());
		// vector<vector<double>> extDataRankVF(3,vector<double>());
		// for(int i=0; i<mesh.cells.size(); ++i){
			// double X = mesh.cells[i].x;
			// double Y = mesh.cells[i].y;
			// double Z = mesh.cells[i].z;

				// // extDataRankX[0].push_back(X);
				// // extDataRankY[0].push_back(Y);
				// // extDataRankVF[0].push_back(mesh.cells[i].var[controls.VF[0]]);
			
			
			
			// if(
			// (X>=0.00027 && X<=0.0045) && (Y>=0.0 && Y<=0.005)
			// ){
				// extDataRankX[0].push_back(X);
				// extDataRankY[0].push_back(Y);
				// extDataRankVF[0].push_back(mesh.cells[i].var[controls.VF[0]]);
			// }
			
			
			
			// // if(
			// // (X>=0.000852096 && X<=0.000865348) && Y<=0.004933
			// // ){
				// // extDataRankX[0].push_back(X);
				// // extDataRankY[0].push_back(Y);
				// // extDataRankVF[0].push_back(mesh.cells[i].var[controls.VF[0]]);
			// // }
			
			// // if(
			// // (X>=0.00246885 && X<=0.0024821) && Y<=0.004933
			// // ){
				// // extDataRankX[1].push_back(X);
				// // extDataRankY[1].push_back(Y);
				// // extDataRankVF[1].push_back(mesh.cells[i].var[controls.VF[0]]);
			// // }
			
			// // if(
			// // (X>=0.00387357 && X<=0.00388682) && Y<=0.004933
			// // ){
				// // extDataRankX[2].push_back(X);
				// // extDataRankY[2].push_back(Y);
				// // extDataRankVF[2].push_back(mesh.cells[i].var[controls.VF[0]]);
			// // }
			
		// }
		

		// vector<vector<double>> extDataX(3,vector<double>());
		// vector<vector<double>> extDataY(3,vector<double>());
		// vector<vector<double>> extDataVF(3,vector<double>());
		
		// // for(int i=0; i<3; ++i){
		// for(int i=0; i<1; ++i){
			// if(size>1){
				
				// if(rank==0){        
					// vector<int> counts(size,0);
					// int sizeData = extDataRankX[i].size();
					// MPI_Gather(&sizeData, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
				
					// vector<int> displacements(size,0);
					// for(int ip=1; ip<size; ++ip){
						// displacements[ip] = displacements[ip-1] + counts[ip-1];
					// }
					
					// extDataX[i].resize(displacements[size-1]+counts[size-1]);
					// MPI_Gatherv(extDataRankX[i].data(), extDataRankX[i].size(), 
						// MPI_DOUBLE, extDataX[i].data(), counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
					// extDataY[i].resize(displacements[size-1]+counts[size-1]);
					// MPI_Gatherv(extDataRankY[i].data(), extDataRankY[i].size(), 
						// MPI_DOUBLE, extDataY[i].data(), counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
					// extDataVF[i].resize(displacements[size-1]+counts[size-1]);
					// MPI_Gatherv(extDataRankVF[i].data(), extDataRankVF[i].size(), 
						// MPI_DOUBLE, extDataVF[i].data(), counts.data(), displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					
				// }
				// else{        
					// int sizeData = extDataRankX[i].size();
					// MPI_Gather(&sizeData, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
					
					// MPI_Gatherv(extDataRankX[i].data(), extDataRankX[i].size(), 
						// MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					// MPI_Gatherv(extDataRankY[i].data(), extDataRankY[i].size(), 
						// MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					
					// MPI_Gatherv(extDataRankVF[i].data(), extDataRankVF[i].size(), 
						// MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				   
				// }
				
				
			// }
			// else{
				
				// // for(auto& j : extDataRankX[i]){
					// // extDataX[i].push_back(j);
				// // }
				// // for(auto& j : extDataRankY[i]){
					// // extDataY[i].push_back(j);
				// // }
				// // for(auto& j : extDataRankVF[i]){
					// // extDataVF[i].push_back(j);
				// // }
				
			// }
		// }
		
		// if(rank==0){
			
			// string fileName;
			// std::ostringstream streamObj;
			// streamObj << fixed;
			// streamObj.precision(8); 
			// streamObj << starttime;
			
			
			// fileName = "regionAll_" + streamObj.str();
			// saveExtractedData("./extract/", fileName, extDataX[0], extDataY[0], extDataVF[0]);
			
			// // fileName = "region0_" + streamObj.str();
			// // saveExtractedData("./extract/", fileName, extDataX[0], extDataY[0], extDataVF[0]);
			
			// // fileName.clear();
			// // fileName = "region1_" + streamObj.str();
			// // saveExtractedData("./extract/", fileName, extDataX[1], extDataY[1], extDataVF[1]);
			
			// // fileName.clear();
			// // fileName = "region2_" + streamObj.str();
			// // saveExtractedData("./extract/", fileName, extDataX[2], extDataY[2], extDataVF[2]);
		// }
		// MPI_Barrier(MPI_COMM_WORLD); 
	// }
	

	// if(rank==0) cout << endl;
	// if(rank==0) cout << "| completed save extract data files" << endl;
	// if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// MPI_Finalize();
	// return 0;
}



// void mkdirs(char *dir_path){
	// char buff[1024];
	// char *p_dir = buff;
	
	// strcpy(buff, dir_path);
	// buff[1024-1] = '\0';
	
	// while(*p_dir){
		// if('/'==*p_dir){
			// *p_dir = '\0';
			// mkdir(buff,0777);
			// *p_dir = '/';
		// }
		// p_dir ++;
	// }
	
// }


// void saveExtractedData(
	// string folder,
	// string fileName,
	// vector<double> X,
	// vector<double> Y,
	// vector<double> var){

	

	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// char folder_name[1000];
	// strcpy(folder_name, folder.c_str());
	// mkdirs(folder_name);

	// ofstream outputFile;
	// string filenamePlot = folder + fileName;
	
	// if(rank==0){
		// // cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute save file (" << folder_name << "extract data...) ... " << endl;
	// }
	
	// outputFile.open(filenamePlot);
	
	// outputFile.precision( 8 );
	
	// if(outputFile.fail()){
		// cerr << "Unable to write file for writing." << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// // MPI_Finalize();
	// }
	
	// // string out_line;
	// for(int i=0; i<var.size(); ++i) 
		// outputFile << scientific << X[i] << " " << Y[i] << " " << var[i] << endl;
	

	// outputFile.close();




// }