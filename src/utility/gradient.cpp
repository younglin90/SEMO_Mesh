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

using namespace std;

#include "../mesh/build.h" 
#include "../load/load.h" 
#include "../mesh/geometric.h" 

#include "../solvers/build.h" 
#include "../controls/build.h" 
#include "../variables/build.h" 

void calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	);
void calcLeastSquare2nd(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	);
void calcLeastSquare1st(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	);
void saveInitialField(string folder);

void initLeastSquare(
	SEMO_Mesh_Builder& mesh);
void initLeastSquare2nd(
	SEMO_Mesh_Builder& mesh);
	
SEMO_Mesh_Builder mesh;
vector<SEMO_Species> species;
SEMO_Controls_Builder controls;
SEMO_MPI_Builder mpi;
SEMO_Mesh_Load load;
SEMO_Utility_Math math;
SEMO_Mesh_Geometric geometric;

double weight;
	
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
			cout << "| -n \"int num.\" : iteration of gradient" << endl;
			cout << "| -s \"int num.\" : gradient schemes, 1-GG, 2-LS1, 3-LS2" << endl;
			cout << "| -w \"double num.\" : weighting value at only LS" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	
	
	int iMax = stoi(mapArgv["-n"]);
	int schemes = stod(mapArgv["-s"]);
	weight = stod(mapArgv["-w"]);
	
	
	
		
	controls.readSpecies(species);
	
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	
	SEMO_Solvers_Builder solvers;
	
	
	
	load.vtu("./save/0/", mesh, controls, species);
	
	
	geometric.init(mesh);
	
	if(schemes==1){
	}
	else if(schemes==2){
		initLeastSquare(mesh);
	}
	else if(schemes==3){
		initLeastSquare2nd(mesh);
	}
	
	
	// UDF value -> pressure
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double x = cell.x;
		double y = cell.y;
		
		// sin & cos function
		// cell.var[controls.P] = sin(x*6.0)*cos(y*6);
		// cell.var[controls.U] = 6.0*cos(6.0*x)*cos(6.0*y);
		// cell.var[controls.V] = -6.0*sin(6.0*x)*sin(6.0*y);
		
		// Sinc function
		// double value = sqrt(pow(x*24.0,2.0)+pow(y*24.0,2.0));
		// cell.var[controls.P] = sin(value)/value;
		// cell.var[controls.U] = x*cos(24.0*sqrt(x*x+y*y))/(x*x+y*y) - x*sin(24.0*sqrt(x*x+y*y))/24.0/pow(x*x+y*y,1.5);
		// cell.var[controls.V] = y*cos(24.0*sqrt(x*x+y*y))/(x*x+y*y) - y*sin(24.0*sqrt(x*x+y*y))/24.0/pow(x*x+y*y,1.5);
		
		// Rosenbrock function
		cell.var[controls.P] = pow(1.0-(x*3.0),2.0)+100.0*pow((y*2.0+0.25)-pow(3.0*x,2.0),2.0);
		cell.var[controls.U] = 32400.0*(-0.000185185+x*x*x+x*(-0.027222-0.222222*y));
		cell.var[controls.V] = -3600.0*(-0.0277778+x*x-0.222222*y);
		
		cell.var[controls.W] = 0.0;
		
	}
	
	
	
	// gradient P
	vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	
	if(schemes==1){
		for(int i=0; i<iMax; ++i){
			calcGaussGreen(mesh, controls.P, controls.fP, gradP);
		}
	}
	else if(schemes==2){
		for(int i=0; i<iMax; ++i){
			calcLeastSquare1st(mesh, controls.P, controls.fP, gradP);
		}
	}
	else if(schemes==3){
		for(int i=0; i<iMax; ++i){
			calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
		}
	}
	

	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.UDV[0]] = gradP[i][0];
		cell.var[controls.UDV[1]] = gradP[i][1];
		cell.var[controls.UDV[2]] = gradP[i][2];
		
	}
	
	
	saveInitialField("./save/0_Gradient/");
	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
}




void initLeastSquare(
	SEMO_Mesh_Builder& mesh) {
	
	vector<vector<double>> vsum(mesh.cells.size(),vector<double>(6,0.0));
	
	
	// double weight = 2.0;
	
	for(auto& face : mesh.faces){
		
		// double wk = 1.0;
		double wk = 1.0 / sqrt(
			pow(face.distCells[0],2.0)+
			pow(face.distCells[1],2.0)+
			pow(face.distCells[2],2.0));
		wk = pow(wk,weight);
		
		vsum[face.owner][0] += wk * face.distCells[0]*face.distCells[0];
		vsum[face.owner][1] += wk * face.distCells[0]*face.distCells[1];
		vsum[face.owner][2] += wk * face.distCells[0]*face.distCells[2];
		vsum[face.owner][3] += wk * face.distCells[1]*face.distCells[1];
		vsum[face.owner][4] += wk * face.distCells[1]*face.distCells[2];
		vsum[face.owner][5] += wk * face.distCells[2]*face.distCells[2];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			vsum[face.neighbour][0] += wk * face.distCells[0]*face.distCells[0];
			vsum[face.neighbour][1] += wk * face.distCells[0]*face.distCells[1];
			vsum[face.neighbour][2] += wk * face.distCells[0]*face.distCells[2];
			vsum[face.neighbour][3] += wk * face.distCells[1]*face.distCells[1];
			vsum[face.neighbour][4] += wk * face.distCells[1]*face.distCells[2];
			vsum[face.neighbour][5] += wk * face.distCells[2]*face.distCells[2];
		}
		
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double detA = 
			vsum[i][0]*vsum[i][3]*vsum[i][5] + vsum[i][1]*vsum[i][4]*vsum[i][2]
		  + vsum[i][2]*vsum[i][1]*vsum[i][4] - vsum[i][0]*vsum[i][4]*vsum[i][4]
		  - vsum[i][2]*vsum[i][3]*vsum[i][2] - vsum[i][1]*vsum[i][1]*vsum[i][5];
			
		if(detA==0.0){
			cerr << "| #Error, detA=0.0 at leat-sqare" << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		cell.coeffLeastSquare1stCellVetexStencil.clear();
		cell.coeffLeastSquare1stCellVetexStencil.resize(6,0.0);
		
		cell.coeffLeastSquare1stCellVetexStencil[0] = 
			(vsum[i][3] * vsum[i][5] - vsum[i][4] * vsum[i][4]) / detA;    // inv_A(1,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[1] = 
			(vsum[i][2] * vsum[i][4] - vsum[i][1] * vsum[i][5]) / detA;    // inv_A(1,2) = (2,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[2] = 
			(vsum[i][1] * vsum[i][4] - vsum[i][2] * vsum[i][3]) / detA;    // inv_A(1,3) = (3,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[3] = 
			(vsum[i][0] * vsum[i][5] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
			
		cell.coeffLeastSquare1stCellVetexStencil[4] = 
			(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][4]) / detA;    // inv_A(2,3) = (3,2)
			
		cell.coeffLeastSquare1stCellVetexStencil[5] = 
			(vsum[i][0] * vsum[i][3] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
		
	}
	
	
}



void initLeastSquare2nd(
	SEMO_Mesh_Builder& mesh) {
		
		
	
	// double weight = 0.5;
	// double weight = 3.0;
	// double weight = 1.5;
	// double weight = 0.1;
	// double weight = 6.0;
	// double weight = 2.0;
	
	
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	
	SEMO_MPI_Builder mpi;
	
	
	vector<vector<double>> vsum(mesh.cells.size(),vector<double>(6,0.0));
	
	// internal cells
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		// cout << cell.stencil.size() << endl;
		
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
			
			double distX = cellSten.x - cell.x;
			double distY = cellSten.y - cell.y;
			double distZ = cellSten.z - cell.z;
			
			// double wk = 1.0;
			double wk = 1.0 / (
				pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			wk = pow(wk,weight);
			
			vsum[i][0] += wk * distX*distX;
			vsum[i][1] += wk * distX*distY;
			vsum[i][2] += wk * distX*distZ;
			vsum[i][3] += wk * distY*distY;
			vsum[i][4] += wk * distY*distZ;
			vsum[i][5] += wk * distZ*distZ;
			
		}
		
	}
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
				wk = pow(wk,weight);
				
				vsum[face.owner][0] += wk * distX*distX;
				vsum[face.owner][1] += wk * distX*distY;
				vsum[face.owner][2] += wk * distX*distZ;
				vsum[face.owner][3] += wk * distY*distY;
				vsum[face.owner][4] += wk * distY*distZ;
				vsum[face.owner][5] += wk * distZ*distZ;
					
	
			}
		}
	}
	
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				// auto& cell = mesh.cells[face.owner];
				
				// for(auto j : face.points){
					
					// auto& point = mesh.points[j]; 
					
					// double distX = point.x - cell.x;
					// double distY = point.y - cell.y;
					// double distZ = point.z - cell.z;
					
					// // double wk = 1.0;
					// double wk = 1.0 / (
						// pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
					// wk = pow(wk,weight);
					
					// vsum[face.owner][0] += wk * distX*distX;
					// vsum[face.owner][1] += wk * distX*distY;
					// vsum[face.owner][2] += wk * distX*distZ;
					// vsum[face.owner][3] += wk * distY*distY;
					// vsum[face.owner][4] += wk * distY*distZ;
					// vsum[face.owner][5] += wk * distZ*distZ;
					
				// }
	
			// }
		// }
	// }
	
	
	
		
	// processor faces
	if(size>1){

		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				++proc_num;
			}
		}
		
		vector<vector<double>> vsum_send(6,vector<double>(proc_num,0.0));
		proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// // double wk = 1.0;
				// double wk = 1.0 / (
					// pow(face.distCells[0],2.0)+
					// pow(face.distCells[1],2.0)+
					// pow(face.distCells[2],2.0));
				// wk = pow(wk,weight);
				
				// vsum[face.owner][0] += wk * face.distCells[0]*face.distCells[0];
				// vsum[face.owner][1] += wk * face.distCells[0]*face.distCells[1];
				// vsum[face.owner][2] += wk * face.distCells[0]*face.distCells[2];
				// vsum[face.owner][3] += wk * face.distCells[1]*face.distCells[1];
				// vsum[face.owner][4] += wk * face.distCells[1]*face.distCells[2];
				// vsum[face.owner][5] += wk * face.distCells[2]*face.distCells[2];
				
				for(auto j : face.stencil){
					auto& cellSten = mesh.cells[j];
						
					double distX = cellSten.x - (mesh.cells[face.owner].x + face.distCells[0]);
					double distY = cellSten.y - (mesh.cells[face.owner].y + face.distCells[1]);
					double distZ = cellSten.z - (mesh.cells[face.owner].z + face.distCells[2]);
					
					// double wk = 1.0;
					double wk = 1.0 / (
						pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
					wk = pow(wk,weight);
					
					vsum_send[0][proc_num] += wk * distX*distX;
					vsum_send[1][proc_num] += wk * distX*distY;
					vsum_send[2][proc_num] += wk * distX*distZ;
					vsum_send[3][proc_num] += wk * distY*distY;
					vsum_send[4][proc_num] += wk * distY*distZ;
					vsum_send[5][proc_num] += wk * distZ*distZ;
				}
				
				++proc_num;
			}
		}
		
		vector<vector<double>> vsum_recv(6,vector<double>(proc_num,0.0));
		for(int j=0; j<6; ++j){
			mpi.setProcsFaceDatas(
						vsum_send[j], vsum_recv[j],
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
		}
		vsum_send.clear();
		
		
		proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				for(int j=0; j<6; ++j){
					vsum[face.owner][j] += vsum_recv[j][proc_num];
				}
				++proc_num;
			}
		}
		
	}
	
	
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double detA = 
			vsum[i][0]*vsum[i][3]*vsum[i][5] + vsum[i][1]*vsum[i][4]*vsum[i][2]
		  + vsum[i][2]*vsum[i][1]*vsum[i][4] - vsum[i][0]*vsum[i][4]*vsum[i][4]
		  - vsum[i][2]*vsum[i][3]*vsum[i][2] - vsum[i][1]*vsum[i][1]*vsum[i][5];
			
		if(detA==0.0){
			cerr << "| #Error, detA=0.0 at leat-sqare" << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		cell.coeffLeastSquare1stCellVetexStencil.clear();
		cell.coeffLeastSquare1stCellVetexStencil.resize(6,0.0);
		
		cell.coeffLeastSquare1stCellVetexStencil[0] = 
			(vsum[i][3] * vsum[i][5] - vsum[i][4] * vsum[i][4]) / detA;    // inv_A(1,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[1] = 
			(vsum[i][2] * vsum[i][4] - vsum[i][1] * vsum[i][5]) / detA;    // inv_A(1,2) = (2,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[2] = 
			(vsum[i][1] * vsum[i][4] - vsum[i][2] * vsum[i][3]) / detA;    // inv_A(1,3) = (3,1)
			
		cell.coeffLeastSquare1stCellVetexStencil[3] = 
			(vsum[i][0] * vsum[i][5] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
			
		cell.coeffLeastSquare1stCellVetexStencil[4] = 
			(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][4]) / detA;    // inv_A(2,3) = (3,2)
			
		cell.coeffLeastSquare1stCellVetexStencil[5] = 
			(vsum[i][0] * vsum[i][3] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
		
	}
	
	
}


void calcLeastSquare1st(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	
	// double weight = 0.5;
	// double weight = 3.0;
	// double weight = 1.5;
	// double weight = 0.1;
	// double weight = 6.0;
	// double weight = 2.0;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
	
	
	for(auto& face : mesh.faces){
		
		// double wk = 1.0;
		double wk = 1.0 / sqrt(
			pow(face.distCells[0],2.0)+
			pow(face.distCells[1],2.0)+
			pow(face.distCells[2],2.0));
		wk = pow(wk,weight);
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double DVar = 
				mesh.cells[face.neighbour].var[cn] - mesh.cells[face.owner].var[cn];
			
			grad_tmp[face.owner][0] += wk * face.distCells[0] * DVar;
			grad_tmp[face.owner][1] += wk * face.distCells[1] * DVar;
			grad_tmp[face.owner][2] += wk * face.distCells[2] * DVar;
			
			grad_tmp[face.neighbour][0] += wk * face.distCells[0] * DVar;
			grad_tmp[face.neighbour][1] += wk * face.distCells[1] * DVar;
			grad_tmp[face.neighbour][2] += wk * face.distCells[2] * DVar;
			
		}
		else{
			
			double DVar = 0.0;
				// face.varR[fn] - mesh.cells[face.owner].var[cn];
			DVar += gradient[face.owner][0]*(face.x-mesh.cells[face.owner].x);
			DVar += gradient[face.owner][1]*(face.y-mesh.cells[face.owner].y);
			DVar += gradient[face.owner][2]*(face.z-mesh.cells[face.owner].z);
			
			grad_tmp[face.owner][0] += wk * face.distCells[0] * DVar;
			grad_tmp[face.owner][1] += wk * face.distCells[1] * DVar;
			grad_tmp[face.owner][2] += wk * face.distCells[2] * DVar;
			
		}
		
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare1stCellVetexStencil[0] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[3] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[5] * grad_tmp[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
}






void calcLeastSquare2nd(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	
	// double weight = 0.5;
	// double weight = 3.0;
	// double weight = 1.5;
	// double weight = 0.1;
	// double weight = 6.0;
	// double weight = 2.0;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
	
	// internal cells
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
				
			double distX = cellSten.x - cell.x;
			double distY = cellSten.y - cell.y;
			double distZ = cellSten.z - cell.z;
			
			// double wk = 1.0;
			double wk = 1.0 / (
				pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			wk = pow(wk,weight);
				
			double DVar = cellSten.var[cn] - cell.var[cn];
			
			grad_tmp[i][0] += wk * distX * DVar;
			grad_tmp[i][1] += wk * distY * DVar;
			grad_tmp[i][2] += wk * distZ * DVar;
			
		}
		
	}
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
				
				double DVar = 0.0;//face.varR[fn] - cell.var[cn];
				// if(boundary.type[cn] == "fixedValue"){
					// varF = boundary.var[cn];
					// // varF = 0.5*varF;
				// }
				// else if(boundary.type[cn] == "zeroGradient"){
					DVar += gradient[face.owner][0]*(face.x-cell.x);
					DVar += gradient[face.owner][1]*(face.y-cell.y);
					DVar += gradient[face.owner][2]*(face.z-cell.z);
				// }
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(distX,2.0)+
					pow(distY,2.0)+
					pow(distZ,2.0));
				wk = pow(wk,weight);
					
				grad_tmp[face.owner][0] += wk * distX * DVar;
				grad_tmp[face.owner][1] += wk * distY * DVar;
				grad_tmp[face.owner][2] += wk * distZ * DVar;
					
			}
			
		}
	}
	

		
	// processor faces
	if(size>1){

		int proc_num=0;
		vector<double> var_send, var_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				var_send.push_back(mesh.cells[face.owner].var[cn]);
				++proc_num;
			}
		}
		mpi.setProcsFaceDatas(
					var_send, var_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		var_send.clear();
		
		vector<double> gradientX_send(proc_num,0.0);
		vector<double> gradientY_send(proc_num,0.0);
		vector<double> gradientZ_send(proc_num,0.0);
		proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// // double wk = 1.0;
				// double wk = 1.0 / (
					// pow(face.distCells[0],2.0)+
					// pow(face.distCells[1],2.0)+
					// pow(face.distCells[2],2.0));
				// wk = pow(wk,weight);
				
				// double DVar = face.varR[fn] - mesh.cells[face.owner].var[cn];
				
				// gradient[face.owner][0] += wk * face.distCells[0] * DVar;
				// gradient[face.owner][1] += wk * face.distCells[1] * DVar;
				// gradient[face.owner][2] += wk * face.distCells[2] * DVar;

				for(auto j : face.stencil){
					auto& cellSten = mesh.cells[j];
						

					double DVar = cellSten.var[cn] - var_recv[proc_num];
						
					double distX = cellSten.x - (mesh.cells[face.owner].x + face.distCells[0]);
					double distY = cellSten.y - (mesh.cells[face.owner].y + face.distCells[1]);
					double distZ = cellSten.z - (mesh.cells[face.owner].z + face.distCells[2]);
					
					// double wk = 1.0;
					double wk = 1.0 / (
						pow(distX,2.0)+
						pow(distY,2.0)+
						pow(distZ,2.0));
					wk = pow(wk,weight);
						
					gradientX_send[proc_num] += wk * distX * DVar;
					gradientY_send[proc_num] += wk * distY * DVar;
					gradientZ_send[proc_num] += wk * distZ * DVar;
				}
				
				++proc_num;
			}
		}
		
		vector<double> gradientX_recv;
		vector<double> gradientY_recv;
		vector<double> gradientZ_recv;
		mpi.setProcsFaceDatas(
					gradientX_send, gradientX_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradientY_send, gradientY_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					gradientZ_send, gradientZ_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		var_send.clear();
		
		
		proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				grad_tmp[face.owner][0] += gradientX_recv[proc_num];
				grad_tmp[face.owner][1] += gradientY_recv[proc_num];
				grad_tmp[face.owner][2] += gradientZ_recv[proc_num];
				
				++proc_num;
			}
		}
		
	}
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare1stCellVetexStencil[0] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[3] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][0] +
			cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][1] +
			cell.coeffLeastSquare1stCellVetexStencil[5] * grad_tmp[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
}





void calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double wCL = face.wC;
			// double wCL = 0.5;
			double wCR = 1.0 - wCL;
			
			double varF = wCL*mesh.cells[face.owner].var[cn] + wCR*mesh.cells[face.neighbour].var[cn];
			
			varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
			varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
			varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
			varF += wCR*gradient[face.neighbour][0]*face.vecSkewness[0];
			varF += wCR*gradient[face.neighbour][1]*face.vecSkewness[1];
			varF += wCR*gradient[face.neighbour][2]*face.vecSkewness[2];
			
			for(int j=0; j<3; ++j){
				grad_tmp[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				grad_tmp[face.neighbour][j] -= 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			}
		}
	}
	
	
	
	
	// boundary face's nodes
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
				
				double varF = cell.var[cn];
				// if(boundary.type[cn] == "fixedValue"){
					// varF = boundary.var[cn];
					// // varF = 0.5*varF;
				// }
				// else if(boundary.type[cn] == "zeroGradient"){
					varF += gradient[face.owner][0]*(face.x-cell.x);
					varF += gradient[face.owner][1]*(face.y-cell.y);
					varF += gradient[face.owner][2]*(face.z-cell.z);
				// }
				
				for(int j=0; j<3; ++j){
					grad_tmp[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
					
			}
			
		}
	}
	
	
	// // processor faces
	// if(size>1){
		// vector<double> phi_send, phi_recv;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// phi_send.push_back(phi[face.owner]);
			// }
		// }
		
		// SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// phi_send, phi_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
					
		// int proc_num=0;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double wCL = face.wC;
			// // double wCL = 0.5;
			// double wCR = 1.0 - wCL;
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double varF = wCL*phi[face.owner]+wCR*phi_recv[proc_num];
				
				// // varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
				// // varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
				// // varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
				// // varF += wCR*gradient[face.neighbour][0]*face.vecSkewness[0];
				// // varF += wCR*gradient[face.neighbour][1]*face.vecSkewness[1];
				// // varF += wCR*gradient[face.neighbour][2]*face.vecSkewness[2];
				
				// for(int j=0; j<3; ++j){
					// grad_tmp[face.owner][j] += 
						// varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// }
				
				
				// ++proc_num;
			// }
		// }
	// }
	
	
	
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		gradient[i][0] = grad_tmp[i][0];
		gradient[i][1] = grad_tmp[i][1];
		gradient[i][2] = grad_tmp[i][2];
		
	}
	
	
}














void saveInitialmkdirs(char *dir_path){
	char buff[1024];
	char *p_dir = buff;
	
	strcpy(buff, dir_path);
	buff[1024-1] = '\0';
	
	while(*p_dir){
		if('/'==*p_dir){
			*p_dir = '\0';
			mkdir(buff,0777);
			*p_dir = '/';
		}
		p_dir ++;
	}
	
}






void saveInitialField(string folder){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	saveInitialmkdirs(folder_name);

	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << folder_name << "plot...) ... ";
	}
	
	outputFile.open(filenamePlot);
	
	outputFile.precision(20);
	
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	

	// Field data
	outputFile << "    <FieldData>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"totalError\" NumberOfTuples=\"1\" format=\"ascii\">" << endl;
	double totalError = 0.0;
	for(auto& cell : mesh.cells){
		totalError += abs( cell.var[controls.U] - cell.var[controls.UDV[0]] )*cell.volume;
		totalError += abs( cell.var[controls.V] - cell.var[controls.UDV[1]] )*cell.volume;
		totalError += abs( cell.var[controls.W] - cell.var[controls.UDV[2]] )*cell.volume;
	}
	outputFile << totalError << endl;
	outputFile << "     </DataArray>" << endl;
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"pointLevels\" format=\"ascii\">" << endl;
	for(auto& point : mesh.points) outputFile << point.level << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "    </PointData>" << endl;
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[0] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[controls.U] << " " << cell.var[controls.V] << " " << cell.var[controls.W] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.T] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.VF[0]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"massFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.MF[0]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"grad-pressure\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[controls.UDV[0]] << " " << cell.var[controls.UDV[1]] << " " << cell.var[controls.UDV[2]] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellLevels\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.level << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellGroups\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.group << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	// mesh errors
	outputFile << "     <DataArray type=\"Float64\" Name=\"skewness\" format=\"ascii\">" << endl;
	vector<double> skewness(mesh.cells.size(),0.0);
	for(auto& face : mesh.faces){
		double value1=0.0;
		double value2=0.0;
		value1 += face.vecSkewness[0]*face.vecSkewness[0];
		value1 += face.vecSkewness[1]*face.vecSkewness[1];
		value1 += face.vecSkewness[2]*face.vecSkewness[2];
		value2 += face.vecPN[0]*face.vecPN[0];
		value2 += face.vecPN[1]*face.vecPN[1];
		value2 += face.vecPN[2]*face.vecPN[2];
		double value = sqrt(value1)/sqrt(value2);
		skewness[face.owner] = max(skewness[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			skewness[face.neighbour] = max(skewness[face.neighbour],value);
		}
	}
	for(int i=0; i<mesh.cells.size(); ++i){
		outputFile << skewness[i] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"nonOrthogonality\" format=\"ascii\">" << endl;
	vector<double> nonOrthogonality(mesh.cells.size(),0.0);
	for(auto& face : mesh.faces){
		double value1=0.0;
		double value2=0.0;
		value1 += face.vecPN[0]*face.unitNormals[0];
		value1 += face.vecPN[1]*face.unitNormals[1];
		value1 += face.vecPN[2]*face.unitNormals[2];
		value2 += face.vecPN[0]*face.vecPN[0];
		value2 += face.vecPN[1]*face.vecPN[1];
		value2 += face.vecPN[2]*face.vecPN[2];
		double value = acos(abs(value1)/sqrt(value2));
		nonOrthogonality[face.owner] = max(nonOrthogonality[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			nonOrthogonality[face.neighbour] = max(nonOrthogonality[face.neighbour],value);
		}
	}
	for(int i=0; i<mesh.cells.size(); ++i){
		outputFile << nonOrthogonality[i] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"uniformity\" format=\"ascii\">" << endl;
	vector<double> uniformity(mesh.cells.size(),1.0);
	for(auto& face : mesh.faces){
		double value1=0.0;
		double value2=0.0;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			value1 += (face.x-mesh.cells[face.owner].x)*(face.x-mesh.cells[face.owner].x);
			value1 += (face.y-mesh.cells[face.owner].y)*(face.y-mesh.cells[face.owner].y);
			value1 += (face.z-mesh.cells[face.owner].z)*(face.z-mesh.cells[face.owner].z);
			value1 = sqrt(value1);
			value2 += (face.x-mesh.cells[face.neighbour].x)*(face.x-mesh.cells[face.neighbour].x);
			value2 += (face.y-mesh.cells[face.neighbour].y)*(face.y-mesh.cells[face.neighbour].y);
			value2 += (face.z-mesh.cells[face.neighbour].z)*(face.z-mesh.cells[face.neighbour].z);
			value2 = sqrt(value2);
			value2 = value1 + value2;
		}
		else{
			value1 += (face.x-mesh.cells[face.owner].x)*(face.x-mesh.cells[face.owner].x);
			value1 += (face.y-mesh.cells[face.owner].y)*(face.y-mesh.cells[face.owner].y);
			value1 += (face.z-mesh.cells[face.owner].z)*(face.z-mesh.cells[face.owner].z);
			value1 = sqrt(value1);
			value2 = value1*2.0;
		}
		double value = value1/value2;
		value = 1.0-2.0*abs(value-0.5);
		uniformity[face.owner] = min(uniformity[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			uniformity[face.neighbour] = min(uniformity[face.neighbour],value);
		}
	}
	for(int i=0; i<mesh.cells.size(); ++i){
		outputFile << uniformity[i] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	

	
	outputFile << "     <DataArray type=\"Float64\" Name=\"grad-error\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[controls.U] - cell.var[controls.UDV[0]] << " " << cell.var[controls.V] - cell.var[controls.UDV[1]] << " " << cell.var[controls.W] - cell.var[controls.UDV[2]] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "    </CellData>" << endl;
	
	
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.points){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto& cell : mesh.cells){
		outputFile << cell.faces.size() << endl;
		for(auto& i : cell.faces){
			outputFile << mesh.faces[i].points.size() << " ";
			for(auto& j : mesh.faces[i].points){
				outputFile << j << " ";
			}
			outputFile << endl;
		}
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.faces.size();
		for(auto& i : cell.faces){
			numbering += mesh.faces[i].points.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	

	// additional informations
	outputFile << " <owner>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.owner << " ";
	}
	outputFile << endl;
	outputFile << " </owner>" << endl;
	
	outputFile << " <neighbour>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.neighbour << " ";
	}
	outputFile << endl;
	outputFile << " </neighbour>" << endl;
	
	// outputFile << " <faceLevels>" << endl;
	// for(auto& face : mesh.faces){
		// outputFile << face.level << " ";
	// }
	// outputFile << endl;
	// outputFile << " </faceLevels>" << endl;
	
	outputFile << " <faceGroups>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.group << " ";
	}
	outputFile << endl;
	outputFile << " </faceGroups>" << endl;
	
	
	outputFile << " <bcName>" << endl;
	for(auto& boundary : mesh.boundary){
		// cout << boundary.name << endl;
		// read.trim;
		string bcName = boundary.name;
		
		bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
		bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		
		
		// outputFile << boundary.name << " ";
		outputFile << bcName << " ";
	}
	outputFile << endl;
	outputFile << " </bcName>" << endl;
	
	outputFile << " <bcStartFace>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.startFace << " ";
	}
	outputFile << endl;
	outputFile << " </bcStartFace>" << endl;
	
	outputFile << " <bcNFaces>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.nFaces << " ";
	}
	outputFile << endl;
	outputFile << " </bcNFaces>" << endl;
	
	outputFile << " <bcNeighbProcNo>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.neighbProcNo << " ";
	}
	outputFile << endl;
	outputFile << " </bcNeighbProcNo>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	

	if(rank==0){
		string filenamePvtu = folder + "../plot.0_Gradient.pvtu";
		
		outputFile.open(filenamePvtu);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// string out_line;
		outputFile << "<?xml version=\"1.0\"?>" << endl;
		outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		outputFile << "  <PUnstructuredGrid>" << endl;

		outputFile << "   <PFieldData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"totalError\"/>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0; ip<size; ++ip){
			string filenamevtus = "./0_Gradient/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"pointLevels\"/>" << endl;
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"massFraction-" << species[0].name << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"grad-pressure\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellLevels\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellGroups\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"skewness\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"nonOrthogonality\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"uniformity\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"grad-error\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}
