#include <algorithm>
#include <cmath>
#include <limits>
#include "math.h"

using namespace std;


//=========================================
// 1st-stencil least-square
void SEMO_Utility_Math::initLeastSquare(
	SEMO_Mesh_Builder& mesh) {
	
	vector<vector<double>> vsum(mesh.cells.size(),vector<double>(6,0.0));
	
	for(auto& face : mesh.faces){
		
		double wk = 1.0;
		// double wk = 1.0 / sqrt(
			// pow(face.distCells[0],2.0)+
			// pow(face.distCells[1],2.0)+
			// pow(face.distCells[2],2.0));
		
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
		
		cell.coeffLeastSquare.clear();
		cell.coeffLeastSquare.resize(6,0.0);
		
		cell.coeffLeastSquare[0] = 
			(vsum[i][3] * vsum[i][5] - vsum[i][4] * vsum[i][4]) / detA;    // inv_A(1,1)
			
		cell.coeffLeastSquare[1] = 
			(vsum[i][2] * vsum[i][4] - vsum[i][1] * vsum[i][5]) / detA;    // inv_A(1,2) = (2,1)
			
		cell.coeffLeastSquare[2] = 
			(vsum[i][1] * vsum[i][4] - vsum[i][2] * vsum[i][3]) / detA;    // inv_A(1,3) = (3,1)
			
		cell.coeffLeastSquare[3] = 
			(vsum[i][0] * vsum[i][5] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
			
		cell.coeffLeastSquare[4] = 
			(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][4]) / detA;    // inv_A(2,3) = (3,2)
			
		cell.coeffLeastSquare[5] = 
			(vsum[i][0] * vsum[i][3] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
		
	}
	
	
}

void SEMO_Utility_Math::calcLeastSquare(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	
	for(auto& face : mesh.faces){
		
		double wk = 1.0;
		// double wk = 1.0 / sqrt(
			// pow(face.distCells[0],2.0)+
			// pow(face.distCells[1],2.0)+
			// pow(face.distCells[2],2.0));
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double DVar = 
				mesh.cells[face.neighbour].var[cn] - mesh.cells[face.owner].var[cn];
			
			gradient[face.owner][0] += wk * face.distCells[0] * DVar;
			gradient[face.owner][1] += wk * face.distCells[1] * DVar;
			gradient[face.owner][2] += wk * face.distCells[2] * DVar;
			
			gradient[face.neighbour][0] += wk * face.distCells[0] * DVar;
			gradient[face.neighbour][1] += wk * face.distCells[1] * DVar;
			gradient[face.neighbour][2] += wk * face.distCells[2] * DVar;
			
		}
		else{
			
			double DVar = 
				face.varR[fn] - mesh.cells[face.owner].var[cn];
			
			gradient[face.owner][0] += wk * face.distCells[0] * DVar;
			gradient[face.owner][1] += wk * face.distCells[1] * DVar;
			gradient[face.owner][2] += wk * face.distCells[2] * DVar;
			
		}
		
	}
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare[0] * gradient[i][0] +
			cell.coeffLeastSquare[1] * gradient[i][1] +
			cell.coeffLeastSquare[2] * gradient[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare[1] * gradient[i][0] +
			cell.coeffLeastSquare[3] * gradient[i][1] +
			cell.coeffLeastSquare[4] * gradient[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare[2] * gradient[i][0] +
			cell.coeffLeastSquare[4] * gradient[i][1] +
			cell.coeffLeastSquare[5] * gradient[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
	
	
	
	
	
}


void SEMO_Utility_Math::calcLeastSquare(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double> phi,
	vector<vector<double>>& gradient
	) {
		
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	
	for(auto& face : mesh.faces){
		
		double wk = 1.0;
		// double wk = 1.0 / sqrt(
			// pow(face.distCells[0],2.0)+
			// pow(face.distCells[1],2.0)+
			// pow(face.distCells[2],2.0));
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double DVar = 
				phi[face.neighbour] - phi[face.owner];
			
			gradient[face.owner][0] += wk * face.distCells[0] * DVar;
			gradient[face.owner][1] += wk * face.distCells[1] * DVar;
			gradient[face.owner][2] += wk * face.distCells[2] * DVar;
			
			gradient[face.neighbour][0] += wk * face.distCells[0] * DVar;
			gradient[face.neighbour][1] += wk * face.distCells[1] * DVar;
			gradient[face.neighbour][2] += wk * face.distCells[2] * DVar;
			
		}
		// else{
			
			// double DVar = 
				// face.varR[fn] - phi[face.owner];
			
			// gradient[face.owner][0] += wk * face.distCells[0] * DVar;
			// gradient[face.owner][1] += wk * face.distCells[1] * DVar;
			// gradient[face.owner][2] += wk * face.distCells[2] * DVar;
			
		// }
		
	}
	
	
	

	// boundary face's nodes
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
				
				double DVar = 0.0;
				if(boundary.type[cn] == "fixedValue"){
					DVar = 0.0 - phi[face.owner];
				}
				
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				
				double wk = 1.0;
				// double wk = 1.0 / sqrt(
					// pow(distX,2.0)+
					// pow(distY,2.0)+
					// pow(distZ,2.0));
				
				gradient[face.owner][0] += wk * distX * DVar;
				gradient[face.owner][1] += wk * distY * DVar;
				gradient[face.owner][2] += wk * distZ * DVar;
					
			}
			
		}
	}
	
	
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(boundary.type[cn] == "fixedValue"){
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					// auto& cell = mesh.cells[face.owner];
					
					// double DVar = 0.0 - phi[face.owner];
					
					// for(auto j : face.points){
						
						// auto& point = mesh.points[j]; 
						
						// double wk = 1.0;
						// // double wk = 1.0 / sqrt(
							// // pow(face.distCells[0],2.0)+
							// // pow(face.distCells[1],2.0)+
							// // pow(face.distCells[2],2.0));
							
						// double distX = point.x - cell.x;
						// double distY = point.y - cell.y;
						// double distZ = point.z - cell.z;
						
						
						// gradient[face.owner][0] += wk * distX * DVar;
						// gradient[face.owner][1] += wk * distY * DVar;
						// gradient[face.owner][2] += wk * distZ * DVar;
						
					// }
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// auto& face = mesh.faces[i];
					// auto& cell = mesh.cells[face.owner];
					
					// double DVar = 0.0;
					
					// for(auto j : face.points){
						
						// auto& point = mesh.points[j]; 
						
						// double wk = 1.0;
						// // double wk = 1.0 / sqrt(
							// // pow(face.distCells[0],2.0)+
							// // pow(face.distCells[1],2.0)+
							// // pow(face.distCells[2],2.0));
							
						// double distX = point.x - cell.x;
						// double distY = point.y - cell.y;
						// double distZ = point.z - cell.z;
						
						
						// gradient[face.owner][0] += wk * distX * DVar;
						// gradient[face.owner][1] += wk * distY * DVar;
						// gradient[face.owner][2] += wk * distZ * DVar;
						
					// }
				// }
	
			// }
			
		// }
	// }
	
	
	// processor faces
	if(size>1){
		vector<double> phi_send, phi_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send.push_back(phi[face.owner]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send, phi_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				double wk = 1.0;
				// double wk = 1.0 / sqrt(
					// pow(face.distCells[0],2.0)+
					// pow(face.distCells[1],2.0)+
					// pow(face.distCells[2],2.0));
				
				double DVar = phi_recv[proc_num] - phi[face.owner];
				
				gradient[face.owner][0] += wk * face.distCells[0] * DVar;
				gradient[face.owner][1] += wk * face.distCells[1] * DVar;
				gradient[face.owner][2] += wk * face.distCells[2] * DVar;
				
				++proc_num;
			}
		}
	}
	
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare[0] * gradient[i][0] +
			cell.coeffLeastSquare[1] * gradient[i][1] +
			cell.coeffLeastSquare[2] * gradient[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare[1] * gradient[i][0] +
			cell.coeffLeastSquare[3] * gradient[i][1] +
			cell.coeffLeastSquare[4] * gradient[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare[2] * gradient[i][0] +
			cell.coeffLeastSquare[4] * gradient[i][1] +
			cell.coeffLeastSquare[5] * gradient[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
	
	
	
	
	
}







//=========================================
// 2nd-stencil least-square
void SEMO_Utility_Math::initLeastSquare2nd(
	SEMO_Mesh_Builder& mesh) {
		
		
	
	// double weight = 0.5;
	// double weight = 3.0;
	double weight = 1.5;
	// double weight = 0.1;
	
	
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
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
					
					// // double wk = 1.0;
					// double wk = 1.0 / (
						// pow(face.distCells[0],2.0)+
						// pow(face.distCells[1],2.0)+
						// pow(face.distCells[2],2.0));
					// wk = pow(wk,weight);
						
					// double distX = point.x - cell.x;
					// double distY = point.y - cell.y;
					// double distZ = point.z - cell.z;
					
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
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double wk = 1.0;
				double wk = 1.0 / (
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
		
		cell.coeffLeastSquare.clear();
		cell.coeffLeastSquare.resize(6,0.0);
		
		cell.coeffLeastSquare[0] = 
			(vsum[i][3] * vsum[i][5] - vsum[i][4] * vsum[i][4]) / detA;    // inv_A(1,1)
			
		cell.coeffLeastSquare[1] = 
			(vsum[i][2] * vsum[i][4] - vsum[i][1] * vsum[i][5]) / detA;    // inv_A(1,2) = (2,1)
			
		cell.coeffLeastSquare[2] = 
			(vsum[i][1] * vsum[i][4] - vsum[i][2] * vsum[i][3]) / detA;    // inv_A(1,3) = (3,1)
			
		cell.coeffLeastSquare[3] = 
			(vsum[i][0] * vsum[i][5] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
			
		cell.coeffLeastSquare[4] = 
			(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][4]) / detA;    // inv_A(2,3) = (3,2)
			
		cell.coeffLeastSquare[5] = 
			(vsum[i][0] * vsum[i][3] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
		
	}
	
	
}

void SEMO_Utility_Math::calcLeastSquare2nd(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	
	// double weight = 0.5;
	// double weight = 3.0;
	double weight = 1.5;
	// double weight = 0.1;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
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
			
			gradient[i][0] += wk * distX * DVar;
			gradient[i][1] += wk * distY * DVar;
			gradient[i][2] += wk * distZ * DVar;
			
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
				
				double DVar = face.varR[fn] - cell.var[cn];
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(distX,2.0)+
					pow(distY,2.0)+
					pow(distZ,2.0));
				wk = pow(wk,weight);
					
				gradient[face.owner][0] += wk * distX * DVar;
				gradient[face.owner][1] += wk * distY * DVar;
				gradient[face.owner][2] += wk * distZ * DVar;
					
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
				
				// double DVar = face.varR[fn] - cell.var[cn];
				
				// for(auto j : face.points){
					
					// auto& point = mesh.points[j]; 
						
					// double distX = point.x - cell.x;
					// double distY = point.y - cell.y;
					// double distZ = point.z - cell.z;
					
					// // double wk = 1.0;
					// double wk = 1.0 / (
						// pow(distX,2.0)+
						// pow(distY,2.0)+
						// pow(distZ,2.0));
					// wk = pow(wk,weight);
					
					
					// gradient[face.owner][0] += wk * distX * DVar;
					// gradient[face.owner][1] += wk * distY * DVar;
					// gradient[face.owner][2] += wk * distZ * DVar;
					
				// }
	
			// }
		// }
	// }
	
	
	// processor faces
	if(size>1){
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(face.distCells[0],2.0)+
					pow(face.distCells[1],2.0)+
					pow(face.distCells[2],2.0));
				wk = pow(wk,weight);
				
				double DVar = face.varR[fn] - mesh.cells[face.owner].var[cn];
				
				gradient[face.owner][0] += wk * face.distCells[0] * DVar;
				gradient[face.owner][1] += wk * face.distCells[1] * DVar;
				gradient[face.owner][2] += wk * face.distCells[2] * DVar;
			}
		}
	}
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare[0] * gradient[i][0] +
			cell.coeffLeastSquare[1] * gradient[i][1] +
			cell.coeffLeastSquare[2] * gradient[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare[1] * gradient[i][0] +
			cell.coeffLeastSquare[3] * gradient[i][1] +
			cell.coeffLeastSquare[4] * gradient[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare[2] * gradient[i][0] +
			cell.coeffLeastSquare[4] * gradient[i][1] +
			cell.coeffLeastSquare[5] * gradient[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
}




void SEMO_Utility_Math::calcLeastSquare2nd(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double> phi,
	vector<vector<double>>& gradient
	) {
	
	
	// double weight = 0.5;
	// double weight = 3.0;
	double weight = 1.5;
	// double weight = 0.1;
	
	
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
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
				pow(distX,2.0)+
				pow(distY,2.0)+
				pow(distZ,2.0));
			wk = pow(wk,weight);
				
			double DVar = phi[j] - phi[i];
			
			gradient[i][0] += wk * distX * DVar;
			gradient[i][1] += wk * distY * DVar;
			gradient[i][2] += wk * distZ * DVar;
			
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
				
				double DVar = 0.0;
				if(boundary.type[cn] == "fixedValue"){
					DVar = 0.0 - phi[face.owner];
				}
				
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(distX,2.0)+
					pow(distY,2.0)+
					pow(distZ,2.0));
				wk = pow(wk,weight);
				
				gradient[face.owner][0] += wk * distX * DVar;
				gradient[face.owner][1] += wk * distY * DVar;
				gradient[face.owner][2] += wk * distZ * DVar;
					
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
				
				// double DVar = 0.0;
				// if(boundary.type[cn] == "fixedValue"){
					// DVar = 0.0 - phi[face.owner];
				// }
				
				// for(auto j : face.points){
					
					// auto& point = mesh.points[j]; 
						
					// double distX = point.x - cell.x;
					// double distY = point.y - cell.y;
					// double distZ = point.z - cell.z;
					
					// // double wk = 1.0;
					// double wk = 1.0 / (
						// pow(distX,2.0)+
						// pow(distY,2.0)+
						// pow(distZ,2.0));
					// wk = pow(wk,weight);
					
					
					// gradient[face.owner][0] += wk * distX * DVar;
					// gradient[face.owner][1] += wk * distY * DVar;
					// gradient[face.owner][2] += wk * distZ * DVar;
					
				// }
			// }
			
		// }
	// }
	
	
	// processor faces
	if(size>1){
		vector<double> phi_send, phi_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send.push_back(phi[face.owner]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send, phi_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double wk = 1.0;
				double wk = 1.0 / (
					pow(face.distCells[0],2.0)+
					pow(face.distCells[1],2.0)+
					pow(face.distCells[2],2.0));
				wk = pow(wk,weight);
				
				double DVar = phi_recv[proc_num] - phi[face.owner];
				
				gradient[face.owner][0] += wk * face.distCells[0] * DVar;
				gradient[face.owner][1] += wk * face.distCells[1] * DVar;
				gradient[face.owner][2] += wk * face.distCells[2] * DVar;
				
				++proc_num;
			}
		}
	}
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double tmp0 = 
			cell.coeffLeastSquare[0] * gradient[i][0] +
			cell.coeffLeastSquare[1] * gradient[i][1] +
			cell.coeffLeastSquare[2] * gradient[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare[1] * gradient[i][0] +
			cell.coeffLeastSquare[3] * gradient[i][1] +
			cell.coeffLeastSquare[4] * gradient[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare[2] * gradient[i][0] +
			cell.coeffLeastSquare[4] * gradient[i][1] +
			cell.coeffLeastSquare[5] * gradient[i][2];
			
		gradient[i][0] = tmp0;
		gradient[i][1] = tmp1;
		gradient[i][2] = tmp2;
		
	}
}



//=========================================
// 1st-stencil gauss-green
void SEMO_Utility_Math::calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double varF = 
				face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*mesh.cells[face.neighbour].var[cn];
			
			for(int j=0; j<3; ++j){
				gradient[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				gradient[face.neighbour][j] -= 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			}
		}
		else{
			
			double varF = face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*face.varR[fn];
			for(int j=0; j<3; ++j){
				gradient[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			}
			
		}
	}
	
}



void SEMO_Utility_Math::calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double> phi,
	vector<vector<double>>& gradient
	) {
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double varF = 
				face.wC*phi[face.owner]+(1.0-face.wC)*phi[face.neighbour];
			
			for(int j=0; j<3; ++j){
				gradient[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				gradient[face.neighbour][j] -= 
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
				
				double varF = phi[face.owner];
				if(boundary.type[cn] == "fixedValue"){
					varF = 0.5*varF;
				}
				
				for(int j=0; j<3; ++j){
					gradient[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
					
			}
			
		}
	}
	
	
	// processor faces
	if(size>1){
		vector<double> phi_send, phi_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi_send.push_back(phi[face.owner]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					phi_send, phi_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				double varF = 
					face.wC*phi[face.owner]+(1.0-face.wC)*phi_recv[proc_num];
				
				for(int j=0; j<3; ++j){
					gradient[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				
				
				++proc_num;
			}
		}
	}
	
	
	
	
	
}




//=========================================
// 1st-stencil gauss-green + least-square
void SEMO_Utility_Math::calcGGLSQ(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
	
	
	
	vector<double> maxXPF(mesh.cells.size(),0.0);
	vector<double> maxS(mesh.cells.size(),0.0);
	vector<vector<double>> gradientGG(mesh.cells.size(),vector<double>(3,0.0));
	vector<vector<double>> gradientLSQ(mesh.cells.size(),vector<double>(3,0.0));
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double varF = 
				face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*mesh.cells[face.neighbour].var[cn];
			
			for(int j=0; j<3; ++j){
				gradientGG[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				gradientGG[face.neighbour][j] -= 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			}
			
			
			
			double xPF = sqrt(
				pow(face.x - mesh.cells[face.owner].x,2.0)+
				pow(face.y - mesh.cells[face.owner].y,2.0)+
				pow(face.z - mesh.cells[face.owner].z,2.0));
				
			maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			double xNF = sqrt(
				pow(face.x - mesh.cells[face.neighbour].x,2.0)+
				pow(face.y - mesh.cells[face.neighbour].y,2.0)+
				pow(face.z - mesh.cells[face.neighbour].z,2.0));
				
			maxXPF[face.neighbour] = max(maxXPF[face.neighbour],xNF);
			
			maxS[face.owner] = max(maxS[face.owner],face.area);
			maxS[face.neighbour] = max(maxS[face.neighbour],face.area);
			
			
		}
		else{
			
			double varF = face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*face.varR[fn];
			for(int j=0; j<3; ++j){
				gradientGG[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			}
			
			
			double xPF = sqrt(
				pow(face.x - mesh.cells[face.owner].x,2.0)+
				pow(face.y - mesh.cells[face.owner].y,2.0)+
				pow(face.z - mesh.cells[face.owner].z,2.0));
				
			maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			maxS[face.owner] = max(maxS[face.owner],face.area);
			
		}
	}
	
	

	for(auto& face : mesh.faces){
		
		double wk = 1.0;
		// double wk = 1.0 / sqrt(
			// pow(face.distCells[0],2.0)+
			// pow(face.distCells[1],2.0)+
			// pow(face.distCells[2],2.0));
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double DVar = 
				mesh.cells[face.neighbour].var[cn] - mesh.cells[face.owner].var[cn];
			
			gradientLSQ[face.owner][0] += wk * face.distCells[0] * DVar;
			gradientLSQ[face.owner][1] += wk * face.distCells[1] * DVar;
			gradientLSQ[face.owner][2] += wk * face.distCells[2] * DVar;
			
			gradientLSQ[face.neighbour][0] += wk * face.distCells[0] * DVar;
			gradientLSQ[face.neighbour][1] += wk * face.distCells[1] * DVar;
			gradientLSQ[face.neighbour][2] += wk * face.distCells[2] * DVar;
			
		}
		else{
			
			double DVar = 
				face.varR[fn] - mesh.cells[face.owner].var[cn];
			
			gradientLSQ[face.owner][0] += face.distCells[0] * DVar;
			gradientLSQ[face.owner][1] += face.distCells[1] * DVar;
			gradientLSQ[face.owner][2] += face.distCells[2] * DVar;
			
		}
		
	}
	

	gradient.clear();
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		double beta0 = min(1.0, 2.0/AR);
		
		double tmp0 = 
			cell.coeffLeastSquare[0] * gradientLSQ[i][0] +
			cell.coeffLeastSquare[1] * gradientLSQ[i][1] +
			cell.coeffLeastSquare[2] * gradientLSQ[i][2];
			
		double tmp1 = 
			cell.coeffLeastSquare[1] * gradientLSQ[i][0] +
			cell.coeffLeastSquare[3] * gradientLSQ[i][1] +
			cell.coeffLeastSquare[4] * gradientLSQ[i][2];
			
		double tmp2 = 
			cell.coeffLeastSquare[2] * gradientLSQ[i][0] +
			cell.coeffLeastSquare[4] * gradientLSQ[i][1] +
			cell.coeffLeastSquare[5] * gradientLSQ[i][2];
			
		vector<double> tmpVec;
		tmpVec.push_back(beta0*tmp0+(1.0-beta0)*gradientGG[i][0]);
		tmpVec.push_back(beta0*tmp1+(1.0-beta0)*gradientGG[i][1]);
		tmpVec.push_back(beta0*tmp2+(1.0-beta0)*gradientGG[i][2]);
		
		gradient.push_back(tmpVec);
		
	}
	
	
	
}







//=========================================
// 1st-stencil modified gauss-green
void SEMO_Utility_Math::calcMGG(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	int iterMax, double toler, 
	vector<vector<double>>& gradient
	) {
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	vector<vector<double>> tmpGradient(mesh.cells.size(),vector<double>(3,0.0));
	vector<double> tmpGradientX_recv, tmpGradientY_recv, tmpGradientZ_recv;

	if(size>1){
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				tmpGradientX_recv.push_back(0.0);
				tmpGradientY_recv.push_back(0.0);
				tmpGradientZ_recv.push_back(0.0);
			}
		}
	}
	
	gradient.clear();
	gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	for(int iter=0; iter<iterMax; ++iter){

		// internal faces
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double area = face.area;
			vector<double> nvec(3,0.0);
			nvec[0] = face.unitNormals[0];
			nvec[1] = face.unitNormals[1];
			nvec[2] = face.unitNormals[2];
			
			vector<double> distanceCells(3,0.0);
			distanceCells[0] = face.distCells[0];
			distanceCells[1] = face.distCells[1];
			distanceCells[2] = face.distCells[2];
			
			double dPN = sqrt(pow(distanceCells[0],2.0) + 
							  pow(distanceCells[1],2.0) + 
							  pow(distanceCells[2],2.0));
			double Ef[3];
			Ef[0] = distanceCells[0]/dPN;
			Ef[1] = distanceCells[1]/dPN;
			Ef[2] = distanceCells[2]/dPN;
			
			double alpha = nvec[0]*Ef[0]+nvec[1]*Ef[1]+nvec[2]*Ef[2];
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double delPhiFAlphaRf = 
					alpha*(mesh.cells[face.neighbour].var[cn]-mesh.cells[face.owner].var[cn])/dPN;
				
				double barDelPhiFNfMAlphaRf = 0.0;
				for(int j=0; j<3; ++j){
					double gradVarF = face.wC*tmpGradient[face.owner][j]+(1.0-face.wC)*tmpGradient[face.neighbour][j];
					// double gradVarF = 0.5*tmpGradient[face.owner][j]+0.5*tmpGradient[face.neighbour][j];
					barDelPhiFNfMAlphaRf += gradVarF*(nvec[j]-alpha*Ef[j]);
				}
				
				double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				double distOwnertoFace[3];
				distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				double distNeighbourtoFace[3];
				distNeighbourtoFace[0] = mesh.cells[face.neighbour].x-face.x;
				distNeighbourtoFace[1] = mesh.cells[face.neighbour].y-face.y;
				distNeighbourtoFace[2] = mesh.cells[face.neighbour].z-face.z;
				
				for(int j=0; j<3; ++j){
					
					gradient[face.owner][j] += 
						dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
					gradient[face.neighbour][j] += 
						dPhidNF*distNeighbourtoFace[j]*area/mesh.cells[face.neighbour].volume;
						
				}
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				double delPhiFAlphaRf = 
					alpha*(face.varR[fn]-mesh.cells[face.owner].var[cn])/dPN;
				
				double barDelPhiFNfMAlphaRf = 0.0;
				
				double gradVarF[3];
				gradVarF[0] = face.wC*tmpGradient[face.owner][0]+(1.0-face.wC)*tmpGradientX_recv[proc_num];
				gradVarF[1] = face.wC*tmpGradient[face.owner][1]+(1.0-face.wC)*tmpGradientY_recv[proc_num];
				gradVarF[2] = face.wC*tmpGradient[face.owner][2]+(1.0-face.wC)*tmpGradientZ_recv[proc_num];
				for(int j=0; j<3; ++j){
					barDelPhiFNfMAlphaRf += gradVarF[j]*(nvec[j]-alpha*Ef[j]);
				}
				
				double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				double distOwnertoFace[3];
				distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				for(int j=0; j<3; ++j){
					
					gradient[face.owner][j] += 
						dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
				}
				
				++proc_num;
			}
			else{
				
				double delPhiFAlphaRf = 
					alpha*(face.varR[fn]-mesh.cells[face.owner].var[cn])/dPN;
				
				double barDelPhiFNfMAlphaRf = 0.0;
				for(int j=0; j<3; ++j){
					double gradVarF = tmpGradient[face.owner][j];
					barDelPhiFNfMAlphaRf += gradVarF*(nvec[j]-alpha*Ef[j]);
				}
				
				double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				double distOwnertoFace[3];
				distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				for(int j=0; j<3; ++j){
					
					gradient[face.owner][j] += 
						dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
				}
				
			}
		}
		
		
		double tolerancelReduced = 0.0;
		double tolerance = 0.0;
		for(int i=0; i<mesh.cells.size(); ++i){
			tolerance += pow(tmpGradient[i][0]-gradient[i][0],2.0);
			tolerance += pow(tmpGradient[i][1]-gradient[i][1],2.0);
			tolerance += pow(tmpGradient[i][2]-gradient[i][2],2.0);
		}
		MPI_Allreduce(&tolerance, &tolerancelReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tolerancelReduced = sqrt(tolerancelReduced);
		
		// cout << iter << iterMax << endl;
		
		if(tolerancelReduced < toler) break;
		
		
		if( iter != iterMax-1 ){
			
			for(int i=0; i<mesh.cells.size(); ++i){
				tmpGradient[i][0] = gradient[i][0];
				tmpGradient[i][1] = gradient[i][1];
				tmpGradient[i][2] = gradient[i][2];
				
				gradient[i][0] = 0.0;
				gradient[i][1] = 0.0;
				gradient[i][2] = 0.0;
			}
			
			if(size>1){
				// processor faces
				vector<double> gradX_send, gradY_send, gradZ_send;
				for(int i=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					
					if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						gradX_send.push_back(tmpGradient[face.owner][0]);
						gradY_send.push_back(tmpGradient[face.owner][1]);
						gradZ_send.push_back(tmpGradient[face.owner][2]);
					}
				}
				SEMO_MPI_Builder mpi;
							
				mpi.setProcsFaceDatas(
							gradX_send, tmpGradientX_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							gradY_send, tmpGradientY_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							gradZ_send, tmpGradientZ_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				gradX_send.clear();
				gradY_send.clear();
				gradZ_send.clear();
			}
		}
	}

	
}







//=========================================
// gauss-green boundary treatment
void SEMO_Utility_Math::gradientGGBoundaryTreatment(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& gradient
	) {
	
	
	// zeroGradient boundary treatment
	// for(int iterUDF=0; iterUDF<3; ++iterUDF){
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				if(
				boundary.type[0] == "zeroGradient" ||
				boundary.type[0] == "noSlip" ||
				boundary.type[0] == "slip"
				){
					
					// gauss-green
					for(int i=str; i<end; ++i){
						auto& face = mesh.faces[i];
						auto& cell = mesh.cells[face.owner];
						
						double PF = 0.0; 
						PF += gradient[face.owner][0]*face.distCells[0];
						PF += gradient[face.owner][1]*face.distCells[1];
						PF += gradient[face.owner][2]*face.distCells[2];
						for(int j=0; j<3; ++j){
							gradient[face.owner][j] += 
								PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
						}
						
					}
					
				}
			}
		}
	// }
	
	
	
}





void SEMO_Utility_Math::gradientLSQBoundaryTreatment(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& gradient
	) {
	
	
	// zeroGradient boundary treatment
	// for(int iterUDF=0; iterUDF<3; ++iterUDF){
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				if(
				boundary.type[0] == "zeroGradient" ||
				boundary.type[0] == "noSlip" ||
				boundary.type[0] == "slip"
				){
					
					
					// least-square
					for(int i=str; i<end; ++i){
						auto& face = mesh.faces[i];
			
						double wk = 1.0;
						// double wk = 1.0 / sqrt(
							// pow(face.distCells[0],2.0)+
							// pow(face.distCells[1],2.0)+
							// pow(face.distCells[2],2.0));
						
						double resiGardVar = 
							gradient[face.owner][0]*face.distCells[0] +
							gradient[face.owner][1]*face.distCells[1] +
							gradient[face.owner][2]*face.distCells[2];
						
						double DVar = 
							0.0 - resiGardVar;
							
						double tmpXvar = wk * face.distCells[0] * DVar;
						double tmpYvar = wk * face.distCells[1] * DVar;
						double tmpZvar = wk * face.distCells[2] * DVar;
						
						gradient[face.owner][0] += (
							mesh.cells[face.owner].coeffLeastSquare[0]*tmpXvar + 
							mesh.cells[face.owner].coeffLeastSquare[1]*tmpYvar + 
							mesh.cells[face.owner].coeffLeastSquare[2]*tmpZvar
							);
						gradient[face.owner][1] += (
							mesh.cells[face.owner].coeffLeastSquare[1]*tmpXvar + 
							mesh.cells[face.owner].coeffLeastSquare[3]*tmpYvar + 
							mesh.cells[face.owner].coeffLeastSquare[4]*tmpZvar
							);
						gradient[face.owner][2] += (
							mesh.cells[face.owner].coeffLeastSquare[2]*tmpXvar + 
							mesh.cells[face.owner].coeffLeastSquare[4]*tmpYvar + 
							mesh.cells[face.owner].coeffLeastSquare[5]*tmpZvar
							);
						
					}
				}
			}
		}
	// }
	
	
	
}






//=========================================
// gradient limiter
void SEMO_Utility_Math::calcLimiterGradient(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient,
	vector<double>& limGrad
	) {
	
	// vector<double> maxPhi(mesh.cells.size(),-1.e100);
	// vector<double> minPhi(mesh.cells.size(),+1.e100);
	
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// double ownngbMax = max(mesh.cells[face.owner].var[cn],
                                   // mesh.cells[face.neighbour].var[cn]);
			// double ownngbMin = min(mesh.cells[face.owner].var[cn],
                                   // mesh.cells[face.neighbour].var[cn]);
			
			// maxPhi[face.owner] = max(maxPhi[face.owner],ownngbMax);
			// minPhi[face.owner] = min(minPhi[face.owner],ownngbMin);
			
			// maxPhi[face.neighbour] = max(maxPhi[face.neighbour],ownngbMax);
			// minPhi[face.neighbour] = min(minPhi[face.neighbour],ownngbMin);
			
		// }
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// maxPhi[face.owner] = max(maxPhi[face.owner],face.varR[fn]);
			// minPhi[face.owner] = min(minPhi[face.owner],face.varR[fn]);
		// }
	// }
	
	vector<double> maxPhi;
	vector<double> minPhi;
	for(auto& cell : mesh.cells){
		maxPhi.push_back(cell.var[cn]);
		minPhi.push_back(cell.var[cn]);
	}
	// internal cells
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
			maxPhi[i] = max(maxPhi[i],cellSten.var[cn]);
			minPhi[i] = min(minPhi[i],cellSten.var[cn]);
		}
	}
	// PROCESSOR FACE
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			maxPhi[face.owner] = max(maxPhi[face.owner],face.varR[fn]);
			minPhi[face.owner] = min(minPhi[face.owner],face.varR[fn]);
		}
	}
	
	
	limGrad.resize(mesh.cells.size(),1.0);
	for(auto& face : mesh.faces){
		
		int own = face.owner;
		int ngb = face.neighbour;
		
		double delPF = 
			(face.x - mesh.cells[own].x)*gradient[own][0]+
			(face.y - mesh.cells[own].y)*gradient[own][1]+
			(face.z - mesh.cells[own].z)*gradient[own][2];
		
		double rPF;
		if( delPF > 0.0 ){
			double maxDelP = maxPhi[own] - mesh.cells[own].var[cn];
			rPF = delPF/(maxDelP+1.e-200);
		}
		else{
			double minDelP = minPhi[own] - mesh.cells[own].var[cn];
			rPF = delPF/(minDelP+1.e-200);
		}
		
		// // venka
		// double alphaPF = (2.0*rPF+1.0)/(rPF*(2.0*rPF+1.0)+1.0);
		
		// min-mod
		// double alphaPF = min(1.0,1.0/(abs(rPF)+1.e-200));
		
		// m-venka
		double alphaPF = (rPF*(2.0*rPF+1.0)+1.0)/(rPF*(rPF*(2.0*rPF+1.0)+1.0)+1.0);
		
		limGrad[own] = min(limGrad[own],alphaPF);
			
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double delNF = 
				(face.x - mesh.cells[ngb].x)*gradient[ngb][0]+
				(face.y - mesh.cells[ngb].y)*gradient[ngb][1]+
				(face.z - mesh.cells[ngb].z)*gradient[ngb][2];
			
			double rNF;
			if( delNF > 0.0 ){
				double maxDelP = maxPhi[ngb] - mesh.cells[ngb].var[cn];
				rNF = delNF/(maxDelP+1.e-200);
			}
			else{
				double minDelP = minPhi[ngb] - mesh.cells[ngb].var[cn];
				rNF = delNF/(minDelP+1.e-200);
			}
			
			// // venka
			// double alphaNF = (2.0*rNF+1.0)/(rNF*(2.0*rNF+1.0)+1.0);
			
			// min-mod
			// double alphaNF = min(1.0,1.0/(abs(rNF)+1.e-200));
			
			// m-venka
			double alphaNF = (rNF*(2.0*rNF+1.0)+1.0)/(rNF*(rNF*(2.0*rNF+1.0)+1.0)+1.0);
			
			limGrad[ngb] = min(limGrad[ngb],alphaNF);
			
		}
		
	}
	
	
	
}










//=========================================
// face gradient 
void SEMO_Utility_Math::calcGradientFace(
	SEMO_Mesh_Builder& mesh,
	vector<vector<double>>& gradient,
	int cn, int fn,
	int inX, int inY, int inZ
	) {
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	// processor faces
	vector<double> gradientX_send;
	vector<double> gradientY_send;
	vector<double> gradientZ_send;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			gradientX_send.push_back(gradient[face.owner][0]);
			gradientY_send.push_back(gradient[face.owner][1]);
			gradientZ_send.push_back(gradient[face.owner][2]);
		}
	}
	
	vector<double> gradientX_recv;
	vector<double> gradientY_recv;
	vector<double> gradientZ_recv;
	if(size>1){
		SEMO_MPI_Builder mpi;
		
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
		gradientX_send.clear();
		gradientY_send.clear();
		gradientZ_send.clear();
	}
	
	
	
	// faces gradient
	// @ref Jalali, Alireza, Mahkame Sharbatdar, and Carl Ollivier-Gooch. 
	// "Accuracy analysis of unstructured finite volume discretization schemes 
	// for diffusive fluxes." Computers & Fluids 101 (2014): 220-232.
	int proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		vector<double> nvec;
		nvec.push_back(face.unitNormals[0]);
		nvec.push_back(face.unitNormals[1]);
		nvec.push_back(face.unitNormals[2]);
		
		vector<double> distanceCells;
		distanceCells.push_back(face.distCells[0]);
		distanceCells.push_back(face.distCells[1]);
		distanceCells.push_back(face.distCells[2]);
		
		double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];

		double dampCoeff = 4.0/3.0;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// face.var[inX] = 0.5*(gradient[face.owner][0]+gradient[face.neighbour][0]);
			// face.var[inY] = 0.5*(gradient[face.owner][1]+gradient[face.neighbour][1]);
			// face.var[inZ] = 0.5*(gradient[face.owner][2]+gradient[face.neighbour][2]);
			
			double phiL = mesh.cells[face.owner].var[cn];
			phiL += gradient[face.owner][0]*(face.x-mesh.cells[face.owner].x);
			phiL += gradient[face.owner][1]*(face.y-mesh.cells[face.owner].y);
			phiL += gradient[face.owner][2]*(face.z-mesh.cells[face.owner].z);
			
			double phiR = mesh.cells[face.neighbour].var[cn];
			phiR += gradient[face.neighbour][0]*(face.x-mesh.cells[face.neighbour].x);
			phiR += gradient[face.neighbour][1]*(face.y-mesh.cells[face.neighbour].y);
			phiR += gradient[face.neighbour][2]*(face.z-mesh.cells[face.neighbour].z);
			
			face.var[inX] = 0.5*( gradient[face.owner][0] + gradient[face.neighbour][0] );
			face.var[inX] += dampCoeff*(phiR-phiL)/dPN_e*nvec[0];
			
			face.var[inY] = 0.5*( gradient[face.owner][1] + gradient[face.neighbour][1] );
			face.var[inY] += dampCoeff*(phiR-phiL)/dPN_e*nvec[1];
			
			face.var[inZ] = 0.5*( gradient[face.owner][2] + gradient[face.neighbour][2] );
			face.var[inZ] += dampCoeff*(phiR-phiL)/dPN_e*nvec[1];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// face.var[inX] = 0.5*(gradient[face.owner][0]+gradientX_recv[proc_num]);
			// face.var[inY] = 0.5*(gradient[face.owner][1]+gradientY_recv[proc_num]);
			// face.var[inZ] = 0.5*(gradient[face.owner][2]+gradientZ_recv[proc_num]);
			
			double phiL = mesh.cells[face.owner].var[cn];
			phiL += gradient[face.owner][0]*(face.x-mesh.cells[face.owner].x);
			phiL += gradient[face.owner][1]*(face.y-mesh.cells[face.owner].y);
			phiL += gradient[face.owner][2]*(face.z-mesh.cells[face.owner].z);
			
			double phiR = face.varR[fn];
			phiR += gradientX_recv[proc_num]*(face.x-mesh.cells[face.owner].x)*face.wC/(1.0-face.wC);
			phiR += gradientY_recv[proc_num]*(face.y-mesh.cells[face.owner].y)*face.wC/(1.0-face.wC);
			phiR += gradientZ_recv[proc_num]*(face.z-mesh.cells[face.owner].z)*face.wC/(1.0-face.wC);
			
			face.var[inX] = 0.5*( gradient[face.owner][0] + gradientX_recv[proc_num] );
			face.var[inX] += dampCoeff*(phiR-phiL)/dPN_e*nvec[0];
			
			face.var[inY] = 0.5*( gradient[face.owner][1] + gradientY_recv[proc_num] );
			face.var[inY] += dampCoeff*(phiR-phiL)/dPN_e*nvec[1];
			
			face.var[inZ] = 0.5*( gradient[face.owner][2] + gradientZ_recv[proc_num] );
			face.var[inZ] += dampCoeff*(phiR-phiL)/dPN_e*nvec[1];
			
			++proc_num;
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			face.var[inX] = gradient[face.owner][0];
			face.var[inY] = gradient[face.owner][1];
			face.var[inZ] = gradient[face.owner][2];
			
			// vector<double> nvec(3,0.0);
			// nvec[0] = face.unitNormals[0];
			// nvec[1] = face.unitNormals[1];
			// nvec[2] = face.unitNormals[2];
		
			// vector<double> distanceCells(3,0.0);
			// distanceCells[0] = face.distCells[0];
			// distanceCells[1] = face.distCells[1];
			// distanceCells[2] = face.distCells[2];
			// double dPN = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
			
			// double varL = face.varL[fn];
			// double varR = face.varR[fn];
			// face.var[inX] = (varR-varL)/dPN_e*nvec[0];
			// face.var[inY] = (varR-varL)/dPN_e*nvec[1];
			// face.var[inZ] = (varR-varL)/dPN_e*nvec[2];
			
			
		}
	}
	

}



