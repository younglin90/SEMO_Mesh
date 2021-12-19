#include <algorithm>
#include <cmath>
#include <limits>
#include "math.h"

using namespace std;


//=========================================
// 1st-stencil least-square
void SEMO_Utility_Math::initLeastSquare(
	SEMO_Mesh_Builder& mesh) {
	
	
	double weight = 1.0;
		
	
	
	// ===================================
	// face stencil, first & second derivative LS
	{
		int rank = MPI::COMM_WORLD.Get_rank(); 
		int size = MPI::COMM_WORLD.Get_size();
		
		vector<vector<double>> vsum(mesh.cells.size(),vector<double>(45,0.0));
		for(auto& face : mesh.faces){
			
			double distX = face.distCells[0];
			double distY = face.distCells[1];
			double distZ = face.distCells[2];
			
			// weighting function
			double wk = 1.0 / sqrt(pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			wk = pow(wk,weight);
			
			vector<double> vari(9,0.0);
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
			vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
			vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
			vari[8] = 0.5*distZ*distZ;
			
			// 대칭행렬
			for(int ii=0, num=0; ii<9; ++ii){
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						vsum[face.owner][num++] += wk * vari[ii]*vari[jj];
					}
				}
			}
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				vari[0] *= (-1.0); vari[1] *= (-1.0); vari[2] *= (-1.0);
				
				for(int ii=0, num=0; ii<9; ++ii){
					for(int jj=0; jj<9; ++jj){
						if(jj>=ii){
							vsum[face.neighbour][num++] += wk * vari[ii]*vari[jj];
						}
					}
				}
			}
			
		}
		
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
		
			// 1차 least-square 저장
			double detA = 
				vsum[i][0]*vsum[i][9]*vsum[i][17] + vsum[i][1]*vsum[i][10]*vsum[i][2]
			  + vsum[i][2]*vsum[i][1]*vsum[i][10] - vsum[i][0]*vsum[i][10]*vsum[i][10]
			  - vsum[i][2]*vsum[i][9]*vsum[i][2] - vsum[i][1]*vsum[i][1]*vsum[i][17];
				
			if(detA==0.0){
				cerr << "| #Error, detA=0.0 at leat-sqare" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			cell.coeffLeastSquare1stFaceStencil.clear();
			cell.coeffLeastSquare1stFaceStencil.resize(6,0.0);
			
			cell.coeffLeastSquare1stFaceStencil[0] = 
				(vsum[i][9] * vsum[i][17] - vsum[i][10] * vsum[i][10]) / detA;    // inv_A(1,1)
				
			cell.coeffLeastSquare1stFaceStencil[1] = 
				(vsum[i][2] * vsum[i][10] - vsum[i][1] * vsum[i][17]) / detA;    // inv_A(1,2) = (2,1)
				
			cell.coeffLeastSquare1stFaceStencil[2] = 
				(vsum[i][1] * vsum[i][10] - vsum[i][2] * vsum[i][9]) / detA;    // inv_A(1,3) = (3,1)
				
			cell.coeffLeastSquare1stFaceStencil[3] = 
				(vsum[i][0] * vsum[i][17] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
				
			cell.coeffLeastSquare1stFaceStencil[4] = 
				(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][10]) / detA;    // inv_A(2,3) = (3,2)
				
			cell.coeffLeastSquare1stFaceStencil[5] = 
				(vsum[i][0] * vsum[i][9] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
			
			
		
		
			// 2차 least-square 저장
			// 대칭행렬의 역행렬도 대칭행렬
			vector<vector<double>> U(9,vector<double>(9,0.0));
			vector<int> str(9,0);
			for(int ii=0, num=0; ii<9; ++ii){
				str[ii] = num;
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						U[ii][jj] = vsum[i][num++];
					}
					else{
						int newNum = str[jj] + ii - jj;
						// cout << ii << " " << jj << " " << newNum << endl;
						U[ii][jj] = vsum[i][newNum];
					}
					// cout << "A(" << ii+1 << "," << jj+1 << ")=" << U[ii][jj] << endl;
				}
			}
			
		
			// SVD 수도 역행렬 구하기
			vector<double> D(U[0].size());
			vector<vector<double>> V(U[0].size(), vector<double>(U[0].size(), 0.0));
			
			SEMO_Utility_Math::svdcmp(U, D, V);
			
			
			vector<vector<double>> invD(D.size(), vector<double>(D.size(), 0.0));
			for (int ii = 0; ii < D.size(); ++ii) {
				if (abs(D[ii]) > 1.e-4 ) {
					invD[ii][ii] = 1.0/D[ii];
				}
			}
	
			vector<vector<double>> transU;
			SEMO_Utility_Math::transpose(U, transU);
			vector<vector<double>> out;
			SEMO_Utility_Math::matmul(invD, transU, out);
			vector<vector<double>> pseudoInvMatrix;
			SEMO_Utility_Math::matmul(V, out, pseudoInvMatrix);
			
			// 9x9 대칭행렬이라 총 45개만 저장하면 됨
			cell.coeffLeastSquare2ndFaceStencil.clear();
			cell.coeffLeastSquare2ndFaceStencil.resize(45,0.0);
			
			// 대칭행렬 저장
			for(int ii=0, num=0; ii<9; ++ii){
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						// cout << pseudoInvMatrix[ii][jj] << endl;
						cell.coeffLeastSquare2ndFaceStencil[num++] = pseudoInvMatrix[ii][jj];
					}
				}
			}
		}
	}
	

	// ===================================
	// cell vertex stencil 2nd-order LS
	{
			
		int rank = MPI::COMM_WORLD.Get_rank(); 
		int size = MPI::COMM_WORLD.Get_size();
		
		SEMO_MPI_Builder mpi;
		
		vector<vector<double>> vsum(mesh.cells.size(),vector<double>(45,0.0));
		
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
				
				vector<double> vari(9,0.0);
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
				vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
				vari[8] = 0.5*distZ*distZ;
				
				// 대칭행렬
				for(int ii=0, num=0; ii<9; ++ii){
					for(int jj=0; jj<9; ++jj){
						if(jj>=ii){
							vsum[i][num++] += wk * vari[ii]*vari[jj];
						}
					}
				}
				
			}
			
		}
		// processor faces
		if(size>1){

			int proc_num=0;
			for(int i=0; i<mesh.faces.size(); ++i){
				if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
			}
			
			vector<vector<double>> vsum_send(45,vector<double>(proc_num,0.0));
			for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				auto& cell = mesh.cells[face.owner];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					
					for(auto j : face.stencil){
					// for(auto j : cell.stencil){
						auto& cellSten = mesh.cells[j];
							
						double distX = cellSten.x - (mesh.cells[face.owner].x + face.distCells[0]);
						double distY = cellSten.y - (mesh.cells[face.owner].y + face.distCells[1]);
						double distZ = cellSten.z - (mesh.cells[face.owner].z + face.distCells[2]);
						
						double wk = 1.0 / (pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
						wk = pow(wk,weight);
						
						vector<double> vari(9,0.0);
						vari[0] = distX; vari[1] = distY; vari[2] = distZ;
						vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
						vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
						vari[8] = 0.5*distZ*distZ;
						
						// 대칭행렬
						for(int ii=0, num=0; ii<9; ++ii){
							for(int jj=0; jj<9; ++jj){
								if(jj>=ii){
									vsum_send[num++][ip] += wk * vari[ii]*vari[jj];
								}
							}
						}
					}
					
					++ip;
				}
			}
			
			vector<vector<double>> vsum_recv(45,vector<double>());
			for(int j=0; j<45; ++j){
				mpi.setProcsFaceDatas(vsum_send[j], vsum_recv[j],
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
			}
			
			for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					for(int j=0; j<45; ++j) vsum[face.owner][j] += vsum_recv[j][ip];
					++ip;
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
						
					double distX = face.x - cell.x;
					double distY = face.y - cell.y;
					double distZ = face.z - cell.z;
					
					double wk = 1.0 / (pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
					wk = pow(wk,weight);
				
					vector<double> vari(9,0.0);
					vari[0] = distX; vari[1] = distY; vari[2] = distZ;
					vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
					vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
					vari[8] = 0.5*distZ*distZ;
					
					// 대칭행렬
					for(int ii=0, num=0; ii<9; ++ii){
						for(int jj=0; jj<9; ++jj){
							if(jj>=ii){
								vsum[face.owner][num++] += wk * vari[ii]*vari[jj];
							}
						}
					}
		
				}
			}
		}
		// 역행렬 계산 뒤 저장
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			
			// 1차 least-square 저장
			double detA = 
				vsum[i][0]*vsum[i][9]*vsum[i][17] + vsum[i][1]*vsum[i][10]*vsum[i][2]
			  + vsum[i][2]*vsum[i][1]*vsum[i][10] - vsum[i][0]*vsum[i][10]*vsum[i][10]
			  - vsum[i][2]*vsum[i][9]*vsum[i][2] - vsum[i][1]*vsum[i][1]*vsum[i][17];
				
			if(detA==0.0){
				cerr << "| #Error, detA=0.0 at leat-sqare" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			cell.coeffLeastSquare1stCellVetexStencil.clear();
			cell.coeffLeastSquare1stCellVetexStencil.resize(6,0.0);
			
			cell.coeffLeastSquare1stCellVetexStencil[0] = 
				(vsum[i][9] * vsum[i][17] - vsum[i][10] * vsum[i][10]) / detA;    // inv_A(1,1)
				
			cell.coeffLeastSquare1stCellVetexStencil[1] = 
				(vsum[i][2] * vsum[i][10] - vsum[i][1] * vsum[i][17]) / detA;    // inv_A(1,2) = (2,1)
				
			cell.coeffLeastSquare1stCellVetexStencil[2] = 
				(vsum[i][1] * vsum[i][10] - vsum[i][2] * vsum[i][9]) / detA;    // inv_A(1,3) = (3,1)
				
			cell.coeffLeastSquare1stCellVetexStencil[3] = 
				(vsum[i][0] * vsum[i][17] - vsum[i][2] * vsum[i][2]) / detA;    // inv_A(2,2)
				
			cell.coeffLeastSquare1stCellVetexStencil[4] = 
				(vsum[i][2] * vsum[i][1] - vsum[i][0] * vsum[i][10]) / detA;    // inv_A(2,3) = (3,2)
				
			cell.coeffLeastSquare1stCellVetexStencil[5] = 
				(vsum[i][0] * vsum[i][9] - vsum[i][1] * vsum[i][1]) / detA;    // inv_A(3,3)
			
			
			
		
			// 2차 least-square 저장
			// 대칭행렬의 역행렬도 대칭행렬
			vector<vector<double>> U(9,vector<double>(9,0.0));
			vector<int> str(9,0);
			for(int ii=0, num=0; ii<9; ++ii){
				str[ii] = num;
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						U[ii][jj] = vsum[i][num++];
					}
					else{
						int newNum = str[jj] + ii - jj;
						// cout << ii << " " << jj << " " << newNum << endl;
						U[ii][jj] = vsum[i][newNum];
					}
					// cout << "A(" << ii+1 << "," << jj+1 << ")=" << U[ii][jj] << endl;
				}
			}
			
		
			// SVD 수도 역행렬 구하기
			vector<double> D(U[0].size());
			vector<vector<double>> V(U[0].size(), vector<double>(U[0].size(), 0.0));
			
			SEMO_Utility_Math::svdcmp(U, D, V);
			
			
			vector<vector<double>> invD(D.size(), vector<double>(D.size(), 0.0));
			for (int ii = 0; ii < D.size(); ++ii) {
				if (abs(D[ii]) > 1.e-8 ) {
					invD[ii][ii] = 1.0/D[ii];
				}
			}
	
			vector<vector<double>> transU;
			SEMO_Utility_Math::transpose(U, transU);
			vector<vector<double>> out;
			SEMO_Utility_Math::matmul(invD, transU, out);
			vector<vector<double>> pseudoInvMatrix;
			SEMO_Utility_Math::matmul(V, out, pseudoInvMatrix);
			
			// 9x9 대칭행렬이라 총 45개만 저장하면 됨
			cell.coeffLeastSquare2ndCellVertexStencil.clear();
			cell.coeffLeastSquare2ndCellVertexStencil.resize(45,0.0);
			
			// 대칭행렬 저장
			for(int ii=0, num=0; ii<9; ++ii){
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						// cout << pseudoInvMatrix[ii][jj] << " " << pseudoInvMatrix[jj][ii] << endl;
						cell.coeffLeastSquare2ndCellVertexStencil[num++] = pseudoInvMatrix[ii][jj];
					}
				}
			}
		}
		
		
	}
		
		
	
	
	
}



void SEMO_Utility_Math::calcLeastSquare(
	SEMO_Mesh_Builder& mesh,
	string sStencil, string sOrder, string sCellorInp,
	int cn, int fn,
	vector<double>& phi,
	vector<vector<double>>& gradient
	) {
	
	
	double weight = 1.0;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>());
	
	if(sOrder=="1st"){
		for(int i=0; i<mesh.cells.size(); ++i){
			grad_tmp[i].resize(3,0.0);
		}
	}
	else if(sOrder=="2nd"){
		for(int i=0; i<mesh.cells.size(); ++i){
			grad_tmp[i].resize(9,0.0);
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| Error, sOrder not defined" << endl;
		cout << endl;
		cout << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	




	if(sStencil=="face"){
		

		// processor faces
		int proc_num = 0;
		vector<double> phi_recv;
		if(size>1){
			vector<double> phi_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					if(sCellorInp=="cell"){
						phi_send.push_back(mesh.cells[face.owner].var[cn]);
					}
					else if(sCellorInp=="input"){
						phi_send.push_back(phi[face.owner]);
					}
					else{
						cout << endl; cout << endl;
						cout << "| Error, sCellorInp not defined" << endl;
						cout << endl; cout << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					++proc_num;
				}
			}
			
			SEMO_MPI_Builder mpi;
			
			mpi.setProcsFaceDatas(
						phi_send, phi_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
		}
		
		// internal cells
		int ip=0;
		for(auto& face : mesh.faces){
			
			if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
	

			double distX = face.distCells[0];
			double distY = face.distCells[1];
			double distZ = face.distCells[2];
			
			// double wk = 1.0;
			double wk = 1.0 / sqrt(pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
			wk = pow(wk,weight);
			
			double DVar = 0.0;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				if(sCellorInp=="cell"){
					DVar = mesh.cells[face.neighbour].var[cn] - mesh.cells[face.owner].var[cn];
				}
				else if(sCellorInp=="input"){
					DVar = phi[face.neighbour] - phi[face.owner];
				}
				else{
					cout << endl; cout << endl;
					cout << "| Error, sCellorInp not defined" << endl; 
					cout << endl; cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				if(sCellorInp=="cell"){
					DVar = phi_recv[ip] - mesh.cells[face.owner].var[cn];
				}
				else if(sCellorInp=="input"){
					DVar = phi_recv[ip] - phi[face.owner];
				}
				else{
					cout << endl; cout << endl;
					cout << "| Error, sCellorInp not defined" << endl; 
					cout << endl; cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
		
			vector<double> vari(9,0.0);
			vari[0] = distX; vari[1] = distY; vari[2] = distZ;
			vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
			vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
			vari[8] = 0.5*distZ*distZ;
	
			if(sOrder=="1st"){
				for(int ii=0; ii<3; ++ii){
					grad_tmp[face.owner][ii] += wk * vari[ii] * DVar;
				}
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					vari[0] *= (-1.0); vari[1] *= (-1.0); vari[2] *= (-1.0);
					DVar *= (-1.0);
					for(int ii=0; ii<3; ++ii){
						grad_tmp[face.neighbour][ii] += wk * vari[ii] * DVar;
					}
				}
			}
			else if(sOrder=="2nd"){
				for(int ii=0; ii<9; ++ii){
					grad_tmp[face.owner][ii] += wk * vari[ii] * DVar;
				}
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					vari[0] *= (-1.0); vari[1] *= (-1.0); vari[2] *= (-1.0);
					DVar *= (-1.0);
					for(int ii=0; ii<9; ++ii){
						grad_tmp[face.neighbour][ii] += wk * vari[ii] * DVar;
					}
				}
			}
			else{
				cout << endl; cout << endl;
				cout << "| Error, sOrder not defined" << endl;
				cout << endl; cout << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
				
			if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
		}
		
		
	
	}
	else if(sStencil=="cellVertex"){

		// internal cells
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			
			for(auto j : cell.stencil){
				auto& cellSten = mesh.cells[j];
					
				double distX = cellSten.x - cell.x;
				double distY = cellSten.y - cell.y;
				double distZ = cellSten.z - cell.z;
				
				// double wk = 1.0;
				double wk = 1.0 / (pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
				wk = pow(wk,weight);
					
				double DVar = 0.0;
				if(sCellorInp=="cell"){
					DVar = cellSten.var[cn] - cell.var[cn];
				}
				else if(sCellorInp=="input"){
					DVar = phi[j] - phi[i];
				}
				else{
					cout << endl; cout << endl;
					cout << "| Error, sCellorInp not defined" << endl;
					cout << endl; cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				vector<double> vari(9,0.0);
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
				vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
				vari[8] = 0.5*distZ*distZ;

				if(sOrder=="1st"){
					for(int jj=0; jj<3; ++jj){
						grad_tmp[i][jj] += wk * vari[jj] * DVar;
					}
				}
				else if(sOrder=="2nd"){
					for(int jj=0; jj<9; ++jj){
						grad_tmp[i][jj] += wk * vari[jj] * DVar;
					}
				}
				else{
					cout << endl; cout << endl;
					cout << "| Error, sOrder not defined" << endl;
					cout << endl; cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
					if(sCellorInp=="cell"){
						var_send.push_back(mesh.cells[face.owner].var[cn]);
					}
					else if(sCellorInp=="input"){
						var_send.push_back(phi[face.owner]);
					}
					else{
						cout << endl; cout << endl;
						cout << "| Error, sCellorInp not defined" << endl;
						cout << endl; cout << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					++proc_num;
				}
			}
			mpi.setProcsFaceDatas(var_send, var_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			var_send.clear();
			

			vector<vector<double>> gradient_send, gradient_recv;
			if(sOrder=="1st"){
				gradient_send.resize(3);
				gradient_recv.resize(3);
				for(int i=0; i<3; ++i){
					gradient_send[i].resize(proc_num,0.0);
				}
			}
			else if(sOrder=="2nd"){
				gradient_send.resize(9);
				gradient_recv.resize(9);
				for(int i=0; i<9; ++i){
					gradient_send[i].resize(proc_num,0.0);
				}
			}
			else{
				cout << endl;
				cout << endl;
				cout << "| Error, sOrder not defined" << endl;
				cout << endl;
				cout << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					
					for(auto j : face.stencil){
						auto& cellSten = mesh.cells[j];
							
						double DVar = 0.0;
						if(sCellorInp=="cell"){
							DVar = cellSten.var[cn] - var_recv[ip];
						}
						else if(sCellorInp=="input"){
							DVar = phi[j] - var_recv[ip];
						}
						else{
							cout << endl;
							cout << endl;
							cout << "| Error, sCellorInp not defined" << endl;
							cout << endl;
							cout << endl;
							MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
						}
							
						double distX = cellSten.x - (mesh.cells[face.owner].x + face.distCells[0]);
						double distY = cellSten.y - (mesh.cells[face.owner].y + face.distCells[1]);
						double distZ = cellSten.z - (mesh.cells[face.owner].z + face.distCells[2]);
						
						// double wk = 1.0;
						double wk = 1.0 / (pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
						wk = pow(wk,weight);
							
						vector<double> vari(9,0.0);
						vari[0] = distX; vari[1] = distY; vari[2] = distZ;
						vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
						vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
						vari[8] = 0.5*distZ*distZ;

						if(sOrder=="1st"){
							for(int jj=0; jj<3; ++jj){
								gradient_send[jj][ip] += wk * vari[jj] * DVar;
							}
						}
						else if(sOrder=="2nd"){
							for(int jj=0; jj<9; ++jj){
								gradient_send[jj][ip] += wk * vari[jj] * DVar;
							}
						}
						else{
							cout << endl; cout << endl;
							cout << "| Error, sOrder not defined" << endl;
							cout << endl; cout << endl;
							MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
						}
				
					}
					
					++ip;
				}
			}
			
			if(sOrder=="1st"){
				for(int i=0; i<3; ++i){
					mpi.setProcsFaceDatas(gradient_send[i], gradient_recv[i],
								mesh.countsProcFaces, mesh.countsProcFaces, 
								mesh.displsProcFaces, mesh.displsProcFaces);
				}
				for(int i=0, ip=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					
					if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						for(int jj=0; jj<3; ++jj){
							grad_tmp[face.owner][jj] += gradient_recv[jj][ip];
						}
						++ip;
					}
				}
			}
			else if(sOrder=="2nd"){
				for(int i=0; i<9; ++i){
					mpi.setProcsFaceDatas(gradient_send[i], gradient_recv[i],
								mesh.countsProcFaces, mesh.countsProcFaces, 
								mesh.displsProcFaces, mesh.displsProcFaces);
				}
				for(int i=0, ip=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					
					if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						for(int jj=0; jj<9; ++jj){
							grad_tmp[face.owner][jj] += gradient_recv[jj][ip];
						}
						
						++ip;
					}
				}
			}
			else{
				cout << endl;
				cout << endl;
				cout << "| Error, sOrder not defined" << endl;
				cout << endl;
				cout << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
		}
		
		
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| Error, sStencil not defined" << endl;
		cout << endl;
		cout << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
				if(sCellorInp=="cell"){
					DVar = face.varR[fn] - cell.var[cn];
					if(boundary.type[cn] == "zeroGradient"){
						DVar += gradient[face.owner][0]*(face.x-cell.x);
						DVar += gradient[face.owner][1]*(face.y-cell.y);
						DVar += gradient[face.owner][2]*(face.z-cell.z);
					}
				}
				else if(sCellorInp=="input"){
					if(cn==-1){
						DVar += gradient[face.owner][0]*(face.x-cell.x);
						DVar += gradient[face.owner][1]*(face.y-cell.y);
						DVar += gradient[face.owner][2]*(face.z-cell.z);
					}
					else{
						if(boundary.type[cn] == "fixedValue"){
							DVar = boundary.var[cn] - phi[face.owner];
						}
						
						if(boundary.type[cn] == "zeroGradient"){
							DVar += gradient[face.owner][0]*(face.x-cell.x);
							DVar += gradient[face.owner][1]*(face.y-cell.y);
							DVar += gradient[face.owner][2]*(face.z-cell.z);
						}
					}
				}
				else{
					cout << endl;
					cout << endl;
					cout << "| Error, sCellorInp not defined" << endl;
					cout << endl;
					cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			
					
				double distX = face.x - cell.x;
				double distY = face.y - cell.y;
				double distZ = face.z - cell.z;
				
				// double wk = 1.0;
				double wk = 1.0 / (pow(distX,2.0)+pow(distY,2.0)+pow(distZ,2.0));
				wk = pow(wk,weight);
				
				vector<double> vari(9,0.0);
				vari[0] = distX; vari[1] = distY; vari[2] = distZ;
				vari[3] = 0.5*distX*distX; vari[4] = distX*distY; vari[5] = distX*distZ;
				vari[6] = 0.5*distY*distY; vari[7] = distY*distZ;
				vari[8] = 0.5*distZ*distZ;

				if(sOrder=="1st"){
					for(int ii=0; ii<3; ++ii){
						grad_tmp[face.owner][ii] += wk * vari[ii] * DVar;
					}
				}
				else if(sOrder=="2nd"){
					for(int ii=0; ii<9; ++ii){
						grad_tmp[face.owner][ii] += wk * vari[ii] * DVar;
					}
				
				}
				else{
					cout << endl;
					cout << endl;
					cout << "| Error, sOrder not defined" << endl;
					cout << endl;
					cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
					
			}
			
		}
	}
	
	

	
	if(sStencil=="face" && sOrder=="1st"){

		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];

			// 대칭행렬
			vector<int> str(3,0);
			for(int ii=0, num=0; ii<3; ++ii){
				str[ii] = num;
				double tmp0 = 0.0;
				for(int jj=0; jj<3; ++jj){
					if(jj>=ii){
						tmp0 += cell.coeffLeastSquare1stFaceStencil[num++] * grad_tmp[i][jj];
					}
					else{
						int newNum = str[jj] + ii - jj;
						tmp0 += cell.coeffLeastSquare1stFaceStencil[newNum] * grad_tmp[i][jj];
					}
				}
				gradient[i][ii] = tmp0;
			}
		}
			
	}
	else if(sStencil=="face" && sOrder=="2nd"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];

			// 대칭행렬
			vector<int> str(9,0);
			for(int ii=0, num=0; ii<9; ++ii){
				str[ii] = num;
				double tmp0 = 0.0;
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						tmp0 += cell.coeffLeastSquare2ndFaceStencil[num++] * grad_tmp[i][jj];
					}
					else{
						int newNum = str[jj] + ii - jj;
						tmp0 += cell.coeffLeastSquare2ndFaceStencil[newNum] * grad_tmp[i][jj];
					}
				}
				gradient[i][ii] = tmp0;
			}
		}
		
		
	}
	else if(sStencil=="cellVertex" && sOrder=="1st"){
		
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];

			// 대칭행렬
			vector<int> str(3,0);
			for(int ii=0, num=0; ii<3; ++ii){
				str[ii] = num;
				double tmp0 = 0.0;
				for(int jj=0; jj<3; ++jj){
					if(jj>=ii){
						tmp0 += cell.coeffLeastSquare1stCellVetexStencil[num++] * grad_tmp[i][jj];
					}
					else{
						int newNum = str[jj] + ii - jj;
						tmp0 += cell.coeffLeastSquare1stCellVetexStencil[newNum] * grad_tmp[i][jj];
					}
				}
				gradient[i][ii] = tmp0;
			}
		}
		
	}
	else if(sStencil=="cellVertex" && sOrder=="2nd"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];

			// 대칭행렬
			vector<int> str(9,0);
			for(int ii=0, num=0; ii<9; ++ii){
				str[ii] = num;
				double tmp0 = 0.0;
				for(int jj=0; jj<9; ++jj){
					if(jj>=ii){
						tmp0 += cell.coeffLeastSquare2ndCellVertexStencil[num++] * grad_tmp[i][jj];
					}
					else{
						int newNum = str[jj] + ii - jj;
						tmp0 += cell.coeffLeastSquare2ndCellVertexStencil[newNum] * grad_tmp[i][jj];
					}
				}
				gradient[i][ii] = tmp0;
			}
			
		}
		
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| Error, sStencil or sOrder not defined" << endl;
		cout << endl;
		cout << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	
	
	
}




// void SEMO_Utility_Math::calcLeastSquare2nd(
	// SEMO_Mesh_Builder& mesh,
	// int cn,
	// vector<double> phi,
	// vector<vector<double>>& gradient
	// ) {
	
	
	// // double weight = 0.5;
	// // double weight = 3.0;
	// // double weight = 1.5;
	// // double weight = 0.1;
	// // double weight = 6.0;
	// // double weight = 1.0;
	// double weight = 1.0;
	
	
	
	
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
		
		
	// // gradient.clear();
	// // gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
	
	
	// // internal cells
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// for(auto j : cell.stencil){
			// auto& cellSten = mesh.cells[j];

			// double distX = cellSten.x - cell.x;
			// double distY = cellSten.y - cell.y;
			// double distZ = cellSten.z - cell.z;
			
			// // double wk = 1.0;
			// double wk = 1.0 / (
				// pow(distX,2.0)+
				// pow(distY,2.0)+
				// pow(distZ,2.0));
			// wk = pow(wk,weight);
				
			// double DVar = phi[j] - phi[i];
			
			// grad_tmp[i][0] += wk * distX * DVar;
			// grad_tmp[i][1] += wk * distY * DVar;
			// grad_tmp[i][2] += wk * distZ * DVar;
			
		// }
		
	// }
	
	// // boundary face's nodes
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				// auto& cell = mesh.cells[face.owner];
				
				// double DVar = 0.0;
				// if(cn==-1){
					// DVar += gradient[face.owner][0]*(face.x-cell.x);
					// DVar += gradient[face.owner][1]*(face.y-cell.y);
					// DVar += gradient[face.owner][2]*(face.z-cell.z);
				// }
				// else{
					// if(boundary.type[cn] == "fixedValue"){
						// DVar = boundary.var[cn] - phi[face.owner];
						// // DVar = 0.0 - phi[face.owner];
					// }
					
					// if(boundary.type[cn] == "zeroGradient"){
						// DVar += gradient[face.owner][0]*(face.x-cell.x);
						// DVar += gradient[face.owner][1]*(face.y-cell.y);
						// DVar += gradient[face.owner][2]*(face.z-cell.z);
					// }
				// }
				
				// double distX = face.x - cell.x;
				// double distY = face.y - cell.y;
				// double distZ = face.z - cell.z;
				
				
				// // double wk = 1.0;
				// double wk = 1.0 / (
					// pow(distX,2.0)+
					// pow(distY,2.0)+
					// pow(distZ,2.0));
				// wk = pow(wk,weight);
				
				// grad_tmp[face.owner][0] += wk * distX * DVar;
				// grad_tmp[face.owner][1] += wk * distY * DVar;
				// grad_tmp[face.owner][2] += wk * distZ * DVar;
					
			// }
			
		// }
	// }
	
	
	// // for(auto& boundary : mesh.boundary){
		
		// // if(boundary.neighbProcNo == -1){
			
			// // int str = boundary.startFace;
			// // int end = str + boundary.nFaces;
			
			// // for(int i=str; i<end; ++i){
				// // auto& face = mesh.faces[i];
				// // auto& cell = mesh.cells[face.owner];
				
				// // double DVar = 0.0;
				// // if(boundary.type[cn] == "fixedValue"){
					// // DVar = 0.0 - phi[face.owner];
				// // }
				
				// // for(auto j : face.points){
					
					// // auto& point = mesh.points[j]; 
						
					// // double distX = point.x - cell.x;
					// // double distY = point.y - cell.y;
					// // double distZ = point.z - cell.z;
					
					// // // double wk = 1.0;
					// // double wk = 1.0 / (
						// // pow(distX,2.0)+
						// // pow(distY,2.0)+
						// // pow(distZ,2.0));
					// // wk = pow(wk,weight);
					
					
					// // gradient[face.owner][0] += wk * distX * DVar;
					// // gradient[face.owner][1] += wk * distY * DVar;
					// // gradient[face.owner][2] += wk * distZ * DVar;
					
				// // }
			// // }
			
		// // }
	// // }
	

	// // processor faces
	// if(size>1){

		// int proc_num=0;
		// vector<double> phi_send, phi_recv;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// phi_send.push_back(phi[face.owner]);
				// ++proc_num;
			// }
		// }
		
		
		// mpi.setProcsFaceDatas(
					// phi_send, phi_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// phi_send.clear();
				
		
		// vector<double> gradientX_send(proc_num,0.0);
		// vector<double> gradientY_send(proc_num,0.0);
		// vector<double> gradientZ_send(proc_num,0.0);
		// proc_num=0;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// auto& cell = mesh.cells[face.owner];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// // // double wk = 1.0;
				// // double wk = 1.0 / (
					// // pow(face.distCells[0],2.0)+
					// // pow(face.distCells[1],2.0)+
					// // pow(face.distCells[2],2.0));
				// // wk = pow(wk,weight);
				
				// // double DVar = phi_recv[proc_num] - phi[face.owner];
				
				// // grad_tmp[face.owner][0] += wk * face.distCells[0] * DVar;
				// // grad_tmp[face.owner][1] += wk * face.distCells[1] * DVar;
				// // grad_tmp[face.owner][2] += wk * face.distCells[2] * DVar;
				

				// for(auto j : face.stencil){
				// // for(auto j : cell.stencil){
					// auto& cellSten = mesh.cells[j];
						

					// double DVar = phi[j] - phi_recv[proc_num];
						
					// double distX = cellSten.x - (mesh.cells[face.owner].x + face.distCells[0]);
					// double distY = cellSten.y - (mesh.cells[face.owner].y + face.distCells[1]);
					// double distZ = cellSten.z - (mesh.cells[face.owner].z + face.distCells[2]);
					
					// // double wk = 1.0;
					// double wk = 1.0 / (
						// pow(distX,2.0)+
						// pow(distY,2.0)+
						// pow(distZ,2.0));
					// wk = pow(wk,weight);
						
					// gradientX_send[proc_num] += wk * distX * DVar;
					// gradientY_send[proc_num] += wk * distY * DVar;
					// gradientZ_send[proc_num] += wk * distZ * DVar;
				// }
				
				// ++proc_num;
			// }
		// }
		
		// vector<double> gradientX_recv;
		// vector<double> gradientY_recv;
		// vector<double> gradientZ_recv;
		// mpi.setProcsFaceDatas(
					// gradientX_send, gradientX_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradientY_send, gradientY_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradientZ_send, gradientZ_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		
		
		// proc_num=0;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// grad_tmp[face.owner][0] += gradientX_recv[proc_num];
				// grad_tmp[face.owner][1] += gradientY_recv[proc_num];
				// grad_tmp[face.owner][2] += gradientZ_recv[proc_num];
				
				// ++proc_num;
			// }
		// }
		
	// }
	
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// double tmp0 = 
			// cell.coeffLeastSquare1stCellVetexStencil[0] * grad_tmp[i][0] +
			// cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][1] +
			// cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][2];
			
		// double tmp1 = 
			// cell.coeffLeastSquare1stCellVetexStencil[1] * grad_tmp[i][0] +
			// cell.coeffLeastSquare1stCellVetexStencil[3] * grad_tmp[i][1] +
			// cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][2];
			
		// double tmp2 = 
			// cell.coeffLeastSquare1stCellVetexStencil[2] * grad_tmp[i][0] +
			// cell.coeffLeastSquare1stCellVetexStencil[4] * grad_tmp[i][1] +
			// cell.coeffLeastSquare1stCellVetexStencil[5] * grad_tmp[i][2];
			
		// gradient[i][0] = tmp0;
		// gradient[i][1] = tmp1;
		// gradient[i][2] = tmp2;
		
	// }
// }



//=========================================
// 1st-stencil gauss-green
void SEMO_Utility_Math::calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn, int fn,
	vector<vector<double>>& gradient
	) {
		
		
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	
	
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));
	
	// // internal faces
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double wCL = face.wC;
			// // double wCL = 0.5;
			// double wCR = 1.0 - wCL;
			
			// double varF = wCL*mesh.cells[face.owner].var[cn] + wCR*mesh.cells[face.neighbour].var[cn];
			
			// for(int j=0; j<3; ++j){
				// gradient[face.owner][j] += 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// gradient[face.neighbour][j] -= 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			// }
		// }
		// else{
			
			// // double varF = face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*face.varR[fn];
			// double varF = 0.5*mesh.cells[face.owner].var[cn]+0.5*face.varR[fn];
			// for(int j=0; j<3; ++j){
				// gradient[face.owner][j] += 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			// }
			
		// }
	// }
	

	// processor faces
	if(size>1){
		vector<double> var_send, var_recv;
		vector<double> gradientX_send, gradientX_recv;
		vector<double> gradientY_send, gradientY_recv;
		vector<double> gradientZ_send, gradientZ_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				var_send.push_back(mesh.cells[face.owner].var[cn]);
				gradientX_send.push_back(gradient[face.owner][0]);
				gradientY_send.push_back(gradient[face.owner][1]);
				gradientZ_send.push_back(gradient[face.owner][2]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					var_send, var_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
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
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double wCL = face.wC;
			double wCR = 1.0 - wCL;
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double varF = wCL*mesh.cells[face.owner].var[cn]+wCR*face.varR[fn];
				double varF = wCL*mesh.cells[face.owner].var[cn]+wCR*var_recv[proc_num];
				
				if(boolSkewnessCorrection){
					varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
					varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
					varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
					varF += wCR*gradientX_recv[proc_num]*face.vecSkewness[0];
					varF += wCR*gradientY_recv[proc_num]*face.vecSkewness[1];
					varF += wCR*gradientZ_recv[proc_num]*face.vecSkewness[2];
				}
				// if(face.unitNormals[2] != 0.0){
					// cout << face.unitNormals[2] << endl;
				// }
				
				for(int j=0; j<3; ++j){
					// if(j==2) varF = mesh.cells[face.owner].var[cn];
					grad_tmp[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				
				
				++proc_num;
			}
		}
	}
	
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double wCL = face.wC;
			// double wCL = 0.5;
			double wCR = 1.0 - wCL;
			
			double varF = wCL*mesh.cells[face.owner].var[cn] + wCR*mesh.cells[face.neighbour].var[cn];
			
			if(boolSkewnessCorrection){
				varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
				varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
				varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
				varF += wCR*gradient[face.neighbour][0]*face.vecSkewness[0];
				varF += wCR*gradient[face.neighbour][1]*face.vecSkewness[1];
				varF += wCR*gradient[face.neighbour][2]*face.vecSkewness[2];
			}
				
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
				if(boundary.type[cn] == "fixedValue"){
					varF = boundary.var[cn];
				}
				else if(boundary.type[cn] == "zeroGradient"){
					varF += gradient[face.owner][0]*(face.x-cell.x);
					varF += gradient[face.owner][1]*(face.y-cell.y);
					varF += gradient[face.owner][2]*(face.z-cell.z);
				}
				else{
					varF = 0.5*(cell.var[cn] + face.varR[fn]);
				}
				
				for(int j=0; j<3; ++j){
					grad_tmp[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
					
			}
			
		}
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		gradient[i][0] = grad_tmp[i][0];
		gradient[i][1] = grad_tmp[i][1];
		gradient[i][2] = grad_tmp[i][2];
	}
	
	
	
	
}



void SEMO_Utility_Math::calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double> phi,
	vector<vector<double>>& gradient
	) {
	
	
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));


	SEMO_MPI_Builder mpi;
	
	
	// processor faces
	if(size>1){
		vector<double> var_send, var_recv;
		vector<double> gradientX_send, gradientX_recv;
		vector<double> gradientY_send, gradientY_recv;
		vector<double> gradientZ_send, gradientZ_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				var_send.push_back(phi[face.owner]);
				gradientX_send.push_back(gradient[face.owner][0]);
				gradientY_send.push_back(gradient[face.owner][1]);
				gradientZ_send.push_back(gradient[face.owner][2]);
			}
		}
		
		
		mpi.setProcsFaceDatas(
					var_send, var_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
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
					
		int proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double wCL = face.wC;
			double wCR = 1.0 - wCL;
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				double varF = wCL*phi[face.owner]+wCR*var_recv[proc_num];
				
				
				if(boolSkewnessCorrection){
					varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
					varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
					varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
					varF += wCR*gradientX_recv[proc_num]*face.vecSkewness[0];
					varF += wCR*gradientY_recv[proc_num]*face.vecSkewness[1];
					varF += wCR*gradientZ_recv[proc_num]*face.vecSkewness[2];
				}
				
				for(int j=0; j<3; ++j){
					grad_tmp[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				
				
				++proc_num;
			}
		}
	}
	
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double wCL = face.wC;
			// double wCL = 0.5;
			double wCR = 1.0 - wCL;
			
			double varF = wCL*phi[face.owner] + wCR*phi[face.neighbour];
			
			if(boolSkewnessCorrection){
				varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
				varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
				varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
				varF += wCR*gradient[face.neighbour][0]*face.vecSkewness[0];
				varF += wCR*gradient[face.neighbour][1]*face.vecSkewness[1];
				varF += wCR*gradient[face.neighbour][2]*face.vecSkewness[2];
			}
			
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
				
				double varF = phi[face.owner];
				if(cn==-1){
					varF += gradient[face.owner][0]*(face.x-cell.x);
					varF += gradient[face.owner][1]*(face.y-cell.y);
					varF += gradient[face.owner][2]*(face.z-cell.z);
				}
				else{
					if(boundary.type[cn] == "fixedValue"){
						varF = boundary.var[cn];
						// varF = 0.5*varF;
					}
					else if(boundary.type[cn] == "zeroGradient"){
						varF += gradient[face.owner][0]*(face.x-cell.x);
						varF += gradient[face.owner][1]*(face.y-cell.y);
						varF += gradient[face.owner][2]*(face.z-cell.z);
					}
				}
				
				for(int j=0; j<3; ++j){
					grad_tmp[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
					
			}
			
		}
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		gradient[i][0] = grad_tmp[i][0];
		gradient[i][1] = grad_tmp[i][1];
		gradient[i][2] = grad_tmp[i][2];
	}
	
	
	
	
	
	
	
	
	
	
}






void SEMO_Utility_Math::calcGaussGreen(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double> phi0, vector<double> phi1, vector<double> phi2,
	vector<vector<double>>& gradient
	) {
	
	
	
	bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	vector<vector<double>> grad_tmp(mesh.cells.size(),vector<double>(3,0.0));


	SEMO_MPI_Builder mpi;
	
	
	// processor faces
	if(size>1){
		vector<double> phi0_send, phi0_recv;
		vector<double> phi1_send, phi1_recv;
		vector<double> phi2_send, phi2_recv;
		// vector<double> gradientX_send, gradientX_recv;
		// vector<double> gradientY_send, gradientY_recv;
		// vector<double> gradientZ_send, gradientZ_recv;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				phi0_send.push_back(phi0[face.owner]);
				phi1_send.push_back(phi1[face.owner]);
				phi2_send.push_back(phi2[face.owner]);
				
				// gradientX_send.push_back(gradient[face.owner][0]);
				// gradientY_send.push_back(gradient[face.owner][1]);
				// gradientZ_send.push_back(gradient[face.owner][2]);
			}
		}
		
		
		mpi.setProcsFaceDatas(
					phi0_send, phi0_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi1_send, phi1_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					phi2_send, phi2_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		// mpi.setProcsFaceDatas(
					// gradientX_send, gradientX_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradientY_send, gradientY_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// gradientZ_send, gradientZ_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
					
		for(int i=0, ip=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double wCL = face.wC;
			double wCR = 1.0 - wCL;
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				double var0F = wCL*phi0[face.owner]+wCR*phi0_recv[ip];
				double var1F = wCL*phi1[face.owner]+wCR*phi1_recv[ip];
				double var2F = wCL*phi2[face.owner]+wCR*phi2_recv[ip];
				
				
				// if(boolSkewnessCorrection){
					// varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
					// varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
					// varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
					// varF += wCR*gradientX_recv[proc_num]*face.vecSkewness[0];
					// varF += wCR*gradientY_recv[proc_num]*face.vecSkewness[1];
					// varF += wCR*gradientZ_recv[proc_num]*face.vecSkewness[2];
				// }
				
				grad_tmp[face.owner][0] += 
					var0F*face.unitNormals[0]*face.area/mesh.cells[face.owner].volume;
				grad_tmp[face.owner][1] += 
					var1F*face.unitNormals[1]*face.area/mesh.cells[face.owner].volume;
				grad_tmp[face.owner][2] += 
					var2F*face.unitNormals[2]*face.area/mesh.cells[face.owner].volume;
				
				++ip;
			}
		}
	}
	
	
	// internal faces
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double wCL = face.wC;
			// double wCL = 0.5;
			double wCR = 1.0 - wCL;
			
			double var0F = wCL*phi0[face.owner] + wCR*phi0[face.neighbour];
			double var1F = wCL*phi1[face.owner] + wCR*phi1[face.neighbour];
			double var2F = wCL*phi2[face.owner] + wCR*phi2[face.neighbour];
			
			// if(boolSkewnessCorrection){
				// varF += wCL*gradient[face.owner][0]*face.vecSkewness[0];
				// varF += wCL*gradient[face.owner][1]*face.vecSkewness[1];
				// varF += wCL*gradient[face.owner][2]*face.vecSkewness[2];
				// varF += wCR*gradient[face.neighbour][0]*face.vecSkewness[0];
				// varF += wCR*gradient[face.neighbour][1]*face.vecSkewness[1];
				// varF += wCR*gradient[face.neighbour][2]*face.vecSkewness[2];
			// }
			
			grad_tmp[face.owner][0] += 
				var0F*face.unitNormals[0]*face.area/mesh.cells[face.owner].volume;
			grad_tmp[face.owner][1] += 
				var1F*face.unitNormals[1]*face.area/mesh.cells[face.owner].volume;
			grad_tmp[face.owner][2] += 
				var2F*face.unitNormals[2]*face.area/mesh.cells[face.owner].volume;
				
			grad_tmp[face.neighbour][0] -= 
				var0F*face.unitNormals[0]*face.area/mesh.cells[face.neighbour].volume;
			grad_tmp[face.neighbour][1] -= 
				var1F*face.unitNormals[1]*face.area/mesh.cells[face.neighbour].volume;
			grad_tmp[face.neighbour][2] -= 
				var2F*face.unitNormals[2]*face.area/mesh.cells[face.neighbour].volume;
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
				
				double var0F = phi0[face.owner];
				double var1F = phi1[face.owner];
				double var2F = phi2[face.owner];
				
				// if(cn==-1){
					// varF += gradient[face.owner][0]*(face.x-cell.x);
					// varF += gradient[face.owner][1]*(face.y-cell.y);
					// varF += gradient[face.owner][2]*(face.z-cell.z);
				// }
				// else{
					// if(boundary.type[cn] == "fixedValue"){
						// varF = boundary.var[cn];
						// // varF = 0.5*varF;
					// }
					// else if(boundary.type[cn] == "zeroGradient"){
						// varF += gradient[face.owner][0]*(face.x-cell.x);
						// varF += gradient[face.owner][1]*(face.y-cell.y);
						// varF += gradient[face.owner][2]*(face.z-cell.z);
					// }
				// }
				
				grad_tmp[face.owner][0] += 
					var0F*face.unitNormals[0]*face.area/mesh.cells[face.owner].volume;
				grad_tmp[face.owner][1] += 
					var1F*face.unitNormals[1]*face.area/mesh.cells[face.owner].volume;
				grad_tmp[face.owner][2] += 
					var2F*face.unitNormals[2]*face.area/mesh.cells[face.owner].volume;
					
			}
			
		}
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		gradient[i][0] = grad_tmp[i][0];
		gradient[i][1] = grad_tmp[i][1];
		gradient[i][2] = grad_tmp[i][2];
	}
	
	
	
	
	
	
	
	
	
	
}




// //=========================================
// // 1st-stencil gauss-green + least-square
// void SEMO_Utility_Math::calcGGLSQ(
	// SEMO_Mesh_Builder& mesh,
	// int cn, int fn,
	// vector<vector<double>>& gradient
	// ) {
	
	
	
	// vector<double> maxXPF(mesh.cells.size(),0.0);
	// vector<double> maxS(mesh.cells.size(),0.0);
	// vector<vector<double>> gradientGG(mesh.cells.size(),vector<double>(3,0.0));
	// vector<vector<double>> gradientLSQ(mesh.cells.size(),vector<double>(3,0.0));
	
	// // internal faces
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double varF = 
				// face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*mesh.cells[face.neighbour].var[cn];
			
			// for(int j=0; j<3; ++j){
				// gradientGG[face.owner][j] += 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// gradientGG[face.neighbour][j] -= 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			// }
			
			
			
			// double xPF = sqrt(
				// pow(face.x - mesh.cells[face.owner].x,2.0)+
				// pow(face.y - mesh.cells[face.owner].y,2.0)+
				// pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// double xNF = sqrt(
				// pow(face.x - mesh.cells[face.neighbour].x,2.0)+
				// pow(face.y - mesh.cells[face.neighbour].y,2.0)+
				// pow(face.z - mesh.cells[face.neighbour].z,2.0));
				
			// maxXPF[face.neighbour] = max(maxXPF[face.neighbour],xNF);
			
			// maxS[face.owner] = max(maxS[face.owner],face.area);
			// maxS[face.neighbour] = max(maxS[face.neighbour],face.area);
			
			
		// }
		// else{
			
			// double varF = face.wC*mesh.cells[face.owner].var[cn]+(1.0-face.wC)*face.varR[fn];
			// for(int j=0; j<3; ++j){
				// gradientGG[face.owner][j] += 
					// varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
			// }
			
			
			// double xPF = sqrt(
				// pow(face.x - mesh.cells[face.owner].x,2.0)+
				// pow(face.y - mesh.cells[face.owner].y,2.0)+
				// pow(face.z - mesh.cells[face.owner].z,2.0));
				
			// maxXPF[face.owner] = max(maxXPF[face.owner],xPF);
			
			// maxS[face.owner] = max(maxS[face.owner],face.area);
			
		// }
	// }
	
	

	// for(auto& face : mesh.faces){
		
		// double wk = 1.0;
		// // double wk = 1.0 / sqrt(
			// // pow(face.distCells[0],2.0)+
			// // pow(face.distCells[1],2.0)+
			// // pow(face.distCells[2],2.0));
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double DVar = 
				// mesh.cells[face.neighbour].var[cn] - mesh.cells[face.owner].var[cn];
			
			// gradientLSQ[face.owner][0] += wk * face.distCells[0] * DVar;
			// gradientLSQ[face.owner][1] += wk * face.distCells[1] * DVar;
			// gradientLSQ[face.owner][2] += wk * face.distCells[2] * DVar;
			
			// gradientLSQ[face.neighbour][0] += wk * face.distCells[0] * DVar;
			// gradientLSQ[face.neighbour][1] += wk * face.distCells[1] * DVar;
			// gradientLSQ[face.neighbour][2] += wk * face.distCells[2] * DVar;
			
		// }
		// else{
			
			// double DVar = 
				// face.varR[fn] - mesh.cells[face.owner].var[cn];
			
			// gradientLSQ[face.owner][0] += face.distCells[0] * DVar;
			// gradientLSQ[face.owner][1] += face.distCells[1] * DVar;
			// gradientLSQ[face.owner][2] += face.distCells[2] * DVar;
			
		// }
		
	// }
	

	// gradient.clear();
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// double AR = 2.0*(maxXPF[i]*maxS[i])/cell.volume;
		// double beta0 = min(1.0, 2.0/AR);
		
		// double tmp0 = 
			// cell.coeffLeastSquare[0] * gradientLSQ[i][0] +
			// cell.coeffLeastSquare[1] * gradientLSQ[i][1] +
			// cell.coeffLeastSquare[2] * gradientLSQ[i][2];
			
		// double tmp1 = 
			// cell.coeffLeastSquare[1] * gradientLSQ[i][0] +
			// cell.coeffLeastSquare[3] * gradientLSQ[i][1] +
			// cell.coeffLeastSquare[4] * gradientLSQ[i][2];
			
		// double tmp2 = 
			// cell.coeffLeastSquare[2] * gradientLSQ[i][0] +
			// cell.coeffLeastSquare[4] * gradientLSQ[i][1] +
			// cell.coeffLeastSquare[5] * gradientLSQ[i][2];
			
		// vector<double> tmpVec;
		// tmpVec.push_back(beta0*tmp0+(1.0-beta0)*gradientGG[i][0]);
		// tmpVec.push_back(beta0*tmp1+(1.0-beta0)*gradientGG[i][1]);
		// tmpVec.push_back(beta0*tmp2+(1.0-beta0)*gradientGG[i][2]);
		
		// gradient.push_back(tmpVec);
		
	// }
	
	
	
// }







// //=========================================
// // 1st-stencil modified gauss-green
// void SEMO_Utility_Math::calcMGG(
	// SEMO_Mesh_Builder& mesh,
	// int cn, int fn,
	// int iterMax, double toler, 
	// vector<vector<double>>& gradient
	// ) {
	
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// vector<vector<double>> tmpGradient(mesh.cells.size(),vector<double>(3,0.0));
	// vector<double> tmpGradientX_recv, tmpGradientY_recv, tmpGradientZ_recv;

	// if(size>1){
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// tmpGradientX_recv.push_back(0.0);
				// tmpGradientY_recv.push_back(0.0);
				// tmpGradientZ_recv.push_back(0.0);
			// }
		// }
	// }
	
	// gradient.clear();
	// gradient.resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int iter=0; iter<iterMax; ++iter){

		// // internal faces
		// int proc_num=0;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// double area = face.area;
			// vector<double> nvec(3,0.0);
			// nvec[0] = face.unitNormals[0];
			// nvec[1] = face.unitNormals[1];
			// nvec[2] = face.unitNormals[2];
			
			// vector<double> distanceCells(3,0.0);
			// distanceCells[0] = face.distCells[0];
			// distanceCells[1] = face.distCells[1];
			// distanceCells[2] = face.distCells[2];
			
			// double dPN = sqrt(pow(distanceCells[0],2.0) + 
							  // pow(distanceCells[1],2.0) + 
							  // pow(distanceCells[2],2.0));
			// double Ef[3];
			// Ef[0] = distanceCells[0]/dPN;
			// Ef[1] = distanceCells[1]/dPN;
			// Ef[2] = distanceCells[2]/dPN;
			
			// double alpha = nvec[0]*Ef[0]+nvec[1]*Ef[1]+nvec[2]*Ef[2];
			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				// double delPhiFAlphaRf = 
					// alpha*(mesh.cells[face.neighbour].var[cn]-mesh.cells[face.owner].var[cn])/dPN;
				
				// double barDelPhiFNfMAlphaRf = 0.0;
				// for(int j=0; j<3; ++j){
					// double gradVarF = face.wC*tmpGradient[face.owner][j]+(1.0-face.wC)*tmpGradient[face.neighbour][j];
					// // double gradVarF = 0.5*tmpGradient[face.owner][j]+0.5*tmpGradient[face.neighbour][j];
					// barDelPhiFNfMAlphaRf += gradVarF*(nvec[j]-alpha*Ef[j]);
				// }
				
				// double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				// double distOwnertoFace[3];
				// distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				// distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				// distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				// double distNeighbourtoFace[3];
				// distNeighbourtoFace[0] = mesh.cells[face.neighbour].x-face.x;
				// distNeighbourtoFace[1] = mesh.cells[face.neighbour].y-face.y;
				// distNeighbourtoFace[2] = mesh.cells[face.neighbour].z-face.z;
				
				// for(int j=0; j<3; ++j){
					
					// gradient[face.owner][j] += 
						// dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
					// gradient[face.neighbour][j] += 
						// dPhidNF*distNeighbourtoFace[j]*area/mesh.cells[face.neighbour].volume;
						
				// }
			// }
			// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				// double delPhiFAlphaRf = 
					// alpha*(face.varR[fn]-mesh.cells[face.owner].var[cn])/dPN;
				
				// double barDelPhiFNfMAlphaRf = 0.0;
				
				// double gradVarF[3];
				// gradVarF[0] = face.wC*tmpGradient[face.owner][0]+(1.0-face.wC)*tmpGradientX_recv[proc_num];
				// gradVarF[1] = face.wC*tmpGradient[face.owner][1]+(1.0-face.wC)*tmpGradientY_recv[proc_num];
				// gradVarF[2] = face.wC*tmpGradient[face.owner][2]+(1.0-face.wC)*tmpGradientZ_recv[proc_num];
				// for(int j=0; j<3; ++j){
					// barDelPhiFNfMAlphaRf += gradVarF[j]*(nvec[j]-alpha*Ef[j]);
				// }
				
				// double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				// double distOwnertoFace[3];
				// distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				// distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				// distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				// for(int j=0; j<3; ++j){
					
					// gradient[face.owner][j] += 
						// dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
				// }
				
				// ++proc_num;
			// }
			// else{
				
				// double delPhiFAlphaRf = 
					// alpha*(face.varR[fn]-mesh.cells[face.owner].var[cn])/dPN;
				
				// double barDelPhiFNfMAlphaRf = 0.0;
				// for(int j=0; j<3; ++j){
					// double gradVarF = tmpGradient[face.owner][j];
					// barDelPhiFNfMAlphaRf += gradVarF*(nvec[j]-alpha*Ef[j]);
				// }
				
				// double dPhidNF = delPhiFAlphaRf + barDelPhiFNfMAlphaRf;
				
				// double distOwnertoFace[3];
				// distOwnertoFace[0] = face.x - mesh.cells[face.owner].x;
				// distOwnertoFace[1] = face.y - mesh.cells[face.owner].y;
				// distOwnertoFace[2] = face.z - mesh.cells[face.owner].z;
				
				// for(int j=0; j<3; ++j){
					
					// gradient[face.owner][j] += 
						// dPhidNF*distOwnertoFace[j]*area/mesh.cells[face.owner].volume;
						
				// }
				
			// }
		// }
		
		
		// double tolerancelReduced = 0.0;
		// double tolerance = 0.0;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// tolerance += pow(tmpGradient[i][0]-gradient[i][0],2.0);
			// tolerance += pow(tmpGradient[i][1]-gradient[i][1],2.0);
			// tolerance += pow(tmpGradient[i][2]-gradient[i][2],2.0);
		// }
		// MPI_Allreduce(&tolerance, &tolerancelReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// tolerancelReduced = sqrt(tolerancelReduced);
		
		// // cout << iter << iterMax << endl;
		
		// if(tolerancelReduced < toler) break;
		
		
		// if( iter != iterMax-1 ){
			
			// for(int i=0; i<mesh.cells.size(); ++i){
				// tmpGradient[i][0] = gradient[i][0];
				// tmpGradient[i][1] = gradient[i][1];
				// tmpGradient[i][2] = gradient[i][2];
				
				// gradient[i][0] = 0.0;
				// gradient[i][1] = 0.0;
				// gradient[i][2] = 0.0;
			// }
			
			// if(size>1){
				// // processor faces
				// vector<double> gradX_send, gradY_send, gradZ_send;
				// for(int i=0; i<mesh.faces.size(); ++i){
					// auto& face = mesh.faces[i];
					
					// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						// gradX_send.push_back(tmpGradient[face.owner][0]);
						// gradY_send.push_back(tmpGradient[face.owner][1]);
						// gradZ_send.push_back(tmpGradient[face.owner][2]);
					// }
				// }
				// SEMO_MPI_Builder mpi;
							
				// mpi.setProcsFaceDatas(
							// gradX_send, tmpGradientX_recv,
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
				// mpi.setProcsFaceDatas(
							// gradY_send, tmpGradientY_recv,
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
				// mpi.setProcsFaceDatas(
							// gradZ_send, tmpGradientZ_recv,
							// mesh.countsProcFaces, mesh.countsProcFaces, 
							// mesh.displsProcFaces, mesh.displsProcFaces);
				// gradX_send.clear();
				// gradY_send.clear();
				// gradZ_send.clear();
			// }
		// }
	// }

	
// }







// //=========================================
// // gauss-green boundary treatment
// void SEMO_Utility_Math::gradientGGBoundaryTreatment(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<vector<double>>& gradient
	// ) {
	
	
	// // zeroGradient boundary treatment
	// // for(int iterUDF=0; iterUDF<3; ++iterUDF){
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// if(
				// boundary.type[0] == "zeroGradient" ||
				// boundary.type[0] == "noSlip" ||
				// boundary.type[0] == "slip"
				// ){
					
					// // gauss-green
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						// auto& cell = mesh.cells[face.owner];
						
						// double PF = 0.0; 
						// PF += gradient[face.owner][0]*face.distCells[0];
						// PF += gradient[face.owner][1]*face.distCells[1];
						// PF += gradient[face.owner][2]*face.distCells[2];
						// for(int j=0; j<3; ++j){
							// gradient[face.owner][j] += 
								// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
						// }
						
					// }
					
				// }
			// }
		// }
	// // }
	
	
	
// }





// void SEMO_Utility_Math::gradientLSQBoundaryTreatment(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<vector<double>>& gradient
	// ) {
	
	
	// // zeroGradient boundary treatment
	// // for(int iterUDF=0; iterUDF<3; ++iterUDF){
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// if(
				// boundary.type[0] == "zeroGradient" ||
				// boundary.type[0] == "noSlip" ||
				// boundary.type[0] == "slip"
				// ){
					
					
					// // least-square
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
			
						// double wk = 1.0;
						// // double wk = 1.0 / sqrt(
							// // pow(face.distCells[0],2.0)+
							// // pow(face.distCells[1],2.0)+
							// // pow(face.distCells[2],2.0));
						
						// double resiGardVar = 
							// gradient[face.owner][0]*face.distCells[0] +
							// gradient[face.owner][1]*face.distCells[1] +
							// gradient[face.owner][2]*face.distCells[2];
						
						// double DVar = 
							// 0.0 - resiGardVar;
							
						// double tmpXvar = wk * face.distCells[0] * DVar;
						// double tmpYvar = wk * face.distCells[1] * DVar;
						// double tmpZvar = wk * face.distCells[2] * DVar;
						
						// gradient[face.owner][0] += (
							// mesh.cells[face.owner].coeffLeastSquare[0]*tmpXvar + 
							// mesh.cells[face.owner].coeffLeastSquare[1]*tmpYvar + 
							// mesh.cells[face.owner].coeffLeastSquare[2]*tmpZvar
							// );
						// gradient[face.owner][1] += (
							// mesh.cells[face.owner].coeffLeastSquare[1]*tmpXvar + 
							// mesh.cells[face.owner].coeffLeastSquare[3]*tmpYvar + 
							// mesh.cells[face.owner].coeffLeastSquare[4]*tmpZvar
							// );
						// gradient[face.owner][2] += (
							// mesh.cells[face.owner].coeffLeastSquare[2]*tmpXvar + 
							// mesh.cells[face.owner].coeffLeastSquare[4]*tmpYvar + 
							// mesh.cells[face.owner].coeffLeastSquare[5]*tmpZvar
							// );
						
					// }
				// }
			// }
		// }
	// // }
	
	
	
// }






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



