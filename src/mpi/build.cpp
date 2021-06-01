#include "build.h"
#include <vector>
#include "mpi.h"
using namespace std;

void SEMO_MPI_Builder::exchangeDatas(
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendValues, vector<int>& recvValues){

	int size = sendCounts.size(); 
	
// cout << size<<endl;
	vector<int> sendDisps(size,0);
	vector<int> recvDisps(size,0);
		
	for(int i=1; i<size; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<size; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
}


void SEMO_MPI_Builder::exchangeDatas(
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<double>& sendValues, vector<double>& recvValues){
	int size = sendCounts.size(); 
	
	vector<int> sendDisps(size,0);
	vector<int> recvDisps(size,0);
		
	for(int i=1; i<size; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<size; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0.0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
	
}


void SEMO_MPI_Builder::setCellDatasToFaceRight(
	SEMO_Mesh_Builder& mesh, 
	int cin, int fin,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps
	){
		
	vector<double> sendValues;
	vector<double> recvValues;
		
	int size = sendCounts.size(); 
	
	sendValues.clear();
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			sendValues.push_back(mesh.cells[face.owner].var[cin]);
		}
	}
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0.0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
				   
	vector<double>::iterator iter = recvValues.begin();
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.varR[fin] = *iter;
			++iter;
		}
	}
	
}



void SEMO_MPI_Builder::setProcsFaceDatas(
	vector<int>& sendValues, vector<int>& recvValues,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	int size = sendCounts.size(); 
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	// int num=0;
	// for(auto& i : recvValues){
		
		// cout << sendValues[num] << " " << i << endl;
		// ++num;
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}



void SEMO_MPI_Builder::setProcsFaceDatasDouble(
	vector<double>& sendValues, vector<double>& recvValues,
	vector<int>& sendCounts, vector<int>& recvCounts, 
	vector<int>& sendDisps, vector<int>& recvDisps){
		
		
	int size = sendCounts.size(); 
	
	recvValues.clear();
	recvValues.resize(recvDisps[size-1] + recvCounts[size-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
	
}


