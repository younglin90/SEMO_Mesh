#pragma once
#include <iostream>
#include <vector>
#include <array>
#include <list>
#include "mpi.h"
using namespace std;

#include "../mesh/build.h"


class SEMO_MPI_Builder{
	public:
		SEMO_MPI_Builder() {
		}
		void init(){
			MPI::Init(); 
		}
		void setRank(){
			rank = MPI::COMM_WORLD.Get_rank(); 
		}
		void setSize(){
			size = MPI::COMM_WORLD.Get_size();
		}
		int getRank() {
			return rank;
		}
		int getSize() {
			return size;
		}
		void setSendCounts(int in){
			sendCounts = new int[in];
		}
		void setRecvCounts(int in){
			recvCounts = new int[in];
		}
		void setSendDispls(int in){
			sendDispls = new int[in];
		}
		void setRecvDispls(int in){
			recvDispls = new int[in];
		}
		
		void exchangeDatas(
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendValues, vector<int>& recvValues);
		void exchangeDatas(
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<double>& sendValues, vector<double>& recvValues);
		void setCellDatasToFaceRight(
				SEMO_Mesh_Builder& mesh, 
				int cin, int fin,
				vector<int>& sendCounts, vector<int>& recvCounts, 
				vector<int>& sendDisps, vector<int>& recvDisps
				);
		
		
		void setProcsFaceDatas(
			vector<int>& sendValues, vector<int>& recvValues,
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendDisps, vector<int>& recvDisps);
		
		void setProcsFaceDatas(
			vector<double>& sendValues, vector<double>& recvValues,
			vector<int>& sendCounts, vector<int>& recvCounts, 
			vector<int>& sendDisps, vector<int>& recvDisps);
		
		// void setProcsFaceDatas(
			// vector<bool>& sendValues, vector<bool>& recvValues,
			// vector<int>& sendCounts, vector<int>& recvCounts, 
			// vector<int>& sendDisps, vector<int>& recvDisps);
		
		
		// void setProcsFaceDatasDouble(
			// vector<double>& sendValues, vector<double>& recvValues,
			// vector<int>& sendCounts, vector<int>& recvCounts, 
			// vector<int>& sendDisps, vector<int>& recvDisps);
		
		int* sendCounts;
		int* recvCounts;
		int* sendDispls;
		int* recvDispls;
		int nSend;
		int nRecv;
		
	private:
		int rank;
		int size;
		MPI::Status status;
		
};

