#pragma once
#include <iostream>
#include "mpi.h"



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

