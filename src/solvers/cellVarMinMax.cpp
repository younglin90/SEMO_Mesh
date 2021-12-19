
#include "build.h"
#include <cmath>
#include <array>
#include <numeric>


void SEMO_Solvers_Builder::setCellVarMinMax(
	SEMO_Mesh_Builder& mesh,
	int cn, int cnMax, int cnMin){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		cell.var[cnMax] = cell.var[cn];
		cell.var[cnMin] = cell.var[cn];
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
			cell.var[cnMax] = max(cell.var[cnMax],cellSten.var[cn]);
			cell.var[cnMin] = min(cell.var[cnMin],cellSten.var[cn]);
		}
	}
	// processor faces
	if(size>1){
		SEMO_MPI_Builder mpi;
		
		vector<double> sendMax, recvMax;
		vector<double> sendMin, recvMin;
		for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				double varMax = mesh.cells[face.owner].var[cn];
				double varMin = mesh.cells[face.owner].var[cn];
				for(auto j : face.stencil){
					auto& cellSten = mesh.cells[j];
					
					varMax = max(varMax,cellSten.var[cn]);
					varMin = min(varMin,cellSten.var[cn]);
				}
				sendMax.push_back(varMax);
				sendMin.push_back(varMin);
				
				++proc_num;
			}
		}
		mpi.setProcsFaceDatas(
					sendMax, recvMax,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					sendMin, recvMin,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		
		for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			auto& cell = mesh.cells[face.owner];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cell.var[cnMax] = max(cell.var[cnMax],recvMax[proc_num]);
				cell.var[cnMin] = min(cell.var[cnMin],recvMax[proc_num]);
				
				++proc_num;
			}
		}
		
	}
	
}