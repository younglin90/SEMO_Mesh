
#include "build.h"
#include <cmath>
#include <array>
#include "mpi.h"

void SEMO_Mesh_Builder::check(){
	
	int owner_max=-1;
	int neighbour_max=-1;
	int checkOwnerNeighbourReverse=0;
	for(auto& face : (*this).faces){
		owner_max = max(owner_max , face.owner);
		neighbour_max = max(neighbour_max , face.neighbour);
	}
	
	// if(owner_max < neighbour_max){
		// cerr << "#error, owner_max(" << owner_max << ") < neighbour_max(" << neighbour_max << ")" << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	
	for(auto& face : (*this).faces){
		if(face.neighbour > owner_max){
			// cerr << "#error face, owner = " << face.owner << " neighbour = " << face.neighbour << endl;
			int tempnum = face.owner;
			face.owner = face.neighbour;
			face.neighbour = tempnum;
			std::reverse(face.points.begin(),face.points.end());
			++checkOwnerNeighbourReverse;
		}
	}
	cout << "  #check face owner < neighbour, executed reverse : " << checkOwnerNeighbourReverse << endl;
	
	
}


void SEMO_Mesh_Builder::loadFile(string filetype){
	SEMO_Mesh_Load meshLoad;
	 
	if(filetype=="OpenFoam"){
		meshLoad.OpenFoam(*this);
	} 
	if(filetype=="vtu"){
		meshLoad.vtu(*this);
	} 
}



void SEMO_Mesh_Builder::saveFile(string filetype){
	SEMO_Mesh_Save meshSave;
	 
	if(filetype=="vtu"){
		meshSave.vtu(*this);
	} 
}



void SEMO_Mesh_Builder::buildCells(){

	for(auto& i : (*this).cells){
		i.points.clear();
		i.faces.clear();
	}
	(*this).cells.clear();

	// add Cells
	int cell_num=0;
	for(auto& face : (*this).faces){
		cell_num = max(cell_num , face.owner);
		cell_num = max(cell_num , face.neighbour);
		
	
	}
	if((*this).faces.size() > 3) {
		for(int i=0; i<cell_num+1; ++i){
			(*this).addCell();
		}
	}
	
	
	
}



void SEMO_Mesh_Builder::buildCells2(){

	// for(auto& i : (*this).cells){
		// i.points.clear();
		// i.faces.clear();
	// }
	// (*this).cells.clear();

	// // add Cells
	// int cell_num=0;
	// for(auto& face : (*this).faces){
		// // cell_num = max(cell_num , face.owner);
		// cell_num = max(cell_num , face.neighbour);
	// }
		// cout << cell_num << endl;
	// for(int i=0; i<cell_num+1; ++i){
		// // (*this).addCell();
	// }
	
}


void SEMO_Mesh_Builder::setFaceTypes(){

	// set types
	int num_face = 0;
	for(auto& face : (*this).faces){
		face.setType(SEMO_Types::INTERNAL_FACE);
	}
	
	for(int ibc=0; ibc<(*this).boundary.size(); ++ibc){
		int startFace = (*this).boundary[ibc].startFace;
		for(int i=startFace; i<startFace+(*this).boundary[ibc].nFaces; ++i){
			if( (*this).boundary[ibc].myProcNo == -1 ){
				(*this).faces[i].setType(SEMO_Types::BOUNDARY_FACE);
				(*this).faces[i].setTypeBC(ibc);
			}
			else{
				(*this).faces[i].setType(SEMO_Types::PROCESSOR_FACE);
			}
		}
	}
	
}


void SEMO_Mesh_Builder::buildLists(){

	(*this).listPoints.clear();
	(*this).listCells.clear();
	(*this).listInternalFaces.clear();
	(*this).listBoundaryFaces.clear();
	(*this).listProcessorFaces.clear();
	
	for(auto iter=(*this).points.begin(); iter!=(*this).points.end(); iter++){
		(*this).listPoints.push_back(&*iter);
	}
	
	for(auto iter=(*this).cells.begin(); iter!=(*this).cells.end(); iter++){
		(*this).listCells.push_back(&*iter);
	}
	
	for(auto iter=(*this).faces.begin(); iter!=(*this).faces.end(); iter++){
		if(iter->getType() == SEMO_Types::INTERNAL_FACE){
			(*this).listInternalFaces.push_back(&*iter);
		}
		else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			(*this).listBoundaryFaces.push_back(&*iter);
		}
		else if(iter->getType() == SEMO_Types::PROCESSOR_FACE){
			(*this).listProcessorFaces.push_back(&*iter);
		}
		else{
			cout << "TO DO : other face types" << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	}
	
	
	  

}


void SEMO_Mesh_Builder::checkLists(){


	cout << "------------------------------------" << endl;
	cout << "point size : " << (*this).listPoints.size() << endl;
	cout << "cell size : " << (*this).listCells.size() << endl;
	cout << "internal face size : " << (*this).listInternalFaces.size() << endl;
	cout << "boundary face size : " << (*this).listBoundaryFaces.size() << endl;
	cout << "processor face size : " << (*this).listProcessorFaces.size() << endl;

}

void SEMO_Mesh_Builder::connectCelltoFaces(){

	// cell connection (cell's face)
	cout << "------------------------------------" << endl;
	cout << "| execute cell's face connecting ... ";
	int faceNum=0;
	for(auto iter=(*this).faces.begin(); iter!=(*this).faces.end(); iter++){
		if(iter->getType() == SEMO_Types::INTERNAL_FACE){
			(*this).cells[iter->owner].faces.push_back(faceNum);
			(*this).cells[iter->neighbour].faces.push_back(faceNum);
		}
		else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			(*this).cells[iter->owner].faces.push_back(faceNum);
		}
		else if(iter->getType() == SEMO_Types::PROCESSOR_FACE){
			(*this).cells[iter->owner].faces.push_back(faceNum);
		}
		++faceNum;
	}
	cout << "-> completed" << endl;
	

}




void SEMO_Mesh_Builder::connectCelltoPoints(){


	// cell connection (cell's points)
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
	cout << "| execute cell's points connecting ... ";
	for(auto& iter : (*this).faces){
		if(iter.getType() == SEMO_Types::INTERNAL_FACE){
			for(int i=0; i<iter.points.size(); ++i){
				int l=0;
				for(int j=0; j<(*this).cells[iter.owner].points.size(); ++j){
					if(
					(*this).cells[iter.owner].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) (*this).cells[iter.owner].points.push_back(iter.points[i]);
				
				l=0;
				for(int j=0; j<(*this).cells[iter.neighbour].points.size(); ++j){
					if(
					(*this).cells[iter.neighbour].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) (*this).cells[iter.neighbour].points.push_back(iter.points[i]);
			}
			
		}
		else if(iter.getType() == SEMO_Types::BOUNDARY_FACE){
			for(int i=0; i<iter.points.size(); ++i){
				int l=0;
				for(int j=0; j<(*this).cells[iter.owner].points.size(); ++j){
					if(
					(*this).cells[iter.owner].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) (*this).cells[iter.owner].points.push_back(iter.points[i]);
			}
		}
		else if(iter.getType() == SEMO_Types::PROCESSOR_FACE){
			for(int i=0; i<iter.points.size(); ++i){
				int l=0;
				for(int j=0; j<(*this).cells[iter.owner].points.size(); ++j){
					if(
					(*this).cells[iter.owner].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) (*this).cells[iter.owner].points.push_back(iter.points[i]);
			}
		}
	}
	cout << "-> completed" << endl;
	
}


;
void SEMO_Mesh_Builder::setCountsProcFaces(){
	
	// // SEMO_Mesh_Builder& mesh = *this;
	// cout << " a " << endl;
	
	// countsProcFaces.clear();
	
	// cout << " b " << endl;
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	// cout << size << endl;
	
	
	// (*this).countsProcFaces.push_back(0);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	
	(*this).countsProcFaces.resize(size,0);
	for(int ip=0; ip<size; ++ip){
		(*this).countsProcFaces[ip] = 0;
		for(auto& bc : (*this).boundary){
			if(ip == bc.neighbProcNo){
				// cout << bc.nFaces << endl;
				(*this).countsProcFaces[ip] = bc.nFaces;
				break;
			}
		}
		
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	
}


void SEMO_Mesh_Builder::setDisplsProcFaces(){
	
	// SEMO_Mesh_Builder& mesh = *this;
	
	// displsProcFaces.clear();
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	displsProcFaces.resize(size,0);
	displsProcFaces[0] = 0;
	for(int ip=1; ip<size; ++ip){
		displsProcFaces[ip] = displsProcFaces[ip-1] + countsProcFaces[ip-1];
		
	}
	
	
}