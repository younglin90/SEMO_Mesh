
#include "build.h"
#include <cmath>
#include <array>
#include "mpi.h"
#include "geometric.h" 

void SEMO_Mesh_Builder::check(){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	int owner_max=-1;
	int neighbour_max=-1;
	int checkOwnerNeighbourReverse=0;
	for(auto& face : (*this).faces){
		owner_max = max(owner_max , face.owner);
		neighbour_max = max(neighbour_max , face.neighbour);
	}
	int cell_id_max = max(owner_max,neighbour_max);
	
	// check, owner < neighbour 
	for(auto& face : (*this).faces){
		if(face.neighbour > owner_max){
			// cerr << "#error face, owner = " << face.owner << " neighbour = " << face.neighbour << " max cell id = " << owner_max << endl;
			int tempnum = face.owner;
			face.owner = face.neighbour;
			face.neighbour = tempnum;
			
			// std::reverse(face.points.begin(),face.points.end());
			std::reverse(face.points.begin()+1,face.points.end());
			
			++checkOwnerNeighbourReverse;
		}
		// else{
			// if(face.owner < 0)
				// cerr << "#error face, owner = " << face.owner << " neighbour = " << face.neighbour << " max cell id = " << owner_max << endl;
		// }
	}
	
	if(size>1){
		if(rank == 0){
			vector<int> buffer(size,0);
			int gatherValue = checkOwnerNeighbourReverse;
			MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
			int sumBuffer = 0;
			for(auto& i : buffer) {
				sumBuffer += i;
				// cout << i << " | ";
			}
			if(sumBuffer>0){
				cout << "| #check face owner < neighbour, executed reverse : ";
				for(auto& i : buffer) {
					cout << i << " | ";
				}
				cout << endl;
			}
		}
		else{
			int gatherValue = checkOwnerNeighbourReverse;
			MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	else{		
		if(checkOwnerNeighbourReverse>0){	
			cout << "| #check face owner < neighbour, executed reverse : ";
			cout << checkOwnerNeighbourReverse << " | ";
			cout << endl;
		}
	}
	// cout << "ddddd" << endl;
	
	for(int i=0; i<(*this).faces.size(); ++i){
		SEMO_Face& face = (*this).faces[i];
		
		if(face.neighbour == -1){
			if(i < (*this).boundary[0].startFace){
				cout << "| #Warning : face i(" << i << ") < boundary start face(" 
				<< (*this).boundary[0].startFace << ")" << endl;
			}
		}
		
	}
	
	
	
	
	
	
}




void SEMO_Mesh_Builder::checkMatchingProcessorFace(){
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	SEMO_Mesh_Builder& mesh = *this;
	
	if(size>1){
		
		vector<int> countsFacePoints(size,0);
		vector<int> displsFacePoints(size,0);
		vector<int> procFacePoints;
		
		vector<int> counts(size,0);
		vector<int> displs(size,0);
		vector<double> procFaceX;
		vector<double> procFaceY;
		vector<double> procFaceZ;
		// cout<<"AAAAAAAA"<<endl;
		int tmpNum0 = 0;
		for(auto& boundary : mesh.boundary){
			if(boundary.neighbProcNo==-1) continue;
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			
			// counts[boundary.neighbProcNo] = boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				++tmpNum0;
				
				procFacePoints.push_back(face.points.size());
				countsFacePoints[boundary.neighbProcNo] += 1;
				
			// cout << rank << " " << boundary.neighbProcNo << " " << str << " " << end << endl;
				
			
				vector<double> tmpX;
				vector<double> tmpY;
				vector<double> tmpZ;
				for(auto& j : face.points){
					auto& point = mesh.points[j];
					tmpX.push_back(point.x);
					tmpY.push_back(point.y);
					tmpZ.push_back(point.z);
				}
		// cout<<"A"<<endl;
				if(rank > boundary.neighbProcNo){
					std::reverse(tmpX.begin()+1, tmpX.end());
					std::reverse(tmpY.begin()+1, tmpY.end());
					std::reverse(tmpZ.begin()+1, tmpZ.end());
				}
		// cout<<"B"<< face.points.size() << " " << tmpX.size() << endl;
				for(int j=0; j<face.points.size(); ++j){
					procFaceX.push_back(tmpX[j]);
					procFaceY.push_back(tmpY[j]);
					procFaceZ.push_back(tmpZ[j]);
				}
		// cout<<"C"<<endl;
				
				counts[boundary.neighbProcNo] += tmpX.size();
			}
		}
		// cout<<"BBBBBB"<<endl;
		
		displsFacePoints[0] = 0;
		displs[0] = 0;
		for(int ip=1; ip<size; ++ip){
			displsFacePoints[ip] = displsFacePoints[ip-1] + countsFacePoints[ip-1];
			displs[ip] = displs[ip-1] + counts[ip-1];
		}
		
		// for(int ip=0; ip<size; ++ip){
			// cout<< rank << " -> " << ip << " : "<< counts[ip] << endl;
		// }
		
		// int tmpNum = 0;
		// for(auto& i : mesh.faces){
			// if(i.getType() == SEMO_Types::PROCESSOR_FACE){
				// for(auto& j : i.points){
					// ++tmpNum;
				// }
			// }
		// }
		
		vector<int> tmpNum2(size,0);
		for(auto& boundary : mesh.boundary){
			if(boundary.neighbProcNo==-1) continue;
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			
			// counts[boundary.neighbProcNo] = boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				for(int j=0; j<face.points.size(); ++j){
					++tmpNum2[boundary.neighbProcNo];
				}
			}
		}
		// for(int ip=0; ip<size; ++ip){
			// cout<< rank << " -> " << ip << " : "<< tmpNum2[ip] << endl;
		// }
		
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		vector<int> procFacePoints_recv(displsFacePoints[size-1] + countsFacePoints[size-1],0.0);
		
		MPI_Alltoallv( procFacePoints.data(), countsFacePoints.data(), displsFacePoints.data(), MPI_INT, 
					   procFacePoints_recv.data(), countsFacePoints.data(), displsFacePoints.data(), MPI_INT, 
					   MPI_COMM_WORLD);
					   
		for(int i=0; i<procFacePoints.size(); ++i){
			if(procFacePoints[i] != procFacePoints_recv[i]){
				if(rank==0) {
					cout << "| #Warning : not matching # of proc-Face-point " << " " << i << " / " << procFacePoints.size() << " " << 
					procFacePoints[i] << " /= " << procFacePoints_recv[i] << endl;
				}
			}
		}
			
					   
		
		
		
		vector<double> procFaceX_recv(displs[size-1] + counts[size-1],0.0);
		vector<double> procFaceY_recv(displs[size-1] + counts[size-1],0.0);
		vector<double> procFaceZ_recv(displs[size-1] + counts[size-1],0.0);
		
		MPI_Alltoallv( procFaceX.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   procFaceX_recv.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( procFaceY.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   procFaceY_recv.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( procFaceZ.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   procFaceZ_recv.data(), counts.data(), displs.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
					
		SEMO_Utility_Math utility;
					   
		tmpNum0 = 0;
		for(auto& boundary : mesh.boundary){
			if(boundary.neighbProcNo==-1) continue;
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				for(auto& j : face.points){
					auto& point = mesh.points[j];
					double sendX = procFaceX[tmpNum0];
					double recvX = procFaceX_recv[tmpNum0];
					double sendY = procFaceY[tmpNum0];
					double recvY = procFaceY_recv[tmpNum0];
					double sendZ = procFaceZ[tmpNum0];
					double recvZ = procFaceZ_recv[tmpNum0];
					
					if( 
					utility.approximatelyEqualAbsRel(sendX, recvX, 1.e-8, 1.e-6) &&
					utility.approximatelyEqualAbsRel(sendY, recvY, 1.e-8, 1.e-6) &&
					utility.approximatelyEqualAbsRel(sendZ, recvZ, 1.e-8, 1.e-6) 
					){
						// if(rank==0) {
							// cout << "OKOKOK" << " " << i << " " << face.level << " " << point.level << endl;
							// cout << sendX << " " << " " << sendY << " " << " " << sendZ << " " << endl;
							// cout << recvX << " " << " " << recvY << " " << " " << recvZ << " " << endl;
						// }
					}
					else{
						// if(rank==0) {
							cout << "| #Warning : not matching x,y,z of proc-Face-point " << " " << i << " " << face.level << " " << point.level << endl;
							cout << sendX << " " << " " << sendY << " " << " " << sendZ << " " << endl;
							cout << recvX << " " << " " << recvY << " " << " " << recvZ << " " << endl;
							
							MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
						// }
					}
					++tmpNum0;
				}
				
			}
		}
		// for(int i=0; i<procFaceX.size(); ++i){
			// double sendX = procFaceX[i];
			// double recvX = procFaceX_recv[i];
			// double sendY = procFaceY[i];
			// double recvY = procFaceY_recv[i];
			// double sendZ = procFaceZ[i];
			// double recvZ = procFaceZ_recv[i];
			
			// if( 
			// utility.approximatelyEqualAbsRel(sendX, recvX, 1.e-8, 1.e-6) &&
			// utility.approximatelyEqualAbsRel(sendY, recvY, 1.e-8, 1.e-6) &&
			// utility.approximatelyEqualAbsRel(sendZ, recvZ, 1.e-8, 1.e-6) 
			// ){
			// }
			// else{
				// if(rank==0) {
					// cout << "NONONONO" << " " << sendX << " " << recvX << " " << 
					// sendY << " " << recvY << " " << sendZ << " " << recvZ << endl;
				// }
			// }
			
			
			
		// }
					   
		
		
	}
	
	
}




void SEMO_Mesh_Builder::checkQualities(){
	
	// check, orthogonality
	for(auto& face : (*this).faces){
		vector<double> faceCenterXYZ(3,0.0);
		for(auto& i : face.points){
			faceCenterXYZ[0] += (*this).points[i].x;
			faceCenterXYZ[1] += (*this).points[i].y;
			faceCenterXYZ[2] += (*this).points[i].z;
		}
		int facePointsSize = face.points.size();
		faceCenterXYZ[0] /= (double)facePointsSize;
		faceCenterXYZ[1] /= (double)facePointsSize;
		faceCenterXYZ[2] /= (double)facePointsSize;
		
		vector<double> cellCenterXYZ(3,0.0);
		for(auto& i : (*this).cells[face.owner].points){
			cellCenterXYZ[0] += (*this).points[i].x;
			cellCenterXYZ[1] += (*this).points[i].y;
			cellCenterXYZ[2] += (*this).points[i].z;
		}
		int cellPointsSize = (*this).cells[face.owner].points.size();
		cellCenterXYZ[0] /= (double)cellPointsSize;
		cellCenterXYZ[1] /= (double)cellPointsSize;
		cellCenterXYZ[2] /= (double)cellPointsSize;
		
		vector<double> cellToFaceNvec(3,0.0); 
		cellToFaceNvec[0] = faceCenterXYZ[0] - cellCenterXYZ[0];
		cellToFaceNvec[1] = faceCenterXYZ[1] - cellCenterXYZ[1];
		cellToFaceNvec[2] = faceCenterXYZ[2] - cellCenterXYZ[2];
		double normalized = sqrt(cellToFaceNvec[0]*cellToFaceNvec[0] 
								+cellToFaceNvec[1]*cellToFaceNvec[1]
								+cellToFaceNvec[2]*cellToFaceNvec[2]);
		cellToFaceNvec[0] /= normalized;
		cellToFaceNvec[1] /= normalized;
		cellToFaceNvec[2] /= normalized;
		
		
		
		SEMO_Mesh_Geometric geometric;
		vector<double> faceNvec(3,0.0); 
		
		vector<double> Vx, Vy, Vz;
		for(auto& i : face.points){
			Vx.push_back((*this).points[i].x);
			Vy.push_back((*this).points[i].y);
			Vz.push_back((*this).points[i].z);
		}
		
		double VSn;
		vector<double> cellCentroid;
		geometric.calcUnitNormals_Area3dPolygon(
			face.points.size(), Vx,Vy,Vz,
			face.unitNormals, face.area,
			face.x, face.y, face.z, VSn, cellCentroid );
			
		vector<double> orthogonality(3,0.0); 
		orthogonality[0] = face.unitNormals[0] - cellToFaceNvec[0];
		orthogonality[1] = face.unitNormals[1] - cellToFaceNvec[1];
		orthogonality[2] = face.unitNormals[2] - cellToFaceNvec[2];

		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // cout << "INTERNAL_FACE" << endl;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// // cout << "BOUNDARY_FACE" << endl;
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // cout << "PROCESSOR_FACE" << endl;
		// }
		// else{
			// cout << "PROCESSOR_FACE" << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
			
		if( (orthogonality[0] > -0.001 && orthogonality[0] < 0.001) &&
		    (orthogonality[1] > -0.001 && orthogonality[1] < 0.001) &&
			(orthogonality[2] > -0.001 && orthogonality[2] < 0.001) ){
			
			
			
		}
		else{
			
			cout << endl;
			cout << " # Warning : orthogonality large" << endl;
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				cout << "INTERNAL_FACE" << endl;
			}
			else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				cout << "BOUNDARY_FACE" << endl;
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cout << "PROCESSOR_FACE" << endl;
			}
			
			// cout << face.owner << " " << face.neighbour << endl;
			// for(auto& j : face.points){
				// cout << j << endl;
			// }
			// cout << face.points[0] << " " << face.points[1] << endl;
			// cout << (*this).cells[face.owner].points.size() << endl;
			// cout << face.points.size() << endl;
			cout << "cell center : " << cellCenterXYZ[0] << " " << cellCenterXYZ[1] << " " << cellCenterXYZ[2] << endl;
			cout << "face center : " << faceCenterXYZ[0] << " " << faceCenterXYZ[1] << " " << faceCenterXYZ[2] << endl;
			cout << "face Normal : " << face.unitNormals[0] << " " << face.unitNormals[1] << " " << face.unitNormals[2] << endl;
			cout << "cell->face Normal : " << cellToFaceNvec[0] << " " << cellToFaceNvec[1] << " " << cellToFaceNvec[2] << endl;
			cout << "orthogonality : " << orthogonality[0] << " " << orthogonality[1] << " " << orthogonality[2] << endl;
		}
			
	}
	
	
	
	// // check, skewness
	
	
	
	
	
	// // check, aspect ratio 
	
	
	
	
}


void SEMO_Mesh_Builder::loadFile(string filetype, string folder){
	SEMO_Mesh_Load meshLoad;
	 
	if(filetype=="OpenFoam"){
		meshLoad.OpenFoam(*this, folder);
	} 
	if(filetype=="vtu"){
		cout << "| NOT defined vtu" << endl;
		// meshLoad.vtu(*this, folder);
	} 
}



void SEMO_Mesh_Builder::saveFile(string filetype, string folder, SEMO_Controls_Builder &controls){
	SEMO_Mesh_Save meshSave;
	 
	// if(filetype=="vtu"){
		// // meshSave.vtu(*this, controls);
		// meshSave.vtu(folder, *this, controls);
	// } 
	// if(filetype=="vtuZlib"){
		// meshSave.vtuZlib(*this, controls);
	// } 
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
		// cout << startFace << " " << (*this).boundary[ibc].nFaces << " " << (*this).boundary.size() << endl;
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
	(*this).listFaces.clear();
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
		
		(*this).listFaces.push_back(&*iter);
		
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

	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	// cell connection (cell's face)
		
	// if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	// if(rank==0) cout << "| execute cell's face connecting ... ";
	for(int i=0; i<(*this).faces.size(); ++i){
		SEMO_Face& face = (*this).faces[i];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
		// cout << face.owner << " " << face.neighbour << " " << (*this).cells.size() << endl;
			(*this).cells[face.owner].faces.push_back(i);
			(*this).cells[face.neighbour].faces.push_back(i);
		}
		else if(
		face.getType() == SEMO_Types::BOUNDARY_FACE ||
		face.getType() == SEMO_Types::PROCESSOR_FACE ){
			(*this).cells[face.owner].faces.push_back(i);
		}
	}
	// if(rank==0) cout << "-> completed" << endl;
	// if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
		
	

}




void SEMO_Mesh_Builder::connectCelltoPoints(){

	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 

	// cell connection (cell's points)
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
	// if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	// if(rank==0) cout << "| execute cell's points connecting ... ";
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
	// if(rank==0) cout << "-> completed" << endl;
	// if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	
}


void SEMO_Mesh_Builder::searchNeighbProcFaces(){
	
	
	SEMO_Mesh_Builder& mesh = *this;
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_Utility_Math utility;
	
	
	// vector<double> fsize;
	vector<double> fx;
	vector<double> fy;
	vector<double> fz;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			double tmpx=0.0;
			double tmpy=0.0;
			double tmpz=0.0;
			for(auto& point : face.points){
				tmpx += mesh.points[point].x;
				tmpy += mesh.points[point].y;
				tmpz += mesh.points[point].z;
			}
			tmpx /= (double)face.points.size();
			tmpy /= (double)face.points.size();
			tmpz /= (double)face.points.size();
			
			fx.push_back(tmpx);
			fy.push_back(tmpy);
			fz.push_back(tmpz);
		}
	}
	
	int fsize = fx.size();
	
	vector<int> buffer(size,-1);
	MPI_Allgather(&fsize, 1, MPI_INT, buffer.data(), 1, MPI_INT, MPI_COMM_WORLD);
	
	vector<int> counts(size,0);
	vector<int> displacements(size,0);
	for(int ip=0; ip<size; ++ip){
		counts[ip] = buffer[ip];
	}
	for(int ip=1; ip<size; ++ip){
		displacements[ip] = displacements[ip-1] + counts[ip-1];
	}
	vector<double> xBuffer(displacements[size-1]+counts[size-1],0.0);
	vector<double> yBuffer(displacements[size-1]+counts[size-1],0.0);
	vector<double> zBuffer(displacements[size-1]+counts[size-1],0.0);
	
    MPI_Allgatherv(fx.data(), fsize, MPI_DOUBLE, xBuffer.data(), counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(fy.data(), fsize, MPI_DOUBLE, yBuffer.data(), counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgatherv(fz.data(), fsize, MPI_DOUBLE, zBuffer.data(), counts.data(), displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
 
	// if(rank==0){
		// for(int ip=0; ip<size; ++ip){
			// int str = displacements[ip];
			// int end = displacements[ip]+counts[ip];
			// for(int i=str; i<end; ++i){
				// cout << ip << " " << fx[i] << " " << fy[i] << " " << fz[i] << endl;
			// }
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	
	// vector<int> nProcFaces(size,0);
	// vector<vector<int>> testMachI(size,vector<int>(0,0));
	vector<vector<int>> test2(size,vector<int>(0,0));
	
	int strBase = displacements[rank];
	int endBase = displacements[rank]+counts[rank];
	int num2=0;
	for(int base=strBase; base<endBase; ++base){
		// int matchI = -1;
		// int matchIp = -1;
		double minX=100000.0;
		double minY=100000.0;
		double minZ=100000.0;
		bool matching = false;
		int saveNum=-1;
		int saveIp=-1;
		for(int ip=0; ip<size; ++ip){
			if(rank==ip) continue;
			int str = displacements[ip];
			int end = displacements[ip]+counts[ip];
			int num=0;
			for(int i=str; i<end; ++i){
				// int cx = utility.CompareDoubleAbsoulteAndUlps(fx[base], fx[i], 1.0e-2, 4);
				// int cy = utility.CompareDoubleAbsoulteAndUlps(fy[base], fy[i], 1.0e-2, 4);
				// int cz = utility.CompareDoubleAbsoulteAndUlps(fz[base], fz[i], 1.0e-2, 4);
				bool cx = utility.approximatelyEqualAbsRel(xBuffer[base], xBuffer[i], 1.e-8, 1.e-6);
				bool cy = utility.approximatelyEqualAbsRel(yBuffer[base], yBuffer[i], 1.e-8, 1.e-6);
				bool cz = utility.approximatelyEqualAbsRel(zBuffer[base], zBuffer[i], 1.e-8, 1.e-6);
				// bool cx = false;
				// bool cy = false;
				// bool cz = false;
				// if(xBuffer[base]-xBuffer[i]>-0.1 && xBuffer[base]-xBuffer[i]<0.1) cx=true;
				// if(yBuffer[base]-yBuffer[i]>-0.1 && yBuffer[base]-yBuffer[i]<0.1) cy=true;
				// if(zBuffer[base]-zBuffer[i]>-0.1 && zBuffer[base]-zBuffer[i]<0.1) cz=true;
				if(
				abs(xBuffer[base]-xBuffer[i])<minX &&
				abs(yBuffer[base]-yBuffer[i])<minY &&
				abs(zBuffer[base]-zBuffer[i])<minZ ){
					saveIp = ip;
					saveNum = num2;
					minX=abs(xBuffer[base]-xBuffer[i]);
					minY=abs(yBuffer[base]-yBuffer[i]);
					minZ=abs(zBuffer[base]-zBuffer[i]);
				}
				
				if(cx && cy && cz){
				// if(cx==0 && cy==0 && cz==0){
					matching = true;
					// matchI = i;
					// matchIp = ip;
					// ++nProcFaces[ip];
					// testMachI[ip].push_back(num);
					test2[ip].push_back(num2);
					break;
				}
				++num;
			}
			if(matching) break;
		}
		if(!matching){
			cout << endl;
			cout << "#Warning : PROC face NOT matching" << endl;
			// cout << endl;
			// cout << "#ERROR : PROC face NOT matching" << endl;
			cout << xBuffer[base] << " " << yBuffer[base] << " " << zBuffer[base] << endl;
			cout << minX << " " << minY << " " << minZ << endl;
			
			test2[saveIp].push_back(saveNum);
			// cout << endl;
			// cout << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		++num2;
	}
	

	
	
	

	vector<SEMO_Face> procFaces;
	
	// int startFace=0;
	// bool boolStartFace=false;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// boolStartFace = true;
			
			SEMO_Face tmpFace;
			
			for(auto& point : face.points){
				tmpFace.points.push_back(point);
			}
			
			tmpFace.owner = face.owner;
			
			procFaces.push_back(tmpFace);
		}
		// if(!boolStartFace) ++startFace;
	}
	
	int num_proc=0;
	vector<int> saveChangeNum;
	vector<bool> boolExecuted;
	vector<int> matOrder2;
	
	for(int ip=0; ip<size; ++ip){
		for(auto& i : test2[ip]){
			matOrder2.push_back(i);
			boolExecuted.push_back(false);
		}
	}
	
	int num_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// if(boolExecuted[num_proc]==false){
				int changeNum = matOrder2[num_proc];
				
				// boolExecuted[num_proc] = true;
				// boolExecuted[changeNum] = true;
				
				mesh.faces[num_face].points.clear();
				for(auto& point : procFaces[changeNum].points){
					mesh.faces[num_face].points.push_back(point);
				}
				mesh.faces[num_face].owner = procFaces[changeNum].owner;
				
				// mesh.faces[num_face+changeNum].points.clear();
				// for(auto& point : procFaces[num_proc].points){
					// mesh.faces[num_face+changeNum].points.push_back(point);
				// }
				// mesh.faces[num_face+changeNum].owner = procFaces[num_face].owner;
				
			// }
			
			++num_proc;
		}
		++num_face;
	}
	
	
		
	
	int startFace=0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			break;
		}
		++startFace;
	}
	
	int nbcs=0;
	for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
		if(mesh.boundary[ibcs].neighbProcNo == -1) ++nbcs;
	}
	
	mesh.boundary.erase(mesh.boundary.begin() + nbcs, mesh.boundary.end());
	
	for(int ip=0; ip<size; ++ip){
		if( test2[ip].size()>0){
			mesh.addBoundary();
			string bcnames = "procBoundary" + to_string(rank) + "to" + to_string(ip);
			mesh.boundary.back().name = bcnames;
			mesh.boundary.back().nFaces = test2[ip].size();
			mesh.boundary.back().startFace = startFace;
			mesh.boundary.back().myProcNo = rank;
			mesh.boundary.back().neighbProcNo = ip;
			
			startFace += test2[ip].size();
		}
		// test2[ip].size();
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
}
	
	
	
	
	
	
	
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
	
	(*this).countsProcFaces.clear();
	(*this).countsProcFaces.resize(size,0);
	
	for(int ip=0; ip<size; ++ip){
		(*this).countsProcFaces[ip] = 0;
		for(auto& bc : (*this).boundary){
			if(ip == bc.neighbProcNo){
				// if(rank==ip && bc.nFaces>0){
					// cout << endl;
					// cout << "| #Warning : rank(" << rank << ")==ip(" << ip << "), nFaces=" << bc.nFaces << endl;
					// cout << endl;
				// }
				// cout << bc.nFaces << endl;
				
				// cout << rank << " " << ip << " " << bc.nFaces << endl;
				
				// cout << rank << " -> " << ip << " " << bc.nFaces << endl;
				
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
	
	displsProcFaces.clear();
	displsProcFaces.resize(size,0);
	displsProcFaces[0] = 0;
	for(int ip=1; ip<size; ++ip){
		displsProcFaces[ip] = displsProcFaces[ip-1] + countsProcFaces[ip-1];
		
	}
	
	
}








void SEMO_Mesh_Builder::informations(){
	
	SEMO_Mesh_Builder& mesh = *this;
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	// cout << "AAAAAAAA" << endl;
	
	int nbcs=0;
	int nFacesBC = 0;
	int nFacesProc = 0;
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			nFacesBC += boundary.nFaces;
			++nbcs;
		}
		else{
			nFacesProc += boundary.nFaces;
		}
	}
	
	
	// faces type
	int nFacesInt = 0;
	int nTriangle = 0;
	int nQuadrangle = 0;
	int nPolygon = 0;
	for(auto& face : mesh.faces){
		
		if(face.neighbour != -1){
			++nFacesInt;
		}
		
		if(face.points.size() == 3){
			++nTriangle;
		}
		else if(face.points.size() == 4){
			++nQuadrangle;
		}
		else{
			++nPolygon;
		}
	}
	
	
	
	// cells type
	int nTetrahedron = 0;
	int nHexahedron = 0;
	int nPrism = 0;
	int nPyramid = 0;
	int nPolyhedron = 0;
	for(auto& cell : mesh.cells){
		if( cell.points.size() == 4 && cell.faces.size() == 4 ){
			++nTetrahedron;
		}
		else if( cell.points.size() == 5 && cell.faces.size() == 5 ){
			++nPyramid;
		}
		else if( cell.points.size() == 6 && cell.faces.size() == 5 ){
			++nPrism;
		}
		else if( cell.points.size() == 8 && cell.faces.size() == 6 ){
			++nHexahedron;
		}
		else {
			++nPolyhedron;
		}
	}
	// cout << "BBBBBBBB" << endl;
	
    if(rank == 0){
        vector<int> buffer(size,0);
		int gatherValue = mesh.points.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| present MPI size : " << size;
		cout << "| load data MPI size : " << size << endl;
		cout << "| points size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = mesh.faces.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = mesh.cells.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| cells size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nFacesInt;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| internal faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nFacesBC;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| boundary faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nFacesProc;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| processor faces size : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "| boundary types : " << nbcs << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nTriangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Triangle faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nQuadrangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Quadrangle faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPolygon;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Polygon faces : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		gatherValue = nTetrahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Tetrahedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nHexahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Hexahedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPrism;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Prism cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPyramid;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Pyramid cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		gatherValue = nPolyhedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, buffer.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "| Polyhedron cells : ";
		for(auto& i : buffer) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
    }
    else{
		int gatherValue = mesh.points.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = mesh.faces.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = mesh.cells.size();
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nFacesInt;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nFacesBC;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nFacesProc;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		gatherValue = nTriangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nQuadrangle;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPolygon;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		
		gatherValue = nTetrahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nHexahedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPrism;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPyramid;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
		gatherValue = nPolyhedron;
        MPI_Gather(&gatherValue, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}












void SEMO_Mesh_Builder::cellsGlobal(){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	
	SEMO_Mesh_Builder& mesh = *this;
	

	int ncells = mesh.cells.size();
	vector<int> procNcells(size,0);
	MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	int myStart=0;
	for(int i=0; i<rank; ++i){
		myStart += procNcells[i];
	}
	
	
	
	int ncellTot;
	MPI_Allreduce(&ncells, &ncellTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	vector<int> procStart(size,0);
	MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	mesh.startCellGlobal = myStart;
	mesh.startProcCellGlobal.resize(size+1,0);
	for(int i=0; i<size; ++i){
		mesh.startProcCellGlobal[i] = procStart[i];
	}
	mesh.ncellsTotal = ncellTot;
	mesh.startProcCellGlobal[size] = ncellTot;
	
	vector<int> neighbProcNo;
	vector<int> ownerNo;
	vector<int> neighbNo;
	
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo != -1){
			int str=boundary.startFace;
			int end=str+boundary.nFaces;
			for(int i=str; i<end; ++i){
				neighbProcNo.push_back(boundary.neighbProcNo);
				ownerNo.push_back(mesh.faces[i].owner);
				// cout << mesh.faces[i].owner << endl;
			}
		}
	}
	
	if(size>1){
		
		neighbNo.clear();
		neighbNo.resize(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		
		MPI_Alltoallv( ownerNo.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   neighbNo.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);

		// SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// ownerNo, neighbNo,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
	
	}
	
	mesh.neighbProcNo.resize(neighbProcNo.size(),0);
	for(int i=0; i<neighbProcNo.size(); ++i){
		mesh.neighbProcNo[i] = neighbProcNo[i];
	}
	
	mesh.procNeighbCellNo.resize(neighbNo.size(),0);
	for(int i=0; i<neighbNo.size(); ++i){
		mesh.procNeighbCellNo[i] = neighbNo[i];
	}
	
	
	
	
	
	
	// temporary COO format
	int strRow = mesh.startCellGlobal;
	vector<int> COO_row(mesh.cells.size(),0);
	vector<int> COO_col(mesh.cells.size(),0);
	vector<int> COO_iface;
	vector<string> COO_iface_LR;
	for(int i=0; i<mesh.cells.size(); ++i){
		COO_row[i] = strRow + i;
		COO_col[i] = strRow + i;
	}
	for(int i=0, proc_num=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			COO_row.push_back(strRow + face.owner);
			COO_col.push_back(strRow + face.neighbour);
			
			COO_iface_LR.push_back("LR");
			COO_iface.push_back(i);
			
			COO_row.push_back(strRow + face.neighbour);
			COO_col.push_back(strRow + face.owner);
			
			COO_iface_LR.push_back("RL");
			COO_iface.push_back(i);
		}
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int strNeigbRow_Glo = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
			int idNeigbCell_Loc = mesh.procNeighbCellNo[proc_num];
			
			COO_row.push_back(strRow + face.owner);
			COO_col.push_back(strNeigbRow_Glo + idNeigbCell_Loc);
			
			COO_iface_LR.push_back("LR");
			COO_iface.push_back(i);
			
			++proc_num;
		}
	}
	
	// save CSR format
    //compute number of non-zero entries per row of A 
	int n_row = mesh.cells.size();
	int nnz = COO_row.size();
	
	mesh.non_zeros = nnz;
	
	mesh.CRS_ptr.clear();
	mesh.CRS_col.clear();
	mesh.CRS_col_ptr_dig.clear();
	mesh.CRS_col_ptr_LR.clear();
	mesh.CRS_col_ptr_RL.clear();
	
	mesh.CRS_col.resize(nnz,0);
	mesh.CRS_ptr.resize(n_row+1,0);
	mesh.CRS_col_ptr_dig.resize(mesh.cells.size(),-1);
	mesh.CRS_col_ptr_LR.resize(mesh.faces.size(),-1);
	mesh.CRS_col_ptr_RL.resize(mesh.faces.size(),-1);

    for (int n = 0; n < nnz; n++){            
        mesh.CRS_ptr[COO_row[n]-strRow]++;
    }

    //cumsum the nnz per row to get ptr[]
    for(int i = 0, cumsum = 0; i < n_row; i++){     
        int temp = mesh.CRS_ptr[i];
        mesh.CRS_ptr[i] = cumsum;
        cumsum += temp;
    }
    mesh.CRS_ptr[n_row] = nnz; 

    //write col,val into col,val
    for(int n = 0; n < nnz; n++){
        int row  = COO_row[n] - strRow;
        int dest = mesh.CRS_ptr[row];

        mesh.CRS_col[dest] = COO_col[n];
        // val[dest] = A_vals[n];
		
		if(n<n_row){
			mesh.CRS_col_ptr_dig[n] = dest;
		}
		else{
			int iface = COO_iface[n-n_row];
			
			if(COO_iface_LR[n-n_row] == "LR"){
				mesh.CRS_col_ptr_LR[iface] = dest;
			}
			else{
				mesh.CRS_col_ptr_RL[iface] = dest;
			}
			
		}

        mesh.CRS_ptr[row]++;
    }

    for(int i = 0, last = 0; i <= n_row; i++){
        int temp = mesh.CRS_ptr[i];
        mesh.CRS_ptr[i]  = last;
        last = temp;
    }
	
	
}




void SEMO_Mesh_Builder::calcSkewness(vector<double>& output){
	
	SEMO_Mesh_Builder& mesh = *this;
	
	output.clear();
	output.resize(mesh.cells.size(),0.0);

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
		output[face.owner] = max(output[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			output[face.neighbour] = max(output[face.neighbour],value);
		}
	}
		
}

void SEMO_Mesh_Builder::calcNonOrthogonality(vector<double>& output){
	
	SEMO_Mesh_Builder& mesh = *this;
	
	output.clear();
	output.resize(mesh.cells.size(),0.0);
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
		// cout << face.owner << endl;
		output[face.owner] = max(output[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			output[face.neighbour] = max(output[face.neighbour],value);
		}
	}
	
}




void SEMO_Mesh_Builder::calcUniformity(vector<double>& output){
	
	SEMO_Mesh_Builder& mesh = *this;
	
	output.clear();
	output.resize(mesh.cells.size(),1.e8);
	for(auto& face : mesh.faces){
		double value1=0.0;
		double value2=0.0;
		for(int ii=0; ii<3; ++ii){
			value1 += face.vecPF[ii]*face.vecPF[ii];
			value2 += face.vecNF[ii]*face.vecNF[ii];
		}
		value1 = sqrt(value1);
		value2 = sqrt(value2);
		value2 = value1 + value2;
		double value = value1/value2;
		value = 1.0-2.0*abs(value-0.5);
		output[face.owner] = min(output[face.owner],value);
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			output[face.neighbour] = min(output[face.neighbour],value);
		}
	}
		
}