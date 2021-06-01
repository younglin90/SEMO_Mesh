#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "build.h"
#include "hexaOctAMR.h"
// #include <cmath>
// #include <array>
#include "mpi.h"

#include <random>

void SEMO_Mesh_Builder::hexaOctAMR(){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_Mesh_Builder& mesh = *this;
	
	SEMO_Hexa_Oct_AMR AMR;
	
	vector<bool> executeAMR(mesh.cells.size(), false);
	
	vector<SEMO_Cell> newCells;
	vector<SEMO_Face> newFaces;
	vector<SEMO_Point> newPoints;
	
	vector<int> idNewCells;
	vector<int> idNewFaces;
	vector<int> idNewPoints;
	
	vector<int> levelNewCells;
	
	vector<int> startFaceBC;
	vector<int> nFacesBC;
	vector<int> iFaceBC(mesh.listFaces.size(),-1);
	
	vector<int> startFace(size,0);
	vector<int> nFaces(size,0);
	vector<int> neighbProcNo(mesh.listFaces.size(),-1);
	
	vector<vector<int>> idNewFaceProc(size,vector<int>(0,0));
		
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	executeAMR[0] = true;
	
	int num=0;
	int nbcs=0;
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1){
			startFaceBC.push_back(boundary.startFace);
			nFacesBC.push_back(boundary.nFaces);
			for(int i=boundary.startFace; i<boundary.startFace+boundary.nFaces; ++i){
				iFaceBC[i] = num;
			}
			++nbcs;
			++num;
		}
		
	}
	
	vector<vector<int>> idNewFaceBC(startFaceBC.size(),vector<int>(0,0));
	
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo != -1){
			startFace[boundary.neighbProcNo] = boundary.startFace;
			nFaces[boundary.neighbProcNo] = boundary.nFaces;
			for(int i=boundary.startFace; i<boundary.startFace+boundary.nFaces; ++i){
				neighbProcNo[i] = boundary.neighbProcNo;
			}
		}
		
	}
	
	
	//================================================
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			if( mesh.cells[face.owner].level > 
				mesh.cells[face.neighbour].level ){
				
				executeAMR[face.owner] = false;

			}
			if( mesh.cells[face.owner].level < 
				mesh.cells[face.neighbour].level ){
				
				executeAMR[face.neighbour] = false;

			}
		}
	}
	
	
	//================================================
	
		
	
	int newNPoints = 0;
	int newNFaces = 0;
	int newNCells = 0;
	
	num=0;
	for(SEMO_Cell* cell : mesh.listCells){
		
		// SEMO_Cell& cell = cell;
		
		if( cell->points.size() == 8 ){
			bool boolExecuteAMR = false;
			if( executeAMR[num] == true ){
				int four=0;
				for(auto& i : cell->faces){
					if(mesh.faces[i].points.size() == 4) ++four;
				}
				if(four==6) boolExecuteAMR = true;
			}
			if(boolExecuteAMR){
				
				newPoints.clear();
				idNewFaces.clear();
				newFaces.clear();
				idNewCells.clear();
				levelNewCells.clear();
				
				int strPoints = newPoints.size();
				int strFaces = newNFaces;
				int strCells = newNCells;
			
				vector<int> faces;
				for(int i=0; i<6; ++i) faces.push_back(cell->faces[i]);
				
				vector<int> points;
				for(int i=0; i<6; ++i){
					if( &mesh.cells[ mesh.faces[ cell->faces[i]].owner ] == cell ){
						for(int j=0; j<4; ++j) points.push_back( mesh.faces[cell->faces[i]].points[j] );
					}
					else{
						for(int j=3; j>=0; --j) points.push_back( mesh.faces[cell->faces[i]].points[j] );
					}
				}
				vector<int> cellPoints;
				for(int i=0; i<8; ++i) cellPoints.push_back(cell->points[i]);
				
				AMR.searchLoacationsCase0( faces, points, cellPoints );
				
				vector<vector<int>> cellPointsSet(19,vector<int>(2,0));
				
				AMR.sortCellPointsSet(cellPoints, cellPointsSet);
				
					
					
				for(int i=0; i<19; ++i){
					SEMO_Point tempPoint;
					tempPoint.x = 0.0;
					tempPoint.y = 0.0;
					tempPoint.z = 0.0;
					
					// cout << i << " " << cellPointsSet[i].size() << endl;
					for(auto& j : cellPointsSet[i]){
						tempPoint.x += mesh.points[j].x / (double)cellPointsSet[i].size();
						tempPoint.y += mesh.points[j].y / (double)cellPointsSet[i].size();
						tempPoint.z += mesh.points[j].z / (double)cellPointsSet[i].size();
					}
					
					newPoints.push_back(tempPoint);
					++newNPoints;
				}
				
				
				
				int idPoint = mesh.listPoints.size() + strPoints;
				
				vector<vector<int>> newFacesPointsSet(36,vector<int>(4,0));
				
				AMR.sortFacesPointsSet(cellPoints, idPoint, newFacesPointsSet);
				
				for(int i=0; i<2; ++i){
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] == cell ){
						for(int j=0; j<4; ++j) 
							reverse( &newFacesPointsSet[i*4+j][0], &newFacesPointsSet[i*4+j][4]);
					}
				}
				for(int i=2; i<6; ++i){
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] != cell ){
						for(int j=0; j<4; ++j)
							reverse( &newFacesPointsSet[i*4+j][0], &newFacesPointsSet[i*4+j][4]);
					}
				}
				
				vector<vector<int>> newFaceOwnNgbSet(6,vector<int>(4,0));
				
				int strOwnNgb = mesh.listCells.size() + strCells;
				
				AMR.sortOwnNgbOuterFaces(num, strOwnNgb, newFaceOwnNgbSet);

				
				for(int i=0; i<6; ++i){
					
					bool directionBasedOwner = false;
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] == cell )
						directionBasedOwner = true;
					
					for(int j=0; j<4; ++j){
						SEMO_Face tempFace;
						for(int k=0; k<4; ++k){
							tempFace.points.push_back( newFacesPointsSet[i*4+j][k] );
						}
						if(directionBasedOwner){
							tempFace.owner = newFaceOwnNgbSet[i][j];
							tempFace.neighbour = mesh.faces[ faces[i] ].neighbour;
						}
						else{
							tempFace.owner = mesh.faces[ faces[i] ].owner;
							tempFace.neighbour = newFaceOwnNgbSet[i][j];
						}
						
						tempFace.setType( mesh.faces[ faces[i] ].getType() );
						
						int idFace = mesh.listFaces.size() + newNFaces;
						if(j==0) idFace = faces[i];
						
						if( mesh.faces[ faces[i] ].getType() == SEMO_Types::BOUNDARY_FACE ){
							idNewFaceBC[ iFaceBC[ faces[i] ] ].push_back( idFace );
						}
						else if( mesh.faces[ faces[i] ].getType() == SEMO_Types::PROCESSOR_FACE ){
							int procNum = neighbProcNo[faces[i]];
							idNewFaceProc[procNum].push_back( idFace );
						}
						
						idNewFaces.push_back( idFace );
						
						newFaces.push_back(tempFace);
						
						if(j!=0) ++newNFaces;
					}
				}
				
				vector<vector<int>> newFaceOwnSet(3,vector<int>(4,0));
				vector<vector<int>> newFaceNgbSet(3,vector<int>(4,0));
				
				AMR.sortOwnNgbInternalFaces(num, strOwnNgb, newFaceOwnSet, newFaceNgbSet);
				
				
				// cell internal faces
				for(int i=6; i<9; ++i){
					for(int j=0; j<4; ++j){
						SEMO_Face tempFace;
						for(int k=0; k<4; ++k)
							tempFace.points.push_back( newFacesPointsSet[i*4+j][k] );
						
						tempFace.owner = newFaceOwnSet[i-6][j];
						tempFace.neighbour = newFaceNgbSet[i-6][j];
						
						tempFace.setType(SEMO_Types::INTERNAL_FACE);
							
						idNewFaces.push_back(mesh.listFaces.size() + newNFaces);
						
						++newNFaces;
						
						newFaces.push_back(tempFace);
						
					} 
				}
				
				for(int i=0; i<8; ++i){
					
					int idCell = mesh.listCells.size() + newNCells;
					if(i==0) idCell = num;
					
					idNewCells.push_back(idCell);
					levelNewCells.push_back(mesh.cells[num].level + 1);
					
					if(i!=0) ++newNCells;
					
				}
				
				
				AMR.addPoints(mesh, newPoints);
				AMR.addFaces(mesh, idNewFaces, newFaces);
				AMR.addCells(mesh, idNewCells, levelNewCells);
				
			}
		}
		else if( cell->points.size() == 13 ){
			bool boolExecuteAMR = false;
			if( executeAMR[num] == true ){
				int four=0;
				int five=0;
				for(auto& i : cell->faces){
					if(mesh.faces[i].points.size() == 4) ++four;
					if(mesh.faces[i].points.size() == 5) ++five;
				}
				if(four==5 && five==4) boolExecuteAMR = true;
			}
			if(boolExecuteAMR){
				
				newPoints.clear();
				idNewFaces.clear();
				newFaces.clear();
				idNewCells.clear();
				levelNewCells.clear();

				int strPoints = newPoints.size();
				int strFaces = newNFaces;
				int strCells = newNCells;
			
				vector<int> faces;
				for(int i=0; i<9; ++i) faces.push_back(cell->faces[i]);
				
				vector<int> points;
				for(int i=0; i<9; ++i){
					if( &mesh.cells[ mesh.faces[ cell->faces[i] ].owner ] == cell ){
						int psize = mesh.faces[cell->faces[i]].points.size();
						for(int j=0; j<psize; ++j) points.push_back( mesh.faces[cell->faces[i]].points[j] );
					}
					else{
						int psize = mesh.faces[cell->faces[i]].points.size();
						for(int j=psize-1; j>=0; --j) points.push_back( mesh.faces[cell->faces[i]].points[j] );
					}
				}
				vector<int> cellPoints;
				for(int i=0; i<13; ++i) cellPoints.push_back(cell->points[i]);
				
				AMR.searchLoacationsCase0( faces, points, cellPoints );
				
				vector<vector<int>> cellPointsSet(19,vector<int>(2,0));
				
				AMR.sortCellPointsSet(cellPoints, cellPointsSet);
				
					
					
				for(int i=0; i<19; ++i){
					SEMO_Point tempPoint;
					tempPoint.x = 0.0;
					tempPoint.y = 0.0;
					tempPoint.z = 0.0;
					
					// cout << i << " " << cellPointsSet[i].size() << endl;
					for(auto& j : cellPointsSet[i]){
						tempPoint.x += mesh.points[j].x / (double)cellPointsSet[i].size();
						tempPoint.y += mesh.points[j].y / (double)cellPointsSet[i].size();
						tempPoint.z += mesh.points[j].z / (double)cellPointsSet[i].size();
					}
					
					newPoints.push_back(tempPoint);
					++newNPoints;
				}
				
				
				
				int idPoint = mesh.listPoints.size() + strPoints;
				
				vector<vector<int>> newFacesPointsSet(36,vector<int>(4,0));
				
				AMR.sortFacesPointsSet(cellPoints, idPoint, newFacesPointsSet);
				
				for(int i=0; i<2; ++i){
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] == cell ){
						for(int j=0; j<4; ++j) 
							reverse( &newFacesPointsSet[i*4+j][0], &newFacesPointsSet[i*4+j][4]);
					}
				}
				for(int i=2; i<6; ++i){
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] != cell ){
						for(int j=0; j<4; ++j)
							reverse( &newFacesPointsSet[i*4+j][0], &newFacesPointsSet[i*4+j][4]);
					}
				}
				
				vector<vector<int>> newFaceOwnNgbSet(6,vector<int>(4,0));
				
				int strOwnNgb = mesh.listCells.size() + strCells;
				
				AMR.sortOwnNgbOuterFaces(num, strOwnNgb, newFaceOwnNgbSet);

				
				for(int i=0; i<6; ++i){
					
					bool directionBasedOwner = false;
					if( &mesh.cells[ mesh.faces[ faces[i] ].owner ] == cell )
						directionBasedOwner = true;
					
					for(int j=0; j<4; ++j){
						SEMO_Face tempFace;
						for(int k=0; k<4; ++k){
							tempFace.points.push_back( newFacesPointsSet[i*4+j][k] );
						}
						if(directionBasedOwner){
							tempFace.owner = newFaceOwnNgbSet[i][j];
							tempFace.neighbour = mesh.faces[ faces[i] ].neighbour;
						}
						else{
							tempFace.owner = mesh.faces[ faces[i] ].owner;
							tempFace.neighbour = newFaceOwnNgbSet[i][j];
						}
						
						tempFace.setType( mesh.faces[ faces[i] ].getType() );
						
						int idFace = mesh.listFaces.size() + newNFaces;
						if(j==0) idFace = faces[i];
						
						if( mesh.faces[ faces[i] ].getType() == SEMO_Types::BOUNDARY_FACE ){
							idNewFaceBC[ iFaceBC[ faces[i] ] ].push_back( idFace );
						}
						else if( mesh.faces[ faces[i] ].getType() == SEMO_Types::PROCESSOR_FACE ){
							int procNum = neighbProcNo[faces[i]];
							idNewFaceProc[procNum].push_back( idFace );
						}
						
						idNewFaces.push_back( idFace );
						
						newFaces.push_back(tempFace);
						
						if(j!=0) ++newNFaces;
					}
				}
				
				vector<vector<int>> newFaceOwnSet(3,vector<int>(4,0));
				vector<vector<int>> newFaceNgbSet(3,vector<int>(4,0));
				
				AMR.sortOwnNgbInternalFaces(num, strOwnNgb, newFaceOwnSet, newFaceNgbSet);
				
				
				// cell internal faces
				for(int i=6; i<9; ++i){
					for(int j=0; j<4; ++j){
						SEMO_Face tempFace;
						for(int k=0; k<4; ++k)
							tempFace.points.push_back( newFacesPointsSet[i*4+j][k] );
						
						tempFace.owner = newFaceOwnSet[i-6][j];
						tempFace.neighbour = newFaceNgbSet[i-6][j];
						
						tempFace.setType(SEMO_Types::INTERNAL_FACE);
							
						idNewFaces.push_back(mesh.listFaces.size() + newNFaces);
						
						++newNFaces;
						
						newFaces.push_back(tempFace);
						
					} 
				}
				
				for(int i=0; i<8; ++i){
					
					int idCell = mesh.listCells.size() + newNCells;
					if(i==0) idCell = num;
					
					idNewCells.push_back(idCell);
					levelNewCells.push_back(mesh.cells[num].level + 1);
					
					if(i!=0) ++newNCells;
					
				}
				
				
				AMR.addPoints(mesh, newPoints);
				AMR.addFaces(mesh, idNewFaces, newFaces);
				AMR.addCells(mesh, idNewCells, levelNewCells);
				
				
			}
		}
		else if( cell->points.size() == 8+5+4 ){
			
		}
		else if( cell->points.size() == 8+5+5 ){
			
		}
		else if( cell->points.size() == 8+5+4+3 ){
			
		}
		else if( cell->points.size() == 8+5+4+4 ){
			
		}
		else if( cell->points.size() == 8+5+4+4+3 ){
			
		}
		else if( cell->points.size() == 8+5+4+3+3 ){
			
		}
		else if( cell->points.size() == 8+5+5+3+2+2 ){
			
		}
		else if( cell->points.size() == 8+5+5+3+2+2+1 ){
			
		}
		
		
		++num;
	}
	
	
	
	
	
	//===============================================
	
	
	
	
	// list setup
	// list points
	mesh.listPoints.clear();
	for(auto& point : mesh.points){
		mesh.listPoints.push_back(&point);
	}
	
	// list faces
	int faceNum=0;
	mesh.listFaces.clear();
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.listFaces.push_back(&face);
			++faceNum;
		}
	}
	for(int ibc=0; ibc<nbcs; ++ibc){
		int nFaces = 0;
		int strFace = mesh.boundary[ibc].startFace;
		int endFace = strFace + mesh.boundary[ibc].nFaces;
		mesh.boundary[ibc].startFace = faceNum;
		for(int i=0; i<idNewFaceBC[ibc].size(); ++i){
			int id = idNewFaceBC[ibc][i];
			mesh.listFaces.push_back(&mesh.faces[id]);
			++nFaces;
			++faceNum;
			
		}
		mesh.boundary[ibc].nFaces = nFaces;
	}
	if(size>1){
		for(int ip=0; ip<size; ++ip){
			int nFaces = 0;
			int ibc = nbcs+ip;
			int strFace = mesh.boundary[ibc].startFace;
			int endFace = strFace + mesh.boundary[ibc].nFaces;
			mesh.boundary[ibc].startFace = faceNum;
			for(int i=0; i<idNewFaceProc[ip].size(); ++i){
				int id = idNewFaceProc[ip][i];
				mesh.listFaces.push_back(&mesh.faces[id]);
				++nFaces;
				++faceNum;
			}
		}
	}
	
	
	
	// list Cells
	mesh.listCells.clear();
	for(auto& cell : mesh.cells){
		mesh.listCells.push_back(&cell);
	}
	
	
	
	// cell connection (cell's face)
	// delete faces
	for(auto cell : mesh.listCells){
		cell->faces.clear();
		cell->points.clear();
	}
	
	
	faceNum=0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
			mesh.cells[face.neighbour].faces.push_back(faceNum);
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
		}
		++faceNum;
	}
	
	
	
	
	// cell connection (cell's points)
	for(auto face : mesh.listFaces){
			// cout << face->points.size() << endl;
		if(face->getType() == SEMO_Types::INTERNAL_FACE){
			// cout << face->points.size() << endl;
			for(int i=0; i<face->points.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face->owner].points.size(); ++j){
					if(
					mesh.cells[face->owner].points[j] == 
					face->points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face->owner].points.push_back(face->points[i]);
				
				l=0;
				for(int j=0; j<mesh.cells[face->neighbour].points.size(); ++j){
					if(
					mesh.cells[face->neighbour].points[j] == 
					face->points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face->neighbour].points.push_back(face->points[i]);
			}
			
		}
		else if(face->getType() == SEMO_Types::BOUNDARY_FACE){
			for(int i=0; i<face->points.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face->owner].points.size(); ++j){
					if(
					mesh.cells[face->owner].points[j] == 
					face->points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face->owner].points.push_back(face->points[i]);
			}
		}
		else if(face->getType() == SEMO_Types::PROCESSOR_FACE){
			for(int i=0; i<face->points.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[face->owner].points.size(); ++j){
					if(
					mesh.cells[face->owner].points[j] == 
					face->points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[face->owner].points.push_back(face->points[i]);
			}
		}
	}
	
	
	
	
	
	
	
	
	
	
	cout << "┌-----------------------------------" << endl;
	cout << "| AMR : # of new child faces : " << idNewFaces.size() << endl;
	cout << "| AMR : # of new child cells : " << idNewCells.size() << endl;
	cout << "| AMR : # of created points : " << newNPoints << endl;
	cout << "| AMR : # of created faces : " << newNFaces << endl;
	cout << "| AMR : # of created cells : " << newNCells << endl;
	cout << "| AMR : # of modified faces : " << idNewFaces.size()-newNFaces << endl;
	cout << "| AMR : # of modified cells : " << idNewCells.size()-newNCells << endl;
	cout << "└-----------------------------------" << endl;
	
	
	
	mesh.saveFile("vtu");
		
		
		
		
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
}




void SEMO_Hexa_Oct_AMR::searchLoacationsCase0(
	vector<int> &faces, vector<int> &points, vector<int> &cellPoints){
	
	vector<int> strPoint;
	
	// for(auto& i : points){
		// cout << i << endl;
	// }
	// cout << endl;
	
	std::reverse( &points[0], &points[4] );
	
	// for(auto& i : points){
		// cout << i << endl;
	// }
	
	strPoint.push_back(0);
	strPoint.push_back(4);
	strPoint.push_back(8);
	strPoint.push_back(12);
	strPoint.push_back(16);
	strPoint.push_back(20);
	strPoint.push_back(24);
	
	// search face locations 
	// vector<int> nFace_Points(6,0);
	// vector<vector<int>> idFaces_Points(8,vector<int>(0,0));
	vector<int> iFirstPoint_Faces(6,-1);
	iFirstPoint_Faces[0] = points[0];
	vector<int> locationFace(6,0);
	locationFace[0] = 0;
	for(int i=1; i<strPoint.size()-1; ++i){
		int l=0;
		vector<int> nPointFirstFace(4,0);
		for(int j=strPoint[i]; j<strPoint[i+1]; ++j){
			for(int k=strPoint[0]; k<strPoint[1]; ++k){
				if( points[j] == points[k] ){
					// ++nFace_Points[i];
					++nPointFirstFace[k];
				}
			}
			if( nPointFirstFace[0] == 1 && nPointFirstFace[1] == 1 ){
				locationFace[i] = 1;
			}
			else if( nPointFirstFace[1] == 1 && nPointFirstFace[2] == 1 ){
				locationFace[i] = 2;
			}
			else if( nPointFirstFace[2] == 1 && nPointFirstFace[3] == 1 ){
				locationFace[i] = 3;
			}
			else if( nPointFirstFace[3] == 1 && nPointFirstFace[0] == 1 ){
				locationFace[i] = 4;
			}
		}
		if(locationFace[i]==0) locationFace[i] = 5;
	}
	
	
	// search point locations
	int orderPointsFace5[4];
	for(int i=1; i<strPoint.size()-1; ++i){
		if( locationFace[i]==1 ){
			int iFirstPoint;
			int iSecondPoint;
			for(int j=strPoint[i]; j<strPoint[i+1]; ++j){
				if( points[j] == points[0] ) iFirstPoint = j-strPoint[i];
				if( points[j] == points[1] ) iSecondPoint = j-strPoint[i];
			}
			// cout << iFirstPoint << " " << iSecondPoint << endl;
			if(
			(iFirstPoint==0 && iSecondPoint==1) ||
			(iFirstPoint==1 && iSecondPoint==2) ||
			(iFirstPoint==2 && iSecondPoint==3) ||
			(iFirstPoint==3 && iSecondPoint==0) 
			){
				int order;
				order = iFirstPoint-1;
				orderPointsFace5[0] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
				order = iFirstPoint-2;
				orderPointsFace5[1] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
			}
			else{
				int order;
				order = iFirstPoint+1;
				orderPointsFace5[0] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
				order = iFirstPoint+2;
				orderPointsFace5[1] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
			}
		}
			
		if( locationFace[i]==3 ){
			int iThirdPoint;
			int iFourthPoint;
			for(int j=strPoint[i]; j<strPoint[i+1]; ++j){
				// cout << points[j] << endl;
				if( points[j] == points[2] ) iThirdPoint = j-strPoint[i];
				if( points[j] == points[3] ) iFourthPoint = j-strPoint[i];
			}
			if(
			(iThirdPoint==0 && iFourthPoint==1) ||
			(iThirdPoint==1 && iFourthPoint==2) ||
			(iThirdPoint==2 && iFourthPoint==3) ||
			(iThirdPoint==3 && iFourthPoint==0) 
			){
				int order;
				order = iThirdPoint-1;
				orderPointsFace5[2] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
				order = iThirdPoint-2;
				orderPointsFace5[3] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
			}
			else{
				int order;
				order = iThirdPoint+1;
				orderPointsFace5[2] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
				order = iThirdPoint+2;
				orderPointsFace5[3] = strPoint[i] + ( order < 4 && order >= 0 ? order : order-4 );
			}
		}
	}
		
	for(int i=1; i<strPoint.size()-1; ++i){
		
		if( locationFace[i]==5 ){
			int sortedPoints[4];
			for(int j=0; j<4; ++j){
				sortedPoints[j] = points[orderPointsFace5[j]];
			}
			for(int j=0; j<4; ++j) points[strPoint[i]+j] = sortedPoints[j];
		}
	}
	
	// cell point locations
	// vector<int> cellPoints(8,0);
	for(int i=0; i<4; ++i) cellPoints[i] = points[i];
	
	for(int i=1; i<strPoint.size()-1; ++i){
		if( locationFace[i]==5 ){
			for(int j=0; j<4; ++j){
				cellPoints[4+j] = points[ strPoint[i]+j ];
			}
		}
	}
	
	// re order face locations
	vector<int> reorderFace(6,0);
	int num=0;
	for(auto& i : locationFace){
		reorderFace[num] = faces[i];
		++num;
	}
	num=0;
	for(auto& i : reorderFace){
		faces[num] = i;
		++num;
	}
	
	
		// cout << points[0] << endl;
		// cout << points[1] << endl;
		// cout << points[2] << endl;
		// cout << points[3] << endl;
	
		// cout << points[orderPointsFace5[0]] << endl;
		// cout << points[orderPointsFace5[1]] << endl;
		// cout << points[orderPointsFace5[2]] << endl;
		// cout << points[orderPointsFace5[3]] << endl;
	
	// for(auto& i : cellPoints){
		// cout << i << endl;
	// }
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
}






void SEMO_Hexa_Oct_AMR::sortCellPointsSet(
	vector<int> &cellPoints, vector<vector<int>> &cellPointsSet){

	cellPointsSet[0][0] = cellPoints[0]; cellPointsSet[0][1] = cellPoints[1];
	cellPointsSet[1][0] = cellPoints[1]; cellPointsSet[1][1] = cellPoints[2];
	cellPointsSet[2][0] = cellPoints[2]; cellPointsSet[2][1] = cellPoints[3];
	cellPointsSet[3][0] = cellPoints[3]; cellPointsSet[3][1] = cellPoints[0];
	
	cellPointsSet[4][0] = cellPoints[0]; cellPointsSet[4][1] = cellPoints[1];
	cellPointsSet[4].push_back(cellPoints[2]);
	cellPointsSet[4].push_back(cellPoints[3]);
	
	cellPointsSet[5][0] = cellPoints[0]; cellPointsSet[5][1] = cellPoints[4];
	
	cellPointsSet[6][0] = cellPoints[0]; cellPointsSet[6][1] = cellPoints[1];
	cellPointsSet[6].push_back(cellPoints[4]);
	cellPointsSet[6].push_back(cellPoints[5]);
	
	cellPointsSet[7][0] = cellPoints[1]; cellPointsSet[7][1] = cellPoints[5];
	
	cellPointsSet[8][0] = cellPoints[1]; cellPointsSet[8][1] = cellPoints[2];
	cellPointsSet[8].push_back(cellPoints[5]);
	cellPointsSet[8].push_back(cellPoints[6]);
	
	cellPointsSet[9][0] = cellPoints[2]; cellPointsSet[9][1] = cellPoints[6];
	
	cellPointsSet[10][0] = cellPoints[2]; cellPointsSet[10][1] = cellPoints[3];
	cellPointsSet[10].push_back(cellPoints[6]);
	cellPointsSet[10].push_back(cellPoints[7]);
	
	cellPointsSet[11][0] = cellPoints[3]; cellPointsSet[11][1] = cellPoints[7];
	
	cellPointsSet[12][0] = cellPoints[3]; cellPointsSet[12][1] = cellPoints[0];
	cellPointsSet[12].push_back(cellPoints[7]);
	cellPointsSet[12].push_back(cellPoints[4]);
	
	cellPointsSet[13][0] = cellPoints[0]; cellPointsSet[13][1] = cellPoints[1];
	cellPointsSet[13].push_back(cellPoints[2]);
	cellPointsSet[13].push_back(cellPoints[3]);
	cellPointsSet[13].push_back(cellPoints[4]);
	cellPointsSet[13].push_back(cellPoints[5]);
	cellPointsSet[13].push_back(cellPoints[6]);
	cellPointsSet[13].push_back(cellPoints[7]);
	
	cellPointsSet[14][0] = cellPoints[4]; cellPointsSet[14][1] = cellPoints[5];
	cellPointsSet[15][0] = cellPoints[5]; cellPointsSet[15][1] = cellPoints[6];
	cellPointsSet[16][0] = cellPoints[6]; cellPointsSet[16][1] = cellPoints[7];
	cellPointsSet[17][0] = cellPoints[7]; cellPointsSet[17][1] = cellPoints[4];
	
	cellPointsSet[18][0] = cellPoints[4]; cellPointsSet[18][1] = cellPoints[5];
	cellPointsSet[18].push_back(cellPoints[6]);
	cellPointsSet[18].push_back(cellPoints[7]);
}








void SEMO_Hexa_Oct_AMR::sortFacesPointsSet(
	vector<int> &cellPoints, int& idPoint, vector<vector<int>> &newFacesPointsSet){
			
	// face 0
	newFacesPointsSet[0][0] = cellPoints[0];
	newFacesPointsSet[0][1] = idPoint + 0;
	newFacesPointsSet[0][2] = idPoint + 4;
	newFacesPointsSet[0][3] = idPoint + 3;
	
	newFacesPointsSet[1][0] = idPoint + 0;
	newFacesPointsSet[1][1] = cellPoints[1];
	newFacesPointsSet[1][2] = idPoint + 1;
	newFacesPointsSet[1][3] = idPoint + 4;
	
	newFacesPointsSet[2][0] = idPoint + 4;
	newFacesPointsSet[2][1] = idPoint + 1;
	newFacesPointsSet[2][2] = cellPoints[2];
	newFacesPointsSet[2][3] = idPoint + 2;
	
	newFacesPointsSet[3][0] = idPoint + 3;
	newFacesPointsSet[3][1] = idPoint + 4;
	newFacesPointsSet[3][2] = idPoint + 2;
	newFacesPointsSet[3][3] = cellPoints[3];
	
	
	// face 1
	newFacesPointsSet[4][0] = cellPoints[4];
	newFacesPointsSet[4][1] = idPoint + 14;
	newFacesPointsSet[4][2] = idPoint + 6;
	newFacesPointsSet[4][3] = idPoint + 5;
	
	newFacesPointsSet[5][0] = idPoint + 14;
	newFacesPointsSet[5][1] = cellPoints[5];
	newFacesPointsSet[5][2] = idPoint + 7;
	newFacesPointsSet[5][3] = idPoint + 6;
	
	newFacesPointsSet[6][0] = idPoint + 6;
	newFacesPointsSet[6][1] = idPoint + 7;
	newFacesPointsSet[6][2] = cellPoints[1];
	newFacesPointsSet[6][3] = idPoint + 0;
	
	newFacesPointsSet[7][0] = idPoint + 5;
	newFacesPointsSet[7][1] = idPoint + 6;
	newFacesPointsSet[7][2] = idPoint + 0;
	newFacesPointsSet[7][3] = cellPoints[0];
	
	
	// face 2
	newFacesPointsSet[8][0] = cellPoints[6];
	newFacesPointsSet[8][1] = idPoint + 15;
	newFacesPointsSet[8][2] = idPoint + 8;
	newFacesPointsSet[8][3] = idPoint + 9;
	
	newFacesPointsSet[9][0] = idPoint + 15;
	newFacesPointsSet[9][1] = cellPoints[5];
	newFacesPointsSet[9][2] = idPoint + 7;
	newFacesPointsSet[9][3] = idPoint + 8;
	
	newFacesPointsSet[10][0] = idPoint + 8;
	newFacesPointsSet[10][1] = idPoint + 7;
	newFacesPointsSet[10][2] = cellPoints[1];
	newFacesPointsSet[10][3] = idPoint + 1;
	
	newFacesPointsSet[11][0] = idPoint + 9;
	newFacesPointsSet[11][1] = idPoint + 8;
	newFacesPointsSet[11][2] = idPoint + 1;
	newFacesPointsSet[11][3] = cellPoints[2];
	
	// face 3
	newFacesPointsSet[12][0] = cellPoints[7];
	newFacesPointsSet[12][1] = idPoint + 16;
	newFacesPointsSet[12][2] = idPoint + 10;
	newFacesPointsSet[12][3] = idPoint + 11;
	
	newFacesPointsSet[13][0] = idPoint + 16;
	newFacesPointsSet[13][1] = cellPoints[6];
	newFacesPointsSet[13][2] = idPoint + 9;
	newFacesPointsSet[13][3] = idPoint + 10;
	
	newFacesPointsSet[14][0] = idPoint + 10;
	newFacesPointsSet[14][1] = idPoint + 9;
	newFacesPointsSet[14][2] = cellPoints[2];
	newFacesPointsSet[14][3] = idPoint + 2;
	
	newFacesPointsSet[15][0] = idPoint + 11;
	newFacesPointsSet[15][1] = idPoint + 10;
	newFacesPointsSet[15][2] = idPoint + 2;
	newFacesPointsSet[15][3] = cellPoints[3];
	
	// face 4
	newFacesPointsSet[16][0] = cellPoints[7];
	newFacesPointsSet[16][1] = idPoint + 17;
	newFacesPointsSet[16][2] = idPoint + 12;
	newFacesPointsSet[16][3] = idPoint + 11;
	
	newFacesPointsSet[17][0] = idPoint + 17;
	newFacesPointsSet[17][1] = cellPoints[4];
	newFacesPointsSet[17][2] = idPoint + 5;
	newFacesPointsSet[17][3] = idPoint + 12;
	
	newFacesPointsSet[18][0] = idPoint + 12;
	newFacesPointsSet[18][1] = idPoint + 5;
	newFacesPointsSet[18][2] = cellPoints[0];
	newFacesPointsSet[18][3] = idPoint + 3;
	
	newFacesPointsSet[19][0] = idPoint + 11;
	newFacesPointsSet[19][1] = idPoint + 12;
	newFacesPointsSet[19][2] = idPoint + 3;
	newFacesPointsSet[19][3] = cellPoints[3];
	
	// face 5
	newFacesPointsSet[20][0] = cellPoints[4];
	newFacesPointsSet[20][1] = idPoint + 14;
	newFacesPointsSet[20][2] = idPoint + 18;
	newFacesPointsSet[20][3] = idPoint + 17;
	
	newFacesPointsSet[21][0] = idPoint + 14;
	newFacesPointsSet[21][1] = cellPoints[5];
	newFacesPointsSet[21][2] = idPoint + 15;
	newFacesPointsSet[21][3] = idPoint + 18;
	
	newFacesPointsSet[22][0] = idPoint + 18;
	newFacesPointsSet[22][1] = idPoint + 15;
	newFacesPointsSet[22][2] = cellPoints[6];
	newFacesPointsSet[22][3] = idPoint + 16;
	
	newFacesPointsSet[23][0] = idPoint + 17;
	newFacesPointsSet[23][1] = idPoint + 18;
	newFacesPointsSet[23][2] = idPoint + 16;
	newFacesPointsSet[23][3] = cellPoints[7];
	
	// new face 6
	newFacesPointsSet[24][0] = idPoint + 5;
	newFacesPointsSet[24][1] = idPoint + 6;
	newFacesPointsSet[24][2] = idPoint + 13;
	newFacesPointsSet[24][3] = idPoint + 12;
	
	newFacesPointsSet[25][0] = idPoint + 6;
	newFacesPointsSet[25][1] = idPoint + 7;
	newFacesPointsSet[25][2] = idPoint + 8;
	newFacesPointsSet[25][3] = idPoint + 13;
	
	newFacesPointsSet[26][0] = idPoint + 13;
	newFacesPointsSet[26][1] = idPoint + 8;
	newFacesPointsSet[26][2] = idPoint + 9;
	newFacesPointsSet[26][3] = idPoint + 10;
	
	newFacesPointsSet[27][0] = idPoint + 12;
	newFacesPointsSet[27][1] = idPoint + 13;
	newFacesPointsSet[27][2] = idPoint + 10;
	newFacesPointsSet[27][3] = idPoint + 11;
	
	// new face 7
	newFacesPointsSet[28][0] = idPoint + 17;
	newFacesPointsSet[28][1] = idPoint + 18;
	newFacesPointsSet[28][2] = idPoint + 13;
	newFacesPointsSet[28][3] = idPoint + 12;
	
	newFacesPointsSet[29][0] = idPoint + 18;
	newFacesPointsSet[29][1] = idPoint + 15;
	newFacesPointsSet[29][2] = idPoint + 8;
	newFacesPointsSet[29][3] = idPoint + 13;
	
	newFacesPointsSet[30][0] = idPoint + 13;
	newFacesPointsSet[30][1] = idPoint + 8;
	newFacesPointsSet[30][2] = idPoint + 1;
	newFacesPointsSet[30][3] = idPoint + 4;
	
	newFacesPointsSet[31][0] = idPoint + 12;
	newFacesPointsSet[31][1] = idPoint + 13;
	newFacesPointsSet[31][2] = idPoint + 4;
	newFacesPointsSet[31][3] = idPoint + 3;
	
	// new face 8
	newFacesPointsSet[32][0] = idPoint + 14;
	newFacesPointsSet[32][1] = idPoint + 6;
	newFacesPointsSet[32][2] = idPoint + 13;
	newFacesPointsSet[32][3] = idPoint + 18;
	
	newFacesPointsSet[33][0] = idPoint + 6;
	newFacesPointsSet[33][1] = idPoint + 0;
	newFacesPointsSet[33][2] = idPoint + 4;
	newFacesPointsSet[33][3] = idPoint + 13;
	
	newFacesPointsSet[34][0] = idPoint + 13;
	newFacesPointsSet[34][1] = idPoint + 4;
	newFacesPointsSet[34][2] = idPoint + 2;
	newFacesPointsSet[34][3] = idPoint + 10;
	
	newFacesPointsSet[35][0] = idPoint + 18;
	newFacesPointsSet[35][1] = idPoint + 13;
	newFacesPointsSet[35][2] = idPoint + 10;
	newFacesPointsSet[35][3] = idPoint + 16;
}






void SEMO_Hexa_Oct_AMR::sortOwnNgbOuterFaces(
	int& num, int& strOwnNgb, vector<vector<int>> &newFaceOwnNgbSet){
			

	newFaceOwnNgbSet[0][0] = num;
	newFaceOwnNgbSet[0][1] = strOwnNgb + 0;
	newFaceOwnNgbSet[0][2] = strOwnNgb + 1;
	newFaceOwnNgbSet[0][3] = strOwnNgb + 2;
	
	newFaceOwnNgbSet[1][0] = strOwnNgb + 3;
	newFaceOwnNgbSet[1][1] = strOwnNgb + 4;
	newFaceOwnNgbSet[1][2] = strOwnNgb + 0;
	newFaceOwnNgbSet[1][3] = num;
	
	newFaceOwnNgbSet[2][0] = strOwnNgb + 5;
	newFaceOwnNgbSet[2][1] = strOwnNgb + 4;
	newFaceOwnNgbSet[2][2] = strOwnNgb + 0;
	newFaceOwnNgbSet[2][3] = strOwnNgb + 1;

	newFaceOwnNgbSet[3][0] = strOwnNgb + 6;
	newFaceOwnNgbSet[3][1] = strOwnNgb + 5;
	newFaceOwnNgbSet[3][2] = strOwnNgb + 1;
	newFaceOwnNgbSet[3][3] = strOwnNgb + 2;
	
	newFaceOwnNgbSet[4][0] = strOwnNgb + 6;
	newFaceOwnNgbSet[4][1] = strOwnNgb + 3;
	newFaceOwnNgbSet[4][2] = num;
	newFaceOwnNgbSet[4][3] = strOwnNgb + 2;
	
	newFaceOwnNgbSet[5][0] = strOwnNgb + 3;
	newFaceOwnNgbSet[5][1] = strOwnNgb + 4;
	newFaceOwnNgbSet[5][2] = strOwnNgb + 5;
	newFaceOwnNgbSet[5][3] = strOwnNgb + 6;
}








void SEMO_Hexa_Oct_AMR::sortOwnNgbInternalFaces( 
	int& num, int& strOwnNgb, 
	vector<vector<int>> &newFaceOwnSet, vector<vector<int>> &newFaceNgbSet){
	
	newFaceOwnSet[0][0] = num;           newFaceNgbSet[0][0] = strOwnNgb + 3; 
	newFaceOwnSet[0][1] = strOwnNgb + 0; newFaceNgbSet[0][1] = strOwnNgb + 4; 
	newFaceOwnSet[0][2] = strOwnNgb + 1; newFaceNgbSet[0][2] = strOwnNgb + 5; 
	newFaceOwnSet[0][3] = strOwnNgb + 2; newFaceNgbSet[0][3] = strOwnNgb + 6; 
	
	newFaceOwnSet[1][0] = strOwnNgb + 3; newFaceNgbSet[1][0] = strOwnNgb + 6; 
	newFaceOwnSet[1][1] = strOwnNgb + 4; newFaceNgbSet[1][1] = strOwnNgb + 5; 
	newFaceOwnSet[1][2] = strOwnNgb + 0; newFaceNgbSet[1][2] = strOwnNgb + 1; 
	newFaceOwnSet[1][3] = num;           newFaceNgbSet[1][3] = strOwnNgb + 2; 
	
	newFaceOwnSet[2][0] = strOwnNgb + 3; newFaceNgbSet[2][0] = strOwnNgb + 4; 
	newFaceOwnSet[2][1] = num;           newFaceNgbSet[2][1] = strOwnNgb + 0; 
	newFaceOwnSet[2][2] = strOwnNgb + 2; newFaceNgbSet[2][2] = strOwnNgb + 1; 
	newFaceOwnSet[2][3] = strOwnNgb + 6; newFaceNgbSet[2][3] = strOwnNgb + 5; 
}






void SEMO_Hexa_Oct_AMR::addPoints( 
	SEMO_Mesh_Builder& mesh, vector<SEMO_Point>& newPoints){

	// add points
	for(auto& point : newPoints){
		mesh.addPoint();
		mesh.points.back().x = point.x;
		mesh.points.back().y = point.y;
		mesh.points.back().z = point.z;
	}

}

void SEMO_Hexa_Oct_AMR::addFaces( 
	SEMO_Mesh_Builder& mesh, vector<int>& idNewFaces, vector<SEMO_Face>& newFaces){
	
	// add faces
	int num=0;
	int nfaces = mesh.faces.size();
	for(auto& i : idNewFaces){
		
		if( i < nfaces ){
			mesh.faces[i].points.clear();
			for(int j=0; j<newFaces[num].points.size(); ++j){
				mesh.faces[i].points.push_back( newFaces[num].points[j] );
			}
			mesh.faces[i].owner = newFaces[num].owner;
			mesh.faces[i].neighbour = newFaces[num].neighbour;
			mesh.faces[i].setType(newFaces[num].getType());
		}
		else{
			mesh.addFace();
			for(int j=0; j<newFaces[num].points.size(); ++j){
				mesh.faces.back().points.push_back( newFaces[num].points[j] );
			}
			mesh.faces.back().owner = newFaces[num].owner;
			mesh.faces.back().neighbour = newFaces[num].neighbour;
			mesh.faces.back().setType(newFaces[num].getType());
		}
		
		++num;
	}
					
}



void SEMO_Hexa_Oct_AMR::addCells( 
	SEMO_Mesh_Builder& mesh, vector<int>& idNewCells, vector<int>& levelNewCells){
	
	// add cells
	int num=0;
	int ncells = mesh.cells.size();
	for(auto& i : idNewCells){
		if( i <ncells ){
			mesh.cells[i].level = levelNewCells[num];
		}
		else{
			mesh.addCell();
			mesh.cells.back().level = levelNewCells[num];
		}
		++num;
	}
		
}