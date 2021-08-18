#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>

#include "build.h"
#include "mpi.h"
#include "polyAMR.h"
#include "geometric.h"




void SEMO_Poly_AMR_Builder::polyRefine(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
	if(rank==0) cout << "----Refine start-----" << endl;
	
	
	SEMO_Mesh_Geometric geometric;
	
	int nTotalFaceLRVar = mesh.faces[0].varL.size();
	int nTotalFaceVar = mesh.faces[0].var.size();
	
	
	//====================================================
	// Refine 되는 셀 & 면 조사
	
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	
	
	SEMO_Poly_Mesh_Refine refineMesh;
	// refine.iFaces.resize(mesh.faces.size(),vector<int>(0,0));
	// refine.iCells.resize(mesh.cells.size(),vector<int>(0,0));
	
	vector<bool> boolCellRefine(mesh.cells.size(),false);
	for(int i=0; i<mesh.cells.size(); ++i){
		if( distr(eng) > 0.1 ){
			boolCellRefine[i] = true;
		}
		
		
		
			boolCellRefine[i] = true;
		
		
		
	} 
	
	vector<int> cLevel_recv;
	vector<int> cRefine_recv;
	this->mpiLevelRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	this->restrictCellRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);

	vector<int> edgesPoint0;
	vector<int> edgesPoint1;
	vector<vector<int>> facesEdges;
	vector<vector<int>> edgesFaces;
	vector<int> edgesLevel;
	
	this->createEdges(mesh, edgesPoint0, edgesPoint1, facesEdges, edgesFaces, edgesLevel);
	
	
	
	//====================================================
	// AMR : 셀 기반 AMR
	
	// vector<int> childCells(mesh.cells.size(),0);
	// vector<int> childFaces(mesh.faces.size(),0);
	vector<bool> boolEdgeRefine(edgesPoint0.size(),false);
	vector<bool> boolFaceRefine(mesh.faces.size(),false);
	vector<int> edgesCenterPointNumber(edgesPoint0.size(),-1);
	
	vector<int> outFacesParentNumber;
	vector<vector<vector<int>>> outFacesPoints;
	vector<vector<int>> outFacesOwner;
	vector<vector<int>> outFacesNeighbour;
	
	// vector<int> intFacesParentNumber;
	// vector<vector<vector<int>>> intFacesPoints;
	// vector<vector<int>> intFacesOwner;
	// vector<vector<int>> intFacesNeighbour;
	
	vector<int> childOutFacesNumber(mesh.faces.size(),-1);
	
	
	
	cout << "1" << endl;
	
	
	int newPointNumber = mesh.points.size();
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		
		if(boolCellRefine[i]==false) continue;
		
		int cellLevel = mesh.cells[i].level;
		refineMesh.intFacesParentCell.push_back(i);
		
		vector<int> cellVertex;
		this->searchOriginalPoints(mesh, mesh.cells[i].points, cellLevel, cellVertex);
		
		// childCells[i] = cellVertex.size();
		
		// int cellCenterPoint = newPointNumber++;
		// refineMesh.addCenterPoint(mesh, cellVertex, cellLevel);
		int cellCenterPoint = mesh.points.size();
		this->addCenterPoint(mesh, cellVertex, cellLevel);
		
		// face vertex points sort by cell point number
		map<int, int> faceVertexCellOrder;
		for(int j=0; j<cell.points.size(); ++j){
			int vertexPoint = cell.points[j];
			if(mesh.points[vertexPoint].level > cellLevel) continue;
			faceVertexCellOrder.insert(pair<int,int>(vertexPoint,j));
		}
		
		
		map<int, vector<int>> intFacesPoints;
		map<int, int> intFacesOwner;
		map<int, int> intFacesNeighbour;
		for(auto& j : mesh.cells[i].faces){
			
			auto& face = mesh.faces[j];
			int faceLevel = face.level;
			
			if( faceLevel == cellLevel && boolFaceRefine[j] == false ){
				boolFaceRefine[j] = true;
				
				vector<int> faceVertex;
				vector<int> edgeCenterPoints;
				vector<vector<int>> subFaceEdgePoints;
				
	// cout << "1" << endl;
				refineMesh.separateEdgesPoints(mesh, faceLevel, face.points, 
					facesEdges[j], boolEdgeRefine, edgesCenterPointNumber, edgesLevel, 
					faceVertex, edgeCenterPoints, subFaceEdgePoints, newPointNumber);

	for(auto& k : faceVertex){
		cout << "faceVertex " << k << endl;
	}
	for(auto& k : edgeCenterPoints){
		cout << "edgeCenterPoints " << k << endl;
	}
	// for(auto& k : faceVertex){
		// cout << "faceVertex " << k << endl;
	// }
				// face center point
				// int faceCenterPoint = newPointNumber++;
				int faceCenterPoint = mesh.points.size();
				this->addCenterPoint(mesh, faceVertex, faceLevel);
		
				int nChildFaces = faceVertex.size();
				
				// create new faces, face's points, face's owner neighbour
				refineMesh.outFacesParentFace.push_back(j);
				vector<SEMO_Face> childOutFaces(j,SEMO_Face());
				for(int k=0; k<nChildFaces; ++k){
					int vertexPoint = faceVertex[k];
					
					// vector<int> tmpPoints;
					for(auto& l : subFaceEdgePoints[k]){
						childOutFaces[k].points.push_back(l);
					}
					childOutFaces[k].points.push_back(edgeCenterPoints[k]);
					childOutFaces[k].points.push_back(faceCenterPoint);
					
					for(auto& l : subFaceEdgePoints[ 
						k-1<0 ? subFaceEdgePoints.size()-1 : k-1 ]){
						childOutFaces[k].points.push_back(l);
					}
					
					if(face.owner == i){
						childOutFaces[k].owner = faceVertexCellOrder[vertexPoint];
						childOutFaces[k].neighbour = face.neighbour;
					}
					else{
						childOutFaces[k].owner = face.owner;
						childOutFaces[k].neighbour = faceVertexCellOrder[vertexPoint];
					}
				}
				childOutFacesNumber[j] = refineMesh.outFaces.size();
				refineMesh.outFaces.push_back(childOutFaces);
				
				
				// create new cell internal faces
				for(int k=0; k<nChildFaces; ++k){
					int vertexPoint0 = faceVertex[k];
					int vertexPoint1 = faceVertex[ k+1 == nChildFaces ? 0 : k+1 ];
					int edgeCentPoint = edgeCenterPoints[k];

					if( intFacesPoints.find(edgeCentPoint) == intFacesPoints.end() ){
						
						intFacesPoints.insert(pair<int,vector<int>>(edgeCentPoint,vector<int>()));
						intFacesPoints[edgeCentPoint].push_back(edgeCentPoint);
						intFacesPoints[edgeCentPoint].push_back(faceCenterPoint);
						intFacesPoints[edgeCentPoint].push_back(cellCenterPoint);
						
						intFacesOwner.insert(pair<int,int>(edgeCentPoint,vertexPoint0));
						intFacesNeighbour.insert(pair<int,int>(edgeCentPoint,vertexPoint1));
						
					}
					else{
						intFacesPoints[edgeCentPoint].push_back(faceCenterPoint);
					}
					
					if(intFacesPoints[edgeCentPoint].size() > 4){
						cout << "| #Error 4" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
				}
			}
			else if( faceLevel == cellLevel && boolFaceRefine[j] == true ){
				
				
				vector<int> faceVertex;
				int tmpNum = childOutFacesNumber[j];
				for(auto& childOutFace : refineMesh.outFaces[tmpNum]){
					faceVertex.push_back(childOutFace.points[0]);
				}
				
				int nChildFaces = faceVertex.size();
				
				// face - cell connection
				for(int k=0; k<refineMesh.outFaces[tmpNum].size(); ++k){
					auto& childOutFace = refineMesh.outFaces[tmpNum][k];
					int vertexPoint = faceVertex[k];
					
					if(face.owner == i){
						childOutFace.owner = faceVertexCellOrder[vertexPoint];
						childOutFace.neighbour = face.neighbour;
					}
					else{
						childOutFace.owner = face.owner;
						childOutFace.neighbour = faceVertexCellOrder[vertexPoint];
					}
				}
				
				// create cell internal face and face - cell connection
				
				vector<vector<int>> eachChildFaceVertex(refineMesh.outFaces.size(),vector<int>());
				vector<int> edgeCenterPoints(refineMesh.outFaces.size(),-1);
				for(int k=0; k<refineMesh.outFaces[tmpNum].size(); ++k){
					for(auto& childPoints : refineMesh.outFaces[tmpNum][k].points){
						if( mesh.points[childPoints].level <= faceLevel + 1 ){
							eachChildFaceVertex[k].push_back(childPoints);
						}
					}
					
					edgeCenterPoints[k] = eachChildFaceVertex[k][1];
				}
				
				int faceCenterPoint = eachChildFaceVertex[0][2];
				
				for(int k=0; k<nChildFaces; ++k){
					int vertexPoint0 = faceVertex[k];
					int vertexPoint1 = faceVertex[ k+1 == nChildFaces ? 0 : k+1 ];
					int edgeCentPoint = edgeCenterPoints[k];

					if( intFacesPoints.find(edgeCentPoint) == intFacesPoints.end() ){
						
						intFacesPoints.insert(pair<int,vector<int>>(edgeCentPoint,vector<int>()));
						intFacesPoints[edgeCentPoint].push_back(edgeCentPoint);
						intFacesPoints[edgeCentPoint].push_back(faceCenterPoint);
						intFacesPoints[edgeCentPoint].push_back(cellCenterPoint);
						
						intFacesOwner.insert(pair<int,int>(edgeCentPoint,vertexPoint0));
						intFacesNeighbour.insert(pair<int,int>(edgeCentPoint,vertexPoint1));
						
					}
					else{
						intFacesPoints[edgeCentPoint].push_back(faceCenterPoint);
					}
					
					if(intFacesPoints[edgeCentPoint].size() > 4){
						cout << "| #Error 4" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
				}
				
				
			}
			else if( faceLevel == cellLevel + 1 && boolFaceRefine[j] == false ){
				
				
				vector<int> faceVertex;
				vector<int> childFaceVertex;
				for(auto& k : face.points){
					if( mesh.points[k].level <= faceLevel-1 ){
						faceVertex.push_back(k);
					}
					if( mesh.points[k].level <= faceLevel ){
						childFaceVertex.push_back(k);
					}
				}
				
				// face vertex points sort by cell point number
				int nChildFaces = faceVertex.size();
				
				int vertexPoint = faceVertex[0];
				
				// face - cell connection
				if(face.owner == i){
					face.owner = faceVertexCellOrder[vertexPoint];
				}
				else{
					face.neighbour = faceVertexCellOrder[vertexPoint];
				}
			
				// create cell internal face and face - cell connection
				int faceCenterPoint = childFaceVertex[2];
				
				int vertexPoint0 = faceVertex[0];
				int edgeCentPoint0 = childFaceVertex[1];
				int edgeCentPoint1 = childFaceVertex[3];

				if( intFacesPoints.find(edgeCentPoint0) == intFacesPoints.end() ){
					
					intFacesPoints.insert(pair<int,vector<int>>(edgeCentPoint0,vector<int>()));
					intFacesPoints[edgeCentPoint0].push_back(edgeCentPoint0);
					intFacesPoints[edgeCentPoint0].push_back(faceCenterPoint);
					intFacesPoints[edgeCentPoint0].push_back(cellCenterPoint);
					
					if( intFacesOwner.find(edgeCentPoint0) == intFacesOwner.end() ){
						intFacesOwner.insert(pair<int,int>(edgeCentPoint0,vertexPoint0));
					}
					else{
						intFacesNeighbour.insert(pair<int,int>(edgeCentPoint0,vertexPoint0));
					}
					
				}
				else{
					if(std::find(intFacesPoints[edgeCentPoint0].begin(), 
						intFacesPoints[edgeCentPoint0].end(), 
						faceCenterPoint) == 
						intFacesPoints[edgeCentPoint0].end() ){
							
						intFacesPoints[edgeCentPoint0].push_back(faceCenterPoint);
						
					}
				}
				if(intFacesPoints[edgeCentPoint0].size() > 4){
					cout << "| #Error 4" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				

				if( intFacesPoints.find(edgeCentPoint1) == intFacesPoints.end() ){
					
					intFacesPoints.insert(pair<int,vector<int>>(edgeCentPoint1,vector<int>()));
					intFacesPoints[edgeCentPoint1].push_back(edgeCentPoint1);
					intFacesPoints[edgeCentPoint1].push_back(faceCenterPoint);
					intFacesPoints[edgeCentPoint1].push_back(cellCenterPoint);
					
					if( intFacesOwner.find(edgeCentPoint1) == intFacesOwner.end() ){
						intFacesOwner.insert(pair<int,int>(edgeCentPoint1,vertexPoint0));
					}
					else{
						intFacesNeighbour.insert(pair<int,int>(edgeCentPoint1,vertexPoint0));
					}
					
				}
				else{
					if(std::find(intFacesPoints[edgeCentPoint1].begin(), 
						intFacesPoints[edgeCentPoint1].end(), 
						faceCenterPoint) == 
						intFacesPoints[edgeCentPoint1].end() ){
							
						intFacesPoints[edgeCentPoint1].push_back(faceCenterPoint);
						
					}
				}
				if(intFacesPoints[edgeCentPoint1].size() > 4){
					cout << "| #Error 4" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
			}
		}
		
		
		refineMesh.addIntFace(intFacesPoints.size());
		int intFaceNum = refineMesh.intFaces.size()-1;
		int tmpNum = 0;
		for(auto iter : intFacesOwner){
			int edgeCent = iter.first;
			for(auto& k : intFacesPoints[edgeCent]){
				refineMesh.intFaces[intFaceNum][tmpNum].points.push_back(k);
			}
			refineMesh.intFaces[intFaceNum][tmpNum].owner = intFacesOwner[edgeCent];
			refineMesh.intFaces[intFaceNum][tmpNum].neighbour = intFacesNeighbour[edgeCent];
			
			++tmpNum;
		}
		
		
	}
	
	
	cout << "2" << endl;
	

	// mesh.check();
	
	mesh.buildCells();
	
	// mesh.setFaceTypes();

	mesh.buildLists();
	
	mesh.connectCelltoFaces();
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// // int minLevel = 999999999;
		// // for(int j=0; j<cell.faces.size(); ++j){
			// // auto& face = mesh.faces[cell.faces[j]];
			// // minLevel = min(minLevel,face.level);
		// // }
		// // cell.level = minLevel;
		
		// cell.level = cellLeveling[i];
		// cell.group = cellGroupping[i];
		
		
		// // cout << "CELL VAR " << cellVariables[i].size() << endl;
		// for(auto& k : cellVariables[i]){
			// cell.var.push_back(k);
		// }
		
		
		// // cout << " G : " << cellGroupping[i] << " " << cellGroupping.size() << " " << mesh.cells.size() << endl;
	// }
	
	mesh.connectCelltoPoints();
	
	// set processor face counts
	mesh.setCountsProcFaces();
	
	// set processor face displacements
	mesh.setDisplsProcFaces(); 
		
	// mesh.informations();
	
	
	
	SEMO_Mesh_Save save;
	string tmpFile = "./Rf" + to_string(iter);
	// string tmpFile = "./";
	save.vtu(tmpFile, rank, mesh);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// cout << "SIZE : " << cellGroupping.size() << " " << mesh.cells.size() << endl;
		
	// cout << "CELL" << endl;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// cout << i << " " << cell.level << " " << cell.group << endl;
	// }
	
	// cout << "FACE" << endl;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// cout << i << " " << face.level << " " << face.group << endl;
	// }
		

	
	
	// //====================================================
	// // AMR : 엣지 중심 점 생성
	// vector<vector<double>> newEdgeCenterPointsXYZ;
	// vector<int> newEdgeCenterPointsLevel;
	// vector<int> newEdgeCenterPointNum(edgesPoints.size(),-1);
	// int newPointNum = mesh.points.size();
	// int strEdgeCenterPointNum = newPointNum;
	// for(int i=0; i<edgesPoints.size(); ++i){
		
		// if(boolEdgeRefine[i] == false) continue;
		
		// int point0 = edgesPoints[i][0];
		// int point1 = edgesPoints[i][1];
		
		// double x0 = mesh.points[point0].x;
		// double y0 = mesh.points[point0].y;
		// double z0 = mesh.points[point0].z;
		
		// double x1 = mesh.points[point1].x;
		// double y1 = mesh.points[point1].y;
		// double z1 = mesh.points[point1].z;
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(0.5*(x0+x1));
		// tmpXYZ.push_back(0.5*(y0+y1));
		// tmpXYZ.push_back(0.5*(z0+z1));
		
		// newEdgeCenterPointsXYZ.push_back(tmpXYZ);
		// // newEdgeCenterPointsLevel.push_back(mesh.points[point0].level+1);
		// newEdgeCenterPointsLevel.push_back(edgeLevel[i]+1);
		
		// newEdgeCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	// //====================================================
	// // AMR : 면 중심 점 생성
	// vector<vector<double>> newFaceCenterPointsXYZ;
	// vector<int> newFaceCenterPointsLevel;
	// vector<int> newFaceCenterPointNum(mesh.faces.size(),-1);
	// int strFaceCenterPointNum = newPointNum;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// if(boolFaceRefine[i] == false) continue;
		
		// auto& face = mesh.faces[i];
		// double xCenter = 0.0;
		// double yCenter = 0.0;
		// double zCenter = 0.0;
		// double tmpNum = 0.0;
		// for(auto& j : face.points){
			// auto& point = mesh.points[j];
			// if(point.level>face.level) continue;
			// xCenter += point.x;
			// yCenter += point.y;
			// zCenter += point.z;
			// tmpNum = tmpNum+1.0;
		// }
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(xCenter/tmpNum);
		// tmpXYZ.push_back(yCenter/tmpNum);
		// tmpXYZ.push_back(zCenter/tmpNum);
		
		// newFaceCenterPointsXYZ.push_back(tmpXYZ);
		// newFaceCenterPointsLevel.push_back(face.level+1);
		
		// newFaceCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 중심 점 생성
	// vector<vector<double>> newCellCenterPointsXYZ;
	// vector<int> newCellCenterPointsLevel;
	// vector<int> newCellCenterPointNum(mesh.cells.size(),-1);
	// int strCellCenterPointNum = newPointNum;
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) continue;
		
		// auto& cell = mesh.cells[i];
		// double xCenter = 0.0;
		// double yCenter = 0.0;
		// double zCenter = 0.0;
		// double tmpNum = 0.0;
		// for(auto& j : cell.points){
			// auto& point = mesh.points[j];
			// if(point.level>cell.level) continue;
			// xCenter += point.x;
			// yCenter += point.y;
			// zCenter += point.z;
			// tmpNum = tmpNum+1.0;
		// }
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(xCenter/tmpNum);
		// tmpXYZ.push_back(yCenter/tmpNum);
		// tmpXYZ.push_back(zCenter/tmpNum);
		
		// newCellCenterPointsXYZ.push_back(tmpXYZ);
		// newCellCenterPointsLevel.push_back(cell.level+1);
		
		// newCellCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	
	// if(rank==0) cout << "point add start" << endl;
	

	// //====================================================
	// // AMR : 포인트 추가 (엣지 중심 -> 면 중심 -> 셀 중심)
	// newPointNum = 0;
	// for(auto& i : newEdgeCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newEdgeCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newEdgeCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newEdgeCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newEdgeCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	// newPointNum = 0;
	// for(auto& i : newFaceCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newFaceCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newFaceCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newFaceCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newFaceCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	// newPointNum = 0;
	// for(auto& i : newCellCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newCellCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newCellCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newCellCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newCellCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	
	
	
	// //====================================================
	// // AMR : 셀 포인트로 넘버링된 "내부"면 - 셀 연결
	// // AMR : 셀 포인트로 넘버링된 "내부"면 - 포인트 연결
	
	// // 내부 면
	
	
	// vector<int> newCellNumbering(mesh.cells.size(),-1);
	// vector<vector<int>> cellVertexPoints;
	// vector<vector<int>> cellEdges(mesh.cells.size(),vector<int>(0,0));
	// int newCellNum = 0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) {
			
			// newCellNumbering[i] = newCellNum;
			
			// newCellNum += 1;
			
			// vector<int> tmpPoints;
			// cellVertexPoints.push_back(tmpPoints);
			
			// continue;
		// }
		
		// auto& cell = mesh.cells[i];
		// int level = cell.level;
		
		// vector<int> tmpPoints;
		// for(auto& j : cell.points){
			// if( mesh.points[j].level <= level ){
				// tmpPoints.push_back(j);
			// }
		// }
		// cellVertexPoints.push_back(tmpPoints);
		
		// vector<int> tmpCellEdges;
		// for(auto& j : cell.faces){
			// for(auto& k : facesEdges[j]){
				// if ( std::find( tmpCellEdges.begin(), tmpCellEdges.end(), k ) 
					// == tmpCellEdges.end() ) {
					// tmpCellEdges.push_back(k);
				// }
			// }
			
		// }
		// cellEdges[i] = tmpCellEdges;
		
		// newCellNumbering[i] = newCellNum;
		
		// newCellNum += tmpPoints.size();
		
	// }
	
	
	// vector<vector<int>> faceVertexPoints;
	// vector<vector<int>> faceEdgeCenterPoints;
	// // vector<int> faceCenterPoints;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];
		// int level = face.level;
		
		// vector<int> tmpPoints0;
		// for(auto& j : face.points){
			// if( mesh.points[j].level <= level ){
				// tmpPoints0.push_back(j);
			// }
		// }
		// faceVertexPoints.push_back(tmpPoints0);
			
		// if(boolFaceRefine[i] == false) {
			
			// vector<int> tmpPoints1;
			// faceEdgeCenterPoints.push_back(tmpPoints1);
			
			// // faceCenterPoints.push_back(newFaceCenterPointNum[i]);
			
			// continue;
		// }
		
		// vector<int> tmpPoints1;
		// int tmpNum = 0;
		// for(auto& j : face.points){
			// if( mesh.points[j].level == level+1 ){
				// tmpPoints1.push_back(j);
			// }
			// int edgeNum = facesEdges[i][tmpNum];
			// if( 
			// boolEdgeRefine[edgeNum] == true && 
			// edgeLevel[edgeNum] == level ){
				// tmpPoints1.push_back( newEdgeCenterPointNum[edgeNum] );
			// }
			// ++tmpNum;
		// }
		// faceEdgeCenterPoints.push_back(tmpPoints1);
		
		// // faceCenterPoints.push_back(newFaceCenterPointNum[i]);
		
	// }
	
	
	
	// vector<int> newIntFaceOwner;
	// vector<int> newIntFaceNeighbour;
	// vector<vector<int>> newIntFacePoints;
	// vector<vector<int>> newCellIntEdgeCenterPoints(mesh.cells.size(),vector<int>(0,0));
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) continue;
		
		// auto& cell = mesh.cells[i];
		// int level = cell.level;
		
		// vector<int> tmpIntEdges;
		// for(auto& j : cellEdges[i]){
			
			// if( edgeLevel[j] != level+1 ) continue;
			
			// int point0 = edgesPoints[j][0];
			// int point1 = edgesPoints[j][1];
			
			// bool boolContinue=true;
			// for(auto& k : cellVertexPoints[i]){
				// if(point0 == k || point1 == k){
					// boolContinue=false;
					// break;
				// }
			// }
			// if(boolContinue){
				// tmpIntEdges.push_back(j);
			// }
		// }
		
		// vector<int> tmpEdgeCenterPoint;
		// vector<int> tmpFaceCenterPoint;
		// for(auto& j : tmpIntEdges){
			// int point0 = edgesPoints[j][0];
			// int point1 = edgesPoints[j][1];
			// int matchNum0 = 0;
			// int matchNum1 = 0;
			// for(auto& k : tmpIntEdges){
				// int point2 = edgesPoints[k][0];
				// int point3 = edgesPoints[k][1];
				// if(point0==point2) ++matchNum0;
				// if(point0==point3) ++matchNum0;
				// if(point1==point2) ++matchNum1;
				// if(point1==point3) ++matchNum1;
			// }
			// if(matchNum0>=3){
				// if ( std::find( tmpFaceCenterPoint.begin(), tmpFaceCenterPoint.end(), point0 ) 
					// == tmpFaceCenterPoint.end() ) {
					// tmpFaceCenterPoint.push_back(point0);
				// }
			// }
			// else{
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), point0 ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(point0);
				// }
			// }
			// if(matchNum1>=3){
				// if ( std::find( tmpFaceCenterPoint.begin(), tmpFaceCenterPoint.end(), point1 ) 
					// == tmpFaceCenterPoint.end() ) {
					// tmpFaceCenterPoint.push_back(point1);
				// }
			// }
			// else{
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), point1 ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(point1);
				// }
			// }
		// }
		
		// for(auto& j : cell.faces){
			// auto& face = mesh.faces[j];
			
			// if( face.level != level ) continue;
			
			// tmpFaceCenterPoint.push_back(newFaceCenterPointNum[j]);
			
			// for(auto& k : faceEdgeCenterPoints[j]){
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), k ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(k);
				// }
			// }
		// }
		
		
		

		
		
		// vector<vector<int>> tmpIntFaceCenterPoints(tmpEdgeCenterPoint.size(),vector<int>(0,0));
		// vector<vector<int>> tmpIntFaceOwnNgb(tmpEdgeCenterPoint.size(),vector<int>(0,0));
		// int tmpNum = 0;
		// for(auto& j : cell.faces){
			// auto& face = mesh.faces[j];
			
			// if( face.level == level ) {
				
				// vector<int> tmpTmpPoints;
				// int tmpNumber=0;
				// for(auto& k : faceEdgeCenterPoints[j]){
					// int saveNum=0;
					// for(auto& l : tmpEdgeCenterPoint){
						// // cout << k << " " << l << endl;
						// if(k==l) break;
						// ++saveNum;
					// }
					
					
					// int firNum = tmpNumber;
					// int lasNum = tmpNumber+1;
					// if( tmpNumber == faceEdgeCenterPoints[j].size()-1 ){
						// lasNum = 0;
					// }
					// int tmpVerP0 = faceVertexPoints[j][firNum];
					// int tmpVerP1 = faceVertexPoints[j][lasNum];
					
					
					// // cout << tmpVerP0 << " " << tmpVerP1 << endl;
					
					
					
	// // cout << tmpIntFaceOwnNgb.size() << " " << saveNum << endl;
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpVerP0 ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						
						// tmpIntFaceOwnNgb[saveNum].push_back( tmpVerP0 );
						
					// }
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpVerP1 ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						
						// tmpIntFaceOwnNgb[saveNum].push_back( tmpVerP1 );
						
					// }
	// // cout << "BBBBBBBB" << endl;
					
					// int tmpCellFaceCenterPoint = newFaceCenterPointNum[j];
					// if ( std::find( tmpIntFaceCenterPoints[saveNum].begin(), 
									// tmpIntFaceCenterPoints[saveNum].end(), 
									// tmpCellFaceCenterPoint ) 
						// == tmpIntFaceCenterPoints[saveNum].end() ) {
						// tmpIntFaceCenterPoints[saveNum].push_back(tmpCellFaceCenterPoint);
					// }
					
					
					// ++tmpNumber;
					
				// }
				
			
			// }
			// else{
				
				// int tmpCellFaceVertexPoint=-1;
				// int tmpCellFaceCenterPoint=-1;
				// vector<int> tmpFaceEdgeCenterPoint;
				
				// bool boolBreak=false;
				// for(auto& k : face.points){
					// for(auto& l : cellVertexPoints[i]){
						// if(k==l){
							// tmpCellFaceVertexPoint = k;
							// boolBreak=true;
							// break;
						// }
					// }
					// if(boolBreak) break;
				// }
				// boolBreak=false;
				// for(auto& k : face.points){
					// for(auto& l : tmpFaceCenterPoint){
						// if(k==l){
							// tmpCellFaceCenterPoint = k;
							// boolBreak=true;
							// break;
						// }
					// }
					// if(boolBreak) break;
				// }
				// for(auto& k : face.points){
					// for(auto& l : tmpEdgeCenterPoint){
						// if(k==l){
							// tmpFaceEdgeCenterPoint.push_back(k);
						// }
						// if(tmpFaceEdgeCenterPoint.size()==2) break;
					// }
					// if(tmpFaceEdgeCenterPoint.size()==2) break;
				// }
				
				
				// for(auto& k : tmpFaceEdgeCenterPoint){
					// int saveNum=0;
					// for(auto& l : tmpEdgeCenterPoint){
						// if(k==l) break;
						// ++saveNum;
					// }
					
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpCellFaceVertexPoint ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						// tmpIntFaceOwnNgb[saveNum].push_back(tmpCellFaceVertexPoint);
					// }
					
					// if ( std::find( tmpIntFaceCenterPoints[saveNum].begin(), 
									// tmpIntFaceCenterPoints[saveNum].end(), 
									// tmpCellFaceCenterPoint ) == tmpIntFaceCenterPoints[saveNum].end() ) {
						// tmpIntFaceCenterPoints[saveNum].push_back(tmpCellFaceCenterPoint);
					// }
					
				// }
				
			// }
			
		// }
		
		
		// newCellIntEdgeCenterPoints[i] = tmpEdgeCenterPoint;
		
		
		
		// // input
		// int cellNumber = newCellNumbering[i];
		// for(int j=0; j<tmpEdgeCenterPoint.size(); ++j){
			
			
			// int ownIPoint = tmpIntFaceOwnNgb[j][0];
			// int ngbIPoint = tmpIntFaceOwnNgb[j][1];
			
			// // cout << tmpEdgeCenterPoint[j] << " " << ownIPoint << " " << ngbIPoint << endl;
			
			// int ownGIPoint = 0;
			// for(auto& k : cellVertexPoints[i]){
				// if(k==ownIPoint) break;
				// ++ownGIPoint;
			// }
			// int ngbGIPoint = 0;
			// for(auto& k : cellVertexPoints[i]){
				// if(k==ngbIPoint) break;
				// ++ngbGIPoint;
			// }
			
			
			// newIntFaceOwner.push_back(cellNumber+ownGIPoint);
			// newIntFaceNeighbour.push_back(cellNumber+ngbGIPoint);
			
			// vector<int> tmpPoints;
			// int edgeCentP = tmpEdgeCenterPoint[j];
			// int faceCentP0 = tmpIntFaceCenterPoints[j][0];
			// int faceCentP1 = tmpIntFaceCenterPoints[j][1];
			// int cellCentP = newCellCenterPointNum[i];
			// tmpPoints.push_back(edgeCentP);
			// tmpPoints.push_back(faceCentP0);
			// tmpPoints.push_back(cellCentP);
			// tmpPoints.push_back(faceCentP1);
			
			
			
			
			// // 포인트 순서 판단
			// vector<double> vectorEdge;
			// vectorEdge.push_back( mesh.points[ngbIPoint].x - mesh.points[edgeCentP].x );
			// vectorEdge.push_back( mesh.points[ngbIPoint].y - mesh.points[edgeCentP].y );
			// vectorEdge.push_back( mesh.points[ngbIPoint].z - mesh.points[edgeCentP].z );
			
			// vector<double> Vx, Vy, Vz;
			// for(auto& k : tmpPoints){
				// Vx.push_back(mesh.points[k].x);
				// Vy.push_back(mesh.points[k].y);
				// Vz.push_back(mesh.points[k].z);
			// }
			
			// vector<double> vectorNormFace;
			// double dummy1;
			// geometric.calcUnitNormals_Area3dPolygon(4, Vx,Vy,Vz, vectorNormFace, dummy1);
			
			// double cosTheta = 
				// vectorNormFace[0]*vectorEdge[0]+
				// vectorNormFace[1]*vectorEdge[1]+
				// vectorNormFace[2]*vectorEdge[2];
				
			// if(cosTheta < 0.0){
				// std::reverse(tmpPoints.begin(), tmpPoints.end());
			// }
			
			
			// newIntFacePoints.push_back(tmpPoints);
			
		// }
		
		
	// }
	
	// // cout << "BBBBBBBB" << endl;
	
	
	

	
	// if(rank==0) cout << "face-point start" << endl;
	
	// //====================================================
	// // AMR : 면 - 포인트 연결 
	// vector<vector<int>> newFacePoints;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// if(boolFaceRefine[i] == false) {
			// vector<int> tmpNewFacePoints(0,0);
			// newFacePoints.push_back(tmpNewFacePoints);
			// continue;
		// }
		
		// auto& face = mesh.faces[i];
		// int level = face.level;
		// // vector<int> tmpPointsRefine;
		// // vector<int> tmpPointsNumOrigin;
		// // for(int j=0; j<face.points.size(); ++j){
			// // int iPoint = face.points[j];
			// // if( mesh.points[iPoint].level == level+1 ){
				// // tmpPointsRefine.push_back(iPoint);
			// // }
			// // if( mesh.points[iPoint].level <= level ){
				// // tmpPointsNumOrigin.push_back(j);
			// // }
		// // }
		// // for(auto& j : tmpPointsRefine){
			// // cout << " pr :" << j << endl;
		// // }
		
		// vector<int> tmpTotalPoints;
		// for(int j=0; j<face.points.size(); ++j){
			// int iPoint = face.points[j];
			// tmpTotalPoints.push_back(iPoint);
			
			// int tmpIEdge = facesEdges[i][j];
			// if(
			// boolEdgeRefine[tmpIEdge] == true
			// ){
				// int tmpEdgeCentP = newEdgeCenterPointNum[tmpIEdge];
				// tmpTotalPoints.push_back(tmpEdgeCentP);
			// }
		// }
		
		// vector<int> tmpRefineVertexPoints;
		// for(int j=0; j<face.points.size(); ++j){
			// int iPoint = face.points[j];
			// if(mesh.points[iPoint].level <= level+1){
				// tmpRefineVertexPoints.push_back(iPoint);
			// }
			
			// int tmpIEdge = facesEdges[i][j];
			// if(
			// boolEdgeRefine[tmpIEdge] == true &&
			// edgeLevel[tmpIEdge] == level
			// ){
				// int tmpEdgeCentP = newEdgeCenterPointNum[tmpIEdge];
				// tmpRefineVertexPoints.push_back(tmpEdgeCentP);
			// }
			// // cout << iPoint << " " << edgeLevel[tmpIEdge] << " " << level << " " << boolEdgeRefine[tmpIEdge] << " " << 
			// // mesh.points[edgesPoints[tmpIEdge][0]].level << " " << mesh.points[edgesPoints[tmpIEdge][1]].level << " " << endl;
		// }
			
		// // for(int j=0; j<tmpRefineVertexPoints.size(); ++j){
			// // cout << j << " " << tmpRefineVertexPoints[j] << endl;
		// // }
		// // for(int j=0; j<tmpTotalPoints.size(); ++j){
			// // cout << j << " " << tmpTotalPoints[j] << endl;
		// // }
		
		// vector<vector<int>> tmpSubEdgePoints;
		// vector<int> tmpSubPoints;
		// int iSubEdges = 1;
		// for(int j=0; j<tmpTotalPoints.size(); ++j){
			// int iPoint = tmpTotalPoints[j];
			// tmpSubPoints.push_back(iPoint);
			
			// if(
			// tmpRefineVertexPoints[iSubEdges] == iPoint
			// ){
				// tmpSubEdgePoints.push_back(tmpSubPoints);
				// tmpSubPoints.clear();
				// tmpSubPoints.push_back(iPoint);
				// ++iSubEdges;
				// if(tmpRefineVertexPoints.size()-1 < iSubEdges){
					// iSubEdges = 0;
				// }
			// }
		// }
		// tmpSubPoints.push_back(tmpTotalPoints[0]);
		// tmpSubEdgePoints.push_back(tmpSubPoints);
		
		
		
		// // for(int j=0; j<tmpSubEdgePoints.size(); ++j){
			// // for(auto& k : tmpSubEdgePoints[j]){
				// // cout << j << " " << k << endl;
			// // }
		// // }
		
		// int faceRefineSize = tmpSubEdgePoints.size()/2;
		
		// // cout << faceRefineSize << " " << tmpRefineVertexPoints.size() << " " << faceVertexPoints[i].size() << endl;
		
		// for(int j=0; j<faceRefineSize; ++j){
			// int fEdgeSet = j*2;
			// int lEdgeSet = j*2-1;
			// if(j==0) lEdgeSet = tmpSubEdgePoints.size()-1;
			
			// vector<int> tmpNewFacePoints;
			
			// for(int k=0; k<tmpSubEdgePoints[fEdgeSet].size(); ++k){
				// int tmpPoi = tmpSubEdgePoints[fEdgeSet][k];
				// tmpNewFacePoints.push_back(tmpPoi);
			// }
			
			// tmpNewFacePoints.push_back( newFaceCenterPointNum[i] );
			
			// for(int k=0; k<tmpSubEdgePoints[lEdgeSet].size()-1; ++k){
				// int tmpPoi = tmpSubEdgePoints[lEdgeSet][k];
				// tmpNewFacePoints.push_back(tmpPoi);
			// }
			
			// newFacePoints.push_back(tmpNewFacePoints);
			
			// // for(auto& k : tmpNewFacePoints){
				// // if(k==2) cout << j << " " << k << endl;
			// // }
			
		// }
		
	// }
	
	
	// int test_newFacesNum = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(boolFaceRefine[i] == true){
			// for(int j=0; j<faceVertexPoints[i].size(); ++j){
				// ++test_newFacesNum;
			// }
		// }
		// else{
			// ++test_newFacesNum;
		// }
		
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 포인트로 넘버링된 "겉" 면 - 셀 연결
	// vector<int> newFaceOwner;
	// vector<int> newFaceNeighbour;
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// int tmpNewCellNumL = newCellNumbering[face.owner];
		
		// vector<int> tmpPointsSave;
		


		// if( 
		// boolFaceRefine[i]==true &&
		// boolCellRefine[face.owner] == true
		// ){
			// for(auto& k : face.points){
				// int tmpNum = 0;
				// int level = mesh.cells[face.owner].level;
				// int pLevel = mesh.points[k].level;
				// if(pLevel > level) continue;
				// // for(auto& l : mesh.cells[face.owner].points){
				// for(auto& l : cellVertexPoints[face.owner]){
					// if(k==l){
						// tmpPointsSave.push_back(tmpNum);
						// break;
					// }
					// ++tmpNum;
				// }
			// }
		// }
		// else if( 
		// boolFaceRefine[i]==false &&
		// boolCellRefine[face.owner] == true
		// ){
			// vector<int> tmpPoint; 
			// for(auto& k : faceVertexPoints[i]){
				// int tmpPNum = 0;
				// for(auto& l : cellVertexPoints[face.owner]){
					// if(k==l){
						// tmpPoint.push_back(tmpPNum);
						// break;
					// }
					// ++tmpPNum;
				// }
				// if(tmpPoint.size() == 1) break;
			// }
			// for(auto& k : tmpPoint){
				// tmpPointsSave.push_back(k);
			// }
		// }
		// else if( 
		// boolFaceRefine[i]==true &&
		// boolCellRefine[face.owner] == false
		// ){
			// for(auto& k : face.points){
				// int level = mesh.cells[face.owner].level;
				// int pLevel = mesh.points[k].level;
				// if(pLevel > level) continue;
				// tmpPointsSave.push_back(0);
			// }
		// }
		// else{
			// tmpPointsSave.push_back(0);
		// }
		
		// for(auto& k : tmpPointsSave){
				// // cout << " L " << tmpNewCellNumL << " " << k << endl;
			// newFaceOwner.push_back( tmpNewCellNumL + k );
		// }
		
		
		
		
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // cout << "INTERNAL_FACE" << endl;
			
			// int tmpNewCellNumR = newCellNumbering[face.neighbour];
			
			// tmpPointsSave.clear();
			
			// if( 
			// boolFaceRefine[i]==true &&
			// boolCellRefine[face.neighbour] == true
			// ){
				// for(auto& k : face.points){
					// int tmpNum = 0;
					// int level = mesh.cells[face.neighbour].level;
					// int pLevel = mesh.points[k].level;
					// if(pLevel > level) continue;
					// // for(auto& l : mesh.cells[face.neighbour].points){
					// for(auto& l : cellVertexPoints[face.neighbour]){
						// if(k==l){
							// tmpPointsSave.push_back(tmpNum);
							// break;
						// }
						// ++tmpNum;
					// }
				// }
			// }
			// else if( 
			// boolFaceRefine[i]==false &&
			// boolCellRefine[face.neighbour] == true
			// ){
				// vector<int> tmpPoint; 
				// for(auto& k : faceVertexPoints[i]){
					// int tmpPNum = 0;
					// for(auto& l : cellVertexPoints[face.neighbour]){
						// if(k==l){
							// tmpPoint.push_back(tmpPNum);
							// break;
						// }
						// ++tmpPNum;
					// }
					// if(tmpPoint.size() == 1) break;
				// }
				// for(auto& k : tmpPoint){
					// tmpPointsSave.push_back(k);
				// }
			// }
			// else if( 
			// boolFaceRefine[i]==true &&
			// boolCellRefine[face.neighbour] == false
			// ){
				// for(auto& k : face.points){
					// int level = mesh.cells[face.neighbour].level;
					// int pLevel = mesh.points[k].level;
					// if(pLevel > level) continue;
					// tmpPointsSave.push_back(0);
				// }
			// }
			// else{
				// tmpPointsSave.push_back(0);
			// }
			
			
			
			// for(auto& k : tmpPointsSave){
				// // cout << " R " << tmpNewCellNumR + k << endl;
				// newFaceNeighbour.push_back( tmpNewCellNumR + k );
			// }
		
			
		// }
		// else{
			// // cout << "BC_FACE" << endl;
			// if( boolFaceRefine[i]==true ){
				// for(auto& k : faceVertexPoints[i]){
					// newFaceNeighbour.push_back( -1 );
				// }
			// }
			// else{
				// newFaceNeighbour.push_back( -1 );
			// }
		// }
		
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 그룹핑
	// vector<int> cellGroupping;
	// vector<int> cellLeveling;
	// vector<vector<double>> cellVariables;
	// newPointNum = 0;
	// // for(auto& i : newCellCenterPointNum){
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// int level = mesh.cells[i].level;
		// int group = mesh.cells[i].group;
		// if(boolCellRefine[i] == true){
			// // cout << cellVertexPoints[i].size() << " AAAAAAAA " << endl;
			// for(auto& j : cellVertexPoints[i]){
				// cellLeveling.push_back(level+1);
				// cellGroupping.push_back(group);
				
				
				// vector<double> tmpVariables;
				// for(auto& k : cell.var){
					// tmpVariables.push_back(k);
				// }
				// cellVariables.push_back(tmpVariables);
			// }
		// }
		// else{
			// vector<double> tmpVariables;
			// cellLeveling.push_back(level);
			// cellGroupping.push_back(group);
			// for(auto& k : cell.var){
				// tmpVariables.push_back(k);
			// }
			// cellVariables.push_back(tmpVariables);
		// }
	// }
	
	
	
	// if(rank==0) cout << "face-add start" << endl;
	
	// //====================================================
	// // AMR : 면 추가 (면 -> new 셀 내부 면 -> proc면 -> B.C.면)
	// // AMR : 면 - 셀 리넘버링
	// // AMR : 면 - 포인트 리넘버링
	// // AMR : 면 타입 셋팅
	// // AMR : 면 바운더리 셋팅
	
	// SEMO_Mesh_Builder newMesh;
	
	// // 원래 내부 면
	// newPointNum = 0;
	// int newFacesNum = 0;
	// // cout << "AAAAAAAAAAAAAAA" << endl;
	// // for(auto& i : newFaceCenterPointNum){
	// for(int i=0; i<mesh.faces.size(); ++i){
		// int level = mesh.faces[i].level;
		// int group = mesh.faces[i].group;
		
		// // cout << level << " " << group << endl;
		
		// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			// if(boolFaceRefine[i] == true){
				// // if(faceVertexPoints.size()-1 < i) cout << "ERROR!!!!!!!!!!" << endl;
				// // if(faceVertexPoints[i].size()>8) cout << faceVertexPoints[i].size() << endl;
				// for(int j=0; j<faceVertexPoints[i].size(); ++j){
					// newMesh.addFace();
					// newMesh.faces.back().group = group;
					// newMesh.faces.back().level = level + 1;
					
					// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
					// newMesh.faces.back().neighbour = newFaceNeighbour[newFacesNum];
					
					// for(auto& k : newFacePoints[newFacesNum]){
						// // if(k==2) cout << k << endl;
						// newMesh.faces.back().points.push_back(k);
					// }
					
					// newMesh.faces.back().setType(mesh.faces[i].getType());
					
					// ++newFacesNum;
				// }
				
			// }
			// else{
				// newMesh.addFace();
				// newMesh.faces.back().group = group;
				// newMesh.faces.back().level = level;
				
				// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
				// newMesh.faces.back().neighbour = newFaceNeighbour[newFacesNum];
				
				// newMesh.faces.back().setType(mesh.faces[i].getType());
				
				// for(int j=0; j<facesEdges[i].size(); ++j){
					// int tmpNum = facesEdges[i][j];
					// if(boolEdgeRefine[tmpNum] == true){
						// int point0 = mesh.faces[i].points[j];
						// newMesh.faces.back().points.push_back(point0);
						// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
					// }
					// else{
						// int point0 = mesh.faces[i].points[j]; 
						// newMesh.faces.back().points.push_back(point0);
					// }
				// }
				
				// ++newFacesNum;
			// }
		// }
	// }
	
	
	// // cout << "BBBBBBBBBBB" << endl;
	// // MPI_Barrier(MPI_COMM_WORLD);
	
	
	// // 셀 내부 면
	// int newIntFaceNum = 0;
	// // for(auto& i : newCellCenterPointNum){
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == true){
			
			// int level = mesh.cells[i].level;
			// int group = mesh.faces[i].group;
			
			// for(int j=0; j<newCellIntEdgeCenterPoints[i].size(); ++j){
				// newMesh.addFace();
				// newMesh.faces.back().group = mesh.faces.size() + newIntFaceNum;
				// newMesh.faces.back().level = level + 1;
				// newMesh.faces.back().owner = newIntFaceOwner[newIntFaceNum];
				// newMesh.faces.back().neighbour = newIntFaceNeighbour[newIntFaceNum];
				// newMesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
				
				// for(auto& k : newIntFacePoints[newIntFaceNum]){
					// newMesh.faces.back().points.push_back(k);
				// }
					
				// ++newIntFaceNum;
			// }
		// }
	// }
	
	
	
	
	// if(rank==0) cout << "proc & bc face start" << endl;
	
	
	// // 프록, 경계면
	// vector<int> newBoundaryStr;
	// int BCnum = 0;
	// int test_num2 = 0;
	
	// // for(auto& i : newFaceCenterPointNum){
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];
		
		
		// int procNeighbNo = mesh.boundary[BCnum].neighbProcNo;
		
		// if(BCnum < mesh.boundary.size()){
			// if( mesh.boundary[BCnum].startFace == i ){
				
				// // if(mesh.boundary.size() <= BCnum) {
					// // cout << "| #Error : mesh.boundary.size() <= BCnum" << endl;
				// // }
				
				// newMesh.addBoundary();
				// newMesh.boundary[BCnum].name = mesh.boundary[BCnum].name;
				// newMesh.boundary[BCnum].type = mesh.boundary[BCnum].type;
				// for(auto& j : mesh.boundary[BCnum].var){
					// newMesh.boundary[BCnum].var.push_back(j);
				// }
				// newMesh.boundary[BCnum].startFace = newMesh.faces.size();
				// if(BCnum != 0) {
					// newMesh.boundary[BCnum-1].nFaces = 
						// newMesh.faces.size()-newMesh.boundary[BCnum-1].startFace;
				// }
				// newMesh.boundary[BCnum].myProcNo = rank;
				// newMesh.boundary[BCnum].neighbProcNo = mesh.boundary[BCnum].neighbProcNo;
				// ++BCnum;
			// }
		// }
		
		// if(face.getType() != SEMO_Types::INTERNAL_FACE){
			// if(boolFaceRefine[i] == true){
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // cout << rank << " yes refine " << test_num2 << endl;
					// // ++test_num2;
				// // }
				
				// int level = face.level;
				// int group = face.group;
				
				// // int nChildCells = mesh.faces[i].points.size();
				// int nChildCells = faceVertexPoints[i].size();
				
				
				// vector<int> vecNewPointNumNum;
				// for(int j=0; j<nChildCells; ++j){
					// vecNewPointNumNum.push_back(newFacesNum);
					// ++newFacesNum;
				// }
				
				// // if(
				// // mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
				// // cout << nChildCells << endl;
				// // }
				// if(
				// face.getType() == SEMO_Types::PROCESSOR_FACE &&
				// rank > procNeighbNo
				// ){
					// std::reverse(vecNewPointNumNum.begin(), vecNewPointNumNum.end());
				// }
				
				
				// for(auto& j : vecNewPointNumNum){
					// newMesh.addFace();
					// newMesh.faces.back().group = group;
					// newMesh.faces.back().level = level + 1;
					
					// newMesh.faces.back().owner = newFaceOwner[j];
					// newMesh.faces.back().neighbour = -1;
					// newMesh.faces.back().setType(face.getType());
					
					// for(auto& k : newFacePoints[j]){
						// newMesh.faces.back().points.push_back(k);
					// }
					
					// if(
					// face.getType() == SEMO_Types::PROCESSOR_FACE &&
					// rank > procNeighbNo
					// ){
						// std::reverse(newMesh.faces.back().points.begin()+1, 
									 // newMesh.faces.back().points.end());
						// std::reverse(newMesh.faces.back().points.begin(), 
									 // newMesh.faces.back().points.end());
					// }
					
				// }
				
				
			// }
			// else{
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // cout << rank << " no refine " << test_num2 << endl;
					// // ++test_num2;
				// // }
				// newMesh.addFace();
				// newMesh.faces.back().group = face.group;
				// newMesh.faces.back().level = face.level;
				// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
				// newMesh.faces.back().neighbour = -1;
				// newMesh.faces.back().setType(face.getType());
				
				// // cout << face.points.size() << " " << facesEdges[i].size() << endl;
				
				// if(
				// face.getType() == SEMO_Types::PROCESSOR_FACE &&
				// rank > procNeighbNo
				// ){
					
					// int lastEdNum = facesEdges[i][face.points.size()-1];
					// if(boolEdgeRefine[lastEdNum] == true){
						// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[lastEdNum]);
						// // ++test_num2;
					// }
						
					// for(int j=0; j<face.points.size(); ++j){
						// int tmpNum = facesEdges[i][j];
						// if(boolEdgeRefine[tmpNum] == true){
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
							// if(j!=face.points.size()-1) 
								// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
							// // ++test_num2;
						// }
						// else{
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
						// }
					// }
				// }
				// else{
					// for(int j=0; j<face.points.size(); ++j){
						// int tmpNum = facesEdges[i][j];
						// if(boolEdgeRefine[tmpNum] == true){
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
							// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
							// // ++test_num2;
						// }
						// else{
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
						// }
					// }
				// }
				
				
				
				// ++newFacesNum;
				
				
			// }
		// }
	// }
	// if(mesh.boundary.size() < BCnum) {
		// cout << "| #Error 2 : mesh.boundary.size() <= BCnum" << endl;
	// }
	// newMesh.boundary[BCnum-1].nFaces = 
		// newMesh.faces.size()-newMesh.boundary[BCnum-1].startFace;
		
		
		
	// // for(int i=0; i<mesh.faces.size(); ++i){
		
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // if(boolFaceRefine[i] == true){
			// // }
			// // else{
				// // for(int j=0; j<face.points.size(); ++j){
					// // int tmpNum = facesEdges[i][j];
					// // int iNewPoint0 = face.points[j];
					// // // if(rank==1) cout << j << " P : " << test_num2 << " " << mesh.points[iNewPoint0].x << " " << mesh.points[iNewPoint0].y << " " << mesh.points[iNewPoint0].z << " " << endl;
					// // if(boolEdgeRefine[tmpNum] == true){
						// // int iNewPoint = newEdgeCenterPointNum[tmpNum];
						// // // if(rank==1) cout << j << " : " << test_num2 << " " << mesh.points[iNewPoint].x << " " << mesh.points[iNewPoint].y << " " << mesh.points[iNewPoint].z << " " << endl;
						// // ++test_num2;
					// // }
				// // }
			// // }
			
		// // }
	// // }
		
		
	// // cout << "AGAGAG " << test_num2 << endl;
	
	
	// if(rank==0) cout << "mesh copy start" << endl;
		
	// //====================================================
	// // 원래 메쉬 클리어 및 복사
	// mesh.cells.clear();
	// mesh.faces.clear();
	// // mesh.points.clear();
	// mesh.boundary.clear();
	
	// for(int i=0; i<newMesh.faces.size(); ++i){
		// auto& newFace = newMesh.faces[i];
		
		// mesh.addFace();
		// mesh.faces.back().group = newFace.group;
		// mesh.faces.back().level = newFace.level;
		// mesh.faces.back().owner = newFace.owner;
		// mesh.faces.back().neighbour = newFace.neighbour;
		// mesh.faces.back().setType(newFace.getType());
		// for(auto& j : newFace.points){
			// mesh.faces.back().points.push_back(j);
		// }
		
		// mesh.faces.back().varL.resize(nTotalFaceLRVar,0.0);
		// mesh.faces.back().varR.resize(nTotalFaceLRVar,0.0);
		// mesh.faces.back().var.resize(nTotalFaceVar,0.0);
		
	// }
	
	// for(int i=0; i<newMesh.boundary.size(); ++i){
		// auto& newBoundary = newMesh.boundary[i];
		// mesh.addBoundary();
		// mesh.boundary.back().name = newBoundary.name;
		// mesh.boundary.back().type = newBoundary.type;
		// for(auto& j : newBoundary.var){
			// mesh.boundary.back().var.push_back(j);
		// }
		// mesh.boundary.back().startFace = newBoundary.startFace;
		// mesh.boundary.back().nFaces = newBoundary.nFaces;
		// mesh.boundary.back().myProcNo = newBoundary.myProcNo;
		// mesh.boundary.back().neighbProcNo = newBoundary.neighbProcNo;
		
	// }
	
	
	
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // cout << face.owner << " " << face.neighbour << endl;
		
		// // // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // // cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
		// // // }
	// // }
	// // cout << "cccccc" << endl;

	// mesh.check();
	
	
	
	
	// // cout << "cccccc" << endl;
	
	// // mesh.checkMatchingProcessorFace();
	
	// // cout << "ddddd" << endl;
	
	
	
	
	// mesh.buildCells();
	
	// // mesh.setFaceTypes();

	// mesh.buildLists();
	
	// mesh.connectCelltoFaces();
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// // int minLevel = 999999999;
		// // for(int j=0; j<cell.faces.size(); ++j){
			// // auto& face = mesh.faces[cell.faces[j]];
			// // minLevel = min(minLevel,face.level);
		// // }
		// // cell.level = minLevel;
		
		// cell.level = cellLeveling[i];
		// cell.group = cellGroupping[i];
		
		
		// // cout << "CELL VAR " << cellVariables[i].size() << endl;
		// for(auto& k : cellVariables[i]){
			// cell.var.push_back(k);
		// }
		
		
		// // cout << " G : " << cellGroupping[i] << " " << cellGroupping.size() << " " << mesh.cells.size() << endl;
	// }
	
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// mesh.connectCelltoPoints();
	
	// // set processor face counts
	// mesh.setCountsProcFaces();
	
	// // set processor face displacements
	// mesh.setDisplsProcFaces(); 
		
	// // mesh.informations();
	
	
	
	// // SEMO_Mesh_Save save;
	// // string tmpFile = "./Rf" + to_string(iter);
	// // // string tmpFile = "./";
	// // save.vtu(tmpFile, rank, mesh);
	
	
	// // cout << "SIZE : " << cellGroupping.size() << " " << mesh.cells.size() << endl;
		
	// // cout << "CELL" << endl;
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // auto& cell = mesh.cells[i];
		// // cout << i << " " << cell.level << " " << cell.group << endl;
	// // }
	
	// // cout << "FACE" << endl;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // cout << i << " " << face.level << " " << face.group << endl;
	// // }
		

	
	// if(rank==0) cout << "----Refine finished-----" << endl;
}









void SEMO_Poly_Mesh_Refine::separateEdgesPoints(
	SEMO_Mesh_Builder& mesh, 
	int faceLevel,
	vector<int>& points, 
	vector<int>& edges,
	vector<bool>& boolEdgeRefine,
	vector<int>& edgesCenterPointNumber,
	vector<int>& edgesLevel,
	vector<int>& faceVertex,
	vector<int>& edgeCenterPoints,
	vector<vector<int>>& subFaceEdgePoints,
	int& newPointNumber){


	faceVertex.clear();
	edgeCenterPoints.clear();
	subFaceEdgePoints.clear();

	vector<int> tmpVector;
	int iEdge;
	for(int k=0; k<points.size(); ++k){
		
		
		
		int iPoint = points[k];
		int pointLevel = mesh.points[iPoint].level;
		iEdge = edges[k];
		bool saveEdgeRefine = boolEdgeRefine[iEdge];
		int edgeLevel = edgesLevel[iEdge];
		
		if( pointLevel <= faceLevel ){
			faceVertex.push_back(iPoint);
			
			if( tmpVector.size() > 0 ){
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
			
	// cout << k << " " << points.size() << endl;
			if(
			faceVertex.size() != 1 &&
			edgeCenterPoints.size() != faceVertex.size()-1 &&
			saveEdgeRefine == false
			) {
				
				int backVertexNum = faceVertex.back();
				
				boolEdgeRefine[iEdge] = true;
				
				// edgesCenterPointNumber[iEdge] = newPointNumber++;
				edgesCenterPointNumber[iEdge] = mesh.points.size();
				// this->addCenterPoint(mesh, cellVertex, cellLevel);
				
				vector<int> vertex;
				vertex.push_back(backVertexNum);
				vertex.push_back(iPoint);
				this->addCenterPoint(mesh, vertex, faceLevel);
				
				edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
				
				tmpVector.push_back(edgesCenterPointNumber[iEdge]);
			}
			
			
		}
		else if( pointLevel == faceLevel+1 ){
			
			edgeCenterPoints.push_back(iPoint);
			
			if( tmpVector.size() > 0 ){
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
		}
		
		tmpVector.push_back(iPoint);
		
		if( saveEdgeRefine == true ){
			if( faceLevel == edgeLevel ){
				edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
				subFaceEdgePoints.push_back(tmpVector);
				tmpVector.clear();
			}
			
			tmpVector.push_back(edgesCenterPointNumber[iEdge]);
			
		}
	}
	
	subFaceEdgePoints.push_back(tmpVector);
	
	if( edgeCenterPoints.size() != faceVertex.size() ){
		boolEdgeRefine[iEdge] = true;
		// edgesCenterPointNumber[iEdge] = newPointNumber++;
		edgesCenterPointNumber[iEdge] = mesh.points.size();
		
		vector<int> vertex;
		vertex.push_back(faceVertex[0]);
		vertex.push_back(faceVertex.back());
		this->addCenterPoint(mesh, vertex, faceLevel);
		
		edgeCenterPoints.push_back(edgesCenterPointNumber[iEdge]);
	}
	
	if( edgeCenterPoints.size() != faceVertex.size() ){
		cout << "| #Error 1" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
}






































// void SEMO_Poly_AMR_Builder::polyRefine(
	// SEMO_Mesh_Builder& mesh, 
	// SEMO_Controls_Builder& controls,
	// int iter){

	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	
	// if(rank==0) cout << "----Refine start-----" << endl;
	
	
	// SEMO_Mesh_Geometric geometric;
	
	// int nTotalFaceLRVar = mesh.faces[0].varL.size();
	// int nTotalFaceVar = mesh.faces[0].var.size();
	
	
	// // mesh.checkMatchingProcessorFace();

	// //====================================================
	// // Refine 되는 셀 & 면 조사
	
	// // random
    // std::random_device rd;
    // std::default_random_engine eng(rd());
    // std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	
	
	// // SEMO_Utility_Math math;
	// // vector<vector<double>> gradVF;
	// // // math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// // math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	
	// vector<bool> boolCellRefine(mesh.cells.size(),false);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// // if( distr(eng) > 0.9 ){
			// // boolCellRefine[i] = true;
		// // }
		
		// if( mesh.cells[i].var[controls.indicatorAMR] > 1.0 ){
			// boolCellRefine[i] = true;
		// }
		
		// if(mesh.cells[i].level >= 2) boolCellRefine[i] = false;
		
		// // if(mesh.cells[i].volume < 1.e-8) boolCellRefine[i] = false;
		
	// } 
	
	
	// //#################################################
	// // prcoface 앞 셀은 Refine 안되게
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// boolCellRefine[face.owner] = false;
		// }
	// }
	
	
// // cout << "AAAAAAAAAA" << endl;
	
	// //====================================================
	// // MPI
	// vector<int> cLevel_send, cLevel_recv;
	// vector<int> cRefine_send, cRefine_recv;
	// if(size>1){
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// cLevel_send.push_back(mesh.cells[face.owner].level);
				
				// if(boolCellRefine[face.owner]){
					// cRefine_send.push_back(1);
				// }
				// else{
					// cRefine_send.push_back(0);
				// }
			// }
		// }
		
		// // SEMO_MPI_Builder mpi;
		
		// cLevel_recv.resize(cLevel_send.size(),0);
		// cRefine_recv.resize(cRefine_send.size(),0);

		// MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // MPI_COMM_WORLD);
					   
		// MPI_Alltoallv( cRefine_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // cRefine_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   // MPI_COMM_WORLD);
					   
		
		// // mpi.setProcsFaceDatas(
					// // cLevel_send, cLevel_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
		
		// // mpi.setProcsFaceDatas(
					// // cRefine_send, cRefine_recv,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
					   
	// }
	
// // cout << "BBBBBBBB" << endl;
	
	// //====================================================
	// // Refine 되는 셀 옆에 레벨보다 크면 안됨
	// int proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// if(mesh.cells[face.owner].level > mesh.cells[face.neighbour].level){
				// boolCellRefine[face.owner] = false;
			// }
			// if(mesh.cells[face.owner].level < mesh.cells[face.neighbour].level){
				// boolCellRefine[face.neighbour] = false;
			// }
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// if(mesh.cells[face.owner].level > cLevel_recv[proc_num]){
				// boolCellRefine[face.owner] = false;
			// }
			// ++proc_num;
		// }
	// }
	
	
	
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // if(boolCellRefine[i]) cout << rank << " " << i << endl;
	// // }
	
	
	
	
	// //====================================================
	// // 면 Refine
	// vector<bool> boolFaceRefine(mesh.faces.size(),false);
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// if(mesh.cells[face.owner].level == mesh.cells[face.neighbour].level){
				// if(
				// boolCellRefine[face.owner] == true ||
				// boolCellRefine[face.neighbour] == true){
					// boolFaceRefine[i] = true;
				// }
			// }
			
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// if(boolCellRefine[face.owner] == true){
				// boolFaceRefine[i] = true;
			// }
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// if(mesh.cells[face.owner].level == cLevel_recv[proc_num]){
				// if(
				// boolCellRefine[face.owner] == true ||
				// cRefine_recv[proc_num] == 1){
					// boolFaceRefine[i] = true;
				// }
			// }
			// ++proc_num;
		// }
		
	// }
	
	
	// cLevel_send.clear();
	// cRefine_send.clear();
	// cLevel_recv.clear();
	// cRefine_recv.clear();
	

	// //====================================================
	// // 엣지 생성
	// vector<vector<int>> edgesPoints;
	// vector<vector<int>> facesEdges(mesh.faces.size(),vector<int>(0,0));
	// vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// int pointSize = face.points.size();
		// for(int j=0; j<pointSize; ++j){
			// int ipoint0 = face.points[j];
			// int ipoint1 = ( j+1 == pointSize ? face.points[0] : face.points[j+1] );
			// vector<int> matchFaces;
			// for(auto& k : pointsFaces[ipoint0]){
				// for(auto& l : pointsFaces[ipoint1]){
					// if( k == l ) {
						// matchFaces.push_back(l);
					// }
				// }
			// }
			
			// if(matchFaces.size()==0){
				// vector<int> tmpPoints;
				// tmpPoints.push_back(ipoint0);
				// tmpPoints.push_back(ipoint1);
				// edgesPoints.push_back(tmpPoints);
				
				// facesEdges[i].push_back(edgesPoints.size()-1);
				
			// }
			// else{
				// int iFace = matchFaces[0];
				// int iEdgeSave = -1;
				// for(auto& iEdge : facesEdges[iFace]){
					// if(
					// (edgesPoints[iEdge][0]==ipoint0 && edgesPoints[iEdge][1]==ipoint1) ||
					// (edgesPoints[iEdge][1]==ipoint0 && edgesPoints[iEdge][0]==ipoint1) 
					// ){
						// iEdgeSave = iEdge;
						// break;
					// }
				// }
				
				// facesEdges[i].push_back(iEdgeSave);
			// }
			// pointsFaces[ipoint0].push_back(i);
			
		// }
	// }
	// pointsFaces.clear();
	
	
	// vector<vector<int>> edgesFaces(edgesPoints.size(),vector<int>(0,0));
	// for(int i=0; i<mesh.faces.size(); ++i){
		// for(auto& j : facesEdges[i]){
			// edgesFaces[j].push_back(i);
		// }
	// }
	
	
	// //====================================================
	// // AMR : 엣지 레벨 (최대 포인트 레벨)
	// vector<int> edgeLevel(edgesPoints.size(),0);
	// for(int i=0; i<edgesPoints.size(); ++i){
		// int point0 = edgesPoints[i][0];
		// int point1 = edgesPoints[i][1];
		
		// edgeLevel[i] = max(
			// mesh.points[point0].level, 
			// mesh.points[point1].level);
	// }
	
	
	
	
	// //====================================================
	// // AMR : 엣지 refine
	// vector<bool> boolEdgeRefine(edgesPoints.size(),false);
	// for(int i=0; i<edgesPoints.size(); ++i){
		// int point0 = edgesPoints[i][0];
		// int point1 = edgesPoints[i][1];
		
		// bool tmpBool = false;
		// for(auto& j : edgesFaces[i]){
			// int level = mesh.faces[j].level;
			// if(
			// boolFaceRefine[j] == true &&
			// mesh.points[point0].level <= level &&
			// mesh.points[point1].level <= level 
			// ){
				// tmpBool = true;
			// }
		// }
		
		// if(tmpBool) boolEdgeRefine[i] = true;
	// }
	
	
	
	// //====================================================
	// // AMR : 엣지 refine MPI
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // cout<<"AAAAAAAAA"<<endl;
	
	
	// // if(size>1){
	 // // // for(int ip0=0; ip0<100; ++ip0){
		// // vector<int> eRefine_send, eRefine_recv;
		// // vector<int> countsEdge(size,0);
		// // vector<int> displsEdge(size,0);
		
		// // for(int ip=0; ip<size; ++ip){
			// // for(auto& boundary : mesh.boundary){
				// // if(ip == boundary.neighbProcNo){
					// // int str = boundary.startFace;
					// // int end = str + boundary.nFaces;
					// // for(int i=str; i<end; ++i){
						// // SEMO_Face& face = mesh.faces[i];
						// // for(auto& j : facesEdges[i]){
							// // if(boolEdgeRefine[j]){
								// // eRefine_send.push_back(1);
							// // }
							// // else{
								// // eRefine_send.push_back(0);
							// // }
						// // }
						// // countsEdge[ip] += facesEdges[i].size();
					// // }
					// // break;
				// // }
			// // }
		// // }
		
		// // displsEdge[0] = 0;
		// // for(int ip=1; ip<size; ++ip){
			// // displsEdge[ip] = displsEdge[ip-1] + countsEdge[ip-1];
			
		// // }
		
		// // SEMO_MPI_Builder mpi;
		
		// // mpi.setProcsFaceDatas(eRefine_send, eRefine_recv,
                              // // countsEdge, countsEdge, 
                              // // displsEdge, displsEdge);
					
		// // // int tmp_proc_num = 0;
		// // // for(int i=0; i<mesh.faces.size(); ++i){
			// // // auto& face = mesh.faces[i];
			// // // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// // // vector<int> tmpEdgeBool;
				// // // for(auto& j : facesEdges[i]){
					// // // tmpEdgeBool.push_back(eRefine_recv[tmp_proc_num]);
					// // // ++tmp_proc_num;
				// // // }
				
				// // // std::reverse(tmpEdgeBool.begin(), tmpEdgeBool.end());
				
				// // // for(int j=0; j<facesEdges[i].size(); ++j){
					// // // if(tmpEdgeBool[j]==1){
						// // // boolEdgeRefine[facesEdges[i][j]] = true;
					// // // }
				// // // }
			// // // }
		// // // }	
		// // int tmp_proc_num = 0;
		// // for(int ip=0; ip<size; ++ip){
			// // for(auto& boundary : mesh.boundary){
				// // if(ip == boundary.neighbProcNo){
					// // int str = boundary.startFace;
					// // int end = str + boundary.nFaces;
					// // for(int i=str; i<end; ++i){
						// // SEMO_Face& face = mesh.faces[i];
						// // vector<int> tmpEdgeBool;
						// // for(auto& j : facesEdges[i]){
							// // tmpEdgeBool.push_back(eRefine_recv[tmp_proc_num]);
							// // ++tmp_proc_num;
						// // }
						// // // std::reverse(tmpEdgeBool.begin(), tmpEdgeBool.end());
						// // for(int j=0; j<facesEdges[i].size(); ++j){
							// // if(tmpEdgeBool[j]==1){
								// // boolEdgeRefine[facesEdges[i][j]] = true;
							// // }
						// // }
					// // }
					// // break;
				// // }
			// // }
		// // }
	 // // // }		
	// // }
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // cout<<"BBBBBBBB"<<endl;
	
	
	
	// //====================================================
	// // AMR : 엣지 중심 점 생성
	// vector<vector<double>> newEdgeCenterPointsXYZ;
	// vector<int> newEdgeCenterPointsLevel;
	// vector<int> newEdgeCenterPointNum(edgesPoints.size(),-1);
	// int newPointNum = mesh.points.size();
	// int strEdgeCenterPointNum = newPointNum;
	// for(int i=0; i<edgesPoints.size(); ++i){
		
		// if(boolEdgeRefine[i] == false) continue;
		
		// int point0 = edgesPoints[i][0];
		// int point1 = edgesPoints[i][1];
		
		// double x0 = mesh.points[point0].x;
		// double y0 = mesh.points[point0].y;
		// double z0 = mesh.points[point0].z;
		
		// double x1 = mesh.points[point1].x;
		// double y1 = mesh.points[point1].y;
		// double z1 = mesh.points[point1].z;
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(0.5*(x0+x1));
		// tmpXYZ.push_back(0.5*(y0+y1));
		// tmpXYZ.push_back(0.5*(z0+z1));
		
		// newEdgeCenterPointsXYZ.push_back(tmpXYZ);
		// // newEdgeCenterPointsLevel.push_back(mesh.points[point0].level+1);
		// newEdgeCenterPointsLevel.push_back(edgeLevel[i]+1);
		
		// newEdgeCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	// //====================================================
	// // AMR : 면 중심 점 생성
	// vector<vector<double>> newFaceCenterPointsXYZ;
	// vector<int> newFaceCenterPointsLevel;
	// vector<int> newFaceCenterPointNum(mesh.faces.size(),-1);
	// int strFaceCenterPointNum = newPointNum;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// if(boolFaceRefine[i] == false) continue;
		
		// auto& face = mesh.faces[i];
		// double xCenter = 0.0;
		// double yCenter = 0.0;
		// double zCenter = 0.0;
		// double tmpNum = 0.0;
		// for(auto& j : face.points){
			// auto& point = mesh.points[j];
			// if(point.level>face.level) continue;
			// xCenter += point.x;
			// yCenter += point.y;
			// zCenter += point.z;
			// tmpNum = tmpNum+1.0;
		// }
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(xCenter/tmpNum);
		// tmpXYZ.push_back(yCenter/tmpNum);
		// tmpXYZ.push_back(zCenter/tmpNum);
		
		// newFaceCenterPointsXYZ.push_back(tmpXYZ);
		// newFaceCenterPointsLevel.push_back(face.level+1);
		
		// newFaceCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 중심 점 생성
	// vector<vector<double>> newCellCenterPointsXYZ;
	// vector<int> newCellCenterPointsLevel;
	// vector<int> newCellCenterPointNum(mesh.cells.size(),-1);
	// int strCellCenterPointNum = newPointNum;
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) continue;
		
		// auto& cell = mesh.cells[i];
		// double xCenter = 0.0;
		// double yCenter = 0.0;
		// double zCenter = 0.0;
		// double tmpNum = 0.0;
		// for(auto& j : cell.points){
			// auto& point = mesh.points[j];
			// if(point.level>cell.level) continue;
			// xCenter += point.x;
			// yCenter += point.y;
			// zCenter += point.z;
			// tmpNum = tmpNum+1.0;
		// }
		
		// vector<double> tmpXYZ;
		// tmpXYZ.push_back(xCenter/tmpNum);
		// tmpXYZ.push_back(yCenter/tmpNum);
		// tmpXYZ.push_back(zCenter/tmpNum);
		
		// newCellCenterPointsXYZ.push_back(tmpXYZ);
		// newCellCenterPointsLevel.push_back(cell.level+1);
		
		// newCellCenterPointNum[i] = newPointNum;
		// ++newPointNum;
	// }
	
	
	
	
	// if(rank==0) cout << "point add start" << endl;
	

	// //====================================================
	// // AMR : 포인트 추가 (엣지 중심 -> 면 중심 -> 셀 중심)
	// newPointNum = 0;
	// for(auto& i : newEdgeCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newEdgeCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newEdgeCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newEdgeCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newEdgeCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	// newPointNum = 0;
	// for(auto& i : newFaceCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newFaceCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newFaceCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newFaceCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newFaceCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	// newPointNum = 0;
	// for(auto& i : newCellCenterPointNum){
		// if(i != -1){
			// mesh.addPoint();
			// mesh.points.back().x = newCellCenterPointsXYZ[newPointNum][0];
			// mesh.points.back().y = newCellCenterPointsXYZ[newPointNum][1];
			// mesh.points.back().z = newCellCenterPointsXYZ[newPointNum][2];
			// mesh.points.back().level = newCellCenterPointsLevel[newPointNum];
			// ++newPointNum;
		// }
	// }
	
	
	
	
	// //====================================================
	// // AMR : 셀 포인트로 넘버링된 "내부"면 - 셀 연결
	// // AMR : 셀 포인트로 넘버링된 "내부"면 - 포인트 연결
	
	// // 내부 면
	
	
	// vector<int> newCellNumbering(mesh.cells.size(),-1);
	// vector<vector<int>> cellVertexPoints;
	// vector<vector<int>> cellEdges(mesh.cells.size(),vector<int>(0,0));
	// int newCellNum = 0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) {
			
			// newCellNumbering[i] = newCellNum;
			
			// newCellNum += 1;
			
			// vector<int> tmpPoints;
			// cellVertexPoints.push_back(tmpPoints);
			
			// continue;
		// }
		
		// auto& cell = mesh.cells[i];
		// int level = cell.level;
		
		// vector<int> tmpPoints;
		// for(auto& j : cell.points){
			// if( mesh.points[j].level <= level ){
				// tmpPoints.push_back(j);
			// }
		// }
		// cellVertexPoints.push_back(tmpPoints);
		
		// vector<int> tmpCellEdges;
		// for(auto& j : cell.faces){
			// for(auto& k : facesEdges[j]){
				// if ( std::find( tmpCellEdges.begin(), tmpCellEdges.end(), k ) 
					// == tmpCellEdges.end() ) {
					// tmpCellEdges.push_back(k);
				// }
			// }
			
		// }
		// cellEdges[i] = tmpCellEdges;
		
		// newCellNumbering[i] = newCellNum;
		
		// newCellNum += tmpPoints.size();
		
	// }
	
	
	// vector<vector<int>> faceVertexPoints;
	// vector<vector<int>> faceEdgeCenterPoints;
	// // vector<int> faceCenterPoints;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];
		// int level = face.level;
		
		// vector<int> tmpPoints0;
		// for(auto& j : face.points){
			// if( mesh.points[j].level <= level ){
				// tmpPoints0.push_back(j);
			// }
		// }
		// faceVertexPoints.push_back(tmpPoints0);
			
		// if(boolFaceRefine[i] == false) {
			
			// vector<int> tmpPoints1;
			// faceEdgeCenterPoints.push_back(tmpPoints1);
			
			// // faceCenterPoints.push_back(newFaceCenterPointNum[i]);
			
			// continue;
		// }
		
		// vector<int> tmpPoints1;
		// int tmpNum = 0;
		// for(auto& j : face.points){
			// if( mesh.points[j].level == level+1 ){
				// tmpPoints1.push_back(j);
			// }
			// int edgeNum = facesEdges[i][tmpNum];
			// if( 
			// boolEdgeRefine[edgeNum] == true && 
			// edgeLevel[edgeNum] == level ){
				// tmpPoints1.push_back( newEdgeCenterPointNum[edgeNum] );
			// }
			// ++tmpNum;
		// }
		// faceEdgeCenterPoints.push_back(tmpPoints1);
		
		// // faceCenterPoints.push_back(newFaceCenterPointNum[i]);
		
	// }
	
	
	
	// vector<int> newIntFaceOwner;
	// vector<int> newIntFaceNeighbour;
	// vector<vector<int>> newIntFacePoints;
	// vector<vector<int>> newCellIntEdgeCenterPoints(mesh.cells.size(),vector<int>(0,0));
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == false) continue;
		
		// auto& cell = mesh.cells[i];
		// int level = cell.level;
		
		// vector<int> tmpIntEdges;
		// for(auto& j : cellEdges[i]){
			
			// if( edgeLevel[j] != level+1 ) continue;
			
			// int point0 = edgesPoints[j][0];
			// int point1 = edgesPoints[j][1];
			
			// bool boolContinue=true;
			// for(auto& k : cellVertexPoints[i]){
				// if(point0 == k || point1 == k){
					// boolContinue=false;
					// break;
				// }
			// }
			// if(boolContinue){
				// tmpIntEdges.push_back(j);
			// }
		// }
		
		// vector<int> tmpEdgeCenterPoint;
		// vector<int> tmpFaceCenterPoint;
		// for(auto& j : tmpIntEdges){
			// int point0 = edgesPoints[j][0];
			// int point1 = edgesPoints[j][1];
			// int matchNum0 = 0;
			// int matchNum1 = 0;
			// for(auto& k : tmpIntEdges){
				// int point2 = edgesPoints[k][0];
				// int point3 = edgesPoints[k][1];
				// if(point0==point2) ++matchNum0;
				// if(point0==point3) ++matchNum0;
				// if(point1==point2) ++matchNum1;
				// if(point1==point3) ++matchNum1;
			// }
			// if(matchNum0>=3){
				// if ( std::find( tmpFaceCenterPoint.begin(), tmpFaceCenterPoint.end(), point0 ) 
					// == tmpFaceCenterPoint.end() ) {
					// tmpFaceCenterPoint.push_back(point0);
				// }
			// }
			// else{
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), point0 ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(point0);
				// }
			// }
			// if(matchNum1>=3){
				// if ( std::find( tmpFaceCenterPoint.begin(), tmpFaceCenterPoint.end(), point1 ) 
					// == tmpFaceCenterPoint.end() ) {
					// tmpFaceCenterPoint.push_back(point1);
				// }
			// }
			// else{
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), point1 ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(point1);
				// }
			// }
		// }
		
		// for(auto& j : cell.faces){
			// auto& face = mesh.faces[j];
			
			// if( face.level != level ) continue;
			
			// tmpFaceCenterPoint.push_back(newFaceCenterPointNum[j]);
			
			// for(auto& k : faceEdgeCenterPoints[j]){
				// if ( std::find( tmpEdgeCenterPoint.begin(), tmpEdgeCenterPoint.end(), k ) 
					// == tmpEdgeCenterPoint.end() ) {
					// tmpEdgeCenterPoint.push_back(k);
				// }
			// }
		// }
		
		
		

		
		
		// vector<vector<int>> tmpIntFaceCenterPoints(tmpEdgeCenterPoint.size(),vector<int>(0,0));
		// vector<vector<int>> tmpIntFaceOwnNgb(tmpEdgeCenterPoint.size(),vector<int>(0,0));
		// int tmpNum = 0;
		// for(auto& j : cell.faces){
			// auto& face = mesh.faces[j];
			
			// if( face.level == level ) {
				
				// vector<int> tmpTmpPoints;
				// int tmpNumber=0;
				// for(auto& k : faceEdgeCenterPoints[j]){
					// int saveNum=0;
					// for(auto& l : tmpEdgeCenterPoint){
						// // cout << k << " " << l << endl;
						// if(k==l) break;
						// ++saveNum;
					// }
					
					
					// int firNum = tmpNumber;
					// int lasNum = tmpNumber+1;
					// if( tmpNumber == faceEdgeCenterPoints[j].size()-1 ){
						// lasNum = 0;
					// }
					// int tmpVerP0 = faceVertexPoints[j][firNum];
					// int tmpVerP1 = faceVertexPoints[j][lasNum];
					
					
					// // cout << tmpVerP0 << " " << tmpVerP1 << endl;
					
					
					
	// // cout << tmpIntFaceOwnNgb.size() << " " << saveNum << endl;
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpVerP0 ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						
						// tmpIntFaceOwnNgb[saveNum].push_back( tmpVerP0 );
						
					// }
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpVerP1 ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						
						// tmpIntFaceOwnNgb[saveNum].push_back( tmpVerP1 );
						
					// }
	// // cout << "BBBBBBBB" << endl;
					
					// int tmpCellFaceCenterPoint = newFaceCenterPointNum[j];
					// if ( std::find( tmpIntFaceCenterPoints[saveNum].begin(), 
									// tmpIntFaceCenterPoints[saveNum].end(), 
									// tmpCellFaceCenterPoint ) 
						// == tmpIntFaceCenterPoints[saveNum].end() ) {
						// tmpIntFaceCenterPoints[saveNum].push_back(tmpCellFaceCenterPoint);
					// }
					
					
					// ++tmpNumber;
					
				// }
				
			
			// }
			// else{
				
				// int tmpCellFaceVertexPoint=-1;
				// int tmpCellFaceCenterPoint=-1;
				// vector<int> tmpFaceEdgeCenterPoint;
				
				// bool boolBreak=false;
				// for(auto& k : face.points){
					// for(auto& l : cellVertexPoints[i]){
						// if(k==l){
							// tmpCellFaceVertexPoint = k;
							// boolBreak=true;
							// break;
						// }
					// }
					// if(boolBreak) break;
				// }
				// boolBreak=false;
				// for(auto& k : face.points){
					// for(auto& l : tmpFaceCenterPoint){
						// if(k==l){
							// tmpCellFaceCenterPoint = k;
							// boolBreak=true;
							// break;
						// }
					// }
					// if(boolBreak) break;
				// }
				// for(auto& k : face.points){
					// for(auto& l : tmpEdgeCenterPoint){
						// if(k==l){
							// tmpFaceEdgeCenterPoint.push_back(k);
						// }
						// if(tmpFaceEdgeCenterPoint.size()==2) break;
					// }
					// if(tmpFaceEdgeCenterPoint.size()==2) break;
				// }
				
				
				// for(auto& k : tmpFaceEdgeCenterPoint){
					// int saveNum=0;
					// for(auto& l : tmpEdgeCenterPoint){
						// if(k==l) break;
						// ++saveNum;
					// }
					
					// if ( std::find( tmpIntFaceOwnNgb[saveNum].begin(), 
									// tmpIntFaceOwnNgb[saveNum].end(), 
									// tmpCellFaceVertexPoint ) == tmpIntFaceOwnNgb[saveNum].end() ) {
						// tmpIntFaceOwnNgb[saveNum].push_back(tmpCellFaceVertexPoint);
					// }
					
					// if ( std::find( tmpIntFaceCenterPoints[saveNum].begin(), 
									// tmpIntFaceCenterPoints[saveNum].end(), 
									// tmpCellFaceCenterPoint ) == tmpIntFaceCenterPoints[saveNum].end() ) {
						// tmpIntFaceCenterPoints[saveNum].push_back(tmpCellFaceCenterPoint);
					// }
					
				// }
				
			// }
			
		// }
		
		
		// newCellIntEdgeCenterPoints[i] = tmpEdgeCenterPoint;
		
		
		
		// // input
		// int cellNumber = newCellNumbering[i];
		// for(int j=0; j<tmpEdgeCenterPoint.size(); ++j){
			
			
			// int ownIPoint = tmpIntFaceOwnNgb[j][0];
			// int ngbIPoint = tmpIntFaceOwnNgb[j][1];
			
			// // cout << tmpEdgeCenterPoint[j] << " " << ownIPoint << " " << ngbIPoint << endl;
			
			// int ownGIPoint = 0;
			// for(auto& k : cellVertexPoints[i]){
				// if(k==ownIPoint) break;
				// ++ownGIPoint;
			// }
			// int ngbGIPoint = 0;
			// for(auto& k : cellVertexPoints[i]){
				// if(k==ngbIPoint) break;
				// ++ngbGIPoint;
			// }
			
			
			// newIntFaceOwner.push_back(cellNumber+ownGIPoint);
			// newIntFaceNeighbour.push_back(cellNumber+ngbGIPoint);
			
			// vector<int> tmpPoints;
			// int edgeCentP = tmpEdgeCenterPoint[j];
			// int faceCentP0 = tmpIntFaceCenterPoints[j][0];
			// int faceCentP1 = tmpIntFaceCenterPoints[j][1];
			// int cellCentP = newCellCenterPointNum[i];
			// tmpPoints.push_back(edgeCentP);
			// tmpPoints.push_back(faceCentP0);
			// tmpPoints.push_back(cellCentP);
			// tmpPoints.push_back(faceCentP1);
			
			
			
			
			// // 포인트 순서 판단
			// vector<double> vectorEdge;
			// vectorEdge.push_back( mesh.points[ngbIPoint].x - mesh.points[edgeCentP].x );
			// vectorEdge.push_back( mesh.points[ngbIPoint].y - mesh.points[edgeCentP].y );
			// vectorEdge.push_back( mesh.points[ngbIPoint].z - mesh.points[edgeCentP].z );
			
			// vector<double> Vx, Vy, Vz;
			// for(auto& k : tmpPoints){
				// Vx.push_back(mesh.points[k].x);
				// Vy.push_back(mesh.points[k].y);
				// Vz.push_back(mesh.points[k].z);
			// }
			
			// vector<double> vectorNormFace;
			// double dummy1;
			// geometric.calcUnitNormals_Area3dPolygon(4, Vx,Vy,Vz, vectorNormFace, dummy1);
			
			// double cosTheta = 
				// vectorNormFace[0]*vectorEdge[0]+
				// vectorNormFace[1]*vectorEdge[1]+
				// vectorNormFace[2]*vectorEdge[2];
				
			// if(cosTheta < 0.0){
				// std::reverse(tmpPoints.begin(), tmpPoints.end());
			// }
			
			
			// newIntFacePoints.push_back(tmpPoints);
			
		// }
		
		
	// }
	
	// // cout << "BBBBBBBB" << endl;
	
	
	

	
	// if(rank==0) cout << "face-point start" << endl;
	
	// //====================================================
	// // AMR : 면 - 포인트 연결 
	// vector<vector<int>> newFacePoints;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// if(boolFaceRefine[i] == false) {
			// vector<int> tmpNewFacePoints(0,0);
			// newFacePoints.push_back(tmpNewFacePoints);
			// continue;
		// }
		
		// auto& face = mesh.faces[i];
		// int level = face.level;
		// // vector<int> tmpPointsRefine;
		// // vector<int> tmpPointsNumOrigin;
		// // for(int j=0; j<face.points.size(); ++j){
			// // int iPoint = face.points[j];
			// // if( mesh.points[iPoint].level == level+1 ){
				// // tmpPointsRefine.push_back(iPoint);
			// // }
			// // if( mesh.points[iPoint].level <= level ){
				// // tmpPointsNumOrigin.push_back(j);
			// // }
		// // }
		// // for(auto& j : tmpPointsRefine){
			// // cout << " pr :" << j << endl;
		// // }
		
		// vector<int> tmpTotalPoints;
		// for(int j=0; j<face.points.size(); ++j){
			// int iPoint = face.points[j];
			// tmpTotalPoints.push_back(iPoint);
			
			// int tmpIEdge = facesEdges[i][j];
			// if(
			// boolEdgeRefine[tmpIEdge] == true
			// ){
				// int tmpEdgeCentP = newEdgeCenterPointNum[tmpIEdge];
				// tmpTotalPoints.push_back(tmpEdgeCentP);
			// }
		// }
		
		// vector<int> tmpRefineVertexPoints;
		// for(int j=0; j<face.points.size(); ++j){
			// int iPoint = face.points[j];
			// if(mesh.points[iPoint].level <= level+1){
				// tmpRefineVertexPoints.push_back(iPoint);
			// }
			
			// int tmpIEdge = facesEdges[i][j];
			// if(
			// boolEdgeRefine[tmpIEdge] == true &&
			// edgeLevel[tmpIEdge] == level
			// ){
				// int tmpEdgeCentP = newEdgeCenterPointNum[tmpIEdge];
				// tmpRefineVertexPoints.push_back(tmpEdgeCentP);
			// }
			// // cout << iPoint << " " << edgeLevel[tmpIEdge] << " " << level << " " << boolEdgeRefine[tmpIEdge] << " " << 
			// // mesh.points[edgesPoints[tmpIEdge][0]].level << " " << mesh.points[edgesPoints[tmpIEdge][1]].level << " " << endl;
		// }
			
		// // for(int j=0; j<tmpRefineVertexPoints.size(); ++j){
			// // cout << j << " " << tmpRefineVertexPoints[j] << endl;
		// // }
		// // for(int j=0; j<tmpTotalPoints.size(); ++j){
			// // cout << j << " " << tmpTotalPoints[j] << endl;
		// // }
		
		// vector<vector<int>> tmpSubEdgePoints;
		// vector<int> tmpSubPoints;
		// int iSubEdges = 1;
		// for(int j=0; j<tmpTotalPoints.size(); ++j){
			// int iPoint = tmpTotalPoints[j];
			// tmpSubPoints.push_back(iPoint);
			
			// if(
			// tmpRefineVertexPoints[iSubEdges] == iPoint
			// ){
				// tmpSubEdgePoints.push_back(tmpSubPoints);
				// tmpSubPoints.clear();
				// tmpSubPoints.push_back(iPoint);
				// ++iSubEdges;
				// if(tmpRefineVertexPoints.size()-1 < iSubEdges){
					// iSubEdges = 0;
				// }
			// }
		// }
		// tmpSubPoints.push_back(tmpTotalPoints[0]);
		// tmpSubEdgePoints.push_back(tmpSubPoints);
		
		
		
		// // for(int j=0; j<tmpSubEdgePoints.size(); ++j){
			// // for(auto& k : tmpSubEdgePoints[j]){
				// // cout << j << " " << k << endl;
			// // }
		// // }
		
		// int faceRefineSize = tmpSubEdgePoints.size()/2;
		
		// // cout << faceRefineSize << " " << tmpRefineVertexPoints.size() << " " << faceVertexPoints[i].size() << endl;
		
		// for(int j=0; j<faceRefineSize; ++j){
			// int fEdgeSet = j*2;
			// int lEdgeSet = j*2-1;
			// if(j==0) lEdgeSet = tmpSubEdgePoints.size()-1;
			
			// vector<int> tmpNewFacePoints;
			
			// for(int k=0; k<tmpSubEdgePoints[fEdgeSet].size(); ++k){
				// int tmpPoi = tmpSubEdgePoints[fEdgeSet][k];
				// tmpNewFacePoints.push_back(tmpPoi);
			// }
			
			// tmpNewFacePoints.push_back( newFaceCenterPointNum[i] );
			
			// for(int k=0; k<tmpSubEdgePoints[lEdgeSet].size()-1; ++k){
				// int tmpPoi = tmpSubEdgePoints[lEdgeSet][k];
				// tmpNewFacePoints.push_back(tmpPoi);
			// }
			
			// newFacePoints.push_back(tmpNewFacePoints);
			
			// // for(auto& k : tmpNewFacePoints){
				// // if(k==2) cout << j << " " << k << endl;
			// // }
			
		// }
		
	// }
	
	
	// int test_newFacesNum = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(boolFaceRefine[i] == true){
			// for(int j=0; j<faceVertexPoints[i].size(); ++j){
				// ++test_newFacesNum;
			// }
		// }
		// else{
			// ++test_newFacesNum;
		// }
		
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 포인트로 넘버링된 "겉" 면 - 셀 연결
	// vector<int> newFaceOwner;
	// vector<int> newFaceNeighbour;
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// int tmpNewCellNumL = newCellNumbering[face.owner];
		
		// vector<int> tmpPointsSave;
		


		// if( 
		// boolFaceRefine[i]==true &&
		// boolCellRefine[face.owner] == true
		// ){
			// for(auto& k : face.points){
				// int tmpNum = 0;
				// int level = mesh.cells[face.owner].level;
				// int pLevel = mesh.points[k].level;
				// if(pLevel > level) continue;
				// // for(auto& l : mesh.cells[face.owner].points){
				// for(auto& l : cellVertexPoints[face.owner]){
					// if(k==l){
						// tmpPointsSave.push_back(tmpNum);
						// break;
					// }
					// ++tmpNum;
				// }
			// }
		// }
		// else if( 
		// boolFaceRefine[i]==false &&
		// boolCellRefine[face.owner] == true
		// ){
			// vector<int> tmpPoint; 
			// for(auto& k : faceVertexPoints[i]){
				// int tmpPNum = 0;
				// for(auto& l : cellVertexPoints[face.owner]){
					// if(k==l){
						// tmpPoint.push_back(tmpPNum);
						// break;
					// }
					// ++tmpPNum;
				// }
				// if(tmpPoint.size() == 1) break;
			// }
			// for(auto& k : tmpPoint){
				// tmpPointsSave.push_back(k);
			// }
		// }
		// else if( 
		// boolFaceRefine[i]==true &&
		// boolCellRefine[face.owner] == false
		// ){
			// for(auto& k : face.points){
				// int level = mesh.cells[face.owner].level;
				// int pLevel = mesh.points[k].level;
				// if(pLevel > level) continue;
				// tmpPointsSave.push_back(0);
			// }
		// }
		// else{
			// tmpPointsSave.push_back(0);
		// }
		
		// for(auto& k : tmpPointsSave){
				// // cout << " L " << tmpNewCellNumL << " " << k << endl;
			// newFaceOwner.push_back( tmpNewCellNumL + k );
		// }
		
		
		
		
			
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // cout << "INTERNAL_FACE" << endl;
			
			// int tmpNewCellNumR = newCellNumbering[face.neighbour];
			
			// tmpPointsSave.clear();
			
			// if( 
			// boolFaceRefine[i]==true &&
			// boolCellRefine[face.neighbour] == true
			// ){
				// for(auto& k : face.points){
					// int tmpNum = 0;
					// int level = mesh.cells[face.neighbour].level;
					// int pLevel = mesh.points[k].level;
					// if(pLevel > level) continue;
					// // for(auto& l : mesh.cells[face.neighbour].points){
					// for(auto& l : cellVertexPoints[face.neighbour]){
						// if(k==l){
							// tmpPointsSave.push_back(tmpNum);
							// break;
						// }
						// ++tmpNum;
					// }
				// }
			// }
			// else if( 
			// boolFaceRefine[i]==false &&
			// boolCellRefine[face.neighbour] == true
			// ){
				// vector<int> tmpPoint; 
				// for(auto& k : faceVertexPoints[i]){
					// int tmpPNum = 0;
					// for(auto& l : cellVertexPoints[face.neighbour]){
						// if(k==l){
							// tmpPoint.push_back(tmpPNum);
							// break;
						// }
						// ++tmpPNum;
					// }
					// if(tmpPoint.size() == 1) break;
				// }
				// for(auto& k : tmpPoint){
					// tmpPointsSave.push_back(k);
				// }
			// }
			// else if( 
			// boolFaceRefine[i]==true &&
			// boolCellRefine[face.neighbour] == false
			// ){
				// for(auto& k : face.points){
					// int level = mesh.cells[face.neighbour].level;
					// int pLevel = mesh.points[k].level;
					// if(pLevel > level) continue;
					// tmpPointsSave.push_back(0);
				// }
			// }
			// else{
				// tmpPointsSave.push_back(0);
			// }
			
			
			
			// for(auto& k : tmpPointsSave){
				// // cout << " R " << tmpNewCellNumR + k << endl;
				// newFaceNeighbour.push_back( tmpNewCellNumR + k );
			// }
		
			
		// }
		// else{
			// // cout << "BC_FACE" << endl;
			// if( boolFaceRefine[i]==true ){
				// for(auto& k : faceVertexPoints[i]){
					// newFaceNeighbour.push_back( -1 );
				// }
			// }
			// else{
				// newFaceNeighbour.push_back( -1 );
			// }
		// }
		
	// }
	
	
	
	// //====================================================
	// // AMR : 셀 그룹핑
	// vector<int> cellGroupping;
	// vector<int> cellLeveling;
	// vector<vector<double>> cellVariables;
	// newPointNum = 0;
	// // for(auto& i : newCellCenterPointNum){
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// int level = mesh.cells[i].level;
		// int group = mesh.cells[i].group;
		// if(boolCellRefine[i] == true){
			// // cout << cellVertexPoints[i].size() << " AAAAAAAA " << endl;
			// for(auto& j : cellVertexPoints[i]){
				// cellLeveling.push_back(level+1);
				// cellGroupping.push_back(group);
				
				
				// vector<double> tmpVariables;
				// for(auto& k : cell.var){
					// tmpVariables.push_back(k);
				// }
				// cellVariables.push_back(tmpVariables);
			// }
		// }
		// else{
			// vector<double> tmpVariables;
			// cellLeveling.push_back(level);
			// cellGroupping.push_back(group);
			// for(auto& k : cell.var){
				// tmpVariables.push_back(k);
			// }
			// cellVariables.push_back(tmpVariables);
		// }
	// }
	
	
	
	// if(rank==0) cout << "face-add start" << endl;
	
	// //====================================================
	// // AMR : 면 추가 (면 -> new 셀 내부 면 -> proc면 -> B.C.면)
	// // AMR : 면 - 셀 리넘버링
	// // AMR : 면 - 포인트 리넘버링
	// // AMR : 면 타입 셋팅
	// // AMR : 면 바운더리 셋팅
	
	// SEMO_Mesh_Builder newMesh;
	
	// // 원래 내부 면
	// newPointNum = 0;
	// int newFacesNum = 0;
	// // cout << "AAAAAAAAAAAAAAA" << endl;
	// // for(auto& i : newFaceCenterPointNum){
	// for(int i=0; i<mesh.faces.size(); ++i){
		// int level = mesh.faces[i].level;
		// int group = mesh.faces[i].group;
		
		// // cout << level << " " << group << endl;
		
		// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			// if(boolFaceRefine[i] == true){
				// // if(faceVertexPoints.size()-1 < i) cout << "ERROR!!!!!!!!!!" << endl;
				// // if(faceVertexPoints[i].size()>8) cout << faceVertexPoints[i].size() << endl;
				// for(int j=0; j<faceVertexPoints[i].size(); ++j){
					// newMesh.addFace();
					// newMesh.faces.back().group = group;
					// newMesh.faces.back().level = level + 1;
					
					// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
					// newMesh.faces.back().neighbour = newFaceNeighbour[newFacesNum];
					
					// for(auto& k : newFacePoints[newFacesNum]){
						// // if(k==2) cout << k << endl;
						// newMesh.faces.back().points.push_back(k);
					// }
					
					// newMesh.faces.back().setType(mesh.faces[i].getType());
					
					// ++newFacesNum;
				// }
				
			// }
			// else{
				// newMesh.addFace();
				// newMesh.faces.back().group = group;
				// newMesh.faces.back().level = level;
				
				// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
				// newMesh.faces.back().neighbour = newFaceNeighbour[newFacesNum];
				
				// newMesh.faces.back().setType(mesh.faces[i].getType());
				
				// for(int j=0; j<facesEdges[i].size(); ++j){
					// int tmpNum = facesEdges[i][j];
					// if(boolEdgeRefine[tmpNum] == true){
						// int point0 = mesh.faces[i].points[j];
						// newMesh.faces.back().points.push_back(point0);
						// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
					// }
					// else{
						// int point0 = mesh.faces[i].points[j]; 
						// newMesh.faces.back().points.push_back(point0);
					// }
				// }
				
				// ++newFacesNum;
			// }
		// }
	// }
	
	
	// // cout << "BBBBBBBBBBB" << endl;
	// // MPI_Barrier(MPI_COMM_WORLD);
	
	
	// // 셀 내부 면
	// int newIntFaceNum = 0;
	// // for(auto& i : newCellCenterPointNum){
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// if(boolCellRefine[i] == true){
			
			// int level = mesh.cells[i].level;
			// int group = mesh.faces[i].group;
			
			// for(int j=0; j<newCellIntEdgeCenterPoints[i].size(); ++j){
				// newMesh.addFace();
				// newMesh.faces.back().group = mesh.faces.size() + newIntFaceNum;
				// newMesh.faces.back().level = level + 1;
				// newMesh.faces.back().owner = newIntFaceOwner[newIntFaceNum];
				// newMesh.faces.back().neighbour = newIntFaceNeighbour[newIntFaceNum];
				// newMesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
				
				// for(auto& k : newIntFacePoints[newIntFaceNum]){
					// newMesh.faces.back().points.push_back(k);
				// }
					
				// ++newIntFaceNum;
			// }
		// }
	// }
	
	
	
	
	// if(rank==0) cout << "proc & bc face start" << endl;
	
	
	// // 프록, 경계면
	// vector<int> newBoundaryStr;
	// int BCnum = 0;
	// int test_num2 = 0;
	
	// // for(auto& i : newFaceCenterPointNum){
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];
		
		
		// int procNeighbNo = mesh.boundary[BCnum].neighbProcNo;
		
		// if(BCnum < mesh.boundary.size()){
			// if( mesh.boundary[BCnum].startFace == i ){
				
				// // if(mesh.boundary.size() <= BCnum) {
					// // cout << "| #Error : mesh.boundary.size() <= BCnum" << endl;
				// // }
				
				// newMesh.addBoundary();
				// newMesh.boundary[BCnum].name = mesh.boundary[BCnum].name;
				// newMesh.boundary[BCnum].type = mesh.boundary[BCnum].type;
				// for(auto& j : mesh.boundary[BCnum].var){
					// newMesh.boundary[BCnum].var.push_back(j);
				// }
				// newMesh.boundary[BCnum].startFace = newMesh.faces.size();
				// if(BCnum != 0) {
					// newMesh.boundary[BCnum-1].nFaces = 
						// newMesh.faces.size()-newMesh.boundary[BCnum-1].startFace;
				// }
				// newMesh.boundary[BCnum].myProcNo = rank;
				// newMesh.boundary[BCnum].neighbProcNo = mesh.boundary[BCnum].neighbProcNo;
				// ++BCnum;
			// }
		// }
		
		// if(face.getType() != SEMO_Types::INTERNAL_FACE){
			// if(boolFaceRefine[i] == true){
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // cout << rank << " yes refine " << test_num2 << endl;
					// // ++test_num2;
				// // }
				
				// int level = face.level;
				// int group = face.group;
				
				// // int nChildCells = mesh.faces[i].points.size();
				// int nChildCells = faceVertexPoints[i].size();
				
				
				// vector<int> vecNewPointNumNum;
				// for(int j=0; j<nChildCells; ++j){
					// vecNewPointNumNum.push_back(newFacesNum);
					// ++newFacesNum;
				// }
				
				// // if(
				// // mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
				// // cout << nChildCells << endl;
				// // }
				// if(
				// face.getType() == SEMO_Types::PROCESSOR_FACE &&
				// rank > procNeighbNo
				// ){
					// std::reverse(vecNewPointNumNum.begin(), vecNewPointNumNum.end());
				// }
				
				
				// for(auto& j : vecNewPointNumNum){
					// newMesh.addFace();
					// newMesh.faces.back().group = group;
					// newMesh.faces.back().level = level + 1;
					
					// newMesh.faces.back().owner = newFaceOwner[j];
					// newMesh.faces.back().neighbour = -1;
					// newMesh.faces.back().setType(face.getType());
					
					// for(auto& k : newFacePoints[j]){
						// newMesh.faces.back().points.push_back(k);
					// }
					
					// if(
					// face.getType() == SEMO_Types::PROCESSOR_FACE &&
					// rank > procNeighbNo
					// ){
						// std::reverse(newMesh.faces.back().points.begin()+1, 
									 // newMesh.faces.back().points.end());
						// std::reverse(newMesh.faces.back().points.begin(), 
									 // newMesh.faces.back().points.end());
					// }
					
				// }
				
				
			// }
			// else{
				// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					// // cout << rank << " no refine " << test_num2 << endl;
					// // ++test_num2;
				// // }
				// newMesh.addFace();
				// newMesh.faces.back().group = face.group;
				// newMesh.faces.back().level = face.level;
				// newMesh.faces.back().owner = newFaceOwner[newFacesNum];
				// newMesh.faces.back().neighbour = -1;
				// newMesh.faces.back().setType(face.getType());
				
				// // cout << face.points.size() << " " << facesEdges[i].size() << endl;
				
				// if(
				// face.getType() == SEMO_Types::PROCESSOR_FACE &&
				// rank > procNeighbNo
				// ){
					
					// int lastEdNum = facesEdges[i][face.points.size()-1];
					// if(boolEdgeRefine[lastEdNum] == true){
						// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[lastEdNum]);
						// // ++test_num2;
					// }
						
					// for(int j=0; j<face.points.size(); ++j){
						// int tmpNum = facesEdges[i][j];
						// if(boolEdgeRefine[tmpNum] == true){
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
							// if(j!=face.points.size()-1) 
								// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
							// // ++test_num2;
						// }
						// else{
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
						// }
					// }
				// }
				// else{
					// for(int j=0; j<face.points.size(); ++j){
						// int tmpNum = facesEdges[i][j];
						// if(boolEdgeRefine[tmpNum] == true){
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
							// newMesh.faces.back().points.push_back(newEdgeCenterPointNum[tmpNum]);
							// // ++test_num2;
						// }
						// else{
							// int point0 = face.points[j]; //edgesPoints[j][0];
							// newMesh.faces.back().points.push_back(point0);
						// }
					// }
				// }
				
				
				
				// ++newFacesNum;
				
				
			// }
		// }
	// }
	// if(mesh.boundary.size() < BCnum) {
		// cout << "| #Error 2 : mesh.boundary.size() <= BCnum" << endl;
	// }
	// newMesh.boundary[BCnum-1].nFaces = 
		// newMesh.faces.size()-newMesh.boundary[BCnum-1].startFace;
		
		
		
	// // for(int i=0; i<mesh.faces.size(); ++i){
		
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // if(boolFaceRefine[i] == true){
			// // }
			// // else{
				// // for(int j=0; j<face.points.size(); ++j){
					// // int tmpNum = facesEdges[i][j];
					// // int iNewPoint0 = face.points[j];
					// // // if(rank==1) cout << j << " P : " << test_num2 << " " << mesh.points[iNewPoint0].x << " " << mesh.points[iNewPoint0].y << " " << mesh.points[iNewPoint0].z << " " << endl;
					// // if(boolEdgeRefine[tmpNum] == true){
						// // int iNewPoint = newEdgeCenterPointNum[tmpNum];
						// // // if(rank==1) cout << j << " : " << test_num2 << " " << mesh.points[iNewPoint].x << " " << mesh.points[iNewPoint].y << " " << mesh.points[iNewPoint].z << " " << endl;
						// // ++test_num2;
					// // }
				// // }
			// // }
			
		// // }
	// // }
		
		
	// // cout << "AGAGAG " << test_num2 << endl;
	
	
	// if(rank==0) cout << "mesh copy start" << endl;
		
	// //====================================================
	// // 원래 메쉬 클리어 및 복사
	// mesh.cells.clear();
	// mesh.faces.clear();
	// // mesh.points.clear();
	// mesh.boundary.clear();
	
	// for(int i=0; i<newMesh.faces.size(); ++i){
		// auto& newFace = newMesh.faces[i];
		
		// mesh.addFace();
		// mesh.faces.back().group = newFace.group;
		// mesh.faces.back().level = newFace.level;
		// mesh.faces.back().owner = newFace.owner;
		// mesh.faces.back().neighbour = newFace.neighbour;
		// mesh.faces.back().setType(newFace.getType());
		// for(auto& j : newFace.points){
			// mesh.faces.back().points.push_back(j);
		// }
		
		// mesh.faces.back().varL.resize(nTotalFaceLRVar,0.0);
		// mesh.faces.back().varR.resize(nTotalFaceLRVar,0.0);
		// mesh.faces.back().var.resize(nTotalFaceVar,0.0);
		
	// }
	
	// for(int i=0; i<newMesh.boundary.size(); ++i){
		// auto& newBoundary = newMesh.boundary[i];
		// mesh.addBoundary();
		// mesh.boundary.back().name = newBoundary.name;
		// mesh.boundary.back().type = newBoundary.type;
		// for(auto& j : newBoundary.var){
			// mesh.boundary.back().var.push_back(j);
		// }
		// mesh.boundary.back().startFace = newBoundary.startFace;
		// mesh.boundary.back().nFaces = newBoundary.nFaces;
		// mesh.boundary.back().myProcNo = newBoundary.myProcNo;
		// mesh.boundary.back().neighbProcNo = newBoundary.neighbProcNo;
		
	// }
	
	
	
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		
		// // cout << face.owner << " " << face.neighbour << endl;
		
		// // // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // // cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
		// // // }
	// // }
	// // cout << "cccccc" << endl;

	// mesh.check();
	
	
	
	
	// // cout << "cccccc" << endl;
	
	// // mesh.checkMatchingProcessorFace();
	
	// // cout << "ddddd" << endl;
	
	
	
	
	// mesh.buildCells();
	
	// // mesh.setFaceTypes();

	// mesh.buildLists();
	
	// mesh.connectCelltoFaces();
	
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// // int minLevel = 999999999;
		// // for(int j=0; j<cell.faces.size(); ++j){
			// // auto& face = mesh.faces[cell.faces[j]];
			// // minLevel = min(minLevel,face.level);
		// // }
		// // cell.level = minLevel;
		
		// cell.level = cellLeveling[i];
		// cell.group = cellGroupping[i];
		
		
		// // cout << "CELL VAR " << cellVariables[i].size() << endl;
		// for(auto& k : cellVariables[i]){
			// cell.var.push_back(k);
		// }
		
		
		// // cout << " G : " << cellGroupping[i] << " " << cellGroupping.size() << " " << mesh.cells.size() << endl;
	// }
	
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// mesh.connectCelltoPoints();
	
	// // set processor face counts
	// mesh.setCountsProcFaces();
	
	// // set processor face displacements
	// mesh.setDisplsProcFaces(); 
		
	// // mesh.informations();
	
	
	
	// // SEMO_Mesh_Save save;
	// // string tmpFile = "./Rf" + to_string(iter);
	// // // string tmpFile = "./";
	// // save.vtu(tmpFile, rank, mesh);
	
	
	// // cout << "SIZE : " << cellGroupping.size() << " " << mesh.cells.size() << endl;
		
	// // cout << "CELL" << endl;
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // auto& cell = mesh.cells[i];
		// // cout << i << " " << cell.level << " " << cell.group << endl;
	// // }
	
	// // cout << "FACE" << endl;
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // cout << i << " " << face.level << " " << face.group << endl;
	// // }
		

	
	// if(rank==0) cout << "----Refine finished-----" << endl;
// }