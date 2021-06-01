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
	vector<int> iFaceBC(mesh.faces.size(),-1);
	
	vector<int> startFace(size,0);
	vector<int> nFaces(size,0);
	vector<int> neighbProcNo(mesh.faces.size(),-1);
	
	vector<vector<int>> idNewFaceProc(size,vector<int>(0,0));
		
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		if( distr(eng) > 0.5 ){
			executeAMR[i] = true;
		}
	}
	
	// executeAMR[0] = true;
	// executeAMR[1] = false;
	
	
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
	
	num=0;
	// // cout << mesh.boundary[0].startFace << " " << mesh.boundary[0].nFaces << endl;
	// for(auto& i : mesh.faces){
		// if(i.getType() == SEMO_Types::BOUNDARY_FACE){
			// if(iFaceBC[num]==-1) {
				// cerr << num << " " << mesh.faces.size() << endl;
				// cerr << "#ERROR : iFaceBC == -1" << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
		// }
		// ++num;
	// }
	// cout << startFaceBC.size() << endl;
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
	
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // cout << mesh.cells[face.owner].level << endl;
			// // cout << mesh.cells[face.neighbour].level << endl;
			// if( mesh.cells[face.owner].level > 
				// mesh.cells[face.neighbour].level ){
				
				// executeAMR[face.owner] = false;

			// }
			// if( mesh.cells[face.owner].level < 
				// mesh.cells[face.neighbour].level ){
				
				// executeAMR[face.neighbour] = false;

			// }
		// }
	// }
	
	
	//================================================
	
		
		
	// points -> faces connection
	vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	num=0;
	for(auto& face : mesh.faces){
		for(auto& i : face.points){
			int l=0;
			for(auto& j : pointsFaces[i]){
				if(num==j) ++l;
			}
			if(l==0) {
				pointsFaces[i].push_back(num);
	// cout << num << endl;
			}
		}
		++num;
	}
	
	vector<vector<int>> faceToFaces(mesh.faces.size(),vector<int>(0,0));
	num=0;
	for(auto& point : mesh.points){
		
		for(auto& i : pointsFaces[num]){
			for(auto& j : pointsFaces[num]){
				
				if(i==j) continue;
				int l=0;
				for(auto& p0 : mesh.faces[i].points){
					for(auto& p1 : mesh.faces[j].points){
						if(p0==p1) ++l;
					}
				}
							// cout << mesh.faces.size() << " " << i << endl;
				if(l>=2){
					if(faceToFaces[i].size() > 0){
						if ( std::find( faceToFaces[i].begin(), faceToFaces[i].end(), j ) 
							== faceToFaces[i].end() ) {
							faceToFaces[i].push_back(j);
						// cout << i << " " << j << endl;
						}
					}
					else{
						faceToFaces[i].push_back(j);
						// cout << i << " " << j << endl;
					}
				}
				
			}
		}
		
		++num;
	}
	
	pointsFaces.clear();
		
	//====================================================
	// face AMR
	
	vector<int> test_num1(20,0);
	
	cout << "| Face AMR start ... ";
	
	int newNPoints = 0;
	int newNFaces = 0;
	int newNCells = 0;
	
	int executedNFaceAMR=0;
	
	vector<int> idExecutedFaceAMR;
	
	num=0;
	int orgNFaces = mesh.faces.size();
	for(int iface=0; iface<orgNFaces; ++iface){
		auto face = mesh.faces[iface];
		
		int cellPointsSize = mesh.cells[face.owner].points.size();
		int cellFacesSize = mesh.cells[face.owner].faces.size();
		
		int facePointsSize = mesh.faces[iface].points.size();
		// cout << cellPointsSize << " " << cellFacesSize << endl;
		
		bool cellFacesHexaShapes = true;
		vector<int> ownerCellShapes;
		vector<int> neighbourCellShapes;
		for(auto& i : mesh.cells[face.owner].faces){
			ownerCellShapes.push_back(mesh.faces[i].points.size());
		}
		

			// cout << iface << endl;
			// cout << cellPointsSize << " " << cellFacesSize << " " << facePointsSize << endl;
		
		// bool kindHexaShapes = true;
		bool canExecutedHexaAMR = true;
		bool boolExecuteAMR = false;
		
		
		
		// if(facePointsSize>8) canExecutedHexaAMR = false;
		
		
		if(mesh.faces[iface].getType() == SEMO_Types::INTERNAL_FACE){
			int owner = mesh.faces[iface].owner;
			for(auto& i : mesh.cells[owner].faces){
				vector<vector<int>> tmp_parPoints;
				AMR.searchParallelPointsPolygon(mesh, mesh.faces[i].points, tmp_parPoints);
				if(tmp_parPoints.size() != 4){
					canExecutedHexaAMR=false;
				}
			}
			
			int neighbour = mesh.faces[iface].neighbour;
			for(auto& i : mesh.cells[neighbour].faces){
				vector<vector<int>> tmp_parPoints;
				AMR.searchParallelPointsPolygon(mesh, mesh.faces[i].points, tmp_parPoints);
				if(tmp_parPoints.size() != 4){
					canExecutedHexaAMR=false;
				}
			}
		
			if(canExecutedHexaAMR){
				if( 
				executeAMR[face.owner] == true ||
				executeAMR[face.neighbour] == true 
				){
					boolExecuteAMR = true;
				}
			}
		}
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int owner = mesh.faces[iface].owner;
			for(auto& i : mesh.cells[owner].faces){
				vector<vector<int>> tmp_parPoints;
				AMR.searchParallelPointsPolygon(mesh, mesh.faces[i].points, tmp_parPoints);
				if(tmp_parPoints.size() != 4){
					canExecutedHexaAMR=false;
				}
			}
			
			if(canExecutedHexaAMR){
				if( 
				executeAMR[face.owner] == true
				){
					boolExecuteAMR = true;
				}
			}
		}
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int owner = mesh.faces[iface].owner;
			for(auto& i : mesh.cells[owner].faces){
				vector<vector<int>> tmp_parPoints;
				AMR.searchParallelPointsPolygon(mesh, mesh.faces[i].points, tmp_parPoints);
				if(tmp_parPoints.size() != 4){
					canExecutedHexaAMR=false;
				}
			}
			
			if(canExecutedHexaAMR){
				if( 
				executeAMR[face.owner] == true
				){
					boolExecuteAMR = true;
				}
			}
		}
		
		
		// B.C. & Proc face save
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			if( iFaceBC[iface] == -1 ){
				cerr << "#ERROR : iFaceBC[iface] == -1" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			
			}
			idNewFaceBC[ iFaceBC[iface] ].push_back( iface );
		}
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int procNum = neighbProcNo[iface];
			idNewFaceProc[procNum].push_back( iface );
		}
		
		
		
		
		if(boolExecuteAMR){
			// cout << "CONPLE1" << endl;
			// cout << "CONPLE1" << endl;
			
			bool boolContinue = false;
			
			// for(int i=0; i<faceToFaces[iface].size(); ++i){
				
			// }
			
			
			vector<int> parentPoints;
			
			int psize = face.points.size();
			for(int i=0; i<psize; ++i){
				// cout << face.points[i] << endl;
				parentPoints.push_back(face.points[i]);
			}
			
			int pointsize = mesh.points.size();
			vector<vector<int>> parPoints;
			AMR.searchParallelPointsPolygon(mesh, parentPoints, parPoints);
			
			
			for(int i=0; i<parPoints.size(); ++i){
				
				if(parPoints[i].size() >3) boolContinue=true;
				
			}
			
			if(boolContinue) continue;
			

			vector<int> vertexPoints(0,0);
			vector<int> edgeCenterPoints(0,0);
			vector<int> idNewPoints(0,0);
			int tmpNum = 0;
			for(int i=0; i<parPoints.size(); ++i){
				if(parPoints[i].size()==2){
					edgeCenterPoints.push_back(pointsize+tmpNum);
					idNewPoints.push_back(pointsize+tmpNum);
					
					// cout << pointsize << " " << tmpNum << endl;
					
					int p0 = parPoints[i][0];
					int p1 = parPoints[i][1];
					
					mesh.addPoint();
					++newNPoints;
					mesh.points.back().x = 0.5*(mesh.points[p0].x + mesh.points[p1].x);
					mesh.points.back().y = 0.5*(mesh.points[p0].y + mesh.points[p1].y);
					mesh.points.back().z = 0.5*(mesh.points[p0].z + mesh.points[p1].z);
					
					int l=0;
					int m=0;
					for(int j=0; j<vertexPoints.size(); ++j){
						if(vertexPoints[j]==p0) ++l;
						if(vertexPoints[j]==p1) ++m;
					}
					if(l==0) vertexPoints.push_back(p0);
					if(m==0) vertexPoints.push_back(p1);
					
					++tmpNum;
				}
				else if(parPoints[i].size()==3){
					edgeCenterPoints.push_back(parPoints[i][1]);
					
					// cout << parPoints[i][1] << endl;
					
					int l=0;
					int m=0;
					for(int j=0; j<vertexPoints.size(); ++j){
						if(vertexPoints[j]==parPoints[i][0]) ++l;
						if(vertexPoints[j]==parPoints[i][2]) ++m;
					}
					if(l==0) vertexPoints.push_back(parPoints[i][0]);
					if(m==0) vertexPoints.push_back(parPoints[i][2]);
					
				}
			}
			
			// cout << "CONPLE2" << endl;
			
			if(std::find(edgeCenterPoints.begin(),edgeCenterPoints.end(),-1)!=edgeCenterPoints.end()){
				cerr << "# ERROR :  edgeCenterPoints == -1" << endl;
				cerr << parPoints[0].size() << endl;
				cerr << parPoints[1].size() << endl;
				cerr << parPoints[2].size() << endl;
				cerr << parPoints[3].size() << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			

			// cout << "CONPLE1" << endl;
			// cout << parentPoints.size() << endl;
			// cout << vertexPoints.size() << endl;
			bool returnFaceAMR = false;
			for(auto& i : faceToFaces[iface]){
				if(i<=iface) continue;
				vector<int> idConnPoints;
				for(auto& j : parentPoints){
					for(auto& point : mesh.faces[i].points){
						if(j==point) idConnPoints.push_back(j);
					}
				}
				if(idConnPoints.size()==2){
					int addCenterPoint;
					if( 
					(idConnPoints[0] == vertexPoints[0] && idConnPoints[1] == vertexPoints[1]) ||
					(idConnPoints[0] == vertexPoints[1] && idConnPoints[1] == vertexPoints[0]) ){
					}
					else if( 
					(idConnPoints[0] == vertexPoints[1] && idConnPoints[1] == vertexPoints[2]) ||
					(idConnPoints[0] == vertexPoints[2] && idConnPoints[1] == vertexPoints[1]) ){
					}
					else if( 
					(idConnPoints[0] == vertexPoints[2] && idConnPoints[1] == vertexPoints[3]) ||
					(idConnPoints[0] == vertexPoints[3] && idConnPoints[1] == vertexPoints[2]) ){
					}
					else if( 
					(idConnPoints[0] == vertexPoints[3] && idConnPoints[1] == vertexPoints[0]) ||
					(idConnPoints[0] == vertexPoints[0] && idConnPoints[1] == vertexPoints[3]) ){
					}
					else{
						returnFaceAMR = true;
					}
				}
			}
			
			if(returnFaceAMR) continue;
			
			
			// cout << "CONPLE3" << endl;
			
			// // cout << endl;
			// for(auto& iii : vertexPoints){
				// cout << iii << " ";
			// }
			// cout << endl;
			
			
			if(vertexPoints.size() != 4) {
				cout << " DELETE " << endl;
				for(int i=0; i<idNewPoints.size(); ++i){
					mesh.points.pop_back();
				}
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// if(face.getType() == SEMO_Types::INTERNAL_FACE){
					// executeAMR[face.owner] = false;
					// executeAMR[face.neighbour] = false;
				// }
				// else{
					// executeAMR[face.owner] = false;
				// }
				continue;
			}
			
			
			++executedNFaceAMR;
			
			if(vertexPoints.size() != 4) {
				cerr << endl;
				cout << "┌-----------------------------------" << endl;
				cerr << "| #ERROR : vertexPoints.size() != 4, = " << vertexPoints.size() << endl;
				cerr << "| face points size : " << facePointsSize <<endl;
				cerr << "| cell points size : " << cellPointsSize << " | cell faces size : " << cellFacesSize <<endl;
				cerr << "| owner face shape : ";
				for(auto& i : ownerCellShapes){
					cerr << i << " ";
				}
				cerr << endl;
				cerr << "| neighbour face shape : ";
				for(auto& i : neighbourCellShapes){
					cerr << i << " ";
				}
				cerr << endl;
				cerr << "| vertexPoints : ";
				for(auto& i : vertexPoints){
					cerr << i << " ";
				}
				cerr << endl;
				cerr << "| parentPoints : ";
				for(auto& i : parentPoints){
					cerr << i << " ";
				}
				cerr << endl;
				for(auto& i : parPoints){
					cerr << "| parPoints : ";
					for(auto& j : i){
						cerr << j << " ";
					}
					cerr << endl;
				}
				cerr << endl;
				for(auto& i : parentPoints){
					cerr << "| parentPoints x y z : ";
					cerr << mesh.points[i].x << " " << mesh.points[i].y << " " << mesh.points[i].z << endl;
				}
				cout << "└-----------------------------------" << endl;
				cerr << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			
			// cout << "CONPLE3" << endl;
			
			// face center point
			int faceCenterPoint = pointsize+tmpNum;
			mesh.addPoint();
			++newNPoints;
			mesh.points.back().x = 0.0;
			mesh.points.back().y = 0.0;
			mesh.points.back().z = 0.0;
			for(int i=0; i<4; ++i){
				mesh.points.back().x += 0.25*mesh.points[vertexPoints[i]].x;
				mesh.points.back().y += 0.25*mesh.points[vertexPoints[i]].y;
				mesh.points.back().z += 0.25*mesh.points[vertexPoints[i]].z;
			}
			
			int parentFace = iface;
			
			idExecutedFaceAMR.push_back(iface);
			
			int nfaces = mesh.faces.size();
			vector<int> idNewFacesTmp(4,0);
			idNewFacesTmp[0] = parentFace;
			idNewFacesTmp[1] = nfaces;
			idNewFacesTmp[2] = nfaces+1;
			idNewFacesTmp[3] = nfaces+2;
			
			// if(vertexPoints[0]==0 || vertexPoints[1]==0 || vertexPoints[2]==0 || vertexPoints[3]==0 ||
			// edgeCenterPoints[0]==0 || edgeCenterPoints[1]==0 || edgeCenterPoints[2]==0 || edgeCenterPoints[3]==0 ||
			// faceCenterPoint==0){
				// cout<<  vertexPoints[0]<<" " << vertexPoints[1]<<" " << vertexPoints[2] << " " << vertexPoints[3] << endl;
				// cout<<  edgeCenterPoints[0]<<" " << edgeCenterPoints[1]<<" " << edgeCenterPoints[2] << " " << edgeCenterPoints[3] << endl;
				// cerr << parPoints[0].size() << endl;
				// cerr << parPoints[1].size() << endl;
				// cerr << parPoints[2].size() << endl;
				// cerr << parPoints[3].size() << endl;
				// cout << faceCenterPoint << endl << endl;
				// // cout << parPoints[i][1] << endl;
			// }
			
			mesh.faces[parentFace].points.clear();
			mesh.faces[parentFace].points.push_back(vertexPoints[0]);
			mesh.faces[parentFace].points.push_back(edgeCenterPoints[0]);
			mesh.faces[parentFace].points.push_back(faceCenterPoint);
			mesh.faces[parentFace].points.push_back(edgeCenterPoints[3]);
			
			// cout << "CONPLE4" << endl;
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(edgeCenterPoints[0]);
			mesh.faces.back().points.push_back(vertexPoints[1]);
			mesh.faces.back().points.push_back(edgeCenterPoints[1]);
			mesh.faces.back().points.push_back(faceCenterPoint);
			mesh.faces.back().owner = mesh.faces[parentFace].owner;
			mesh.faces.back().neighbour = mesh.faces[parentFace].neighbour;
			mesh.faces.back().setType(mesh.faces[parentFace].getType());
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(faceCenterPoint);
			mesh.faces.back().points.push_back(edgeCenterPoints[1]);
			mesh.faces.back().points.push_back(vertexPoints[2]);
			mesh.faces.back().points.push_back(edgeCenterPoints[2]);
			mesh.faces.back().owner = mesh.faces[parentFace].owner;
			mesh.faces.back().neighbour = mesh.faces[parentFace].neighbour;
			mesh.faces.back().setType(mesh.faces[parentFace].getType());
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(edgeCenterPoints[3]);
			mesh.faces.back().points.push_back(faceCenterPoint);
			mesh.faces.back().points.push_back(edgeCenterPoints[2]);
			mesh.faces.back().points.push_back(vertexPoints[3]);
			mesh.faces.back().owner = mesh.faces[parentFace].owner;
			mesh.faces.back().neighbour = mesh.faces[parentFace].neighbour;
			mesh.faces.back().setType(mesh.faces[parentFace].getType());
			
			// cout << "CONPLE5" << endl;

			if(mesh.faces[parentFace].getType() == SEMO_Types::INTERNAL_FACE){
				mesh.cells[ mesh.faces[parentFace].owner ].points.push_back(faceCenterPoint);
				mesh.cells[ mesh.faces[parentFace].neighbour ].points.push_back(faceCenterPoint);
				
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[1]);
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[2]);
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[3]);
				mesh.cells[ mesh.faces[parentFace].neighbour ].faces.push_back(idNewFacesTmp[1]);
				mesh.cells[ mesh.faces[parentFace].neighbour ].faces.push_back(idNewFacesTmp[2]);
				mesh.cells[ mesh.faces[parentFace].neighbour ].faces.push_back(idNewFacesTmp[3]);
			}
			else{
				mesh.cells[ mesh.faces[parentFace].owner ].points.push_back(faceCenterPoint);
				
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[1]);
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[2]);
				mesh.cells[ mesh.faces[parentFace].owner ].faces.push_back(idNewFacesTmp[3]);
			}
			
			// cout << "CONPLE6" << endl;
			
			if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				if( iFaceBC[parentFace] == -1 ){
					
					cerr << "#ERROR : iFaceBC[iface] == -1" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				
				}
			// cout << "CONPLE6.5" << endl;
			// cout << iFaceBC[parentFace] << " " << parentFace << " " << idNewFacesTmp[0] << " " << idNewFacesTmp[1] << " " << idNewFacesTmp[2] << " " << idNewFacesTmp[3] << endl;
				// idNewFaceBC[ iFaceBC[parentFace] ].push_back( idNewFacesTmp[0] );
				idNewFaceBC[ iFaceBC[parentFace] ].push_back( idNewFacesTmp[1] );
				idNewFaceBC[ iFaceBC[parentFace] ].push_back( idNewFacesTmp[2] );
				idNewFaceBC[ iFaceBC[parentFace] ].push_back( idNewFacesTmp[3] );
			}
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				int procNum = neighbProcNo[parentFace];
				// idNewFaceProc[procNum].push_back( idNewFacesTmp[0] );
				idNewFaceProc[procNum].push_back( idNewFacesTmp[1] );
				idNewFaceProc[procNum].push_back( idNewFacesTmp[2] );
				idNewFaceProc[procNum].push_back( idNewFacesTmp[3] );
			}
			
			
			// cout << "CONPLE7" << endl;
			
			
			// // add points at other faces
			vector<vector<int>> tempPointsFaces(4,vector<int>(0,0));
			int nnn = 0;
			for(auto& i : faceToFaces[iface]){
				// cout << "AA" << " " << nnn << " " << faceToFaces[iface].size() << endl;
				++nnn;
				if(i<=iface) continue;
				// cout << "BB" << endl;
				vector<int> idConnPoints;
				vector<int> sortConnFace;
				int l=0;
				vector<int> orders;
				for(auto& j : parentPoints){
				// cout << "BA" << endl;
					int n=0;
					for(auto& point : mesh.faces[i].points){
						if(j==point) idConnPoints.push_back(j);
						if(j==point) sortConnFace.push_back(n);
						if(j==point) orders.push_back(l);
						++n;
					}
					++l;
				// cout << "BAA" << endl;
				}
				// cout << idConnPoints.size()<< endl;
				if(idConnPoints.size()==2){
					int addCenterPoint;
					if( 
					(idConnPoints[0] == vertexPoints[0] && idConnPoints[1] == vertexPoints[1]) ||
					(idConnPoints[0] == vertexPoints[1] && idConnPoints[1] == vertexPoints[0]) ){
						addCenterPoint = edgeCenterPoints[0];
					}
					else if( 
					(idConnPoints[0] == vertexPoints[1] && idConnPoints[1] == vertexPoints[2]) ||
					(idConnPoints[0] == vertexPoints[2] && idConnPoints[1] == vertexPoints[1]) ){
						addCenterPoint = edgeCenterPoints[1];
					}
					else if( 
					(idConnPoints[0] == vertexPoints[2] && idConnPoints[1] == vertexPoints[3]) ||
					(idConnPoints[0] == vertexPoints[3] && idConnPoints[1] == vertexPoints[2]) ){
						addCenterPoint = edgeCenterPoints[2];
					}
					else if( 
					(idConnPoints[0] == vertexPoints[3] && idConnPoints[1] == vertexPoints[0]) ||
					(idConnPoints[0] == vertexPoints[0] && idConnPoints[1] == vertexPoints[3]) ){
						addCenterPoint = edgeCenterPoints[3];
					}
					else{
						continue;
						cerr << "#ERROR " << endl;
						cerr << idConnPoints[0] << " " << idConnPoints[1] << " " << endl;
						cerr << vertexPoints[0] << " " << vertexPoints[1] << " " << vertexPoints[2] << " " << vertexPoints[3] << " " << endl;
						cerr << "| parentPoints : ";
						for(auto& i : parentPoints){
							cerr << i << " ";
						}
						cerr << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
							// cout << " 1 " << addCenterPoint << " " << edgeCenterPoints[0] << endl;
							// cout << " 1 " << addCenterPoint << " " << edgeCenterPoints[1] << endl;
							// cout << " 1 " << addCenterPoint << " " << edgeCenterPoints[2] << endl;
							// cout << " 1 " << addCenterPoint << " " << edgeCenterPoints[3] << endl;
					
					int order = sortConnFace[0] < sortConnFace[1] ? sortConnFace[1] : sortConnFace[0];
					if(
					(sortConnFace[0]==0 && sortConnFace[1]==mesh.faces[i].points.size()-1) ||
					(sortConnFace[0]==mesh.faces[i].points.size()-1 && sortConnFace[1]==0) ){
						mesh.faces[i].points.push_back(addCenterPoint);
					}
					else{
						// vector<int>::iterator iter;
						// iter = mesh.faces[i].points.begin();
						// cout << order << " " << mesh.faces[i].points.size() << " " << addCenterPoint << endl;
						// mesh.faces[i].points.insert(iter+order, addCenterPoint);
						
						vector<int> tmp_points;
						int number = 0;
						for(auto& point : mesh.faces[i].points){
							if(number == order) tmp_points.push_back(addCenterPoint);
							tmp_points.push_back(point);
							++number;
						}
						mesh.faces[i].points.clear();
						for(auto& point : tmp_points){
							mesh.faces[i].points.push_back(point);
						}
					}
					
					// add cell's points
					if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
						if ( std::find( mesh.cells[ mesh.faces[i].owner ].points.begin(), 
						                mesh.cells[ mesh.faces[i].owner ].points.end(), addCenterPoint ) 
						      == mesh.cells[ mesh.faces[i].owner ].points.end() ) {
								  
							mesh.cells[ mesh.faces[i].owner ].points.push_back(addCenterPoint);
							
						}
						if ( std::find( mesh.cells[ mesh.faces[i].neighbour ].points.begin(), 
						                mesh.cells[ mesh.faces[i].neighbour ].points.end(), addCenterPoint ) 
						      == mesh.cells[ mesh.faces[i].neighbour ].points.end() ) {
								  
							mesh.cells[ mesh.faces[i].neighbour ].points.push_back(addCenterPoint);
							
						}
					}
					else{
						if ( std::find( mesh.cells[ mesh.faces[i].owner ].points.begin(), 
						mesh.cells[ mesh.faces[i].owner ].points.end(), addCenterPoint ) 
						== mesh.cells[ mesh.faces[i].owner ].points.end() ) {
							mesh.cells[ mesh.faces[i].owner ].points.push_back(addCenterPoint);
						}
					}
				// cout << "BB" << endl;
				// cout << endl;
				// cout << addCenterPoint << " " << faceCenterPoint << " " << mesh.cells[ mesh.faces[i].neighbour ].points.size() << endl;
					
					
				}
			}
			
			// cout << "CONPLE8" << endl;
			
		}
		
		++num;
		
	}
	
	
	
	//====================================================
	// face re-ordering : internal face -> 
	//                    physical boundary face (0 1 ... nbcs-1) -> 
	//                    processor face (0 1 2.... size-1)
	
	cout << " -> completed" << endl;
	
	cout << "| Face Re-ordering ...";
	
	vector<SEMO_Face> temp_Faces;
	num=0;
	int num_n = orgNFaces;
	for(int i=0; i<orgNFaces; ++i){
		SEMO_Face temp_Faces_serial;
		
		for(int j=0; j<mesh.faces[i].points.size(); ++j){
			temp_Faces_serial.points.push_back(mesh.faces[i].points[j]);
			// cout << mesh.faces[i].points[j] << endl;
		}
			// cout << endl;
		temp_Faces_serial.owner = mesh.faces[i].owner;
		temp_Faces_serial.neighbour = mesh.faces[i].neighbour;
		temp_Faces_serial.setType(mesh.faces[i].getType());
		
		temp_Faces.push_back(temp_Faces_serial);
		
		if(idExecutedFaceAMR.size() > 0){
			if( i == idExecutedFaceAMR[num] ){
				
				for(int ic=num_n; ic<num_n+3; ++ic){
					
					SEMO_Face temp_Faces_serial2;
					
					for(int j=0; j<mesh.faces[ic].points.size(); ++j){
						temp_Faces_serial2.points.push_back(mesh.faces[ic].points[j]);
						// cout << mesh.faces[ic].points[j] << endl;
					}
			// cout << endl;
					temp_Faces_serial2.owner = mesh.faces[ic].owner;
					temp_Faces_serial2.neighbour = mesh.faces[ic].neighbour;
					temp_Faces_serial2.setType(mesh.faces[ic].getType());
					
					temp_Faces.push_back(temp_Faces_serial2);
					
				}
				
				num_n += 3;
				++num;
			}
		}
	}
	
	
	
	mesh.faces.clear();
	
	num=0;
	for(auto& face : temp_Faces){
		mesh.addFace();
		for(auto& point : face.points){
			mesh.faces.back().points.push_back(point);
			// cout << num << " " << point << endl;
		}
		mesh.faces.back().owner = face.owner;
		mesh.faces.back().neighbour = face.neighbour;
		mesh.faces.back().setType(face.getType());
		++num;
	}
	
	temp_Faces.clear();
	
	cout << " -> completed" << endl;
	
	
	// cout << idNewFaceBC[1].size() << endl;

	// num=0;
	// for(auto& face : mesh.faces){
		// for(auto& point : face.points){
			// cout << num << " " << point << endl;
		// }
		// ++num;
	// }
	
	
	//====================================================
	// cell AMR
	
	
	int orgNCells = mesh.cells.size();
	for(int icell=0; icell<orgNCells; ++icell){
		mesh.cells[icell].faces.clear();
		mesh.cells[icell].points.clear();
	}
	
	num=0;
	int orgNFaces2 = mesh.faces.size();
	for(int iface=0; iface<orgNFaces2; ++iface){
		if(mesh.faces[iface].getType() == SEMO_Types::INTERNAL_FACE){
			int owner = mesh.faces[iface].owner;
			int neighbour = mesh.faces[iface].neighbour;
			mesh.cells[owner].faces.push_back(iface);
			mesh.cells[neighbour].faces.push_back(iface);
			// cout << iface << endl;
		}
		else if(
		mesh.faces[iface].getType() == SEMO_Types::BOUNDARY_FACE ||
		mesh.faces[iface].getType() == SEMO_Types::PROCESSOR_FACE ){
			int owner = mesh.faces[iface].owner;
			mesh.cells[owner].faces.push_back(iface);
			// cout << iface << endl;
			
		}
		
		// if( iface > 10000000 ) cout << "AAAAAAAAAA " << endl;
		
	}
	
	

	// cell connection (cell's points)
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() ) {
					mesh.cells[face.owner].points.push_back(point);
				}
				vector<int> point_vec2 = mesh.cells[face.neighbour].points; 
				if ( std::find( point_vec2.begin(), point_vec2.end(), point ) == point_vec2.end() ) {
					mesh.cells[face.neighbour].points.push_back(point);
				}
			}
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() ) 
					mesh.cells[face.owner].points.push_back(point);
			}
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() )
					mesh.cells[face.owner].points.push_back(point);
			}
		}
	}
	
	
	
	cout << "| Cell AMR start ... ";
	
	vector<int> accumulateCells(orgNCells,0);
	
	
	orgNFaces2 = mesh.faces.size();
	
	int executedNCellAMR=0;
	int addNumCells = 0;
	num=0;
	for(int icell=0; icell<orgNCells; ++icell){
		// continue
		// cout << mesh.cells[icell].points.size() << endl;
		// cout << mesh.cells[icell].faces.size() << endl;
		
		bool canExecutedHexaAMR = true;
		
		if( executeAMR[icell] == false ) canExecutedHexaAMR = false;
		if( mesh.cells[icell].faces.size() != 24 ) canExecutedHexaAMR = false;
		if( mesh.cells[icell].points.size() != 26 ) canExecutedHexaAMR = false;
		for(auto& i : mesh.cells[icell].faces){
			if(mesh.faces[i].points.size() != 4) canExecutedHexaAMR = false;
		}
		
		if( canExecutedHexaAMR ){
			// cout << "AAAAAAAAAA" << endl;
			
			++executedNCellAMR;
			
			vector<int> sortedCellPoints;
			
			// vector<int> points;
			// for(auto& i : mesh.cells[icell].points){
				// points.push_back(i);
			// }
			
			// vector<int> ordFace;
			// int arr1[] = {0,6,7,8,1,9,10,11,2,12,13,14,3,15,16,17,4,18,19,20,5,21,22,23};
			// ordFace.insert(ordFace.begin(), arr1, arr1+24);
			// vector<int> ordPoint;
			// int arr2[] = {0,1,1,2,2,3,3,0,1};
			// ordPoint.insert(ordPoint.begin(), arr2, arr2+24);
			
			// vector<int> idPoints(54,0);
			// vector<int> strPoints(7,0);
			
			// for(int i=0; i<7; ++i) strPoints[i] = i*9;
			
			int iface;
			
			
			// outer points
			vector<vector<int>> idPointsOuter(6,vector<int>(8,0));
			
			// face center points
			vector<int> idPointsCenter(6,0);
			
			for(int i=0; i<6; ++i){
				iface = mesh.cells[icell].faces[i*4];
				idPointsOuter[i][0] = mesh.faces[iface].points[0];
				idPointsOuter[i][1] = mesh.faces[iface].points[1];
				iface = mesh.cells[icell].faces[i*4+1];
				idPointsOuter[i][2] = mesh.faces[iface].points[1];
				idPointsOuter[i][3] = mesh.faces[iface].points[2];
				iface = mesh.cells[icell].faces[i*4+2];
				idPointsOuter[i][4] = mesh.faces[iface].points[2];
				idPointsOuter[i][5] = mesh.faces[iface].points[3];
				iface = mesh.cells[icell].faces[i*4+3];
				idPointsOuter[i][6] = mesh.faces[iface].points[3];
				idPointsOuter[i][7] = mesh.faces[iface].points[0];
				
				iface = mesh.cells[icell].faces[i*4];
				idPointsCenter[i] = mesh.faces[iface].points[2];
			}
			
			
			
			// search vertex points
			vector<vector<int>> idPointsVertex(6,vector<int>(4,0));
			for(int i=0; i<6; ++i){
				for(int j=0; j<4; ++j){
					iface = mesh.cells[icell].faces[i*4+j];
					idPointsVertex[i][j] = mesh.faces[iface].points[j];
					// cout << i << " " << mesh.faces[iface].points[j] << endl;
				}
			}
			
			// ordering faces
			vector<int> ordRightFace(2,-1);
			vector<int> ordBottomFace(2,-1);
			vector<int> ordLeftFace(2,-1);
			vector<int> ordTopFace(2,-1);
			vector<int> ordBackFace(2,-1);
			
			int rightFace=-1;
			int bottomFace=-1;
			int leftFace=-1;
			int topFace=-1;
			int backFace=-1;
			for(int i=1; i<6; ++i){
				if( rightFace == -1 ){
					int l=0;
					for(int j=0; j<4; ++j){
						if( idPointsVertex[i][j] ==  idPointsVertex[0][0] ) {
							ordRightFace[0] = j;
							++l;
						}
						if( idPointsVertex[i][j] ==  idPointsVertex[0][1] ) {
							ordRightFace[1] = j;
							++l;
						}
					}
					if(l==2) {
						rightFace = i;
						continue;
					}
				}
				if( bottomFace == -1 ){
					int l=0;
					for(int j=0; j<4; ++j){
						if( idPointsVertex[i][j] ==  idPointsVertex[0][1] ) {
							ordBottomFace[0] = j;
							++l;
						}
						if( idPointsVertex[i][j] ==  idPointsVertex[0][2] ) {
							ordBottomFace[1] = j;
							++l;
						}
					}
					if(l==2) {
						bottomFace = i;
						continue;
					}
				}
				if( leftFace == -1 ){
					int l=0;
					for(int j=0; j<4; ++j){
						if( idPointsVertex[i][j] ==  idPointsVertex[0][2] ) {
							ordLeftFace[0] = j;
							++l;
						}
						if( idPointsVertex[i][j] ==  idPointsVertex[0][3] ) {
							ordLeftFace[1] = j;
							++l;
						}
					}
					if(l==2) {
						leftFace = i;
						continue;
					}
				}
				if( topFace == -1 ){
					int l=0;
					for(int j=0; j<4; ++j){
						if( idPointsVertex[i][j] ==  idPointsVertex[0][3] ) {
							ordTopFace[0] = j;
							++l;
						}
						if( idPointsVertex[i][j] ==  idPointsVertex[0][0] ) {
							ordTopFace[1] = j;
							++l;
						}
					}
					if(l==2) {
						topFace = i;
						continue;
					}
				}
				backFace=i;
			}
			
			// back face ordering
			if(backFace != -1){
				int l=0;
				int ord0 = ordRightFace[0];
				int ord1 = ordRightFace[1];
				int sort0, sort1;
				
				if(     ord0==0 && ord1==1){
					sort0=3; sort1=2;
				}
				else if(ord0==1 && ord1==2){
					sort0=0; sort1=3;
				}
				else if(ord0==2 && ord1==3){
					sort0=1; sort1=0;
				}
				else if(ord0==3 && ord1==0){
					sort0=2; sort1=1;
				}
				else if(ord0==0 && ord1==3){
					sort0=1; sort1=2;
				}
				else if(ord0==3 && ord1==2){
					sort0=0; sort1=1;
				}
				else if(ord0==2 && ord1==1){
					sort0=3; sort1=0;
				}
				else if(ord0==1 && ord1==0){
					sort0=2; sort1=3;
				}
				
				for(int j=0; j<4; ++j){
					if( idPointsVertex[backFace][j] ==  idPointsVertex[rightFace][sort0] ) 
						ordBackFace[0] = j;
					if( idPointsVertex[backFace][j] ==  idPointsVertex[rightFace][sort1] ) 
						ordBackFace[1] = j;
				}
			}
			
			
			// sort cell faces
			vector<int> idFaces(24,0);
			
			// cell internal faces
			vector<int> idFacesInternal(12,0);
			
			// sort cell points
			vector<int> idPoints(27,0);
			iface = 0;
			idPoints[0] = idPointsOuter[0][0];
			idPoints[1] = idPointsOuter[0][1];
			idPoints[2] = idPointsOuter[0][2];
			idPoints[3] = idPointsOuter[0][3];
			idPoints[4] = idPointsOuter[0][4];
			idPoints[5] = idPointsOuter[0][5];
			idPoints[6] = idPointsOuter[0][6];
			idPoints[7] = idPointsOuter[0][7];
			idPoints[8] = idPointsCenter[0];
			
			idFaces[0] = mesh.cells[icell].faces[0];
			idFaces[1] = mesh.cells[icell].faces[1];
			idFaces[2] = mesh.cells[icell].faces[2];
			idFaces[3] = mesh.cells[icell].faces[3];

			
			int ord0, ord1, ord2, ord3;
			int ordP0,ordP1,ordP2,ordP3,ordP4,ordP5,ordP6,ordP7;
			if(ordRightFace[0]==0 && ordRightFace[1]==1){
				ord0=3; ord1=2; ord2=1; ord3=0;
				ordP0=6;ordP1=5;ordP2=4;ordP3=3;ordP4=2;ordP5=1;ordP6=0;ordP7=7;
			}
			else if(ordRightFace[0]==1 && ordRightFace[1]==2){
				ord0=0; ord1=3; ord2=2; ord3=1;
				ordP0=0;ordP1=7;ordP2=6;ordP3=5;ordP4=4;ordP5=3;ordP6=2;ordP7=1;
			}
			else if(ordRightFace[0]==2 && ordRightFace[1]==3){
				ord0=1; ord1=0; ord2=3; ord3=2;
				ordP0=2;ordP1=1;ordP2=0;ordP3=7;ordP4=6;ordP5=5;ordP6=4;ordP7=3;
			}
			else if(ordRightFace[0]==3 && ordRightFace[1]==0){
				ord0=2; ord1=1; ord2=0; ord3=3;
				ordP0=4;ordP1=3;ordP2=2;ordP3=1;ordP4=0;ordP5=7;ordP6=6;ordP7=5;
			}
			else if(ordRightFace[0]==0 && ordRightFace[1]==3){
				ord0=1; ord1=2; ord2=3; ord3=0;
				ordP0=2;ordP1=3;ordP2=4;ordP3=5;ordP4=6;ordP5=7;ordP6=0;ordP7=1;
			}
			else if(ordRightFace[0]==3 && ordRightFace[1]==2){
				ord0=0; ord1=1; ord2=2; ord3=3;
				ordP0=0;ordP1=1;ordP2=2;ordP3=3;ordP4=4;ordP5=5;ordP6=6;ordP7=7;
			}
			else if(ordRightFace[0]==2 && ordRightFace[1]==1){
				ord0=3; ord1=0; ord2=1; ord3=2;
				ordP0=6;ordP1=7;ordP2=0;ordP3=1;ordP4=2;ordP5=3;ordP6=4;ordP7=5;
			}
			else if(ordRightFace[0]==1 && ordRightFace[1]==0){
				ord0=2; ord1=3; ord2=0; ord3=1;
				ordP0=4;ordP1=5;ordP2=6;ordP3=7;ordP4=0;ordP5=1;ordP6=2;ordP7=3;
			}
			
			idPoints[9] = idPointsOuter[rightFace][ordP7];
			idPoints[10] = idPointsCenter[rightFace];
			idPoints[11] = idPointsOuter[rightFace][ordP3];
			
			idPoints[12] = idPointsCenter[bottomFace];
			
			idFaces[4] = mesh.cells[icell].faces[rightFace*4+ord0];
			idFaces[5] = mesh.cells[icell].faces[rightFace*4+ord1];
			idFaces[6] = mesh.cells[icell].faces[rightFace*4+ord2];
			idFaces[7] = mesh.cells[icell].faces[rightFace*4+ord3];
			
			

			if(ordBottomFace[0]==0 && ordBottomFace[1]==1){
				ord0=2; ord1=3; ord2=0; ord3=1;
			}
			else if(ordBottomFace[0]==1 && ordBottomFace[1]==2){
				ord0=3; ord1=0; ord2=1; ord3=2;
			}
			else if(ordBottomFace[0]==2 && ordBottomFace[1]==3){
				ord0=0; ord1=1; ord2=2; ord3=3;
			}
			else if(ordBottomFace[0]==3 && ordBottomFace[1]==0){
				ord0=1; ord1=2; ord2=3; ord3=0;
			}
			else if(ordBottomFace[0]==0 && ordBottomFace[1]==3){
				ord0=2; ord1=1; ord2=0; ord3=3;
			}
			else if(ordBottomFace[0]==3 && ordBottomFace[1]==2){
				ord0=1; ord1=0; ord2=3; ord3=2;
			}
			else if(ordBottomFace[0]==2 && ordBottomFace[1]==1){
				ord0=0; ord1=3; ord2=2; ord3=1;
			}
			else if(ordBottomFace[0]==1 && ordBottomFace[1]==0){
				ord0=3; ord1=2; ord2=1; ord3=0;
			}
			
			// cout << ordBottomFace[0] << " " << ordBottomFace[1] << endl;
			
			idFaces[8] = mesh.cells[icell].faces[bottomFace*4+ord0];
			idFaces[9] = mesh.cells[icell].faces[bottomFace*4+ord1];
			idFaces[10] = mesh.cells[icell].faces[bottomFace*4+ord2];
			idFaces[11] = mesh.cells[icell].faces[bottomFace*4+ord3];

			// cout << bottomFace << " " <<  ord1 << " " << idFaces[9] << endl;


			if(ordLeftFace[0]==0 && ordLeftFace[1]==1){
				ord0=2; ord1=3; ord2=0; ord3=1;
				ordP0=4;ordP1=5;ordP2=6;ordP3=7;ordP4=0;ordP5=1;ordP6=2;ordP7=3;
			}
			else if(ordLeftFace[0]==1 && ordLeftFace[1]==2){
				ord0=3; ord1=0; ord2=1; ord3=2;
				ordP0=6;ordP1=7;ordP2=0;ordP3=1;ordP4=2;ordP5=3;ordP6=4;ordP7=5;
			}
			else if(ordLeftFace[0]==2 && ordLeftFace[1]==3){
				ord0=0; ord1=1; ord2=2; ord3=3;
				ordP0=0;ordP1=1;ordP2=2;ordP3=3;ordP4=4;ordP5=5;ordP6=6;ordP7=7;
			}
			else if(ordLeftFace[0]==3 && ordLeftFace[1]==0){
				ord0=1; ord1=2; ord2=3; ord3=0;
				ordP0=2;ordP1=3;ordP2=4;ordP3=5;ordP4=6;ordP5=7;ordP6=0;ordP7=1;
			}
			else if(ordLeftFace[0]==0 && ordLeftFace[1]==3){
				ord0=2; ord1=1; ord2=0; ord3=3;
				ordP0=4;ordP1=3;ordP2=2;ordP3=1;ordP4=0;ordP5=7;ordP6=6;ordP7=5;
			}
			else if(ordLeftFace[0]==3 && ordLeftFace[1]==2){
				ord0=1; ord1=0; ord2=3; ord3=2;
				ordP0=2;ordP1=1;ordP2=0;ordP3=7;ordP4=6;ordP5=5;ordP6=4;ordP7=3;
			}
			else if(ordLeftFace[0]==2 && ordLeftFace[1]==1){
				ord0=0; ord1=3; ord2=2; ord3=1;
				ordP0=0;ordP1=7;ordP2=6;ordP3=5;ordP4=4;ordP5=3;ordP6=2;ordP7=1;
			}
			else if(ordLeftFace[0]==1 && ordLeftFace[1]==0){
				ord0=3; ord1=2; ord2=1; ord3=0;
				ordP0=6;ordP1=5;ordP2=4;ordP3=3;ordP4=2;ordP5=1;ordP6=0;ordP7=7;
			}
			
			
			idPoints[13] = idPointsOuter[leftFace][ordP3];
			idPoints[14] = idPointsCenter[leftFace];
			idPoints[15] = idPointsOuter[leftFace][ordP7];
			
			idPoints[16] = idPointsCenter[topFace];
			
			
			idFaces[12] = mesh.cells[icell].faces[leftFace*4+ord0];
			idFaces[13] = mesh.cells[icell].faces[leftFace*4+ord1];
			idFaces[14] = mesh.cells[icell].faces[leftFace*4+ord2];
			idFaces[15] = mesh.cells[icell].faces[leftFace*4+ord3];
			
			


			// if(icell==0) {
				// // for(auto& ii : mesh.faces[idFaces[7]].points){
					// // cout << ii << endl;
				// // }
					// // cout << endl;
					
				// for(int ii=13; ii<17; ++ii){
					// cout << idPoints[ii] << endl;
				// }
					// cout << endl;
			// }




			if(ordTopFace[0]==0 && ordTopFace[1]==1){
				ord0=3; ord1=2; ord2=1; ord3=0;
			}
			else if(ordTopFace[0]==1 && ordTopFace[1]==2){
				ord0=4; ord1=3; ord2=2; ord3=1;
			}
			else if(ordTopFace[0]==2 && ordTopFace[1]==3){
				ord0=1; ord1=0; ord2=3; ord3=2;
			}
			else if(ordTopFace[0]==3 && ordTopFace[1]==0){
				ord0=2; ord1=1; ord2=0; ord3=3;
			}
			else if(ordTopFace[0]==0 && ordTopFace[1]==3){
				ord0=1; ord1=2; ord2=3; ord3=0;
			}
			else if(ordTopFace[0]==3 && ordTopFace[1]==2){
				ord0=0; ord1=1; ord2=2; ord3=3;
			}
			else if(ordTopFace[0]==2 && ordTopFace[1]==1){
				ord0=3; ord1=0; ord2=1; ord3=2;
			}
			else if(ordTopFace[0]==1 && ordTopFace[1]==0){
				ord0=2; ord1=3; ord2=0; ord3=1;
			}
			
			idFaces[16] = mesh.cells[icell].faces[topFace*4+ord0];
			idFaces[17] = mesh.cells[icell].faces[topFace*4+ord1];
			idFaces[18] = mesh.cells[icell].faces[topFace*4+ord2];
			idFaces[19] = mesh.cells[icell].faces[topFace*4+ord3];





			
			

			if(ordBackFace[0]==0 && ordBackFace[1]==1){
				ord0=0; ord1=1; ord2=2; ord3=3;
				ordP0=0;ordP1=1;ordP2=2;ordP3=3;ordP4=4;ordP5=5;ordP6=6;ordP7=7;
			}
			else if(ordBackFace[0]==1 && ordBackFace[1]==2){
				ord0=1; ord1=2; ord2=3; ord3=0;
				ordP0=2;ordP1=3;ordP2=4;ordP3=5;ordP4=6;ordP5=7;ordP6=0;ordP7=1;
			}
			else if(ordBackFace[0]==2 && ordBackFace[1]==3){
				ord0=2; ord1=3; ord2=0; ord3=1;
				ordP0=4;ordP1=5;ordP2=6;ordP3=7;ordP4=0;ordP5=1;ordP6=2;ordP7=3;
			}
			else if(ordBackFace[0]==3 && ordBackFace[1]==0){
				ord0=3; ord1=0; ord2=1; ord3=2;
				ordP0=6;ordP1=7;ordP2=0;ordP3=1;ordP4=2;ordP5=3;ordP6=4;ordP7=5;
			}
			else if(ordBackFace[0]==0 && ordBackFace[1]==3){
				ord0=0; ord1=3; ord2=2; ord3=1;
				ordP0=0;ordP1=7;ordP2=6;ordP3=5;ordP4=4;ordP5=3;ordP6=2;ordP7=1;
			}
			else if(ordBackFace[0]==3 && ordBackFace[1]==2){
				ord0=3; ord1=2; ord2=1; ord3=0;
				ordP0=6;ordP1=5;ordP2=4;ordP3=3;ordP4=2;ordP5=1;ordP6=0;ordP7=7;
			}
			else if(ordBackFace[0]==2 && ordBackFace[1]==1){
				ord0=2; ord1=1; ord2=0; ord3=3;
				ordP0=4;ordP1=3;ordP2=2;ordP3=1;ordP4=0;ordP5=7;ordP6=6;ordP7=5;
			}
			else if(ordBackFace[0]==1 && ordBackFace[1]==0){
				ord0=1; ord1=0; ord2=3; ord3=2;
				ordP0=2;ordP1=1;ordP2=0;ordP3=7;ordP4=6;ordP5=5;ordP6=4;ordP7=3;
			}
			
			
			idPoints[18] = idPointsOuter[backFace][ordP0];
			idPoints[19] = idPointsOuter[backFace][ordP1];
			idPoints[20] = idPointsOuter[backFace][ordP2];
			idPoints[21] = idPointsOuter[backFace][ordP3];
			idPoints[22] = idPointsOuter[backFace][ordP4];
			idPoints[23] = idPointsOuter[backFace][ordP5];
			idPoints[24] = idPointsOuter[backFace][ordP6];
			idPoints[25] = idPointsOuter[backFace][ordP7];
			idPoints[26] = idPointsCenter[backFace];
			
			
			idFaces[20] = mesh.cells[icell].faces[backFace*4+ord0];
			idFaces[21] = mesh.cells[icell].faces[backFace*4+ord1];
			idFaces[22] = mesh.cells[icell].faces[backFace*4+ord2];
			idFaces[23] = mesh.cells[icell].faces[backFace*4+ord3];
			

			// cell center point
			int strPoint = mesh.points.size();
			mesh.addPoint();
			++newNPoints;
			mesh.points.back().x = 0.0;
			mesh.points.back().y = 0.0;
			mesh.points.back().z = 0.0;
			for(int i=0; i<4; ++i){
				mesh.points.back().x += 0.125*mesh.points[idPointsVertex[0][i]].x;
				mesh.points.back().y += 0.125*mesh.points[idPointsVertex[0][i]].y;
				mesh.points.back().z += 0.125*mesh.points[idPointsVertex[0][i]].z;
			}
			for(int i=0; i<4; ++i){
				mesh.points.back().x += 0.125*mesh.points[idPointsVertex[backFace][i]].x;
				mesh.points.back().y += 0.125*mesh.points[idPointsVertex[backFace][i]].y;
				mesh.points.back().z += 0.125*mesh.points[idPointsVertex[backFace][i]].z;
			}
			
			idPoints[17] = strPoint;
			
			
			
			
			
			// cell internal faces
			int strFace = mesh.faces.size();
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[9]);
			mesh.faces.back().points.push_back(idPoints[10]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[16]);
			mesh.faces.back().owner = icell + addNumCells + 0;
			mesh.faces.back().neighbour = icell + addNumCells + 4;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[10]);
			mesh.faces.back().points.push_back(idPoints[11]);
			mesh.faces.back().points.push_back(idPoints[12]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().owner = icell + addNumCells + 1;
			mesh.faces.back().neighbour = icell + addNumCells + 5;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[12]);
			mesh.faces.back().points.push_back(idPoints[13]);
			mesh.faces.back().points.push_back(idPoints[14]);
			mesh.faces.back().owner = icell + addNumCells + 2;
			mesh.faces.back().neighbour = icell + addNumCells + 6;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[16]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[14]);
			mesh.faces.back().points.push_back(idPoints[15]);
			mesh.faces.back().owner = icell + addNumCells + 3;
			mesh.faces.back().neighbour = icell + addNumCells + 7;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[25]);
			mesh.faces.back().points.push_back(idPoints[26]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[16]);
			mesh.faces.back().owner = icell + addNumCells + 4;
			mesh.faces.back().neighbour = icell + addNumCells + 7;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[26]);
			mesh.faces.back().points.push_back(idPoints[21]);
			mesh.faces.back().points.push_back(idPoints[12]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().owner = icell + addNumCells + 5;
			mesh.faces.back().neighbour = icell + addNumCells + 6;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[12]);
			mesh.faces.back().points.push_back(idPoints[3]);
			mesh.faces.back().points.push_back(idPoints[8]);
			mesh.faces.back().owner = icell + addNumCells + 1;
			mesh.faces.back().neighbour = icell + addNumCells + 2;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[16]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[8]);
			mesh.faces.back().points.push_back(idPoints[7]);
			mesh.faces.back().owner = icell + addNumCells + 0;
			mesh.faces.back().neighbour = icell + addNumCells + 3;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			
			
			
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[19]);
			mesh.faces.back().points.push_back(idPoints[10]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[26]);
			mesh.faces.back().owner = icell + addNumCells + 4;
			mesh.faces.back().neighbour = icell + addNumCells + 5;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[10]);
			mesh.faces.back().points.push_back(idPoints[1]);
			mesh.faces.back().points.push_back(idPoints[8]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().owner = icell + addNumCells + 0;
			mesh.faces.back().neighbour = icell + addNumCells + 1;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[8]);
			mesh.faces.back().points.push_back(idPoints[5]);
			mesh.faces.back().points.push_back(idPoints[14]);
			mesh.faces.back().owner = icell + addNumCells + 3;
			mesh.faces.back().neighbour = icell + addNumCells + 2;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			mesh.addFace();
			++newNFaces;
			mesh.faces.back().points.push_back(idPoints[26]);
			mesh.faces.back().points.push_back(idPoints[17]);
			mesh.faces.back().points.push_back(idPoints[14]);
			mesh.faces.back().points.push_back(idPoints[23]);
			mesh.faces.back().owner = icell + addNumCells + 7;
			mesh.faces.back().neighbour = icell + addNumCells + 6;
			mesh.faces.back().setType(SEMO_Types::INTERNAL_FACE);
			
			
			// outer faces
			if(mesh.faces[idFaces[0]].owner == icell){
				mesh.faces[idFaces[0]].owner = icell + addNumCells + 0;
				mesh.faces[idFaces[1]].owner = icell + addNumCells + 1;
				mesh.faces[idFaces[2]].owner = icell + addNumCells + 2;
				mesh.faces[idFaces[3]].owner = icell + addNumCells + 3;
			}
			else{
				mesh.faces[idFaces[0]].neighbour = icell + addNumCells + 0;
				mesh.faces[idFaces[1]].neighbour = icell + addNumCells + 1;
				mesh.faces[idFaces[2]].neighbour = icell + addNumCells + 2;
				mesh.faces[idFaces[3]].neighbour = icell + addNumCells + 3;
			}
			
			
			if(mesh.faces[idFaces[4]].owner == icell){
				mesh.faces[idFaces[4]].owner = icell + addNumCells + 4;
				mesh.faces[idFaces[5]].owner = icell + addNumCells + 5;
				mesh.faces[idFaces[6]].owner = icell + addNumCells + 1;
				mesh.faces[idFaces[7]].owner = icell + addNumCells + 0;
			}
			else{
				mesh.faces[idFaces[4]].neighbour = icell + addNumCells + 4;
				mesh.faces[idFaces[5]].neighbour = icell + addNumCells + 5;
				mesh.faces[idFaces[6]].neighbour = icell + addNumCells + 1;
				mesh.faces[idFaces[7]].neighbour = icell + addNumCells + 0;
			}
			
			// cout << mesh.faces.size() << " " << mesh.cells[icell].faces.size() << " " << bottomFace*4+ord1 << " " << idFaces[8] << " " << idFaces[9] << " " << idFaces[10] << " " << idFaces[11] << endl;
			
			if(mesh.faces[idFaces[8]].owner == icell){
				mesh.faces[idFaces[8]].owner = icell + addNumCells + 6;
				mesh.faces[idFaces[9]].owner = icell + addNumCells + 5;
				mesh.faces[idFaces[10]].owner = icell + addNumCells + 1;
				mesh.faces[idFaces[11]].owner = icell + addNumCells + 2;
			}
			else{
				mesh.faces[idFaces[8]].neighbour = icell + addNumCells + 6;
				mesh.faces[idFaces[9]].neighbour = icell + addNumCells + 5;
				mesh.faces[idFaces[10]].neighbour = icell + addNumCells + 1;
				mesh.faces[idFaces[11]].neighbour = icell + addNumCells + 2;
			}
			
			// cout << mesh.faces[idFaces[10]].points[0] << " " << mesh.faces[idFaces[10]].points[1] << " " << mesh.faces[idFaces[10]].points[2] << " " << mesh.faces[idFaces[10]].points[3] << endl;
			
			if(mesh.faces[idFaces[12]].owner == icell){
				mesh.faces[idFaces[12]].owner = icell + addNumCells + 7;
				mesh.faces[idFaces[13]].owner = icell + addNumCells + 6;
				mesh.faces[idFaces[14]].owner = icell + addNumCells + 2;
				mesh.faces[idFaces[15]].owner = icell + addNumCells + 3;
			}
			else{
				mesh.faces[idFaces[12]].neighbour = icell + addNumCells + 7;
				mesh.faces[idFaces[13]].neighbour = icell + addNumCells + 6;
				mesh.faces[idFaces[14]].neighbour = icell + addNumCells + 2;
				mesh.faces[idFaces[15]].neighbour = icell + addNumCells + 3;
			}
			
			
			if(mesh.faces[idFaces[16]].owner == icell){
				mesh.faces[idFaces[16]].owner = icell + addNumCells + 7;
				mesh.faces[idFaces[17]].owner = icell + addNumCells + 4;
				mesh.faces[idFaces[18]].owner = icell + addNumCells + 0;
				mesh.faces[idFaces[19]].owner = icell + addNumCells + 3;
			}
			else{
				mesh.faces[idFaces[16]].neighbour = icell + addNumCells + 7;
				mesh.faces[idFaces[17]].neighbour = icell + addNumCells + 4;
				mesh.faces[idFaces[18]].neighbour = icell + addNumCells + 0;
				mesh.faces[idFaces[19]].neighbour = icell + addNumCells + 3;
			}
			
			
			if(mesh.faces[idFaces[20]].owner == icell){
				mesh.faces[idFaces[20]].owner = icell + addNumCells + 4;
				mesh.faces[idFaces[21]].owner = icell + addNumCells + 5;
				mesh.faces[idFaces[22]].owner = icell + addNumCells + 6;
				mesh.faces[idFaces[23]].owner = icell + addNumCells + 7;
			}
			else{
				mesh.faces[idFaces[20]].neighbour = icell + addNumCells + 4;
				mesh.faces[idFaces[21]].neighbour = icell + addNumCells + 5;
				mesh.faces[idFaces[22]].neighbour = icell + addNumCells + 6;
				mesh.faces[idFaces[23]].neighbour = icell + addNumCells + 7;
			}
			
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			mesh.addCell();
			mesh.cells.back().level = mesh.cells[icell].level+1;
			
			++newNCells;
			++newNCells;
			++newNCells;
			++newNCells;
			++newNCells;
			++newNCells;
			++newNCells;
			
			addNumCells += 7;
			
			
			// addNumCells
			
			
			
			
		}
		else{
			
			
			// outer faces
			for(int i=0; i<mesh.cells[icell].faces.size(); ++i){
				int iface = mesh.cells[icell].faces[i];
				if(mesh.faces[iface].owner == icell){
					mesh.faces[iface].owner = icell + addNumCells;
				}
				else{
					mesh.faces[iface].neighbour = icell + addNumCells;
				}
			}
			
		}
		
		
		
		// // AMR.sortCellPoints(mesh.cell[icell].points, sortedCellPoints);
		
	}
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	cout << " -> completed" << endl;
	
	
	
	//====================================================
	// face re-ordering : internal face -> 
	//                    physical boundary face (0 1 ... nbcs-1) -> 
	//                    processor face (0 1 2.... size-1)
	num=0;
	for(int i=0; i<orgNFaces2; ++i){
		if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			SEMO_Face temp_Faces_serial;
			
			for(int j=0; j<mesh.faces[i].points.size(); ++j){
				temp_Faces_serial.points.push_back(mesh.faces[i].points[j]);
			}
			temp_Faces_serial.owner = mesh.faces[i].owner;
			temp_Faces_serial.neighbour = mesh.faces[i].neighbour;
			temp_Faces_serial.setType(mesh.faces[i].getType());
			
			temp_Faces.push_back(temp_Faces_serial);
		}	
	}
	
	for(int i=orgNFaces2; i<mesh.faces.size(); ++i){
		SEMO_Face temp_Faces_serial;
		
		for(int j=0; j<mesh.faces[i].points.size(); ++j){
			temp_Faces_serial.points.push_back(mesh.faces[i].points[j]);
		}
		temp_Faces_serial.owner = mesh.faces[i].owner;
		temp_Faces_serial.neighbour = mesh.faces[i].neighbour;
		temp_Faces_serial.setType(mesh.faces[i].getType());
		
		temp_Faces.push_back(temp_Faces_serial);
	}
	for(int i=0; i<orgNFaces2; ++i){
		if(
		mesh.faces[i].getType() == SEMO_Types::BOUNDARY_FACE ||
		mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE
		){
			SEMO_Face temp_Faces_serial;
			
			for(int j=0; j<mesh.faces[i].points.size(); ++j){
				temp_Faces_serial.points.push_back(mesh.faces[i].points[j]);
			}
			temp_Faces_serial.owner = mesh.faces[i].owner;
			temp_Faces_serial.neighbour = mesh.faces[i].neighbour;
			temp_Faces_serial.setType(mesh.faces[i].getType());
			
			temp_Faces.push_back(temp_Faces_serial);
		}	
	}
	
	mesh.faces.clear();
	
	num=0;
	for(auto& face : temp_Faces){
		mesh.addFace();
		for(auto& point : face.points){
			mesh.faces.back().points.push_back(point);
		}
		mesh.faces.back().owner = face.owner;
		mesh.faces.back().neighbour = face.neighbour;
		mesh.faces.back().setType(face.getType());
		++num;
	}
	
	temp_Faces.clear();
	
	
	
	
	// boundary & processor faces
	int faceNum=0;
	int nnFaceBC = 0;
	int nnFacePR = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			++faceNum;
		}
		else if(mesh.faces[i].getType() == SEMO_Types::BOUNDARY_FACE){
			++nnFaceBC;
		}
		else if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			++nnFacePR;
		}
	}
	
	// cout << faceNum << " " << nnFaceBC << " " << nbcs << endl;
	
	int nnFaceBC2 = 0;
	
	for(int ibc=0; ibc<nbcs; ++ibc){
		mesh.boundary[ibc].startFace = faceNum;
		mesh.boundary[ibc].nFaces = idNewFaceBC[ibc].size();
		faceNum += mesh.boundary[ibc].nFaces;
		nnFaceBC2 += mesh.boundary[ibc].nFaces;
	}
	
	if( nnFaceBC2 != nnFaceBC ){
		cerr << "#Error : nnFaceBC2 != nnFaceBC" << endl;
		cerr << nnFaceBC2 << " " << nnFaceBC << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// cout << mesh.boundary[0].nFaces << " " <<  mesh.faces.size() << endl;
	
	
	if(size>1){
		for(int ip=0; ip<size; ++ip){
			int ibc = nbcs+ip;
			mesh.boundary[ibc].startFace = faceNum;
			mesh.boundary[ibc].nFaces = idNewFaceProc[ip].size();
			faceNum += mesh.boundary[ibc].nFaces;
		}
	}
	
	
	
	//====================================================
	
	// list setup
	// list points
	mesh.listPoints.clear();
	for(auto& point : mesh.points){
		mesh.listPoints.push_back(&point);
	}
	
	// list faces
	mesh.listFaces.clear();
	for(auto& face : mesh.faces){
		mesh.listFaces.push_back(&face);
	}
	
	
	// list Cells
	mesh.listCells.clear();
	for(auto& cell : mesh.cells){
		mesh.listCells.push_back(&cell);
	}
	
	
	
	// cell connection (cell's face)
	// delete faces
	for(auto& cell : mesh.cells){
		// cout << "aa" << endl;
		cell.faces.clear();
		cell.points.clear();
	}
	
	
	// cell face connection
	faceNum=0;
	int faceNumIn = 0;
	int faceNumBC = 0;
	int faceNumProc = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
			mesh.cells[face.neighbour].faces.push_back(faceNum);
			++faceNumIn;
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
			++faceNumBC;
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			mesh.cells[face.owner].faces.push_back(faceNum);
			++faceNumProc;
		}
		++faceNum;
	}
	
	
	// cell point connection (cell's points)
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// cout << "INTERNAL_FACE"<< endl;
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() ) {
					mesh.cells[face.owner].points.push_back(point);
				}
				// else{
					// int l=0;
					// for(auto& i : mesh.cells[face.owner].points){
						// if(i == point) ++l;
					// }
					// cout << l << endl;
				// }
				vector<int> point_vec2 = mesh.cells[face.neighbour].points; 
				if ( std::find( point_vec2.begin(), point_vec2.end(), point ) == point_vec2.end() ) {
					mesh.cells[face.neighbour].points.push_back(point);
				}
				// else{
					// int l=0;
					// for(auto& i : mesh.cells[face.neighbour].points){
						// if(i == point) ++l;
					// }
					// cout << l << endl;
				// }
			}
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// cout << "BOUNDARY_FACE"<< endl;
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() ) {
					int l=0;
					// for(auto& i : mesh.cells[face.owner].points){
						// if(i == point) ++l;
					// }
					// cout << l << endl;
					mesh.cells[face.owner].points.push_back(point);
				}
				// else{
					// int l=0;
					// for(auto& i : mesh.cells[face.owner].points){
						// if(i == point) ++l;
					// }
					// cout << l << endl;
				// }
			}
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// cout << "PROCESSOR_FACE"<< endl;
			for(auto& point : face.points){
				vector<int> point_vec1 = mesh.cells[face.owner].points; 
				if ( std::find( point_vec1.begin(), point_vec1.end(), point ) == point_vec1.end() ) {
					mesh.cells[face.owner].points.push_back(point);
				}
			}
		}
	}
	// for(auto& face : mesh.faces){
		// cout << face.points.size() << endl;
		// for(auto& i : face.points){
			// cout << i << " ";
		// }
		// cout << endl;
	// }
	
	// for(auto& cell : mesh.cells){
		// cout << cell.faces.size() << endl;
		// cout << cell.points.size() << endl;
		// // for(auto& i : cell.points){
			// // cout << i << " ";
		// // }
		// // cout << endl;
	// }
		// mesh.saveFile("vtu");
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	cout << "┌-----------------------------------" << endl;
	cout << "| AMR : # of faces : " << orgNFaces << " -> " << mesh.faces.size() << endl;
	cout << "| AMR : # of cells : " << orgNCells << " -> " << mesh.cells.size() << endl;
	cout << "| AMR : # of created points : " << newNPoints << endl;
	cout << "| AMR : # of executed + created faces : " << executedNFaceAMR << " + " << newNFaces << endl;
	cout << "| AMR : # of executed + created cells : " << executedNCellAMR << " + " << newNCells << endl;
	cout << "└-----------------------------------" << endl;
	
	

	
	
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



void SEMO_Hexa_Oct_AMR::searchParallelPointsPolygon( 
	SEMO_Mesh_Builder& mesh, vector<int>& points, vector<vector<int>>& parPoints){
	
	int psize = points.size();
	
	vector<double> beforeNvec(3,0.0);
	vector<int> temp_parPoints;
	for(int i=0; i<psize; ++i){
		int p0=i;
		int p1=i+1;
		if(i==psize-1) p1=0;
		
		vector<double> nvec(3,0.0);
		
		nvec[0] = mesh.points[points[p1]].x - mesh.points[points[p0]].x;
		nvec[1] = mesh.points[points[p1]].y - mesh.points[points[p0]].y;
		nvec[2] = mesh.points[points[p1]].z - mesh.points[points[p0]].z;
		
		double nomalSize = sqrt( nvec[0]*nvec[0] + nvec[1]*nvec[1] + nvec[2]*nvec[2] );
		nvec[0] /= nomalSize;
		nvec[1] /= nomalSize;
		nvec[2] /= nomalSize;
		
		// cout << psize << endl;
		// double a1 = abs(nvec[0]-beforeNvec[0]);
		// double a2 = abs(nvec[1]-beforeNvec[1]);
		// double a3 = abs(nvec[2]-beforeNvec[2]);
		// cout << nvec[0] << " " << beforeNvec[0] << endl;
		// cout << nvec[1] << " " << beforeNvec[1] << endl;
		// cout << nvec[2] << " " << beforeNvec[2] << endl;
		
		if(i==0){
			temp_parPoints.push_back(points[p0]);
			temp_parPoints.push_back(points[p1]);
		}
		else{
			if(
			(nvec[0]-beforeNvec[0]) < 0.001 && (nvec[0]-beforeNvec[0]) > -0.001 &&
			(nvec[1]-beforeNvec[1]) < 0.001 && (nvec[1]-beforeNvec[1]) > -0.001 &&
			(nvec[2]-beforeNvec[2]) < 0.001 && (nvec[2]-beforeNvec[2]) > -0.001 ){
				temp_parPoints.push_back(points[p1]);
			}
			else {
				parPoints.push_back(temp_parPoints);
				temp_parPoints.clear();
				temp_parPoints.push_back(points[p0]);
				temp_parPoints.push_back(points[p1]);
			}
		}
		
		beforeNvec[0] = nvec[0];
		beforeNvec[1] = nvec[1];
		beforeNvec[2] = nvec[2];
		
		if(i==psize-1) parPoints.push_back(temp_parPoints);
		
	}
	
	// cout << "a : " << psize << endl;
	// for(auto& i : parPoints){
		// for(auto& j : i){
			// cout << "a : " << j << endl;
		// }
		// cout << endl;
	// }
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}








// void SEMO_Hexa_Oct_AMR::sortCellPoints( 
	// vector<int>& points, vector<int>& sortedPoints){
		
		
// }