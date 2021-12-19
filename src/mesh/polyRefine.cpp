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

class faces_refind{

private:

public:
	int id;
	vector<int> points;
	int owner;
	int neighbour;
	bool boolOwner = false;
	bool boolNeighbour = false;
};



class groupMesh_Refine{

private:

public:
	int parent;
	vector<faces_refind> faces;
	vector<int> vertexPoints;
	vector<int> vertexCenterPoints;
	int centerPoint;
	int type;
	vector<vector<int>> subOutEdgePoints;
	vector<vector<int>> subIntEdgePoints;
};





void extractVertexPoints(
	int faceLevel,
	vector<int>& facePoints,
	vector<int>& facePointLevels,
	vector<int>& vertexPoints
) {
	vertexPoints.clear();

	int facePointsSize = (int)facePoints.size();
	for (int i = 0; i < facePointsSize; ++i) {
		if(facePointLevels[i] <= faceLevel){
			vertexPoints.push_back(facePoints[i]);
		}
	}
};


void extractVertexCenterPoints(
	int faceLevel,
	vector<int>& facePoints,
	vector<int>& facePointLevels,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints,
	vector<bool>& addVertexCeterPoints
) {
	vertexCenterPoints.clear();
	addVertexCeterPoints.clear();
	// edgeOrders.clear();

	int vertexPointsSize = (int)vertexPoints.size();
	int facePointsSize = (int)facePoints.size();
	
	// vector<int> facePoints_copy = facePoints;

	int num = 0;
	// auto it = facePoints_copy.begin();
	for (int i = 0; i < facePointsSize; ++i) {
		
		// cout << facePoints[i] << " " << facePointLevels[i] << endl;
		
		
		// cout << i << " " << facePoints.size() << " " << num << " " << vertexPoints.size() << endl;
		// ++it;
		if (facePoints[i] == vertexPoints[num]) {
			int iNext = (i + 1 < facePointsSize) ? i + 1 : 0;
			int numNext = (num + 1 < vertexPointsSize) ? num + 1 : 0;
			
			// cout << iNext << " " << facePoints.size() << " " << numNext << " " << vertexPoints.size() << endl;
			
			
			if (facePoints[iNext] == vertexPoints[numNext]) {
				// // create vertexCenterPoints
				// int newPoint = (int)facePoints_copy.size();
				// it = ++facePoints_copy.insert(it, newPoint);
				int newPoint = -1;
				vertexCenterPoints.push_back(newPoint);
				addVertexCeterPoints.push_back(true);
				// edgeOrders.push_back(i);
			}
			++num;
			if(num >= vertexPointsSize-1) num = vertexPointsSize-1;
		}
		else {
			if (facePointLevels[i] == faceLevel + 1) {
				vertexCenterPoints.push_back(facePoints[i]);
				addVertexCeterPoints.push_back(false);
				// edgeOrders.push_back(i);
			}
		}
	}

	// facePoints = facePoints_copy;
};



void extractSubEdgePoints(
	vector<int>& facePoints,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints,
	vector<vector<int>>& subEdgePoints
) {

	int vertexCeterPointsSize = (int)vertexCenterPoints.size();
	int subEdgeSize = vertexCeterPointsSize * 2;
	subEdgePoints.resize(subEdgeSize,vector<int>());

	int i = 0;
	auto iter0 = find(facePoints.begin(), facePoints.end(), vertexPoints[0]);
	for (i = 0; i < vertexCeterPointsSize-1; ++i) {
		auto iter1 = find(iter0, facePoints.end(), vertexCenterPoints[i]);
		auto iter2 = find(iter1, facePoints.end(), vertexPoints[i+1]);

		int dist = (int)distance(iter0, iter1);
		subEdgePoints[i * 2].insert(subEdgePoints[i * 2].begin(), iter0, iter1);
		//cout << dist << endl;

		dist = (int)distance(iter1, iter2);
		subEdgePoints[i * 2+1].insert(subEdgePoints[i * 2+1].begin(), iter1, iter2);
		//cout << dist << endl;

		iter0 = iter2;
	}
	//cout << "end" << endl;
	auto iter1 = find(iter0, facePoints.end(), vertexCenterPoints[i]);
	auto iter2 = facePoints.end();
	int dist = (int)distance(iter0, iter1);
	subEdgePoints[i * 2].insert(subEdgePoints[i * 2].begin(), iter0, iter1);
	//cout << dist << endl;

	dist = (int)distance(iter1, iter2);
	subEdgePoints[i * 2+1].insert(subEdgePoints[i * 2+1].begin(), iter1, iter2);
	//cout << dist << endl;

};



void addVertexCenterPoint(
	SEMO_Mesh_Builder& mesh, 
	int faceLevel,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints
) {

	int verNum = vertexCenterPoints.size();
	for(int j=0; j<verNum; ++j){
		// if(addVertexCeterPoints[j]==true){
		if(vertexCenterPoints[j] == -1){
			vertexCenterPoints[j] = mesh.points.size();
			mesh.addPoint();
			int pointNum0 = j;
			int pointNum1 = (j+1<verNum) ? j+1 : 0;
			int point0 = vertexPoints[pointNum0];
			int point1 = vertexPoints[pointNum1];
			
			if( point0 >= mesh.points.size() ) cout << "| WARNING : point0 >= point.size" << endl;
			if( point1 >= mesh.points.size() ) cout << "| WARNING : point0 >= point.size" << endl;
			// cout << mesh.points.size() << endl;
			// cout << pointNum0 << " " << pointNum1 << " " << point0 << " " << point1 << endl;
			
			mesh.points.back().x = 0.5*(mesh.points[point0].x + mesh.points[point1].x);
			mesh.points.back().y = 0.5*(mesh.points[point0].y + mesh.points[point1].y);
			mesh.points.back().z = 0.5*(mesh.points[point0].z + mesh.points[point1].z);
			mesh.points.back().level = faceLevel+1;
		}
	}
				

};



void addFaceCenterPoint(
	SEMO_Mesh_Builder& mesh, 
	int faceLevel,
	vector<int>& vertexPoints,
	int& centerPoint
) {

	centerPoint = mesh.points.size();
	mesh.addPoint();
	mesh.points.back().x = 0.0;
	mesh.points.back().y = 0.0;
	mesh.points.back().z = 0.0;
	int verNum = vertexPoints.size();
	double dVerNum = 1.0/(double)verNum;
	for(auto& j : vertexPoints){
		mesh.points.back().x += mesh.points[j].x*dVerNum;
		mesh.points.back().y += mesh.points[j].y*dVerNum;
		mesh.points.back().z += mesh.points[j].z*dVerNum;
	}
	mesh.points.back().level = faceLevel+1;
	
};


void addEdgeFacesPointsToVertexCenterPoint(
	SEMO_Mesh_Builder& mesh,
	int i,
	vector<vector<int>>& facesEdges,
	vector<vector<int>>& edgesFaces,
	vector<bool>& boolEdgeRefine,
	vector<bool>& addVertexCeterPoints,
	vector<int>& canRefineEdgeOrders,
	vector<int>& vertexPoints,
	vector<int>& vertexCenterPoints
) {
	int vertexSize = vertexPoints.size();
	int tmpNewNum = 0;
	for(int j=0; j<vertexSize; ++j){
		if(addVertexCeterPoints[j]==true){
			// int ordEdge = canRefineEdgeOrders[tmpNewNum];
			// ++tmpNewNum;
			// cout << ordEdge << " " << j << endl;
			// int iEdge = facesEdges[i][ordEdge];
			
			int iEdge = canRefineEdgeOrders[tmpNewNum];
			++tmpNewNum;
			
			// if(boolEdgeRefine[iEdge]==false){
				boolEdgeRefine[iEdge] = true;
				int isertP = vertexCenterPoints[j];
				for(auto& k : edgesFaces[iEdge]){
					
					int dist = 0;
					for(int l=0; l<facesEdges[k].size(); ++l){
						int iEdge2 = facesEdges[k][l];
						++dist;
						if(iEdge2 == iEdge) break;
						if(boolEdgeRefine[iEdge2]==true){
							++dist;
						}
					}
					mesh.faces[k].points.insert(
						mesh.faces[k].points.begin()+dist,
						isertP);
					
				}
			// }
		}
	}
			
};



void addCellCenterPoint(
	SEMO_Mesh_Builder& mesh, 
	int cellLevel,
	vector<int>& vertexPoints,
	int& centerPoint
) {

	centerPoint = mesh.points.size();
	mesh.addPoint();
	mesh.points.back().x = 0.0;
	mesh.points.back().y = 0.0;
	mesh.points.back().z = 0.0;
	int verNum = vertexPoints.size();
	double dVerNum = 1.0/(double)verNum;
	for(auto& j : vertexPoints){
		mesh.points.back().x += mesh.points[j].x*dVerNum;
		mesh.points.back().y += mesh.points[j].y*dVerNum;
		mesh.points.back().z += mesh.points[j].z*dVerNum;
	}
	mesh.points.back().level = cellLevel+1;
	
};




void addSubOuterFacesPoints(
	vector<int>& vertexPoints,
	int& centerPoint,
	vector<vector<int>>& subOutEdgePoints,
	vector<faces_refind>& faces
) {

	int vertexSize = (int)vertexPoints.size();
	for(int j=0; j<vertexSize; ++j){
		int leftEdge = ( j==0 ? vertexSize*2-1 : j*2-1 );
		int centerEdge = j*2; 
		int rightEdge = j*2+1;
		
		faces.push_back(faces_refind());
		auto& tmpFace = faces.back();
		for(auto& k : subOutEdgePoints[centerEdge]){
			tmpFace.points.push_back(k);
		}
		tmpFace.points.push_back(subOutEdgePoints[rightEdge][0]);
		tmpFace.points.push_back(centerPoint);
		for(auto& k : subOutEdgePoints[leftEdge]){
			tmpFace.points.push_back(k);
		}
		
	}
	
};



void extractSubInternalEdgesPoints(
	vector<int>& vertexCenterPoints,
	int& centerPoint,
	vector<vector<int>>& subIntEdgePoints
) {

	int verCentPSize = vertexCenterPoints.size();
	// cout << verCentPSize << endl;
	subIntEdgePoints.resize(verCentPSize);
	for(int j=0; j<verCentPSize; ++j){
		subIntEdgePoints[j].push_back(vertexCenterPoints[j]);
		subIntEdgePoints[j].push_back(centerPoint);
	}

};




void extractCellVertexOrder(
	vector<int>& cellVertexPoints,
	map<int,int>& cellVertexOrder
) {

	for(int j=0; j<cellVertexPoints.size(); ++j){
		cellVertexOrder.insert(make_pair(cellVertexPoints[j],j));
	}
	// int tmpNum = 0;
	// for(int j=0; j<cell.points.size(); ++j){
		// int vertex = cell.points[j];
		// if(mesh.points[vertex].level > cellLevel) continue;
		// cellVertexOrder.insert(make_pair(vertex,tmpNum));
		// ++tmpNum;
	// }
};


void extractInternalFacesVertexs(
	SEMO_Cell& cell,
	map<int,int>& cellInternalFaces,
	vector<int>& groupChildFaces_id,
	vector<groupMesh_Refine>& groupChildFaces,
	groupMesh_Refine& groupChildFace,
	int& cellCenterPoint,
	vector<vector<int>>& intFacesVertexs
) {

	intFacesVertexs.resize(cellInternalFaces.size());
	
	for(auto& j : cell.faces){
		auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
		
		int faceCenterPoint = groupOuterFace.centerPoint;
		
		int verCentPSize = groupOuterFace.vertexCenterPoints.size();
		for(int k=0; k<verCentPSize; ++k){
			int verCentP = groupOuterFace.vertexCenterPoints[k];
			int iFace = cellInternalFaces[verCentP];
			auto& face = groupChildFace.faces[iFace];
			// int facePointsSize = face.points.size();
			
			if(
			(std::find(face.points.begin(),face.points.end(),cellCenterPoint) 
				== face.points.end()) 
			){
				
				// cout << "POINTS0 : " << face.points.size() << endl;
				// cout << groupOuterFace.subIntEdgePoints.size() << endl;
				for(auto& l : groupOuterFace.subIntEdgePoints[k]){
					face.points.push_back(l);
				}
				face.points.push_back(cellCenterPoint);
				// cout << "POINTS1 : " << face.points.size() << endl;
			}
			else if(
			(std::find(face.points.begin(),face.points.end(),faceCenterPoint) 
				== face.points.end())
			){
				// cout << "POINTS0 : " << face.points.size() << endl;
				int subEdgeSize = groupOuterFace.subIntEdgePoints[k].size();
				for(int l=subEdgeSize-1; l>=1; --l){
					face.points.push_back(groupOuterFace.subIntEdgePoints[k][l]);
				}
				
				// cout << "POINTS1 : " << face.points.size() << endl;
				
			}
			
			
			if(verCentPSize==2){
				int verPoint = groupOuterFace.vertexPoints[0];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint);
				}
				
			}
			else{
				int verPoint0 = groupOuterFace.vertexPoints[k];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint0) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint0);
				}
				
				int verPoint1 = groupOuterFace.vertexPoints[ ( (k==verCentPSize-1) ? 0 : k+1 ) ];
				if(
				(std::find(intFacesVertexs[iFace].begin(),intFacesVertexs[iFace].end(),verPoint1) 
				== intFacesVertexs[iFace].end())
				){
					intFacesVertexs[iFace].push_back(verPoint1);
				}
				
			}
		}
	}
};



void addInternalFacesOwnerNeighbour(
	SEMO_Mesh_Builder& mesh, 
	int& cellLevel,
	int& totalCellNum,
	map<int,int>& cellVertexOrder,
	vector<vector<int>>& intFacesVertexs,
	vector<faces_refind>& faces
) {

	
	for(int j=0; j<faces.size(); ++j){
		// cout << intFacesVertexs[j].size() << endl;
		int iOwner = intFacesVertexs[j][0];
		int iNeighbour = intFacesVertexs[j][1];
		
		vector<int> tmpVertexPoints;
		for(auto& j : faces[j].points){
			if(mesh.points[j].level > cellLevel+1) continue;
			tmpVertexPoints.push_back(j);
		}
		
		vector<double> unitNormals(3,0.0);
		double x1 = mesh.points[tmpVertexPoints[0]].x;
		double x2 = mesh.points[tmpVertexPoints[1]].x;
		double x3 = mesh.points[tmpVertexPoints[2]].x;
		double y1 = mesh.points[tmpVertexPoints[0]].y;
		double y2 = mesh.points[tmpVertexPoints[1]].y;
		double y3 = mesh.points[tmpVertexPoints[2]].y;
		double z1 = mesh.points[tmpVertexPoints[0]].z;
		double z2 = mesh.points[tmpVertexPoints[1]].z;
		double z3 = mesh.points[tmpVertexPoints[2]].z;
		double v1[3] = {x2-x1,y2-y1,z2-z1};
		double v2[3] = {x3-x1,y3-y1,z3-z1};
		
		unitNormals[0] = v1[1] * v2[2] - v1[2] * v2[1];
		unitNormals[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
		unitNormals[2] = v1[0] * v2[1] - v1[1] * v2[0];
		
		double ownVecX = mesh.points[iOwner].x-mesh.points[iNeighbour].x;
		double ownVecY = mesh.points[iOwner].y-mesh.points[iNeighbour].y;
		double ownVecZ = mesh.points[iOwner].z-mesh.points[iNeighbour].z;
		
		double cosTheta = 
			unitNormals[0]*ownVecX + 
			unitNormals[1]*ownVecY + 
			unitNormals[2]*ownVecZ;
		
		// cout << cosTheta << endl;
		
		if(cosTheta>=0.0){
			faces[j].neighbour = totalCellNum + cellVertexOrder[iOwner];
			faces[j].owner = totalCellNum + cellVertexOrder[iNeighbour];
		}
		else{
			faces[j].owner = totalCellNum + cellVertexOrder[iOwner];
			faces[j].neighbour = totalCellNum + cellVertexOrder[iNeighbour];
		}
		
		
	}
	
			
};


void addOuterFacesOwnerNeighbour(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Cell& cell,
	int& i,
	int& totalCellNum,
	map<int,int>& cellVertexOrder,
	vector<int>& groupChildFaces_id,
	vector<groupMesh_Refine>& groupChildFaces
) {
	for(auto& j : cell.faces){
		auto& face = mesh.faces[j];
		if(groupChildFaces_id[j]==-1) cout << "NONONONONON" << endl;
		auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
		
		int verPSize = groupOuterFace.vertexPoints.size();
		for(int k=0; k<verPSize; ++k){
			int verP = groupOuterFace.vertexPoints[k];
			auto& subFace = groupOuterFace.faces[k];
			if(face.owner == i){
				subFace.owner = totalCellNum + cellVertexOrder[verP];
				cout << " owner = " << subFace.owner << endl;
			}
			else{
				subFace.neighbour = totalCellNum + cellVertexOrder[verP];
				cout << " neighbour = " << subFace.neighbour << endl;
			}
		}
	}

};


void reorderOuterFacesOwnerNeighbour(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Cell& cell,
	int& i,
	int& totalCellNum,
	vector<int>& groupChildFaces_id,
	vector<bool>& boolInputFacesOwner,
	vector<bool>& boolInputFacesNeighbour, 
	vector<groupMesh_Refine>& groupChildFaces
) {

	for(auto& j : cell.faces){
		auto& face = mesh.faces[j];
		
		if(groupChildFaces_id[j]==-1){
			if(face.owner == i && boolInputFacesOwner[j]==false){
				face.owner = totalCellNum;
				boolInputFacesOwner[j] = true;
				// cout << " owner = " << face.owner << " " << face.level << endl;
			}
			else if(face.neighbour == i && boolInputFacesNeighbour[j]==false){
				face.neighbour = totalCellNum;
				boolInputFacesNeighbour[j] = true;
				// cout << " neighbour = " << face.neighbour << " " << face.level << endl;
			}
			else{
				cout << "ERROR 25" << endl;
				cout << face.owner << " " <<  face.neighbour << " " << i << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}
		else{
			auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
			for(auto& subFace : groupOuterFace.faces){
				if(face.owner == i && subFace.boolOwner==false){
					subFace.owner = totalCellNum;
					subFace.boolOwner = true;
					// cout << " owner = " << subFace.owner << endl;
				}
				else if(face.neighbour == i && subFace.boolNeighbour==false){
					subFace.neighbour = totalCellNum;
					subFace.boolNeighbour = true;
					// cout << " neighbour = " << subFace.neighbour << endl;
				}
				else{
					cout << "ERROR 35" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
		}
	}
};























void SEMO_Poly_AMR_Builder::polyRefine(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	
	
	int beforeCellSize = mesh.cells.size();
	int beforeFaceSize = mesh.faces.size();
	int beforePointSize = mesh.points.size();
	int afterCellSize = 0;
	int afterFaceSize = 0;
	int afterPointSize = 0;
	
	
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// for(auto& j : mesh.faces[i].points){
				// if(rank==1) cout << rank << " " << i << " " << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << " " << endl;
			// }
		// }
	// }
	
	
	// cout << mesh.cells[53001].points.size() << endl;
	// cout << mesh.cells[53001].faces.size() << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[0]].points.size() << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[1]].points.size() << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[2]].points.size() << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[3]].points.size() << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[4]].points.size() << endl;
	
	// cout << mesh.faces[mesh.cells[53001].faces[0]].owner << " " << mesh.faces[mesh.cells[53001].faces[0]].neighbour << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[1]].owner << " " << mesh.faces[mesh.cells[53001].faces[1]].neighbour << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[2]].owner << " " << mesh.faces[mesh.cells[53001].faces[2]].neighbour << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[3]].owner << " " << mesh.faces[mesh.cells[53001].faces[3]].neighbour << endl;
	// cout << mesh.faces[mesh.cells[53001].faces[4]].owner << " " << mesh.faces[mesh.cells[53001].faces[4]].neighbour << endl;
	
	// cout << mesh.points[mesh.cells[53001].points[0]].x << " " << mesh.points[mesh.cells[53001].points[0]].y << " " << mesh.points[mesh.cells[53001].points[0]].z << endl;
	
	
	
	if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute AMR - Refine " << endl;
	}
	
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
	
	

	// vector<bool> boolCanNotRefineCells(mesh.cells.size(),false);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(mesh.cells[i].level < 0) boolCanNotRefineCells[i] = true;
	// }
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].level < 0) mesh.faces[i].level = 0;
	// }
	
	
	int maxLevel_AMR = controls.maxLevelRefine;
	// double indicatorRefine_AMR = controls.indicatorRefine;
	double minVolume_AMR = controls.minVolumeRefine;
	int maxCells_AMR = controls.maxCellsRefine;
	
	vector<bool> boolCellRefine(mesh.cells.size(),false);
	if(mesh.cells.size()<maxCells_AMR){
		for(int i=0; i<mesh.cells.size(); ++i){
			
			// if( distr(eng) > 0.8 ){
				// // cout << "CELL REFINE : " << i << endl;
				// boolCellRefine[i] = true;
			// }

			// if(mesh.cells[i].var[controls.indicatorAMR[0]] > indicatorRefine_AMR) 
				// boolCellRefine[i] = true;
			
			for(int level=0; level<controls.indicatorCriterion.size(); ++level){
				double indicatorRefine_AMR = controls.indicatorCriterion[level];
				if( mesh.cells[i].level < level+1 ){
					if( mesh.cells[i].var[controls.indicatorAMR[0]] > indicatorRefine_AMR ){
						boolCellRefine[i] = true;
					}
				}					
			}
			
			
			if(mesh.cells[i].volume < minVolume_AMR) boolCellRefine[i] = false;
			if(mesh.cells[i].level >= maxLevel_AMR) boolCellRefine[i] = false;
			if(mesh.cells[i].level < 0) boolCellRefine[i] = false;
			
			// if(boolCanNotRefineCells[i]==true) {
				// boolCellRefine[i] = false;
				// // mesh.cells[i].level = 0;
			// }
			
			
		} 
		
		
		// vector<bool> boolCellRefine_trm = boolCellRefine;
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			// int cellLevel = cell.level;
			// if(boolCellRefine[i] == false){
				
				// bool boolNgbRefine = false;
				// bool boolNgbNoRefine = false;
				// for(auto& j : cell.faces){
					// auto& face = mesh.faces[i];
					// int ngb = face.neighbour;
					// if(face.owner != i) ngb=face.owner;
					
					// if(ngb!=-1){
						// auto& cellNgb = mesh.cells[ngb];
						// int cellLevelNgb = cellNgb.level;
						
						// if(boolCellRefine[ngb] == true) boolNgbRefine=true;
						// if(boolCellRefine[ngb] == false) boolNgbNoRefine=true;
					// }
				// }
				
				// if(boolNgbRefine==true && boolNgbNoRefine==true){
					// boolCellRefine_trm[i] = true;
				// }
				
			// }
		// } 
		// boolCellRefine = boolCellRefine_trm;
		
		
		
	}
	
	
	// searching same level neighbour
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// if(
			// boolCellRefine[face.owner]==true &&
			// boolCellRefine[face.neighbour]==false
			// ){
				
				
				
			// }
		
			
		// }
		
		// if(mesh.cells[i].var[controls.indicatorAMR[0]] > indicatorRefine_AMR) 
			// boolCellRefine[i] = true;
		
		
		// if(mesh.cells[i].volume < minVolume_AMR) boolCellRefine[i] = false;
		// if(mesh.cells[i].level >= maxLevel_AMR) boolCellRefine[i] = false;
		// if(mesh.cells[i].level < 0) boolCellRefine[i] = false;
		
		// // if(boolCanNotRefineCells[i]==true) {
			// // boolCellRefine[i] = false;
			// // // mesh.cells[i].level = 0;
		// // }
		
		
	// } 

	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// auto& cell = mesh.cells[i];
		
		// if(boolCellRefine[i]==true){
			
			// for(auto k : cell.stencil){
				// auto& cellSten = mesh.cells[k];
				
				// // if(cellSten.level != cell.level){
					
					// boolCellRefine[k] = true;
					
				// // }
				// if(cellSten.volume < minVolume_AMR) boolCellRefine[k] = false;
				// if(cellSten.level >= maxLevel_AMR) boolCellRefine[k] = false;
				// if(cellSten.level < 0) boolCellRefine[k] = false;
				
				// if(boolCellRefine[k]){
					// mesh.cells[k].var[controls.indicatorAMR[0]] = 1.e15;
				// }
				
			// }
			// mesh.cells[i].var[controls.indicatorAMR[0]] = 1.e15;
		// }
		
	// }
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// auto& cell = mesh.cells[i];
		
		// int level = cell.level;
		
		// if(
		// boolCellRefine[i]==false && 
		// cell.level>=1 && 
		// cell.level<maxLevel_AMR &&
		// mesh.cells[i].volume >= minVolume_AMR
		// ){
			
			// bool levelzero_appear=false;
			// bool levelsame_appear=false;
			
			// for(auto k : cell.stencil){
				// auto& cellSten = mesh.cells[k];
				// if(cellSten.level==level-1){
					// levelzero_appear = true;
				// }
				// if(cellSten.level==level){
					// levelsame_appear = true;
				// }
			// }
			
			// if(levelzero_appear==true && levelsame_appear==false){
				// boolCellRefine[i] = true;
				// mesh.cells[i].var[controls.indicatorAMR[0]] = 1.e15;
			// }
		// }
		
	// }
	
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// auto& cell = mesh.cells[i];
		
		// if(boolCellRefine[i]==true){
			
			// for(auto k : cell.stencil){
				// auto& cellSten = mesh.cells[k];
				
				// // if(cellSten.level != cell.level){
					
					// boolCellRefine[k] = true;
					
				// // }
				// if(cellSten.volume < minVolume_AMR) boolCellRefine[k] = false;
				// if(cellSten.level >= maxLevel_AMR) boolCellRefine[k] = false;
				// if(cellSten.level < 0) boolCellRefine[k] = false;
				
				// if(boolCellRefine[k]){
					// mesh.cells[k].var[controls.indicatorAMR[0]] = 1.e15;
				// }
				
			// }
		// }
		
	// }
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 1" << endl;
	
	vector<int> cLevel_recv;
	vector<int> cRefine_recv;
	this->mpiLevelRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	
	
	this->restrictCellRefine(mesh, boolCellRefine, cLevel_recv, cRefine_recv);
	
	
	this->mpiRefines(mesh, boolCellRefine, cRefine_recv);
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if( boolCellRefine[i] ){
			// cout << rank << " CELL REFINE : " << i << endl;
		// }
	// } 
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 2" << endl;

	

	// if(rank==1){
		// for(int ip=0; ip<10000; ++ip){
			// for(int i=0; i<mesh.faces.size(); ++i){
				// // cout << rank << " " << i << " " << mesh.faces[i].level << endl;
				// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
					// for(auto& j : mesh.faces[i].points){
					// }
				// }
			// }
		// }
	// }
	// int proc_num2 = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// cout << rank << " " << i << " " << mesh.cells[mesh.faces[i].owner].level << " " << cLevel_recv[proc_num2] << " " << boolCellRefine[mesh.faces[i].owner] << " " << cRefine_recv[proc_num2] << endl;
			// ++proc_num2;
			// // cout << rank << " " << i << " " << mesh.faces[i].level << endl;
			// // ++proc_num;
			// // for(auto& j : mesh.faces[i].points){
				// // cout << rank << " " << i << " " << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << " " << endl;
			// // }
		// }
	// }
	
	
	
	
	
	
	

	vector<int> edgesPoint0;
	vector<int> edgesPoint1;
	vector<vector<int>> facesEdges;
	vector<vector<int>> edgesFaces;
	vector<int> edgesLevel;
	
	this->createEdges(mesh, edgesPoint0, edgesPoint1, facesEdges, edgesFaces, edgesLevel);
	
	
	vector<bool> boolEdgeRefine(edgesPoint0.size(),false);
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 3" << endl;
	
	// cout << "0000 : " << rank << endl;
	
	
	//====================================================
	vector<int> groupChildFaces_id(mesh.faces.size(),-1);
	vector<bool> groupChildFaces_HighLevel(mesh.faces.size(),false);
	vector<groupMesh_Refine> groupChildFaces;
	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];
		int faceLevel = face.level;
		
		bool ownRefine = boolCellRefine[face.owner];
		bool ngbRefine;
		int ownLevel = mesh.cells[face.owner].level;
		int ngbLevel;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			ngbRefine = boolCellRefine[face.neighbour];
			ngbLevel = mesh.cells[face.neighbour].level;
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			ngbRefine = false;
			if(cRefine_recv[proc_num]==1) ngbRefine = true;
			ngbLevel = cLevel_recv[proc_num];
			++proc_num;
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			ngbRefine = false;
			ngbLevel = faceLevel;
		}
			
			
		if(
		(ownRefine==true && ownLevel==faceLevel) ||
		(ngbRefine==true && ngbLevel==faceLevel) 
		){
			
			// cout << "1-1" << endl;
			
			groupChildFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			
			groupChildFace.parent = i;
			groupChildFace.type = 0;
			
			vector<int> tmpFacePoints;
			vector<int> tmpFacePointLevels;
			for(auto& j : face.points){
				// cout << i << " " << j << endl;
				tmpFacePoints.push_back(j);
				tmpFacePointLevels.push_back(mesh.points[j].level);
			}
			
			extractVertexPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				groupChildFace.vertexPoints);

			
			// cout << groupChildFace.vertexPoints.size() << endl;
			// cout << faceLevel << endl;
			// cout << ownRefine << " " << ownLevel << endl;
			// cout << ngbRefine << " " << ngbLevel << endl;
			
			vector<bool> addVertexCeterPoints;
			extractVertexCenterPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints,
				addVertexCeterPoints);
				
			vector<int> canRefineEdgeOrders;
			for(int j=0; j<facesEdges[i].size(); ++j){
				int iEdge = facesEdges[i][j];
				if( 
				edgesLevel[iEdge]==faceLevel && 
				boolEdgeRefine[iEdge]==false
				){
					canRefineEdgeOrders.push_back(iEdge);
				}
			}
			
			
			// for(int j=0; j<tmpFacePoints.size(); ++j){
				// cout << j << " vertex " << tmpFacePoints[j] << " " << tmpFacePointLevels[j] << endl;
			// }
			// for(int j=0; j<groupChildFace.vertexPoints.size(); ++j){
				// cout << j << " extVertex " << groupChildFace.vertexPoints[j] << endl;
			// }
			// for(int j=0; j<groupChildFace.vertexCenterPoints.size(); ++j){
				// cout << j << " extVertexCent " << groupChildFace.vertexCenterPoints[j] << endl;
			// }
			// for(int j=0; j<addVertexCeterPoints.size(); ++j){
				// cout << j << " extBoolVertexCent " << addVertexCeterPoints[j] << endl;
			// }
			
			
			if( groupChildFace.vertexPoints.size() != groupChildFace.vertexCenterPoints.size() ){
				cout << "| WARNING : vertexPoints.size != vertexCenterPoints.size " << groupChildFace.vertexPoints.size() << " " << groupChildFace.vertexCenterPoints.size() << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
				
			// cout << "1-2" << endl;
			
			addVertexCenterPoint(mesh, faceLevel, 
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints);
				
			// cout << "1-3" << endl;
			
			addEdgeFacesPointsToVertexCenterPoint(
				mesh,
				i,
				facesEdges,
				edgesFaces,
				boolEdgeRefine,
				addVertexCeterPoints,
				canRefineEdgeOrders,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints);
				
			// cout << "1-4" << endl;
					
			
		}
		else if(
		(ownRefine==true && ownLevel+1==faceLevel) ||
		(ngbRefine==true && ngbLevel+1==faceLevel) 
		){
			// cout << "2-1" << endl;
			
			groupChildFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces_HighLevel[i] = true;
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			
			groupChildFace.parent = i;
			groupChildFace.type = 0;
			
			vector<int> tmpFacePoints;
			vector<int> tmpFacePointLevels;
			for(auto& j : face.points){
				tmpFacePoints.push_back(j);
				tmpFacePointLevels.push_back(mesh.points[j].level);
			}
			
			vector<int> vertexPoints;
			extractVertexPoints(faceLevel, 
				tmpFacePoints, tmpFacePointLevels,
				vertexPoints);
			
			groupChildFace.vertexPoints.push_back(vertexPoints[0]);
			groupChildFace.vertexPoints.push_back(vertexPoints[1]);
			groupChildFace.vertexPoints.push_back(vertexPoints[2]);
			groupChildFace.vertexPoints.push_back(vertexPoints[3]);
			
			groupChildFace.vertexCenterPoints.push_back(vertexPoints[1]);
			groupChildFace.centerPoint = vertexPoints[2];
			groupChildFace.vertexCenterPoints.push_back(vertexPoints[3]);
			
			
			// cout << "2-2" << endl;
			
		}
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 4" << endl;
	// cout << "1111 : " << rank << endl;
	
	
	//====================================================
	
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];
		int faceLevel = face.level;
		
		if(
		groupChildFaces_id[i] != -1 &&
		groupChildFaces_HighLevel[i] == false){
			
			// cout << "1 " << faceLevel << endl;
			
			auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
			addFaceCenterPoint(mesh, faceLevel, 
				groupChildFace.vertexPoints,
				groupChildFace.centerPoint);
				
			vector<int> tmpFacePoints;
			for(auto& j : face.points){
				tmpFacePoints.push_back(j);
			}
			
			extractSubEdgePoints(
				tmpFacePoints,
				groupChildFace.vertexPoints,
				groupChildFace.vertexCenterPoints,
				groupChildFace.subOutEdgePoints);
				
			addSubOuterFacesPoints(
				groupChildFace.vertexPoints,
				groupChildFace.centerPoint,
				groupChildFace.subOutEdgePoints,
				groupChildFace.faces);
				
			extractSubInternalEdgesPoints(
				groupChildFace.vertexCenterPoints,
				groupChildFace.centerPoint,
				groupChildFace.subIntEdgePoints);
				
			// cout << groupChildFace.subIntEdgePoints[0].size() << " " << groupChildFace.subIntEdgePoints[1].size() << endl;
				
			// cout << "2" << endl;
			
			for(auto& subFace : groupChildFace.faces){
				if(subFace.points.size()>100){
					cout << "| WARNING 1 : face.points.size > 10" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
				
				
		}
		else if(
		groupChildFaces_id[i] != -1 &&
		groupChildFaces_HighLevel[i] == true){
			
			// cout << "2-1 " << faceLevel << endl;
			
			auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
			
			vector<int> tmpFacePoints;
			for(auto& j : face.points){
				tmpFacePoints.push_back(j);
			}
			
			
			
			// extractSubEdgePoints(
				// tmpFacePoints,
				// groupChildFace.vertexPoints,
				// groupChildFace.vertexCenterPoints,
				// // groupChildFace.vertexPoints,
				// groupChildFace.subOutEdgePoints);
			

			

			// for(int j=0; j<tmpFacePoints.size(); ++j){
				// cout << j << " vertex " << tmpFacePoints[j] << endl;
			// }
			// for(int j=0; j<groupChildFace.vertexPoints.size(); ++j){
				// cout << j << " extVertex " << groupChildFace.vertexPoints[j] << endl;
			// }
			// for(int j=0; j<groupChildFace.vertexCenterPoints.size(); ++j){
				// cout << j << " extVertexCent " << groupChildFace.vertexCenterPoints[j] << endl;
			// }
			// // for(int j=0; j<groupChildFace.subOutEdgePoints.size(); ++j){
				// // for(int k=0; k<groupChildFace.subOutEdgePoints[j].size(); ++k){
					// // cout << j << " subOutEdgePoints " << k << " " << groupChildFace.subOutEdgePoints[j][k] << endl;
				// // }
			// // }
			
			
			groupChildFace.faces.push_back(faces_refind());
			for(auto& j : face.points){
				groupChildFace.faces.back().points.push_back(j);
			}
			
			// cout << "2-3 " << faceLevel << endl;
			// cout << "2-3 " << groupChildFace.vertexPoints.size() << " " << groupChildFace.vertexCenterPoints.size() << " " << groupChildFace.subOutEdgePoints.size() << endl;
			
			
			
			// groupChildFace.subIntEdgePoints.resize(2);
			// // for(auto& j : groupChildFace.subOutEdgePoints[1]){
			// for(auto& j : groupChildFace.subOutEdgePoints[2]){
				// groupChildFace.subIntEdgePoints[0].push_back(j);
			// }
			// // cout << groupChildFace.subIntEdgePoints[0].size() << endl;
			// // groupChildFace.subIntEdgePoints[0].push_back(groupChildFace.subOutEdgePoints[2][0]);
			// // cout << groupChildFace.subIntEdgePoints[0].size() << endl;
			
			// // cout << "2-4 " << faceLevel << endl;
			
			// groupChildFace.subIntEdgePoints[1].push_back(groupChildFace.subOutEdgePoints[3][0]);
			// for(auto& j : groupChildFace.subOutEdgePoints[2]){
				// groupChildFace.subIntEdgePoints[1].push_back(j);
			// }
			
			// // cout << "2-5 " << faceLevel << endl;
			
			// std::reverse(
				// groupChildFace.subIntEdgePoints[1].begin(),
				// groupChildFace.subIntEdgePoints[1].end());
				
				
			auto iter0 = std::find(face.points.begin(), face.points.end(), groupChildFace.vertexPoints[1]);
			auto iter1 = std::find(face.points.begin(), face.points.end(), groupChildFace.vertexPoints[2]);
			auto iter2 = std::find(face.points.begin(), face.points.end(), groupChildFace.vertexPoints[3]);
				
			groupChildFace.subIntEdgePoints.resize(2);
			// cout << "AAAAAAA" << endl;
			std::copy(iter0, iter1+1, std::back_inserter(groupChildFace.subIntEdgePoints[0]));
			// cout << groupChildFace.subIntEdgePoints[0][0] << " " << groupChildFace.subIntEdgePoints[0][1] << endl;
			// cout << "BBBBBBB" << endl;
			std::copy(iter1, iter2+1, std::back_inserter(groupChildFace.subIntEdgePoints[1]));
			std::reverse(
				groupChildFace.subIntEdgePoints[1].begin(),
				groupChildFace.subIntEdgePoints[1].end());
			
			
			// cout << groupChildFace.subIntEdgePoints[0].size() << " " << groupChildFace.subIntEdgePoints[1].size() << endl;
			
			// cout << "2-6" << endl;
			
			if(face.points.size()>100){
				cout << "| WARNING 2 : face.points.size > 10" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
		}
		
			
	}
	
	
	
	
	
	// cout << "2222 : " << rank << endl;
	

	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			// if(groupChildFaces_id[i] != -1){
				// auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
				// for(auto& face : groupChildFace.faces){
					// if(face.points.size()>100){
						// cout << mesh.faces[i].level << " " << mesh.cells[mesh.faces[i].owner].level << " " << mesh.cells[mesh.faces[i].neighbour].level << " " << face.points.size() << endl;
						// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					// }
				// }
			// }
		// }
	// }
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 5" << endl;
	
	//====================================================
	// createCellInternalFaces
	
	vector<int> groupCellInternalFaces_id(mesh.cells.size(),-1);
	
	vector<bool> boolInputFacesOwner(mesh.faces.size(),false);
	vector<bool> boolInputFacesNeighbour(mesh.faces.size(),false);
	
	// vector<vector<int>> groupCell_Levels(mesh.cells.size(),vector<int>());
	// vector<vector<int>> groupCell_Groups(mesh.cells.size(),vector<int>());
	vector<vector<int>> groupCell_Levels;
	vector<vector<int>> groupCell_Groups;
	int totalCellNum = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		int cellLevel = cell.level;
		int cellGroup = cell.group;
		
		if(boolCellRefine[i]==true){
			
			groupCellInternalFaces_id[i] = groupChildFaces.size();
			
			groupChildFaces.push_back(groupMesh_Refine());
			auto& groupChildFace = groupChildFaces.back();
			

			vector<int> cellVertexPoints;
			for(auto& j : cell.points){
				if(mesh.points[j].level > cellLevel) continue;
				cellVertexPoints.push_back(j);
			}
			
			int cellCenterPoint;
			addCellCenterPoint(mesh, cellLevel, 
				cellVertexPoints, cellCenterPoint);
			
			
			map<int,int> cellInternalFaces;
			for(auto& j : cell.faces){
				auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
				for(auto& k : groupOuterFace.vertexCenterPoints){
					if(cellInternalFaces.find(k) == cellInternalFaces.end()){
						int tmpSize = cellInternalFaces.size();
						cellInternalFaces.insert(make_pair(k,tmpSize));
						groupChildFace.faces.push_back(faces_refind());
					}
				}
			}
			
			
			map<int,int> cellVertexOrder;
			extractCellVertexOrder(cellVertexPoints, cellVertexOrder);
			
			
			vector<vector<int>> intFacesVertexs;
			extractInternalFacesVertexs(cell, cellInternalFaces,
				groupChildFaces_id, groupChildFaces, groupChildFace,
				cellCenterPoint, intFacesVertexs);
				
			for(auto& test : intFacesVertexs){
				if( test.size() != 2 ){
					cout << "| WARNING : no own ngb" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
				
			// cout << i << endl;
			addInternalFacesOwnerNeighbour(mesh, cellLevel, totalCellNum,
				cellVertexOrder, intFacesVertexs, groupChildFace.faces);
				
				
			// if(iter==0){
			// addOuterFacesOwnerNeighbour(mesh, cell, i, totalCellNum,
				// cellVertexOrder, groupChildFaces_id, groupChildFaces);
				
			// for(auto& j : cellVertexOrder){
				// cout << j.second << endl;
			// }
				

			// if(iter==0){
			for(auto& j : cell.faces){
				auto& face = mesh.faces[j];
				// if(groupChildFaces_id[j]==-1) cout << "NONONONONON" << endl;
				auto& groupOuterFace = groupChildFaces[groupChildFaces_id[j]];
				
				// int verPSize = groupOuterFace.vertexPoints.size();
				
				int verPSize = groupOuterFace.faces.size();
				
				// if( verPSize != groupOuterFace.faces.size() ){
					// cout << "| WARNING : not matching 1 : " << verPSize << " " << groupOuterFace.faces.size() << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				for(int k=0; k<verPSize; ++k){
					int verP = groupOuterFace.vertexPoints[k];
					auto& subFace = groupOuterFace.faces[k];
					// cout << subFace.points.size() << endl;
					
					if(cellVertexOrder.find(verP) == cellVertexOrder.end()){
						cout << "| WARNING : no find cellVertexOrder.find(verP)" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					
					// cout << k << " " << verP << " " << cellVertexOrder[verP] << endl;
					if(face.owner == i && subFace.boolOwner==false){
						subFace.owner = totalCellNum + cellVertexOrder[verP];
						subFace.boolOwner = true;
						// cout << " owner = " << subFace.owner << endl;
					}
					else if(face.neighbour == i && subFace.boolNeighbour==false){
						subFace.neighbour = totalCellNum + cellVertexOrder[verP];
						subFace.boolNeighbour = true;
						// cout << " neighbour = " << subFace.neighbour << endl;
					}
					else{
						cout << "NONONONONONONONONO 33" << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
				}
			}
			// }
				
				
			// }
			
			groupCell_Levels.push_back(vector<int>());
			for(int jj=0; jj<cellVertexOrder.size(); ++jj){
				groupCell_Levels.back().push_back(cellLevel+1);
			}
			
			groupCell_Groups.push_back(vector<int>());
			for(int jj=0; jj<cellVertexOrder.size(); ++jj){
				groupCell_Groups.back().push_back(cellGroup);
			}
			
			totalCellNum += cellVertexOrder.size();
			// cout << "1-2" << endl;
			
			// cout << cellVertexOrder.size() << endl;
			// for(auto& j : groupChildFace.faces){
				// cout << i << " " << j.points.size() << endl;
			// }
			
			
			// for(auto& j : cell.faces){
				// if(groupChildFaces_id[j] != -1){
					// auto& groupChildFace = groupChildFaces[groupChildFaces_id[j]];
					// for(auto& face : groupChildFace.faces){
						// if(face.points.size()>100){
							// cout << mesh.faces[j].level << " " << mesh.cells[mesh.faces[j].owner].level << " " << mesh.cells[mesh.faces[j].neighbour].level << " " << face.points.size() << endl;
							// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
						// }
					// }
				// }
			// }
			// }
			
		}
		else{
			
			// cout << "2-1" << endl;
			
			reorderOuterFacesOwnerNeighbour(mesh, cell, i, totalCellNum,
				groupChildFaces_id, boolInputFacesOwner, boolInputFacesNeighbour, 
				groupChildFaces);
				
			// groupCell_Levels[i].push_back(cellLevel);
			// groupCell_Groups[i].push_back(cellGroup);
			
			
			groupCell_Levels.push_back(vector<int>());
			// if(boolCanNotRefineCells[i]==true) {
				// groupCell_Levels.back().push_back(-1);
			// }
			// else{
				groupCell_Levels.back().push_back(cellLevel);
			// }
			
			groupCell_Groups.push_back(vector<int>());
			groupCell_Groups.back().push_back(cellGroup);
			
			++totalCellNum;
			
			// cout << "2-2" << endl;
		}
		
	}
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 6" << endl;
	
	// cout << "3333 : " << rank << endl;
	
	
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			// if(groupChildFaces_id[i] != -1){
				// auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
				// for(auto& face : groupChildFace.faces){
					// if(face.points.size()==0){
						// cout << mesh.faces[i].level << " " << mesh.cells[mesh.faces[i].owner].level << " " << mesh.cells[mesh.faces[i].neighbour].level << " " << face.points.size() << endl;
						// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					// }
				// }
			// }
		// }
	// }
	
	
	
	//====================================================
	
	
	// cout << totalCellNum << " " << mesh.cells.size() << endl;
	
	
	int faceResizeNum = 0;
	vector<int> startFaces(mesh.faces.size(),-1);
	
	// cell outer faces
	for(int i=0; i<mesh.faces.size(); ++i){
		if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			if(groupChildFaces_id[i] != -1){
				auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
				startFaces[i] = faceResizeNum;
				faceResizeNum += groupChildFace.faces.size();
			}
			else{
				startFaces[i] = faceResizeNum;
				++faceResizeNum;
			}
		}
	}
	
	// cell internal faces
	vector<int> startIntFaces(mesh.cells.size(),-1);
	for(int i=0; i<mesh.cells.size(); ++i){
		if(groupCellInternalFaces_id[i] != -1){
			auto& groupChildFace = groupChildFaces[groupCellInternalFaces_id[i]];
			startIntFaces[i] = faceResizeNum;
			faceResizeNum += groupChildFace.faces.size();
		}
	}
	
	// Proc & B.C. faces
	for(int i=0; i<mesh.faces.size(); ++i){
		if(mesh.faces[i].getType() != SEMO_Types::INTERNAL_FACE){
			if(groupChildFaces_id[i] != -1){
				auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
				startFaces[i] = faceResizeNum;
				faceResizeNum += groupChildFace.faces.size();
			}
			else{
				startFaces[i] = faceResizeNum;
				++faceResizeNum;
			}
		}
	}
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 7" << endl;
	// cout << "4444 : " << rank << endl;
	
	
	// cout << faceResizeNum << " " << mesh.faces.size() << endl;
	
	int orgFaceSize = mesh.faces.size();
	int orgCellSize = mesh.cells.size();
	
	mesh.faces.resize(faceResizeNum);
	
	// cout << faceResizeNum << " " << mesh.faces.size() << endl;
	
	
	// faceResizeNum = 0;
	
	
	
	// cout << "5555 : " << rank << endl;
	
	
	// Proc & B.C. faces	
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1) continue;
		// cout << "AAAAAAAAAAAAAA :::: " << rank << " " << boundary.neighbProcNo << endl;
		if(boundary.neighbProcNo <= rank) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			if(groupChildFaces_id[i] != -1){
				auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];
				std::reverse(groupChildFace.faces.begin()+1,groupChildFace.faces.end());
			}
		}
	}
	
	int saveI = 0;
	for(int i=orgFaceSize-1; i>=0; --i){
		
		if(mesh.faces[i].getType() != SEMO_Types::INTERNAL_FACE){
		
			int str = startFaces[i];
			
			// cout << str<< endl;
			
			if(groupChildFaces_id[i] != -1){
				auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];

				int tmpNum = 0;
				for(auto& face : groupChildFace.faces){
					mesh.faces[str+tmpNum].points.clear();
					for(auto& j : face.points){
						mesh.faces[str+tmpNum].points.push_back(j);
					}
					mesh.faces[str+tmpNum].owner = face.owner;
					mesh.faces[str+tmpNum].neighbour = -1;
					mesh.faces[str+tmpNum].setType(mesh.faces[i].getType());
					
					
					++tmpNum;
				}
			}
			else{
				mesh.faces[str] = mesh.faces[i];
			}
			
			saveI = i;
		
		}
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 8" << endl;
	
	// cout << "6666 : " << rank << endl;
	
	
	
	// cell internal faces
	for(int i=orgCellSize-1; i>=0; --i){
		if(groupCellInternalFaces_id[i] != -1){
			int str = startIntFaces[i];
		
			// cout << str<< endl;
			
			auto& groupChildFace = groupChildFaces[groupCellInternalFaces_id[i]];

			int tmpNum = 0;
			for(auto& face : groupChildFace.faces){
				mesh.faces[str+tmpNum].points.clear();
				for(auto& j : face.points){
					mesh.faces[str+tmpNum].points.push_back(j);
				}
				mesh.faces[str+tmpNum].owner = face.owner;
				mesh.faces[str+tmpNum].neighbour = face.neighbour;
				mesh.faces[str+tmpNum].setType(SEMO_Types::INTERNAL_FACE);
				
				
				++tmpNum;
			}
		}
	}
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 9" << endl;
	
	
	// cout << "7777 : " << rank << endl;
	// cout << saveI << endl;
	
	
	// cell outer faces
	for(int i=saveI-1; i>=0; --i){
		
		if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
		
			int str = startFaces[i];
			// cout << "AAAAAA" << str<< endl;
			
			if(groupChildFaces_id[i] != -1){
				
				auto& groupChildFace = groupChildFaces[groupChildFaces_id[i]];

				int tmpNum = 0;
				for(auto& face : groupChildFace.faces){
				// cout << "1-1 " << groupChildFace.faces.size() << " " << face.points.size() << endl;
					mesh.faces[str+tmpNum].points.clear();
					for(auto& j : face.points){
				// cout << "1-2 " << face.owner << endl;
						mesh.faces[str+tmpNum].points.push_back(j);
				// cout << "1-3 " << endl;
					}
					mesh.faces[str+tmpNum].owner = face.owner;
					mesh.faces[str+tmpNum].neighbour = face.neighbour;
					mesh.faces[str+tmpNum].setType(SEMO_Types::INTERNAL_FACE);
					
					
					++tmpNum;
				// cout << "1-4" << endl;
				}
				
				
			}
			else{
				mesh.faces[str] = mesh.faces[i];
			}
		
			// cout << "BBBBBB" << endl;
		}
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 10" << endl;
	
	

	// boundary setting
	for (int i=0; i<mesh.boundary.size(); ++i) {
		mesh.boundary[i].startFace = startFaces[ mesh.boundary[i].startFace ];
	}
	
	for (int i=0; i<mesh.boundary.size()-1; ++i) {
		// if(rank==0) cout << " AA : " << mesh.boundary[i+1].startFace << " " << mesh.boundary[i].startFace << endl;
		// mesh.boundary[i+1].startFace = max(mesh.boundary[i+1].startFace,mesh.boundary[i].startFace);
		// if(rank==0) cout << mesh.boundary[i].nFaces << " " << mesh.boundary[i+1].startFace-mesh.boundary[i].startFace << endl;
		mesh.boundary[i].nFaces = mesh.boundary[i+1].startFace-mesh.boundary[i].startFace;
	}
	int maxBDsize = mesh.boundary.size()-1;
	// if(rank==0) cout << mesh.boundary[maxBDsize].nFaces << " " << mesh.faces.size()-mesh.boundary[maxBDsize].startFace << endl;
	mesh.boundary[maxBDsize].nFaces = mesh.faces.size()-mesh.boundary[maxBDsize].startFace;
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 11" << endl;
	
	
	
	// cout << "8888 : " << rank << endl;
	
	
	
	
	// mesh.buildCells();
	
	int tmpCellNum = totalCellNum-1;
	mesh.cells.resize(totalCellNum,SEMO_Cell());
	for(int i=orgCellSize-1; i>=0; --i){
		int subCellSize = groupCell_Levels[i].size();
		for(int j=0; j<subCellSize; ++j){
			// int varSize = mesh.cells[i].var.size();
			// mesh.cells[tmpCellNum].var.resize(controls.nTotalCellVar,0.0);
			mesh.cells[tmpCellNum].var.assign(mesh.cells[i].var.begin(),mesh.cells[i].var.end());
			mesh.cells[tmpCellNum].points.clear();
			mesh.cells[tmpCellNum].faces.clear();
			--tmpCellNum;
		}
	}
	
	
	
	
	// mesh.setFaceTypes();
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "9 : " << rank << endl;
	
	mesh.buildLists();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "10 : " << rank << endl;
	
	mesh.connectCelltoFaces();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "11 : " << rank << endl;
	
	mesh.connectCelltoPoints();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "12 : " << rank << endl;
	
	
	mesh.setCountsProcFaces();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "13 : " << rank << endl;
	
	mesh.setDisplsProcFaces(); 
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "14 : " << rank << endl;
	
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// ++proc_num;
		// }
	// }
	// cout << "11111 : " << rank << " " << proc_num << " " << mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1] << endl;
	
	// if(rank==1){
		// for(int ip=0; ip<10000; ++ip){
			// for(int i=0; i<mesh.faces.size(); ++i){
				// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
					// for(auto& j : mesh.faces[i].points){
					// }
				// }
			// }
		// }
	// }
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// for(auto& j : mesh.faces[i].points){
				// cout << rank << " " << i << " " << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << " " << endl;
			// }
		// }
	// }
	

	// level setting
	tmpCellNum = 0;
	for(auto& groupCell : groupCell_Levels){
		for(auto& level : groupCell){
			mesh.cells[tmpCellNum].level = level;
			++tmpCellNum;
		}
	}
	tmpCellNum = 0;
	for(auto& groupCell : groupCell_Groups){
		for(auto& group : groupCell){
			mesh.cells[tmpCellNum].group = group;
			++tmpCellNum;
		}
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "15 : " << rank << endl;
	
	// cout << tmpCellNum << " " << mesh.cells.size() << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "15 : " << rank << endl;
	
	this->mpiLevels(mesh, cLevel_recv);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "16 : " << rank << " " << controls.nTotalFaceVar << " " << controls.nTotalFaceLRVar << endl;
	
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		face.var.clear();
		face.varL.clear();
		face.varR.clear();
		
		// face.var.resize(controls.nTotalFaceVar,0.0);
		// face.varL.resize(controls.nTotalFaceLRVar,0.0);
		// face.varR.resize(controls.nTotalFaceLRVar,0.0);
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "17 : " << rank << endl;
	
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		face.var.resize(controls.nTotalFaceVar,0.0);
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			if(face.owner > mesh.cells.size()-1){
				cout << " face.owner > mesh.cells.size()-1 " << face.owner << endl;
			}
			if(face.neighbour < 0 || face.neighbour > mesh.cells.size()-1){
				cout << " face.neighbour < 0 || face.neighbour > mesh.cells.size()-1 " << face.neighbour << endl;
			}
			
			int maxLevel = 
				max(mesh.cells[face.owner].level,
					mesh.cells[face.neighbour].level);
			face.level = maxLevel;
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			if(proc_num > cLevel_recv.size()-1){
				cout << " proc_num > cLevel_recv.size()-1 " << " " << cLevel_recv.size() << " " << proc_num << endl;
			}
			
			int maxLevel = 
				max(mesh.cells[face.owner].level,
					cLevel_recv[proc_num]);
			face.level = maxLevel;
			
			++proc_num;
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			if(face.owner > mesh.cells.size()-1){
				cout << " face.owner > mesh.cells.size()-1 " << face.owner << endl;
			}
			
			face.level = mesh.cells[face.owner].level;
		}
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "18 : " << rank << endl;
	
	
	
	afterCellSize = mesh.cells.size();
	afterFaceSize = mesh.faces.size();
	afterPointSize = mesh.points.size();
	
	int addedCellSize = afterCellSize - beforeCellSize;
	int addedFaceSize = afterFaceSize - beforeFaceSize;
	int addedPointSize = afterPointSize - beforePointSize;
	
	int addedCellSize_glb, addedFaceSize_glb, addedPointSize_glb;
	MPI_Allreduce(&addedCellSize, &addedCellSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&addedFaceSize, &addedFaceSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&addedPointSize, &addedPointSize_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	
	
	if(rank==0){
		
		cout << 
		" | refined added cell size = " << addedCellSize_glb <<
		" | refined added face size = " << addedFaceSize_glb <<
		" | refined added point size = " << addedPointSize_glb << endl;
		
		// cout << "| AMR - Refine completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
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





