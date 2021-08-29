#pragma once
#include <iostream>

#include "../controls/build.h" 


// #include "build.h"
class SEMO_Mesh_Builder;

class SEMO_Poly_AMR_Builder {
private:

public:
    SEMO_Poly_AMR_Builder() {
    }
    ~SEMO_Poly_AMR_Builder() {
    }
	
	void polyAMR(
		SEMO_Mesh_Builder& mesh, 
		SEMO_Controls_Builder& controls,
		vector<SEMO_Species>& species,
		int iter);
		
	void polyRefine(
		SEMO_Mesh_Builder& mesh, 
		SEMO_Controls_Builder& controls,
		int iter);
		
	void polyUnrefine(
		SEMO_Mesh_Builder& mesh, 
		SEMO_Controls_Builder& controls,
		int iter);
		
	// void sortCellCanUnrefine(
		// vector<int>& vLevel, 
		// int saveLevel, 
		// int& num, 
		// vector<vector<int>>& vGroupLevel, 
		// vector<vector<int>>& vGroupNumber);
		
	// void sortFaceCanUnrefine(
		// vector<int>& vLevel, 
		// int saveLevel, 
		// int& num, 
		// vector<vector<int>>& vGroupLevel, 
		// vector<vector<int>>& vGroupNumber);
		
	void mpiLevelRefine(
		SEMO_Mesh_Builder& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cLevel_recv, 
		vector<int>& cRefine_recv);
		
	void mpiRefines(
		SEMO_Mesh_Builder& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cRefine_recv);
		
	void mpiLevels(
		SEMO_Mesh_Builder& mesh, 
		vector<int>& cLevel_recv);
		
	void restrictCellRefine(
		SEMO_Mesh_Builder& mesh, 
		vector<bool>& boolCellRefine,
		vector<int>& cLevel_recv, 
		vector<int>& cRefine_recv);
		
	void createEdges(
		SEMO_Mesh_Builder& mesh, 
		vector<int>& edgesPoint0,
		vector<int>& edgesPoint1, 
		vector<vector<int>>& facesEdges,
		vector<vector<int>>& edgesFaces,
		vector<int>& edgeLevel);
		
	void searchOriginalPoints(
		SEMO_Mesh_Builder& mesh, 
		vector<int>& points,
		int targetLevel, 
		vector<int>& originPoints);
		

	void addCenterPoint(
		SEMO_Mesh_Builder& mesh, 
		vector<int> vertex, 
		int level
		){
		mesh.addPoint();
		mesh.points.back().x = 0.0;
		mesh.points.back().y = 0.0;
		mesh.points.back().z = 0.0;
		double n = (double)vertex.size();
		for(auto& i : vertex){
			mesh.points.back().x += mesh.points[i].x/n;
			mesh.points.back().y += mesh.points[i].y/n;
			mesh.points.back().z += mesh.points[i].z/n;
		}
		mesh.points.back().level = level + 1;
	}
		
		
};




class SEMO_Poly_Edge_Refine {
private:

public:
	int level;
	int point0;
	int point1;
	int centerPoint;
	vector<int> faces;
	
	
};




class SEMO_Poly_Mesh_Refine {
private:

public:
    SEMO_Poly_Mesh_Refine() {
    }
    ~SEMO_Poly_Mesh_Refine() {
    }
	
	SEMO_Poly_Mesh_Refine& addPoint(){
		SEMO_Point e;
		points.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addEdge(){
		SEMO_Poly_Edge_Refine e;
		edges.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addOutFace(int n){
		vector<SEMO_Face> e;
		for(int i=0; i<n; ++n){
			SEMO_Face a;
			e.push_back(a);
		}
		outFaces.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addIntFace(int n){
		vector<SEMO_Face> e;
		for(int i=0; i<n; ++n){
			SEMO_Face a;
			e.push_back(a);
		}
		intFaces.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addBdrFace(int n){
		vector<SEMO_Face> e;
		for(int i=0; i<n; ++n){
			SEMO_Face a;
			e.push_back(a);
		}
		bdrFaces.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addPrcFace(int n){
		vector<SEMO_Face> e;
		for(int i=0; i<n; ++n){
			SEMO_Face a;
			e.push_back(a);
		}
		prcFaces.push_back(e);
		return *this;
	}
	SEMO_Poly_Mesh_Refine& addCell(){
		SEMO_Cell e;
		cells.push_back(e);
		return *this;
	}
	
	void addCenterPoint(
		SEMO_Mesh_Builder& mesh, 
		vector<int> vertex, 
		int level
		){
		this->addPoint();
		this->points.back().x = 0.0;
		this->points.back().y = 0.0;
		this->points.back().z = 0.0;
		double n = (double)vertex.size();
		for(auto& i : vertex){
			this->points.back().x += mesh.points[i].x/n;
			this->points.back().y += mesh.points[i].y/n;
			this->points.back().z += mesh.points[i].z/n;
		}
		this->points.back().level = level + 1;
	}
	
	vector<SEMO_Point> points;
	vector<SEMO_Poly_Edge_Refine> edges;
	vector<SEMO_Cell> cells;
	
	vector<int> outFacesParentFace;
	vector<int> intFacesParentCell;
	vector<int> bdrFacesParentFace;
	vector<int> prcFacesParentFace;
	vector<vector<SEMO_Face>> outFaces;
	vector<vector<SEMO_Face>> intFaces;
	vector<vector<SEMO_Face>> bdrFaces;
	vector<vector<SEMO_Face>> prcFaces;
	
	
	void separateEdgesPoints(
		SEMO_Mesh_Builder& mesh, 
		int faceLevel,
		vector<int>& points, 
		vector<int>& edges,
		vector<bool>& boolEdgeRefine,
		vector<int>& edgesCenter,
		vector<int>& edgesLevel,
		vector<int>& faceVertex,
		vector<int>& edgeCenterPoints,
		vector<vector<int>>& subFaceEdgePoints,
		int& newPointNumber);
	
	
};





class SEMO_Poly_Mesh_Unrefine {
private:

public:
    SEMO_Poly_Mesh_Unrefine() {
    }
    ~SEMO_Poly_Mesh_Unrefine() {
    }
	
};
