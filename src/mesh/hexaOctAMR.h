#pragma once
#include <iostream>

// #include "build.h"
class SEMO_Mesh_Builder;

class SEMO_Hexa_Oct_AMR {
private:

public:
    SEMO_Hexa_Oct_AMR() {
    }
    ~SEMO_Hexa_Oct_AMR() {
    }
	
	void searchLoacationsCase0(
		vector<int> &faces, vector<int> &points, vector<int> &cellPoints);
		
	void sortCellPointsSet(
		vector<int> &cellPoints, vector<vector<int>> &cellPointsSet);

	void sortFacesPointsSet(
		vector<int> &cellPoints, int& idPoint, vector<vector<int>> &newFacesPointsSet);
		

	void sortOwnNgbOuterFaces(
		int& num, int& strOwnNgb, vector<vector<int>> &newFaceOwnNgbSet);

	void sortOwnNgbInternalFaces(
		int& num, int& strOwnNgb, 
		vector<vector<int>> &newFaceOwnSet, vector<vector<int>> &newFaceNgbSet);
	

	void addPoints( 
		SEMO_Mesh_Builder& mesh, vector<SEMO_Point>& newPoints);

	void addFaces( 
		SEMO_Mesh_Builder& mesh, vector<int>& idNewFaces, vector<SEMO_Face>& newFaces);
		

	void addCells( 
		SEMO_Mesh_Builder& mesh, vector<int>& idNewCells, vector<int>& levelNewCells);
		

	void searchParallelPointsPolygon( 
		SEMO_Mesh_Builder& mesh, vector<int>& points, vector<vector<int>>& parPoints);
		
};
