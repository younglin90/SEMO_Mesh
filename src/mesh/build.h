#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

#include "load.h"
#include "save.h"
#include "../mpi/build.h"
// #include "geometric.h"

// class SEMO_Mesh_Load;

enum class SEMO_Types{
	INTERNAL_FACE,
	BOUNDARY_FACE,
	PROCESSOR_FACE,
	TO_BE_INTERNAL_FACE,
	TO_BE_PROCESSOR_FACE,
	TO_BE_DELETE_FACE
};


class SEMO_Point{
	public:
		int level;
		double x, y, z;
};

class SEMO_Face{
	public:
		SEMO_Face(){
			neighbour=-1;
		}
		
		void setType(SEMO_Types in){
			type = in;
		}
		SEMO_Types getType(){
			return type;
		}
		void setTypeBC(int in){
			typeBC = in;
		}
		int getTypeBC(){
			return typeBC;
		}
	
	public:
		int level;
		double unitNormals[3];
		double area;
		int owner, neighbour;
		vector<int> points;
		double x, y, z;
		
	private:
		SEMO_Types type;
		int typeBC;
};

class SEMO_Cell{
	public:
		int level;
		vector<int> points;
		vector<int> faces;
		double volume;
};

class SEMO_Boundary{
	public:
		SEMO_Boundary(){
			myProcNo=-1;
			neighbProcNo=-1;
		}
		
		string name;
		int nFaces;
		int startFace;
		int myProcNo;
		int neighbProcNo;
};




class SEMO_Mesh_Builder{
	public:
		SEMO_Mesh_Builder() {}
		SEMO_Mesh_Builder& addPoint(){
			SEMO_Point e;
			points.push_back(e);
			return *this;
		}
		SEMO_Mesh_Builder& addFace(){
			SEMO_Face e;
			faces.push_back(e);
			return *this;
		}
		SEMO_Mesh_Builder& addCell(){
			SEMO_Cell e;
			cells.push_back(e);
			return *this;
		}
		SEMO_Mesh_Builder& addBoundary(){
			SEMO_Boundary e;
			boundary.push_back(e);
			return *this;
		}
		
		void loadFile(string filetype);
		void saveFile(string filetype);
		void check();
		void buildCells();
		void buildCells2();
		void setFaceTypes();
		void buildLists();
		void checkLists();
		
		void connectCelltoFaces();
		void connectCelltoPoints();
		
		void distributeOneToAll(string type);
		void distributeEvenlyOneToAll();
		
		void partitionOneToAll();
		void partitionInit(string type);
		void partition(string type);
		void parMETIS(int nBlocks, int idBlockCell[]);
		void distribute(int nBlocks, int idBlockCell[]);
		
		// processor faces
		void setCountsProcFaces();
		void setDisplsProcFaces();
		
		// mesh datas
		vector<SEMO_Point> points;
		vector<SEMO_Face> faces;
		vector<SEMO_Cell> cells;
		vector<SEMO_Boundary> boundary;
		
		list<SEMO_Point*> listPoints;
		list<SEMO_Cell*> listCells;
		list<SEMO_Face*> listInternalFaces;
		list<SEMO_Face*> listBoundaryFaces;
		list<SEMO_Face*> listProcessorFaces;
		
		vector<int> countsProcFaces;
		vector<int> displsProcFaces;
		
		SEMO_MPI_Builder mpi;
		
		
	private:
		// SEMO_Mesh_Load meshLoad;
		
};

