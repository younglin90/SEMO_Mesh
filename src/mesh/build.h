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
#include "../math/math.h"

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
		
		vector<double> var;
		
		vector<int> stencil;
};

class SEMO_Face{
	public:
		SEMO_Face(){
			unitNormals.resize(3,0.0);
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
		
		// double var(int in){
			// return variables[in];
		// }
	
	public:
		int level;
		vector<double> unitNormals;
		double area;
		double wC;
		double wVC;
		int owner, neighbour;
		vector<int> points;
		double x, y, z;
		vector<double> distCells;
		
		vector<double> var;
		vector<double> varL;
		vector<double> varR;
		
		
	private:
		SEMO_Types type;
		int typeBC;
};

class SEMO_Cell{
	public:
		SEMO_Cell(){
			level=0;
		}
		
		// double var(int in){
			// return variables[in];
		// }
		
		int level;
		vector<int> points;
		vector<int> faces;
		
		double x, y, z;
		double volume;
		
		vector<double> var;
	
		vector<double> coeffLeastSquare;
		
		vector<int> stencil;
		
		int RCM;
		int invRCM;
		
	private:
		
};

class SEMO_Boundary{
	public:
		SEMO_Boundary(){
			myProcNo=-1;
			neighbProcNo=-1;
		}
		
		string name;
		vector<string> type;
		vector<double> var;
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
		
		void loadFile(string filetype, string folder);
		void saveFile(string filetype, string folder, SEMO_Controls_Builder &controls);
		void check();
		void checkQualities();
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
		
		void informations();
		
		// AMR
		void hexaOctAMR();
		
		// processor faces
		void searchNeighbProcFaces();
		void setCountsProcFaces();
		void setDisplsProcFaces();
		
		// mesh datas
		vector<SEMO_Point> points;
		vector<SEMO_Face> faces;
		vector<SEMO_Cell> cells;
		vector<SEMO_Boundary> boundary;
		
		list<SEMO_Point*> listPoints;
		list<SEMO_Cell*> listCells;
		list<SEMO_Face*> listFaces;
		list<SEMO_Face*> listInternalFaces;
		list<SEMO_Face*> listBoundaryFaces;
		list<SEMO_Face*> listProcessorFaces;
		
		vector<int> countsProcFaces;
		vector<int> displsProcFaces;
		
		// SEMO_MPI_Builder mpi;
		
		
	private:
		// SEMO_Mesh_Load meshLoad;
		
};

