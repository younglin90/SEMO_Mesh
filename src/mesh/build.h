#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
// #include <dlfnc.h>
using namespace std;

#include "../load/load.h"
#include "../save/save.h"
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
		SEMO_Point(){
			level=0;
		}
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
			level=0;
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
		int group;
		vector<double> unitNormals;
		double area;
		
		double wC;
		double wVC;
		double magPN;
		double alphaF;
		vector<double> vecPF;
		vector<double> vecNF;
		vector<double> unitNomalsPN;
		vector<double> vecSkewness;
		vector<double> vecPdP;
		vector<double> vecNdN;
		vector<double> vecPN;
		
		int owner, neighbour;
		vector<int> points;
		double x, y, z;
		vector<double> distCells;
		
		vector<double> var;
		vector<double> varL;
		vector<double> varR;
		
		vector<int> stencil;
		
		
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
		
	public:
		int level;
		int group;
		vector<int> points;
		vector<int> faces;
		
		double x, y, z;
		double volume;
		
		vector<double> var;
	
		vector<double> coeffLeastSquare1stFaceStencil;
		vector<double> coeffLeastSquare2ndFaceStencil;
		vector<double> coeffLeastSquare1stCellVetexStencil;
		vector<double> coeffLeastSquare2ndCellVertexStencil;
		
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
		
		// typedef void (*func_t)();
		// func_t setFunctionP;
		typedef void (*setFunc_t)(double, double, double, double, double&);
		vector<setFunc_t> setFunctionVariables;
		
		// void (*setFunctionPressure)(double, double, double, double, double&);
		// void (*setFunctionVelocities)(double, double, double, double, double&, double&, double&);
		// void (*setFunctionTemperature)(double, double, double, double, double&);
		// void (*setFunctionVolumeFractions)(double, double, double, double, vector<double>&);
		// // void (*setFunctionMassFractions)(double, double, double, double, vector<double>&);
		// // void (*setFunctionMassFractions)(double, double, double, double, vector<double>&);
		// vector<void> (*setFunctionTest)(double, double, double, double, double&);
		
		// func_t setFunctionP;
		// int *setFunctionP();
		// void *setFunctionP;
		// void setFunctionP(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionU(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionV(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionW(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionT(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionVF(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		// void setFunctionMF(SEMO_Mesh_Builder& mesh, SEMO_Face& face, int fPhi, int cPhi);
		
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
		void checkMatchingProcessorFace();
		void checkQualities();
		void buildCells();
		void buildCells2();
		void setFaceTypes();
		void buildLists();
		void checkLists();
		
		void cellsGlobal(); 
		
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
		
		
		void reorder();
		
		// 메쉬 에러
		void calcSkewness(vector<double>& output);
		void calcNonOrthogonality(vector<double>& output);
		void calcUniformity(vector<double>& output);
		
		
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
		
		// global
		int startCellGlobal;
		vector<int> neighbProcNo;
		vector<int> startProcCellGlobal;
		vector<int> procNeighbCellNo;
		int ncellsTotal;
		
		// CRS format
		int non_zeros;
		vector<int> CRS_ptr;
		vector<int> CRS_col;
		vector<int> CRS_col_ptr_dig;
		vector<int> CRS_col_ptr_LR;
		vector<int> CRS_col_ptr_RL;
		
		
		// 값들
		vector<vector<double>> cellsProcVar;
		vector<vector<vector<double>>> cellsProcGradientVar;
		vector<vector<vector<double>>> cellsGradientVar;
		
		// SEMO_MPI_Builder mpi;
		
		
		
		
	private:
		// SEMO_Mesh_Load meshLoad;
		
};

