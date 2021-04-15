// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <vector>
// #include <algorithm>
// #include <cmath>
// #include <list>


// using namespace std;

// // builder pattern

// // Part class
// class SEMO_Point{
	// public:
		// int level;
		// double x, y, z;
// };
// class SEMO_Face{
	// public:
		// int level;
		// double unitNormals[3];
		// double area;
		// int owner, neighbour;
		// vector<int> points;
// };
// class SEMO_Cell{
	// public:
		// int level;
		// vector<int> points;
		// vector<int> faces;
// };


// // Product class
// class SEMO_Mesh{
	// public:
		// // SEMO_Mesh() : mPoint(NULL), mFace(NULL), mCell(NULL) {}
		// // ~SEMO_Mesh() {
			// // if(mPoint) delete mPoint;
			// // if(mFace) delete mFace;
			// // if(mCell) delete mCell;
		// // }
		// SEMO_Mesh() {
			
		// }
		// ~SEMO_Mesh() {
			
		// }
		// void setPoint(SEMO_Point& point) { 
			// points.push_back(point); 
		// }
		// void setFace(SEMO_Face& face) { 
			// faces.push_back(face); 
		// }
		// void setCell(SEMO_Cell& cell) { 
			// cells.push_back(cell); 
		// }
		
	// private:
		// vector<SEMO_Point> points;
		// vector<SEMO_Face> faces;
		// vector<SEMO_Cell> cells;
		// // list<int> point;
		// // list<int> face;
		// // list<int> cell;
// };

// // Builder interface , Abstract builder
// class SEMO_Mesh_Builder{
	// public:
		// virtual void setPoint() = 0;
		// virtual void setFace() = 0;
		// virtual void setCell() = 0;
		// virtual SEMO_Mesh* getMesh() = 0;
// };

// // Concrete builder , inheritance class
// class SEMO_Mesh_Point_Builder : public SEMO_Mesh_Builder {
	// public:
		// SEMO_Mesh_Point_Builder(){
			// mesh = new SEMO_Mesh();
		// }
		// virtual void setPoint(SEMO_Point& point) {
			// mesh->setPoint(point);
		// }
		// virtual void setFace(SEMO_Face& face) {
			// mesh->setFace(face);
		// }
		// virtual void setCell(SEMO_Cell& cell) {
			// mesh->setCell(cell);
		// }
		// virtual SEMO_Mesh* getMesh() {
			// return mesh;
		// }
		
	// private:
		// SEMO_Mesh* mesh;
// };

// // Director class
// class SEMO_Mesh_Director {
	// public:
		// void setMeshBuilder(SEMO_Mesh_Builder* meshBuilder){
			// this->meshBuilder = meshBuilder;
		// }
		// SEMO_Mesh* makeMesh(){
			// meshBuilder->setPoint();
			// meshBuilder->setFace();
			// meshBuilder->setCell();
			
			// return meshBuilder->getMesh();
		// }
		
	// private:
		// SEMO_Mesh_Builder* meshBuilder;
// };