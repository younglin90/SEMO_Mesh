#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>

#include "parmetis.h" 

#include "scotch.h" 

using namespace std;
//앞에 있는 개행 문자 제거 
static inline std::string &ltrim(std::string &s) { 
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
	return s; 
};

//뒤에 있는 개행 문자 제거 
static inline std::string &rtrim(std::string &s) { 
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end()); 
	return s; 
};

//양쪽 끝의 개행 문자 제거 
static inline std::string &trim(std::string &s) { 
	return ltrim(rtrim(s)); 
};



enum class SEMO_Types{
	INTERNAL_FACE,
	BOUNDARY_FACE,
	PROCESSOR_FACE
};
// namespace StudentNames { 
// enum StudentNames { 
// KENNY, // 0 KYLE, // 1 STAN, // 2 BUTTERS, // 3 CARTMAN, // 4 WENDY, // 5 MAX_STUDENTS // 6 }; 
// }


// builder pattern

// Part class
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
		string name;
		int nFaces;
		int startFace;
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
		
	// private:
		vector<SEMO_Point> points;
		vector<SEMO_Face> faces;
		vector<SEMO_Cell> cells;
		vector<SEMO_Boundary> boundary;
		
};


class SEMO_Utility_Math{
	public:
		double *crossProduct(double v_A[], double v_B[]) {
			static double c_P[3];
			c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
			c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
			c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
			return c_P;
		}
};


class SEMO_Mesh_Geometric{
	public:
		void calcUnitNormals_Area3dPolygon(
		int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
		double unitNormals[], double& area ){
			double Nx, Ny, Nz;
			int a,b,c,s;
			
			b=n-2; c=n-1;
			Nx=0.0; Ny=0.0; Nz=0.0;
			for(int i = 0; i < n; ++i ) {
			  a = b; b = c; c = i;
		 
			  Nx += Vy[b] * ( Vz[c] - Vz[a] );
			  Ny += Vz[b] * ( Vx[c] - Vx[a] );
			  Nz += Vx[b] * ( Vy[c] - Vy[a] );
			}
			double length = sqrt(pow(Nx, 2) + pow(Ny, 2) + pow(Nz, 2));
			Nx /= length; Ny /= length; Nz /= length;
			
			area = 0.5*length;
			if(area < std::numeric_limits<double>::min()) {
				cerr << endl;
				cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
				cerr << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			unitNormals[0] = Nx; unitNormals[1] = Ny; unitNormals[2] = Nz;
			
		}
		
		
		void calcUnitNormals(
		double unitNormals[],
		double x1,double y1,double z1,
		double x2,double y2,double z2,
		double x3,double y3,double z3){
			
			double v1[3] = {x2-x1,y2-y1,z2-z1};
			double v2[3] = {x3-x1,y3-y1,z3-z1};
			
			unitNormals[0] = v1[1] * v2[2] - v1[2] * v2[1];
			unitNormals[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
			unitNormals[2] = v1[0] * v2[1] - v1[1] * v2[0];
			double length = sqrt(pow(unitNormals[0], 2) + pow(unitNormals[1], 2) + pow(unitNormals[2], 2));
			if(length < std::numeric_limits<double>::min()) {
				cerr << "from calcUnitNormals, length ~= 0.0" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
				

			unitNormals[0] = unitNormals[0]/length;
			unitNormals[1] = unitNormals[1]/length;
			unitNormals[2] = unitNormals[2]/length;
		}
		
		
		double calcArea3dPolygon( 
		int n, 
		vector<double> Vx, vector<double> Vy, vector<double> Vz,
		double Nx, double Ny, double Nz ) {
			double area = 0;
			double an, ax, ay, az; // abs value of normal and its coords
			int  coord;           // coord to ignore: 1=x, 2=y, 3=z
			int  i, j, k;         // loop indices
			
			printf("%d\n",n);

			for (auto item : Vx) {
				printf("%lf\n",item);
			}
			for (auto item : Vy) {
				printf("%lf\n",item);
			}
			for (auto item : Vz) {
				printf("%lf\n",item);
			}
			printf("%lf %lf %lf\n",Vx[0],Vy[0],Vz[0]);
			printf("%lf %lf %lf\n",Nx,Ny,Nz);

			if (n < 3) return 0;  // a degenerate polygon

			// select largest abs coordinate to ignore for projection
			ax = (Nx>0.0 ? Nx : -Nx);    // abs x-coord
			ay = (Ny>0.0 ? Ny : -Ny);    // abs y-coord
			az = (Nz>0.0 ? Nz : -Nz);    // abs z-coord

			coord = 3;                    // ignore z-coord
			if (ax > ay) {
				if (ax > az) coord = 1;   // ignore x-coord
			}
			else if (ay > az) coord = 2;  // ignore y-coord

			// compute area of the 2D projection
			switch (coord) {
			  case 1:
				for (i=1, j=2, k=0; i<n; i++, j++, k++)
					area += (Vy[i] * (Vz[j] - Vz[k]));
				break;
			  case 2:
				for (i=1, j=2, k=0; i<n; i++, j++, k++)
					area += (Vz[i] * (Vx[j] - Vx[k]));
				break;
			  case 3:
				for (i=1, j=2, k=0; i<n; i++, j++, k++)
					area += (Vx[i] * (Vy[j] - Vy[k]));
				break;
			}
			switch (coord) {    // wrap-around term
			  case 1:
				area += (Vy[n] * (Vz[1] - Vz[n-1]));
				break;
			  case 2:
				area += (Vz[n] * (Vx[1] - Vx[n-1]));
				break;
			  case 3:
				area += (Vx[n] * (Vy[1] - Vy[n-1]));
				break;
			}

			// scale to get area before projection
			an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
			switch (coord) {
			  case 1:
				area *= (an / (2.0 * Nx));
				break;
			  case 2:
				area *= (an / (2.0 * Ny));
				break;
			  case 3:
				area *= (an / (2.0 * Nz));
			}
			
			if(area < std::numeric_limits<double>::min()) {
				cerr << endl;
				cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
				cerr << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			return area;
		}
		
		double calcVolumePolyhedron(){
			
		}
		
};


  
int main(int argc, char* argv[]) {

	// printf("%d\n",argc);
	// printf("%s\n",argv[1]);
	
	// MPI initializing
	// int rank, size;
	// MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    // MPI_Comm_size(MPI_COMM_WORLD, &size); 
	
	MPI::Status status;
	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_Mesh_Builder mesh;
	
	
	
	string gridFolderName = "./grid";
	string pointsName = "points";
	string facesName = "faces";
	string ownerName = "owner";
	string neighbourName = "neighbour";
	string boundaryName = "boundary";
	
	bool boolPartitioning = true;
	bool boolPlotting = true;
	
	
	
	ifstream inputFile;
	string openFileName;
	
	// points
	openFileName = gridFolderName + "/" + pointsName;
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		return 1;
	}
	double xyz[3];
	string nextToken;
	bool startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				nextToken.erase(nextToken.find("("),1); 
				nextToken.erase(nextToken.find(")"),1); 
				stringstream sstream(nextToken);
				string word;
				char del = ' ';
				int num=0;
				while (getline(sstream, word, del)){
					xyz[num] = stod(word);
					++num;
				}
				mesh.addPoint();
				mesh.points.back().x = xyz[0];
				mesh.points.back().y = xyz[1];
				mesh.points.back().z = xyz[2];
				
			}
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}			
	}
	cout << "------------------------------------" << endl;
	cout << "point x,y,z size : " << mesh.points.size() << endl;
	inputFile.close();
	
	
	
	
	
	// faces
	openFileName = gridFolderName + "/" + facesName;
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		return EXIT_FAILURE;
	}
	startInput=false;
	bool continueInput=false;
	string saveToken;
	// saveToken.append(" ");
	// saveToken.append("asdfadsf ");
	// saveToken.append("asdfadsf ");
	// saveToken.append("asdfadsf ");
	// saveToken.append("asdfadsf ");
	// cout<<saveToken<< endl;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" && !continueInput ){
				break;
			}
			else{
				if(nextToken.size()==1) continue;
				
				if(nextToken.find(")") == string::npos){
					saveToken.append(" ");
					rtrim(nextToken);
					saveToken.append(nextToken);
					// cout<<saveToken<< endl;
					continueInput = true;
					continue;
				}
				saveToken.append(" ");
				rtrim(nextToken);
				saveToken.append(nextToken);
				
				saveToken.replace(saveToken.find("("),1," ");
				saveToken.replace(saveToken.find(")"),1," ");
				// cout << saveToken << endl;
				istringstream iss(saveToken);
				int tempint;
				iss >> tempint;
				
				
				mesh.addFace();
				
				while(iss >> tempint){
					mesh.faces.back().points.push_back(tempint);
				}
				saveToken.clear();
				continueInput = false;
			}
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}			
	}
	
	
	
	
	cout << "face size : " << mesh.faces.size() << endl;
	inputFile.close();
	
	
	
	// owner
	openFileName = gridFolderName + "/" + ownerName;
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		return EXIT_FAILURE;
	}
	int temp_num = 0;
	startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				
				
				while(iss >> tempint){
					mesh.faces[temp_num].owner = tempint;
					++temp_num;
				}
			}
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}			
	}
	cout << "owner size : " << temp_num << endl;
	inputFile.close();
	
	
	
	
	// neighbour
	openFileName = gridFolderName + "/" + neighbourName;
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		return EXIT_FAILURE;
	}
	temp_num = 0;
	startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					if(tempint<0) break;
					mesh.faces[temp_num].neighbour = tempint;
					
					++temp_num;
				}
			}
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}			
	}
	cout << "neighbour size : " << temp_num << endl;
	inputFile.close();
	
	
	
	// boundary
	openFileName = gridFolderName + "/" + boundaryName;
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		return EXIT_FAILURE;
	}
	vector<string> boundary_name;
	vector<char> boundary_type;
	vector<int> boundary_nFaces;
	vector<int> boundary_startFace;
	
	string backToken;
	startInput=false;
	vector<string> setToken;
	while(getline(inputFile, nextToken)){
		string asignToken;
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			setToken.push_back(nextToken.c_str());
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}
		backToken = nextToken;
	}
	
	string names;
	vector<string> setToken2;
	startInput=false;
	for (auto item : setToken) {
		string asignToken;
		if(startInput){
			if( item.find("}") != string::npos ){
				names.erase(std::remove(names.begin(), names.end(), ' '), names.end());
				boundary_name.push_back(names);
				for (auto item2 : setToken2) {
					if( item2.find("nFaces") != string::npos ){
						istringstream iss(item2);
						string temptemp;
						int temptempint;
						iss >> temptemp >> temptempint;
						boundary_nFaces.push_back(temptempint);
					}
					if( item2.find("startFace") != string::npos ){
						istringstream iss(item2);
						string temptemp;
						int temptempint;
						iss >> temptemp >> temptempint;
						boundary_startFace.push_back(temptempint);
						
					}
				}
				startInput=false;
				setToken2.clear();
			}
			setToken2.push_back(item.c_str());
		}
		else{ 
			if( item.find("{") != string::npos ){
				names = backToken;
				startInput=true;
			}
		}
		backToken = item;
    }
	
	inputFile.close();
	
	temp_num=0;
	for (auto item : boundary_name) {
		mesh.addBoundary();
		mesh.boundary.back().name = item;
		mesh.boundary.back().nFaces = boundary_nFaces[temp_num];
		mesh.boundary.back().startFace = boundary_startFace[temp_num];
		++temp_num;
	}
	
	
	
	//=================== check loading file ================
	
	int owner_max=-1;
	int neighbour_max=-1;
	int checkOwnerNeighbourReverse=0;
	for(auto& face : mesh.faces){
		owner_max = max(owner_max , face.owner);
		neighbour_max = max(neighbour_max , face.neighbour);
	}
	
	// if(owner_max < neighbour_max){
		// cerr << "#error, owner_max(" << owner_max << ") < neighbour_max(" << neighbour_max << ")" << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	
	for(auto& face : mesh.faces){
		if(face.neighbour > owner_max){
			// cerr << "#error face, owner = " << face.owner << " neighbour = " << face.neighbour << endl;
			int tempnum = face.owner;
			face.owner = face.neighbour;
			face.neighbour = tempnum;
			std::reverse(face.points.begin(),face.points.end());
			++checkOwnerNeighbourReverse;
		}
	}
	cout << "  #check face owner < neighbour, executed reverse : " << checkOwnerNeighbourReverse << endl;
	
	
	
	
	
	
	

	
	//=================== Loading End ==========================
	
	
	
	// add Cells
	int cell_num=0;
	for(auto& face : mesh.faces){
		cell_num = max(cell_num , face.owner);
		cell_num = max(cell_num , face.neighbour);
	}
	for(int i=0; i<cell_num+1; ++i){
		mesh.addCell();
	}
	
	
	// set types
	// for (auto item : mesh.faces) {
	for(auto& face : mesh.faces){
		if(face.neighbour == -1){
			face.setType(SEMO_Types::BOUNDARY_FACE);
		}
		else{
			face.setType(SEMO_Types::INTERNAL_FACE);
		}
	}
	
	for(int ibc=0; ibc<mesh.boundary.size(); ++ibc){
		for(int i=mesh.boundary[ibc].startFace; i<mesh.boundary[ibc].startFace + mesh.boundary[ibc].nFaces; ++i){
			mesh.faces[i].setTypeBC(ibc);
		}
	}


	// // owner neighbour order reverse
	// for(auto& i : mesh.faces){
		
		// if(i.getType() == SEMO_Types::BOUNDARY_FACE){
			// if(i.owner
		// }
	// } 

	
	
	// create list
	list<SEMO_Point*> listPoints;
	list<SEMO_Cell*> listCells;
	list<SEMO_Face*> listInternalFaces;
	list<SEMO_Face*> listBoundaryFaces;
	list<SEMO_Face*> listConnectedFaces;

	// listCells.push_back(&mesh.cells[0]);
	for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
		listPoints.push_back(&*iter);
	}
	
	for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		listCells.push_back(&*iter);
	}
	
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		if(iter->getType() == SEMO_Types::INTERNAL_FACE){
			listInternalFaces.push_back(&*iter);
		}
		else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			listBoundaryFaces.push_back(&*iter);
		}
		else if(iter->getType() == SEMO_Types::PROCESSOR_FACE){
			listConnectedFaces.push_back(&*iter);
		}
		else{
			cout << "TO DO : other face types" << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	}
	
	
	cout << "------------------------------------" << endl;
	cout << "point size : " << listPoints.size() << endl;
	cout << "cell size : " << listCells.size() << endl;
	cout << "internal face size : " << listInternalFaces.size() << endl;
	cout << "boundary face size : " << listBoundaryFaces.size() << endl;
	cout << "connected face size : " << listConnectedFaces.size() << endl;
	  


	
	// cell connection (cell's face)
	cout << "------------------------------------" << endl;
	cout << "| execute cell's face connecting ... ";
	int faceNum=0;
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		if(iter->getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[iter->owner].faces.push_back(faceNum);
			mesh.cells[iter->neighbour].faces.push_back(faceNum);
		}
		else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			mesh.cells[iter->owner].faces.push_back(faceNum);
		}
		++faceNum;
	}
	cout << "-> completed" << endl;
	
	// cell connection (cell's points)
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
	cout << "| execute cell's points connecting ... ";
	for(auto& iter : mesh.faces){
		if(iter.getType() == SEMO_Types::INTERNAL_FACE){
			for(int i=0; i<iter.points.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[iter.owner].points.size(); ++j){
					if(
					mesh.cells[iter.owner].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[iter.owner].points.push_back(iter.points[i]);
				
				l=0;
				for(int j=0; j<mesh.cells[iter.neighbour].points.size(); ++j){
					if(
					mesh.cells[iter.neighbour].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[iter.neighbour].points.push_back(iter.points[i]);
			}
			
		}
		else if(iter.getType() == SEMO_Types::BOUNDARY_FACE){
			for(int i=0; i<iter.points.size(); ++i){
				int l=0;
				for(int j=0; j<mesh.cells[iter.owner].points.size(); ++j){
					if(
					mesh.cells[iter.owner].points[j] == 
					iter.points[i] ) ++l;
				}
					
				if(l==0) mesh.cells[iter.owner].points.push_back(iter.points[i]);
			}
		}
	}
	cout << "-> completed" << endl;
	
	
	
	
	// // partitioning




	if(boolPartitioning){




	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	
	int ncommon=3;
	
	
	int nBlocks = 2;
	
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_DBGLVL]=1;
	
	int idBlockCell[ncells];
	int dummy1[npoints];






	// int xadj[ncells+1];
	// xadj[0] = 0;
	// int numt = 0;
	// for(auto& cell : mesh.cells){
		// ++numt;
		// xadj[numt] = xadj[numt-1] + cell.points.size();
	// }
	
	// int adjncy[ xadj[ncells] ];
	// numt = 0;
	// for(auto& cell : mesh.cells){
		// for(auto& i : cell.points){
			// adjncy[numt] = i;
			// ++numt;
		// }
	// }

	// METIS_PartMeshDual(
		// &ncells, &npoints, xadj, adjncy, NULL, NULL,  
		// &ncommon, &nBlocks, NULL, options, &objval, idBlockCell, dummy1); 
	
	
	
	
	

	int xadj[ncells+1];
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		for(auto& face : cell.faces){
			if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE)
				++numt2;
		}
		++numt;
		xadj[numt] = xadj[numt-1] + numt2;
	}
	
	int adjncy[ xadj[ncells] ];
	numt = 0;
	for(auto& cell : mesh.cells){
		for(auto& face : cell.faces){
			if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE) {
				if( &cell == &mesh.cells[ mesh.faces[face].owner ] ){
					adjncy[numt] = mesh.faces[face].neighbour;
				}
				else if( &cell == &mesh.cells[ mesh.faces[face].neighbour ] ){
					adjncy[numt] = mesh.faces[face].owner;
				}
				else {
					cerr << "#error, not matching cell != face owner or neighbour" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				++numt;
			}
		}
	}
	
	// METIS_PartGraphKway(
		// &ncells, &ncon, xadj, adjncy, 
		// NULL, NULL, NULL, &nBlocks, NULL, 
		// NULL, options, &objval, idBlockCell); 
	
	

	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	// tpwgts[0]=0.4;
	// tpwgts[1]=0.9;
	
	real_t ubvec[ncon];
	std::fill_n(ubvec, ncon, 1.05);
	
	int vtxdist[nBlocks+1];
	
	vtxdist[0] = 0;
	vtxdist[1] = ncells;


	ParMETIS_V3_PartKway(
		vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
		&ncon, &nBlocks, tpwgts, ubvec,
		options, &objval, idBlockCell, &comm);


	// //===================================
	
	
	
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0));
	int nCellsLocal[nBlocks];
	int idCellLocal[mesh.cells.size()];
	std::fill_n(nCellsLocal, nBlocks, 0);
	for(int i=0; i<mesh.cells.size(); ++i) {
		for(auto& j : mesh.cells[i].points){
			int l=0;
			for(int k=0; k<idBlockPoint[j].size(); ++k){
				if(idBlockPoint[j][k] == idBlockCell[i]) ++l;
			}
			if(l==0) idBlockPoint[j].push_back( idBlockCell[i] );
		}
		idCellLocal[i] = nCellsLocal[ idBlockCell[i] ];
		++nCellsLocal[ idBlockCell[i] ];
	}
	
	int nTotalLocalPointSize = 0;
	for(auto& i : idBlockPoint) 
		nTotalLocalPointSize += i.size();

	int nPointsLocal[nBlocks];
	std::fill_n(nPointsLocal, nBlocks, 0);
	int idPointLocal[nTotalLocalPointSize];
	int idBlockPointLocal[nTotalLocalPointSize];
	int strPoints[mesh.points.size()+1];
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(int j=0; j<idBlockPoint[i].size(); ++j){
			int idBlock = idBlockPoint[i][j];
			idPointLocal[nIndex] = nPointsLocal[ idBlock ];
			++nPointsLocal[ idBlock ];
			// cout << idPointLocal[nIndex] << endl;
			idBlockPointLocal[nIndex] = idBlock;
			++nIndex;
		}
	}
	strPoints[mesh.points.size()] = nIndex;
	
	
	std::fill_n(nPointsLocal, nBlocks, 0);
	for(int i=0; i<nTotalLocalPointSize; ++i)
		++nPointsLocal[ idBlockPointLocal[i] ];
	
	// open & write local npoints
	ofstream decomPointsName[nBlocks];
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/points." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nPointsLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// wirte local point's x y z
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			decomPointsName[idBlock] << "(" << mesh.points[i].x;
			decomPointsName[idBlock] << " " << mesh.points[i].y;
			decomPointsName[idBlock] << " " << mesh.points[i].z << ")" << endl;
		}
	}
	
	// close points files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	// face setting
	int nFacesLocal[nBlocks];
	int nProcsFace = 0;
	std::fill_n(nFacesLocal, nBlocks, 0);
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal[idBlockOwner];
		
		int idBlockNeighbour;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				face.setType(SEMO_Types::PROCESSOR_FACE);
				++nFacesLocal[idBlockNeighbour];
				++nProcsFace;
			}
		}
	}
	
	// write # of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/faces." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	
	int idFaceLocal[nBlocks];
	std::fill_n(idFaceLocal, nBlocks, 0);
	
	// internal faces
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlock = idBlockCell[face.owner];
			++idFaceLocal[idBlock];
			
			vector<int> idFacePoint;
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
		
			decomPointsName[idBlock] << face.points.size() << "(";
			decomPointsName[idBlock] << idFacePoint[0];
			for(int i=1; i<idFacePoint.size(); ++i){
				decomPointsName[idBlock] << " " << idFacePoint[i];
			}
			decomPointsName[idBlock] << ")" << endl;
		}
	}
	
	// boundary faces
	int nFacesBoundaryLocal[nBlocks];
	int nFacesEachBoundaryLocal[nBlocks][mesh.boundary.size()];
	int nStartFaceEachBoundaryLocal[nBlocks][mesh.boundary.size()];
	for(int i=0; i<nBlocks; ++i){
		nFacesBoundaryLocal[i] = 0;
		for(int j=0; j<mesh.boundary.size(); ++j){
			nFacesEachBoundaryLocal[i][j]=0;
			nStartFaceEachBoundaryLocal[i][j]=0;
		}
	}
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int idBlock = idBlockCell[face.owner];
			++nFacesBoundaryLocal[ idBlock ];
			
			// mesh.getBoundaryNumber();
			++nFacesEachBoundaryLocal[ idBlock ][ face.getTypeBC() ];
			if(nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] == 0)
				nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = idFaceLocal[ idBlock ];
			++idFaceLocal[ idBlock ];
			
			// wirte bc
			vector<int> idFacePoint;
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
		
			decomPointsName[idBlock] << face.points.size() << "(";
			decomPointsName[idBlock] << idFacePoint[0];
			for(int i=1; i<idFacePoint.size(); ++i){
				decomPointsName[idBlock] << " " << idFacePoint[i];
			}
			decomPointsName[idBlock] << ")" << endl;
		}
		
	}
	
	
	
	// processor faces
	int nFacesProcessorLocal[nBlocks];
	int nStartProcessorLocal[nBlocks][nBlocks];
	int sendCounts[nBlocks][nBlocks];
	for(int i=0; i<nBlocks; ++i){
		nFacesProcessorLocal[i] = 0;
		for(int j=0; j<nBlocks; ++j){
			nStartProcessorLocal[i][j]=0;
			sendCounts[i][j]=0;
		}
	}
	vector<int> idFacesProcessor;
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			// cout << mesh.faces[neighbour << endl;
			// cout << face.neighbour << endl;
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			// if(nStartProcessorLocal[ idBlockOwner ][ idBlockNeighbour ] == 0)
				// nStartProcessorLocal[ idBlockOwner ][ idBlockNeighbour ] = idFaceLocal[ idBlockOwner ];
			// if(nStartProcessorLocal[ idBlockNeighbour ][ idBlockOwner ] == 0)
				// nStartProcessorLocal[ idBlockNeighbour ][ idBlockOwner ] = idFaceLocal[ idBlockNeighbour ];
			++sendCounts[idBlockOwner][idBlockNeighbour];
			++sendCounts[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			int k = idFacesProcessor[jp];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			// cout<<k<<endl;
			if(n==ip) {
				int idBlock = m;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].points){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				decomPointsName[idBlock] << mesh.faces[k].points.size() << "(";
				decomPointsName[idBlock] << idFacePoint[0];
				for(int i=1; i<idFacePoint.size(); ++i){
					decomPointsName[idBlock] << " " << idFacePoint[i];
				}
				decomPointsName[idBlock] << ")" << endl;
			}
			else if(m==ip) {
				int idBlock = n;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].points){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				decomPointsName[idBlock] << mesh.faces[k].points.size() << "(";
				decomPointsName[idBlock] << idFacePoint[0];
				for(int i=1; i<idFacePoint.size(); ++i){
					decomPointsName[idBlock] << " " << idFacePoint[i];
				}
				decomPointsName[idBlock] << ")" << endl;
			}
		}
	}
	
	// close faces files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	// write # of owner, same as number of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/owner." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	for(auto& face : mesh.faces){
		if(
		face.getType() == SEMO_Types::INTERNAL_FACE ||
		face.getType() == SEMO_Types::BOUNDARY_FACE
		){
			int j = face.owner;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	for(int i=0; i<nBlocks; ++i){
		for(int j=0; j<idFacesProcessor.size(); ++j){
			int k = idFacesProcessor[j];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			
			if(n==i) {
				decomPointsName[m] << idCellLocal[ mesh.faces[k].owner ] << endl;
			}
			else if(m==i) {
				decomPointsName[n] << idCellLocal[ mesh.faces[k].neighbour ] << endl;
			}
		}
	}
	
	
	// close owner files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	// write # of neighbour
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/neighbour." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] - nFacesBoundaryLocal[i] - nFacesProcessorLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	for(auto& face : mesh.faces){
		if(
		face.getType() == SEMO_Types::INTERNAL_FACE
		){
			int j = face.neighbour;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// close neighbour files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	// write # of boundary
	for(int i=0; i<nBlocks; ++i){
		int n=0;
		for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
			if( nStartFaceEachBoundaryLocal[i][ibcs] > 0) ++n;
		}
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts[i][j] > 0) ++n;
		}
		
		string sFilename = "./grid/boundary." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << n << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	for(int i=0; i<nBlocks; ++i){
		for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
			if( nStartFaceEachBoundaryLocal[i][ibcs] > 0) {
				decomPointsName[i] << "   " << mesh.boundary[ibcs].name << endl;
				decomPointsName[i] << "   {" << endl;
				decomPointsName[i] << "      type Unspecified;" << endl;
				decomPointsName[i] << "      nFaces " << nFacesEachBoundaryLocal[i][ibcs] << endl;
				decomPointsName[i] << "      startFace " << nStartFaceEachBoundaryLocal[i][ibcs] << endl;
				decomPointsName[i] << "   }" << endl;
				decomPointsName[i] << endl;
				
			}
		}
	}
	
	
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts[i][j];
		}
		// ++n;
		// cout << nFacesLocal[i] << endl;
		
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts[i][j] > 0) {
				string bcnames = "   procBoundary" + to_string(i) + "to" + to_string(j);
				decomPointsName[i] << bcnames << endl;
				decomPointsName[i] << "   {" << endl;
				decomPointsName[i] << "      type processor;" << endl;
				decomPointsName[i] << "      nFaces " << sendCounts[i][j] << endl;
				decomPointsName[i] << "      startFace " << n << endl;
				decomPointsName[i] << "      myProcNo " << i << endl;
				decomPointsName[i] << "      neighbProcNo " << j << endl;
				decomPointsName[i] << "   }" << endl;
				decomPointsName[i] << endl;
				
				n += sendCounts[i][j];
			}
		}
	}
	
	
	// close boundary files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	}
	
	
	
	// // geometric
	// cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	// SEMO_Mesh_Geometric geometric;

	
	// // polygon face normal vectors & polygon face area
	// // polygon face center x,y,z
	// // 3D Polygon Areas : https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		// vector<double> Vx, Vy, Vz;
		// for(auto iPoint : iter->points){
			// Vx.push_back(mesh.points[iPoint].x);
			// Vy.push_back(mesh.points[iPoint].y);
			// Vz.push_back(mesh.points[iPoint].z);
		// }
		
		// geometric.calcUnitNormals_Area3dPolygon(
		// iter->points.size(), Vx,Vy,Vz,
		// iter->unitNormals, iter->area );
		
		// iter->x = accumulate(Vx.begin(), Vx.end(), 0.0) / iter->points.size();
		// iter->y = accumulate(Vy.begin(), Vy.end(), 0.0) / iter->points.size();
		// iter->z = accumulate(Vz.begin(), Vz.end(), 0.0) / iter->points.size();
	// }
	
	
	// // polyhedron cell volume (Green-Gauss Theorem.)
	// for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		// iter->volume = 0.0;
	// }
	
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		
		// // cout << iter->owner << endl;
		
		// if(iter->getType() == SEMO_Types::INTERNAL_FACE){
			// mesh.cells[iter->owner].volume += 
				// (iter->x*iter->unitNormals[0]+
				 // iter->y*iter->unitNormals[1]+
				 // iter->z*iter->unitNormals[2])*iter->area;
			
			// mesh.cells[iter->neighbour].volume -= 
				// (iter->x*iter->unitNormals[0]+
				 // iter->y*iter->unitNormals[1]+
				 // iter->z*iter->unitNormals[2])*iter->area;
		// }
		// else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			// mesh.cells[iter->owner].volume += 
				// (iter->x*iter->unitNormals[0]+
				 // iter->y*iter->unitNormals[1]+
				 // iter->z*iter->unitNormals[2])*iter->area;
		// }
		
	// }
	// for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		// iter->volume /= 3.0;
		
		// if(iter->volume < std::numeric_limits<double>::min()) {
			// cerr << endl;
			// cerr << "#error, from calc cell volume, cell volume = " << iter->volume << " < cpu_min_val " << endl;
			// cerr << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		// // cout << iter->volume << endl;
	// }
	// cout << "-> completed" << endl;
	
	
	
	
	// plotting
	
	
	if(boolPlotting){
	
	cout << "| execute save file (.vtu format) ... ";
	ofstream outputFile;
	outputFile.open("./save/plot.vtu");
	// outputFile.open("./save/plot.vtu",ios::binary);
	// outputFile.open("./save/plot.vtu",ios::out | ios::binary);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		return 1;
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << listPoints.size() << "\" NumberOfCells=\"" << listCells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	outputFile << "    </PointData>" << endl;
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float32\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;
	// for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	for(auto iter=listPoints.begin(); iter!=listPoints.end(); iter++){
		outputFile << scientific << (*iter)->x << " " << (*iter)->y << " " << (*iter)->z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		for(auto i : (*iter)->points){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		cellFaceOffset += (*iter)->points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		outputFile << (*iter)->faces.size() << endl;
		for(auto& i : (*iter)->faces){
			outputFile << mesh.faces[i].points.size() << " ";
			for(auto& j : mesh.faces[i].points){
				outputFile << j << " ";
			}
			outputFile << endl;
		}
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		int numbering = 1 + (*iter)->faces.size();
		for(auto& i : (*iter)->faces){
			numbering += mesh.faces[i].points.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	cout << "-> completed" << endl;
	cout << "------------------------------------" << endl;
	

	}
	

	return EXIT_SUCCESS;
}














