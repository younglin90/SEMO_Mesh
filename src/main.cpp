#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>

using namespace std;




enum class SEMO_Types{
	INNER_FACE,
	BOUNDARY_FACE,
	CONNECTION_FACE
};





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
		// void setType(int in){
			// type = in;
		// }
		// int getType(){
			// return type;
		// }
	
	public:
		int level;
		double unitNormals[3];
		double area;
		int owner, neighbour;
		vector<int> points;
		double x, y, z;
		
	private:
		SEMO_Types type;
		// int type;
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
				cerr << "from calcArea3dPolygon, area ~= 0.0" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			return area;
		}
		
		double calcVolumePolyhedron(){
			
		}
		
};



// class SEMO_Mesh_Builder;

// // Product class
// class SEMO_Mesh_Product{
	// public:
	// // private:
		// vector<SEMO_Point> points;
		// vector<SEMO_Face> faces;
		// vector<SEMO_Cell> cells;
		
		// std::string street_address, post_code, city;
		// // list<int> point;
		// // list<int> face;
		// // list<int> cell;
// };

// // Builder interface , Abstract builder
// class SEMO_Mesh_Point_Builder;

// class SEMO_Mesh_Builder_Base{
	// public:
		// operator SEMO_Mesh_Product(){
			// return std::move(product);
		// }
		// SEMO_Mesh_Point_Builder lives() const;
		// // SEMO_Mesh_Point_Builder points() const;
		
		
	// protected:
		// SEMO_Mesh_Product& product;
		// explicit SEMO_Mesh_Builder_Base(SEMO_Mesh_Product& product)
		// : product{ product } {}
		
	// private:
		// SEMO_Point point;
		// SEMO_Face face;
		// SEMO_Cell cell;
// };

// class SEMO_Mesh_Builder : public SEMO_Mesh_Builder_Base{
	// public:
		// SEMO_Mesh_Product p;
		// SEMO_Mesh_Builder() : SEMO_Mesh_Builder_Base{p} {}
	
	// private:
	
// };

// class SEMO_Mesh_Point_Builder : public SEMO_Mesh_Builder_Base{
	// public:
		// typedef SEMO_Mesh_Point_Builder self;
		
		// explicit SEMO_Mesh_Point_Builder(SEMO_Mesh_Product& product)
		// : SEMO_Mesh_Builder_Base{ product } {}
		
		// self& at(std::string street_address){
			// product.street_address = street_address;
			// return *this;
		// }
		
		// self& with_postcode(std::string post_code){}
		// self& in(std::string city){}
		
	// private:
// };






// class PersonBuilder;

// class Person{
	// public:
		// std::string street_address, post_code, city;
		// std::string company_name, position;
		// int annual_income = 0;
		// PersonBuilder create();
		
		// Person(){}
		
// };

// // Builder interface , Abstract builder
// class PersonAddressBuilder;

// class PersonBuilderBase{
	// protected:
		// Person& person;
		// explicit PersonBuilderBase(Person& person)
		// : person{ person } {}
	// public:
		// operator Person(){
			// return std::move(person);
		// }
		// PersonAddressBuilder lives() const;
		// // PersonJobBuilder works() const; 
// };

// class PersonBuilder : public PersonBuilderBase{
	// public:
		// Person p;
		// PersonBuilder() : PersonBuilderBase{p} {}
	
	// private:
	
// };

// class PersonAddressBuilder : public PersonBuilderBase{
		// typedef PersonAddressBuilder self;
		
	// public:
		// explicit PersonAddressBuilder(Person& person)
		// : PersonBuilderBase{ person } {}
		
		// self& at(std::string street_address){
			// person.street_address = street_address;
			// return *this;
		// }
		
		// self& with_postcode(std::string post_code){}
		// self& in(std::string city){}
		
	// private:
// };


// PersonBuilder Person::create() { 
	// return PersonBuilder{}; 
// }

// PersonAddressBuilder PersonBuilderBase::lives() const { 
	// return PersonAddressBuilder{ person }; 
// };

// // PersonJobBuilder PersonBuilderBase::works() const { 
	// // return PersonJobBuilder{ person }; 
// // };










// SEMO_Mesh_Product* SEMO_Mesh_Builder::Build() { 
	// return new SEMO_Mesh_Product(this); 
// }


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








  
int main(int argc, char* argv[]) {

	// printf("%d\n",argc);
	// printf("%s\n",argv[1]);
	
	
	SEMO_Mesh_Builder mesh;
	
	

	ifstream inputFile("./grid/points");
	if(inputFile.fail()){
		cerr << "Unable to open file for reading." << endl;
		return 1;
	}
	
	// mesh.addPoint();
	// vector<double> xyz;
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
				// char del = '(';
				// nextToken.erase(0,1); 
				nextToken.erase(nextToken.find("("),1); 
				nextToken.erase(nextToken.find(")"),1); 
				// nextToken.erase(nextToken.length()-2,1); 
				stringstream sstream(nextToken);
				string word;
				char del = ' ';
				int num=0;
				while (getline(sstream, word, del)){
					// cout << stod(word) << endl;
					// xyz.push_back(stod(word));
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
	
	// for (auto item : xyz) {
        // cout << item << endl;
    // }
	cout << "point x,y,z size : " << mesh.points.size() << endl;
	inputFile.close();
	
	
	
	
	
	
	inputFile.open("./grid/faces");
	if(inputFile.fail()){
		cerr << "Unable to open file for reading." << endl;
		return EXIT_FAILURE;
	}
	// vector<vector<int>> face_node;
	
	// string nextToken;
	startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				nextToken.replace(nextToken.find("("),1," ");
				nextToken.replace(nextToken.find(")"),1," ");
				istringstream iss(nextToken);
				int tempint;
				iss >> tempint;
				// cout << tempint << " ";
				// vector<int> face_node2;
				
				mesh.addFace();
				
				while(iss >> tempint){
					mesh.faces.back().points.push_back(tempint);
					// face_node2.push_back(tempint);
					// cout << tempint << " ";	
				}
				// face_node.push_back(face_node2);
				// cout << endl;
			}
		}
		else{
			if( asignToken.assign(nextToken, 0, 1) == "(" ){
				startInput=true;
			}
		}			
	}
	cout << "face size : " << mesh.faces.size() << endl;
	// cout << "face size : " << face_node[1].size() << endl;
	// for (auto i : face_node) {
		// for(auto j : i) {
			// cout << j << endl;
		// }
	// }
	inputFile.close();
	
	
	
	
	
	
	
	inputFile.open("./grid/owner");
	if(inputFile.fail()){
		cerr << "Unable to open file for reading." << endl;
		return EXIT_FAILURE;
	}
	// vector<int> face_owner;
	
	// string nextToken;
	int temp_num = 0;
	startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				// nextToken.replace(nextToken.find("("),1," ");
				// nextToken.replace(nextToken.find(")"),1," ");
				istringstream iss(nextToken);
				int tempint;
				// iss >> tempint;
				// cout << tempint << " ";	
				while(iss >> tempint){
					// face_owner.push_back(tempint);
					mesh.faces[temp_num].owner = tempint;
					// cout << tempint << " ";	
					++temp_num;
				}
				// cout << endl;
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
	
	
	
	
	
	
	
	inputFile.open("./grid/neighbour");
	if(inputFile.fail()){
		cerr << "Unable to open file for reading." << endl;
		return EXIT_FAILURE;
	}
	// vector<int> face_neighbour;
	
	// string nextToken;
	temp_num = 0;
	startInput=false;
	while(getline(inputFile, nextToken)){
		string asignToken;
		
		if(startInput){
			if( asignToken.assign(nextToken, 0, 1) == ")" ){
				break;
			}
			else{
				// nextToken.replace(nextToken.find("("),1," ");
				// nextToken.replace(nextToken.find(")"),1," ");
				istringstream iss(nextToken);
				int tempint;
				// iss >> tempint;
				// cout << tempint << " ";	
				while(iss >> tempint){
					// face_neighbour.push_back(tempint);
					mesh.faces[temp_num].neighbour = tempint;
					
					// cout << tempint << " ";	
					++temp_num;
				}
				// cout << endl;
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
	
	
	
	
	inputFile.open("./grid/boundary");
	if(inputFile.fail()){
		cerr << "Unable to open file for reading." << endl;
		return EXIT_FAILURE;
	}
	vector<string> boundary_name;
	vector<char> boundary_type;
	vector<int> boundary_nFaces;
	vector<int> boundary_startFace;
	
	// string nextToken;
	string names;
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
				names = backToken;
				startInput=true;
			}
		}
		backToken = nextToken;
	}
	
	
	// cout << names << endl;
	
	// setToken.clear();
	
	vector<string> setToken2;
	startInput=false;
	for (auto item : setToken) {
		string asignToken;
		if(startInput){
			// item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
			// if( asignToken.assign(item, 0, 1) == "}" ){
			if( item.find("}") != string::npos ){
				names.erase(std::remove(names.begin(), names.end(), ' '), names.end());
				boundary_name.push_back(names);
				for (auto item2 : setToken2) {
				// cout << item2 << endl;
					if( item2.find("nFaces") != string::npos ){
						istringstream iss(item2);
						string temptemp;
						int temptempint;
						iss >> temptemp >> temptempint;
				// cout << temptemp << endl;
						boundary_nFaces.push_back(temptempint);
					}
					if( item2.find("startFace") != string::npos ){
						istringstream iss(item2);
						string temptemp;
						int temptempint;
						iss >> temptemp >> temptempint;
						boundary_startFace.push_back(temptempint);
						
					}
					// cout << item2 << endl;
				}
				startInput=false;
				setToken2.clear();
				// break;
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
        // cout << item << endl;
    }
	
	inputFile.close();
	
	temp_num=0;
	for (auto item : boundary_name) {
		mesh.addBoundary();
		mesh.boundary.back().name = item;
		mesh.boundary.back().nFaces = boundary_nFaces[temp_num];
		mesh.boundary.back().startFace = boundary_startFace[temp_num];
		// cout << mesh.boundary.back().name << endl;
		++temp_num;
	}
	
	
	// for (auto item : boundary_name) {
		// cout << item << endl;
	// }
	
	
	//=================== Loading End ==========================
	
	
	
	// add Cells
	int cell_num=0;
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		cell_num = max(cell_num , iter->owner);
	}
	for(int i=0; i<cell_num+1; ++i){
		mesh.addCell();
	}
	
	
	// set types
	// for (auto item : mesh.faces) {
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		
		if(iter->neighbour == -1){
			iter->setType(SEMO_Types::BOUNDARY_FACE);
		}
		else{
			iter->setType(SEMO_Types::INNER_FACE);
		}
	}
	// // for (auto item : mesh.faces) {
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		// if(iter->getType()==SEMO_Types::BOUNDARY_FACE){
			// printf("BOUNDARY_FACE\n");
		// }
	// }
	
	
	
	// connection
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// geometric
	SEMO_Mesh_Geometric geometric;
	
	// for (auto face : mesh.faces) {
		
		
		// // printf("%lf %lf %lf\n",face.unitNormals[0],face.unitNormals[1],face.unitNormals[2]);
		// // cout << face.unitNormals[0] << " " << face.unitNormals[1] << " " << face.unitNormals[2] << endl;
		
	// }
	
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		// iter->unitNormals[0] = -1000;
	// }
	// for(int i=0; i<mesh.faces.size(); ++i){
		// mesh.faces[i].unitNormals[0] = -1000;
	// }	
	
	
	// polygon face normal vectors
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : iter->points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		
		geometric.calcUnitNormals(
		iter->unitNormals,
		mesh.points[iter->points[0]].x,
		mesh.points[iter->points[0]].y,
		mesh.points[iter->points[0]].z,
		mesh.points[iter->points[1]].x,
		mesh.points[iter->points[1]].y,
		mesh.points[iter->points[1]].z,
		mesh.points[iter->points[2]].x,
		mesh.points[iter->points[2]].y,
		mesh.points[iter->points[2]].z
		);
		
	}
	
	// polygon face area
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : iter->points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
			// cout << mesh.points[iPoint].x << endl;
		}
		Vx.push_back(mesh.points[iter->points[0]].x);
		Vy.push_back(mesh.points[iter->points[0]].y);
		Vz.push_back(mesh.points[iter->points[0]].z);
		
		iter->area = 
		geometric.calcArea3dPolygon(
		iter->points.size(),
		Vx,Vy,Vz,
		iter->unitNormals[0],iter->unitNormals[1],iter->unitNormals[2]
		);
		
	}
	
	// polygon face center x,y,z
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		double fx=0, fy=0, fz=0;
		for(auto iPoint : iter->points){
			fx += mesh.points[iPoint].x;
			fy += mesh.points[iPoint].y;
			fz += mesh.points[iPoint].z;
		}
		
		fx /= iter->points.size();
		fy /= iter->points.size();
		fz /= iter->points.size();
		
		iter->x = fx;
		iter->y = fy;
		iter->z = fz;
		
	}
	
	
	// polyhedron cell volume
	for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		iter->volume = 0.0;
	}
	
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		
		// cout << iter->owner << endl;
		
		if(iter->getType() == SEMO_Types::INNER_FACE){
			mesh.cells[iter->owner].volume += 
				(iter->x*iter->unitNormals[0]+
				 iter->y*iter->unitNormals[1]+
				 iter->z*iter->unitNormals[2])*iter->area;
			
			mesh.cells[iter->neighbour].volume -= 
				(iter->x*iter->unitNormals[0]+
				 iter->y*iter->unitNormals[1]+
				 iter->z*iter->unitNormals[2])*iter->area;
		}
		else if(iter->getType() == SEMO_Types::BOUNDARY_FACE){
			mesh.cells[iter->owner].volume += 
				(iter->x*iter->unitNormals[0]+
				 iter->y*iter->unitNormals[1]+
				 iter->z*iter->unitNormals[2])*iter->area;
		}
		
	}
	for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		iter->volume /= 3.0;
		// cout << iter->volume << endl;
	}
	

	return EXIT_SUCCESS;
}














