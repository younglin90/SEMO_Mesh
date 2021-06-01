
/*
* Base64 encoding/decoding (RFC1341)
* Copyright (c) 2005-2011, Jouni Malinen <j@w1.fi>
*
* This software may be distributed under the terms of the BSD license.
* See README for more details.
*/

// 2016-12-12 - Gaspard Petit : Slightly modified to return a std::string 
// instead of a buffer allocated with malloc.

// #include <string>

static const unsigned char base64_table[65] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
* base64_encode - Base64 encode
* @src: Data to be encoded
* @len: Length of the data to be encoded
* @out_len: Pointer to output length variable, or %NULL if not used
* Returns: Allocated buffer of out_len bytes of encoded data,
* or empty string on failure
*/
std::string base64_encode(const unsigned char *src, size_t len)
{
    unsigned char *out, *pos;
    const unsigned char *end, *in;

    size_t olen;

    olen = 4*((len + 2) / 3); /* 3-byte blocks to 4-byte */

    if (olen < len)
        return std::string(); /* integer overflow */

    std::string outStr;
    outStr.resize(olen);
    out = (unsigned char*)&outStr[0];

    end = src + len;
    in = src;
    pos = out;
    while (end - in >= 3) {
        *pos++ = base64_table[in[0] >> 2];
        *pos++ = base64_table[((in[0] & 0x03) << 4) | (in[1] >> 4)];
        *pos++ = base64_table[((in[1] & 0x0f) << 2) | (in[2] >> 6)];
        *pos++ = base64_table[in[2] & 0x3f];
        in += 3;
    }

    if (end - in) {
        *pos++ = base64_table[in[0] >> 2];
        if (end - in == 1) {
            *pos++ = base64_table[(in[0] & 0x03) << 4];
            *pos++ = '=';
        }
        else {
            *pos++ = base64_table[((in[0] & 0x03) << 4) |
                (in[1] >> 4)];
            *pos++ = base64_table[(in[1] & 0x0f) << 2];
        }
        *pos++ = '=';
    }

    return outStr;
}








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






// // 파일 ,쓰기

	// cout << "| execute save file (.vtu format) ... ";
	// ofstream outputFile;
	// outputFile.open("./save/plot.vtu");
	// // outputFile.open("./save/plot.vtu",ios::binary);
	// // outputFile.open("./save/plot.vtu",ios::out | ios::binary);
	// if(outputFile.fail()){
		// cerr << "Unable to write file for writing." << endl;
		// return 1;
	// }
	
	// // string out_line;
	// outputFile << "<?xml version=\"1.0\"?>" << endl;
	// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	// outputFile << "  <UnstructuredGrid>" << endl;
	// outputFile << "   <Piece NumberOfPoints=\"" << listPoints.size() << "\" NumberOfCells=\"" << listCells.size() << "\">" << endl;
	
	// // Points data
	// outputFile << "    <PointData>" << endl;
	// outputFile << "    </PointData>" << endl;
	// // Cells data
	// outputFile << "    <CellData>" << endl;
	// outputFile << "    </CellData>" << endl;
	// // Points
	// outputFile << "    <Points>" << endl;
	// // }
	// outputFile << "     <DataArray type=\"Float32\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	// stringstream streamXYZ;
	// // for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	// for(auto iter=listPoints.begin(); iter!=listPoints.end(); iter++){
		// outputFile << scientific << (*iter)->x << " " << (*iter)->y << " " << (*iter)->z << endl;

	// }
	
	// outputFile << "    </DataArray>" << endl;
	// outputFile << "   </Points>" << endl;
	
	// // cells
	// outputFile << "   <Cells>" << endl; 
	// // connectivity (cell's points)
	// outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	// for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		// for(auto i : (*iter)->points){
			// outputFile << i << " ";
		// }
		// outputFile << endl;
	// }
	
	// outputFile << "    </DataArray>" << endl;
	
	// // offsets (cell's points offset)
	// int cellFaceOffset = 0;
	// outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	// cellFaceOffset = 0;
	// for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		// cellFaceOffset += (*iter)->points.size();
		// outputFile << cellFaceOffset << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // types (cell's type, 42 = polyhedron)
	// outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	// for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		// outputFile << "42" << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // faces (cell's faces number, each face's point number, cell's faces's points)
	// outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// // outputFile << mesh.faces.size() << endl;
	// for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		// outputFile << (*iter)->faces.size() << endl;
		// for(auto& i : (*iter)->faces){
			// outputFile << mesh.faces[i].points.size() << " ";
			// for(auto& j : mesh.faces[i].points){
				// outputFile << j << " ";
			// }
			// outputFile << endl;
		// }
	// }
	
	// outputFile << "    </DataArray>" << endl;
	
	// // faceoffsets (cell's face offset)
	// int cellFacePointOffset = 0;
	
	// outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	// cellFacePointOffset = 0;
	// for(auto iter=listCells.begin(); iter!=listCells.end(); iter++){
		// int numbering = 1 + (*iter)->faces.size();
		// for(auto& i : (*iter)->faces){
			// numbering += mesh.faces[i].points.size();
		// }
		// cellFacePointOffset += numbering;
		// outputFile << cellFacePointOffset << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	// outputFile << "   </Cells>" << endl;
	
	// outputFile << "  </Piece>" << endl;
	// outputFile << " </UnstructuredGrid>" << endl;
	// outputFile << "</VTKFile>" << endl;
	
	// outputFile.close();
	// cout << "-> completed" << endl;
	// cout << "------------------------------------" << endl;
	
















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

using namespace std;



enum class SEMO_Types{
	INTERNAL_FACE,
	BOUNDARY_FACE,
	CONNECTED_FACE
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
	cout << "------------------------------------" << endl;
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
			iter->setType(SEMO_Types::INTERNAL_FACE);
		}
	}

	
	
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
		else if(iter->getType() == SEMO_Types::CONNECTED_FACE){
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
	
	
	
	
	
	
	
	
	
	// geometric
	cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	SEMO_Mesh_Geometric geometric;

	
	// polygon face normal vectors & polygon face area
	// polygon face center x,y,z
	// 3D Polygon Areas : https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : iter->points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		
		geometric.calcUnitNormals_Area3dPolygon(
		iter->points.size(), Vx,Vy,Vz,
		iter->unitNormals, iter->area );
		
		iter->x = accumulate(Vx.begin(), Vx.end(), 0.0) / iter->points.size();
		iter->y = accumulate(Vy.begin(), Vy.end(), 0.0) / iter->points.size();
		iter->z = accumulate(Vz.begin(), Vz.end(), 0.0) / iter->points.size();
	}
	
	
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		iter->volume = 0.0;
	}
	
	for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		
		// cout << iter->owner << endl;
		
		if(iter->getType() == SEMO_Types::INTERNAL_FACE){
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
		
		if(iter->volume < std::numeric_limits<double>::min()) {
			cerr << endl;
			cerr << "#error, from calc cell volume, cell volume = " << iter->volume << " < cpu_min_val " << endl;
			cerr << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// cout << iter->volume << endl;
	}
	cout << "-> completed" << endl;
	
	
	
	
	// plotting
	
	// // // // MPI file reading
	
  // // // // Build the filename for the given process
  // // // std::string filename = "input_" + myrank + ".txt";

  // // // // Open the file stream and read or write
  // // // std::ifstream in(filename.c_str());
  // // // read_file(in);
  // // // in.close();
	
	
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
	

	

	return EXIT_SUCCESS;
}















	
	int nInFacesProcs[nparts]; std::fill_n(nInFacesProcs, nparts, 0);
	int nBdFacesProcs[nparts]; std::fill_n(nBdFacesProcs, nparts, 0);
	int nPrFacesProcs[nparts]; std::fill_n(nPrFacesProcs, nparts, 0);
	vector<int> iInFaceProcs, iBdFaceProcs, iPrFaceProcs;
	for(int ip=0; ip<nparts; ++ip){
		bool boolFace[mesh.faces.size()];
		std::fill_n(boolFace, mesh.faces.size(), false);
		// int orderFaces=0;
		int icell=0;
		for(auto& cell : mesh.cells){
			if(ip==npart[icell]){
				for(int i=0; i<cell.faces.size(); ++i){
					int iface = cell.faces[i];
					if(boolFace[iface] == false){
						if(mesh.faces[iface].getType() == SEMO_Types::INTERNAL_FACE){
							iInFaceProcs.push_back(iface);
							++nInFacesProcs[ip];
						}
						else if(mesh.faces[iface].getType() == SEMO_Types::BOUNDARY_FACE){
							iBdFaceProcs.push_back(iface);
							++nBdFacesProcs[ip];
						}
						else if(mesh.faces[iface].getType() == SEMO_Types::PROCESSOR_FACE){
							iPrFaceProcs.push_back(iface);
							++nPrFacesProcs[ip];
						}
						boolFace[iface] = true;
					}
				}
			}
			++icell;
		}
	}
	int strInFacesProcs[nparts+1]; strInFacesProcs[0]=0;
	int strBdFacesProcs[nparts+1]; strBdFacesProcs[0]=0;
	int strPrFacesProcs[nparts+1]; strPrFacesProcs[0]=0;
	for(int ip=1; ip<nparts+1; ++ip) {
		strInFacesProcs[ip]=strInFacesProcs[ip-1]+nInFacesProcs[ip-1];
		strBdFacesProcs[ip]=strBdFacesProcs[ip-1]+nBdFacesProcs[ip-1];
		strPrFacesProcs[ip]=strPrFacesProcs[ip-1]+nPrFacesProcs[ip-1];
	}
	
	
	vector<vector<int>> realPointsProcs;
	vector<vector<int>> nfacePointsProcs;
	vector<vector<int>> InfacePointsProcs;
	vector<vector<int>> BdfacePointsProcs;
	vector<vector<int>> PrfacePointsProcs;
	vector<vector<int>> InownerProcs;
	vector<vector<int>> BdownerProcs;
	vector<vector<int>> PrownerProcs;
	vector<vector<int>> InneighbourProcs;
	vector<vector<int>> InnfacePointsProcs;
	vector<vector<int>> BdnfacePointsProcs;
	vector<vector<int>> PrnfacePointsProcs;
	
	for(int ip=0; ip<nparts; ++ip){
		bool boolPoints[mesh.points.size()];
		bool boolFaces[mesh.faces.size()];
		std::fill_n(boolPoints, mesh.points.size(), false);
		std::fill_n(boolFaces, mesh.faces.size(), false);
		int orderCells[mesh.cells.size()];
		int orderPoints[mesh.points.size()];
		int orderFaces[mesh.faces.size()];
		int icell=0;
		int ipoint=0;
		int iface=0;
		int jcell=0;
		vector<int> myrealPointsProcs;
		for(auto& cell : mesh.cells){
			if(ip==npart[icell]){
				
				orderCells[icell] =jcell;
				++jcell;
				
				for(auto& point : cell.points){
					if(boolPoints[point] == false){
						orderPoints[point] = ipoint;
						myrealPointsProcs.push_back(point);
						++ipoint;
						boolPoints[point] = true;
					}
				}
				
				for(auto& face : cell.faces){
					if(boolFaces[face] == false){
						boolFaces[face] = true;
					}
				}
			}
			++icell;
		}
		
		realPointsProcs.push_back(myrealPointsProcs);
		
		
		ipoint=0;
		for(auto& point : mesh.points){
			if(boolPoints[ipoint]==true){
			}
			++ipoint;
		}
		
		vector<int> myInnfacePointsProcs;
		vector<int> myBdnfacePointsProcs;
		vector<int> myPrnfacePointsProcs;
		vector<int> myInfacePointsProcs;
		vector<int> myBdfacePointsProcs;
		vector<int> myPrfacePointsProcs;
		vector<int> myInownerProcs;
		vector<int> myBdownerProcs;
		vector<int> myPrownerProcs;
		vector<int> myInneighbourProcs;
		// vector<int> myBdneighbourProcs;
		// vector<int> myPrneighbourProcs;
		iface=0;
		for(auto& face : mesh.faces){
			if(boolFaces[iface]==true){
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					myInnfacePointsProcs.push_back(face.points.size());
					for(auto& point : face.points){
						myInfacePointsProcs.push_back( orderPoints[point] );
					}
					
					if(ip==npart[face.owner] && ip==npart[face.neighbour]){
						myInownerProcs.push_back(orderCells[face.owner]);
						myInneighbourProcs.push_back(orderCells[face.neighbour]);
					}
					else if(ip!=npart[face.owner] && ip==npart[face.neighbour]){
						myInownerProcs.push_back(orderCells[face.neighbour]);
						// myInneighbourProcs.push_back(orderCells[face.neighbour]);
					}
					else if(ip==npart[face.owner] && ip!=npart[face.neighbour]){
						myInownerProcs.push_back(orderCells[face.owner]);
						// myInneighbourProcs.push_back(orderCells[face.neighbour]);
					}
					
				}
				else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
					myBdnfacePointsProcs.push_back(face.points.size());
					for(auto& point : face.points){
						myBdfacePointsProcs.push_back( orderPoints[point] );
					}
					
					myBdownerProcs.push_back(orderCells[face.owner]);
					
				}
				else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					myPrnfacePointsProcs.push_back(face.points.size());
					for(auto& point : face.points){
						myPrfacePointsProcs.push_back( orderPoints[point] );
					}
					
					myPrownerProcs.push_back(orderCells[face.owner]);
					
				}
				
				
				
				
			}
			++iface;
		}
		
		InnfacePointsProcs.push_back(myInnfacePointsProcs);
		BdnfacePointsProcs.push_back(myBdnfacePointsProcs);
		PrnfacePointsProcs.push_back(myPrnfacePointsProcs);
		
		InfacePointsProcs.push_back(myInfacePointsProcs);
		BdfacePointsProcs.push_back(myBdfacePointsProcs);
		PrfacePointsProcs.push_back(myPrfacePointsProcs);
		
		InownerProcs.push_back(myInownerProcs);
		BdownerProcs.push_back(myBdownerProcs);
		PrownerProcs.push_back(myPrownerProcs);
		InneighbourProcs.push_back(myInneighbourProcs);
		
		
	}
	
	
	
	// save
	for(int ip=0; ip<nparts; ++ip){
		ofstream decomName;
		string sFilename = "./grid/decomposit/points." + to_string(ip);
		decomName.open(sFilename);
		decomName << realPointsProcs[ip].size() << endl;
		decomName << "(" << endl;
		for(auto& i : realPointsProcs[ip]){
			decomName << "(";
			decomName << mesh.points[i].x << " ";
			decomName << mesh.points[i].y << " ";
			decomName << mesh.points[i].z << ")" << endl;
		}
		decomName << ")";
		decomName.close();
	}
	
	for(int ip=0; ip<nparts; ++ip){
		ofstream decomName;
		string sFilename = "./grid/decomposit/faces." + to_string(ip);
		decomName.open(sFilename);
		decomName << InnfacePointsProcs[ip].size() + BdnfacePointsProcs[ip].size() + PrnfacePointsProcs[ip].size() << endl;
		decomName << "(" << endl;
		
		int ooooo=0;
		for(auto& nface : InnfacePointsProcs[ip]){
			decomName << nface << "(";
			for(int point=0; point<nface-1; ++point){
				decomName << InfacePointsProcs[ip][ooooo] << " ";
				++ooooo;
			}
			decomName << InfacePointsProcs[ip][ooooo] << ")" << endl;
			++ooooo;
		}
		
		ooooo=0;
		for(auto& nface : BdnfacePointsProcs[ip]){
			decomName << nface << "(";
			for(int point=0; point<nface-1; ++point){
				decomName << BdfacePointsProcs[ip][ooooo] << " ";
				++ooooo;
			}
			decomName << BdfacePointsProcs[ip][ooooo] << ")" << endl;
			++ooooo;
		}
		
		ooooo=0;
		for(auto& nface : PrnfacePointsProcs[ip]){
			decomName << nface << "(";
			for(int point=0; point<nface-1; ++point){
				decomName << PrfacePointsProcs[ip][ooooo] << " ";
				++ooooo;
			}
			decomName << PrfacePointsProcs[ip][ooooo] << ")" << endl;
			++ooooo;
		}
		
		decomName << ")";
		decomName.close();
	}
	
	for(int ip=0; ip<nparts; ++ip){
		ofstream decomName;
		string sFilename = "./grid/decomposit/owner." + to_string(ip);
		decomName.open(sFilename);
		decomName << InownerProcs[ip].size() + BdownerProcs[ip].size() + PrownerProcs[ip].size() << endl;
		decomName << "(" << endl;
		
		for(auto& i : InownerProcs[ip]){
			decomName << i << endl;
		}
		for(auto& i : BdownerProcs[ip]){
			decomName << i << endl;
		}
		for(auto& i : PrownerProcs[ip]){
			decomName << i << endl;
		}
		decomName << ")";
		decomName.close();
	}
	
	for(int ip=0; ip<nparts; ++ip){
		ofstream decomName;
		string sFilename = "./grid/decomposit/neighbour." + to_string(ip);
		decomName.open(sFilename);
		decomName << InneighbourProcs[ip].size() << endl;
		decomName << "(" << endl;
		for(auto& i : InneighbourProcs[ip]){
			decomName << i << endl;
		}
		decomName << ")";
		decomName.close();
	}
	
	
	// for(int ip=0; ip<nparts; ++ip){
		// ofstream decomName;
		// string sFilename = "./grid/decomposit/boundary." + to_string(ip);
		// decomName.open(sFilename);
		// decomName << "(" << endl;
		// for(int i=0; i<nparts; ++i){
			// if(ip==i) continue;
			// if(myadjncyBoundaryNFaces[np][i]==0) continue;
			// string bcnames = "   procBoundary" + to_string(np) + "to" + to_string(i);
			// decomName << bcnames << endl;
			// decomName << "   {" << endl;
			// decomName << "      type processor" << endl;
			// decomName << "      nFaces " << myadjncyBoundaryNFaces[np][i] << endl;
			// decomName << "      startFace " << myadjncyBoundaryStartFace[np][i] << endl;
			// decomName << "      myProcNo " << np << endl;
			// decomName << "      neighbProcNo " << i << endl;
			// decomName << "   }" << endl;
		// }
		// decomName << ")";
		// decomName.close();
	// }
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// int nPointsProcs[nparts];
	// std::fill_n(nPointsProcs, nparts, 0);
	// vector<int> iPointsProcs;
	// vector<int> iPointsAtProcs;
	// vector<int> atFacesFile;
	// for(int ip=0; ip<nparts; ++ip){
		// bool boolPoint[mesh.points.size()];
		// std::fill_n(boolPoint, mesh.points.size(), false);
		// int orderProcs[mesh.points.size()];
		// int orPoint=0;
		// for(int i=strInFacesProcs[ip]; i<strInFacesProcs[ip+1]; ++i){
			// for(auto& ipoint : mesh.faces[i].points){
				// if(boolPoint[ipoint]==false){
					// iPointsProcs.push_back(ipoint);
					// orderProcs[ipoint] = orPoint;
					// ++orPoint;
					// boolPoint[ipoint]=true;
				// }
			// }
		// }
		
		// for(int i=strBdFacesProcs[ip]; i<strBdFacesProcs[ip+1]; ++i){
			// for(auto& ipoint : mesh.faces[i].points){
				// if(boolPoint[ipoint]==false){
					// iPointsProcs.push_back(ipoint);
					// orderProcs[ipoint] = orPoint;
					// ++orPoint;
					// boolPoint[ipoint]=true;
				// }
			// }
		// }
		
		// for(int i=strPrFacesProcs[ip]; i<strPrFacesProcs[ip+1]; ++i){
			// for(auto& ipoint : mesh.faces[i].points){
				// if(boolPoint[ipoint]==false){
					// iPointsProcs.push_back(ipoint);
					// orderProcs[ipoint] = orPoint;
					// ++orPoint;
					// boolPoint[ipoint]=true;
				// }
			// }
		// }
		
		
		
		// atFacesFile.push_back(
		
	// }
	
	
	
	
	
	// vector<vector<bool>> tempPoints(mesh.points.size(),vector<bool>(nparts,false));
	// vector<vector<bool>> tempFaces(mesh.faces.size(),vector<bool>(nparts,false));
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(auto& point : mesh.cells[i].points){
			// tempPoints[point][npart[i]] = true;
		// }
		// for(auto& face : mesh.cells[i].faces){
			// tempFaces[face][npart[i]] = true;
		// }
	// }
	
	
	// int nadjPoints[nparts];
	// std::fill_n(nadjPoints, nparts, 0);
	// vector<int> adjncyPoints;
	// for(int np=0; np<nparts; ++np){
		// for(int i=0; i<mesh.points.size(); ++i){
			// if(tempPoints[i][np] == true){
				// adjncyPoints.push_back(i);
				// ++nadjPoints[np];
			// }
		// }
	// }
	// int xadjPoints[nparts+1];
	// xadjPoints[0]=0;
	// for(int np=1; np<nparts+1; ++np){
		// xadjPoints[np] = xadjPoints[np-1] + nadjPoints[np-1];
	// }
	
	
	
	// // int nadjFaces[nparts];
	// // std::fill_n(nadjFaces, nparts, 0);
	// // vector<int> adjncyFaces;
	// // for(int np=0; np<nparts; ++np){
		// // for(int i=0; i<mesh.faces.size(); ++i){
			// // for(int j=0; j<tempFaces[i].size(); ++j){
				// // if(np==tempFaces[i][j]){
					// // adjncyFaces.push_back(i);
					// // ++nadjFaces[np];
					
				// // }
			// // }
		// // }
	// // }
	// // int xadjFaces[nparts+1];
	// // xadjFaces[0]=0;
	// // for(int np=1; np<nparts+1; ++np){
		// // xadjFaces[np] = xadjFaces[np-1] + nadjFaces[np-1];
	// // }
	
	// //
	// int nadjFaces[nparts];
	// std::fill_n(nadjFaces, nparts, 0);
	// vector<vector<int>> adjncyFaces;
	// for(int np=0; np<nparts; ++np){
		// int temporder = 0;
		// int forFacePoints[mesh.points.size()];
		// bool boolFacePoints[mesh.points.size()];
		// std::fill_n(boolFacePoints, mesh.points.size(), false);
		// for(int i=0; i<mesh.faces.size(); ++i){
			// for(int j=0; j<mesh.faces[i].points.size(); ++j){
				// if(tempPoints[ mesh.faces[i].points[j] ][np] == true){
					// if(boolFacePoints[ mesh.faces[i].points[j] ] == false){
						// forFacePoints[ mesh.faces[i].points[j] ] = temporder;
						// ++temporder;
						// boolFacePoints[ mesh.faces[i].points[j] ] = true;
					// }
				// }
			// }
		// }
		
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// if(tempFaces[i][np]==true){
				// vector<int> temppush;
				// temppush.push_back(mesh.faces[i].points.size());
				// for(int j=0; j<mesh.faces[i].points.size(); ++j){
					// temppush.push_back( forFacePoints[ mesh.faces[i].points[j] ] );
				// }
				
				// adjncyFaces.push_back(temppush);
				// ++nadjFaces[np];
			// }
		// }
	// }
	// int xadjFaces[nparts+1];
	// xadjFaces[0]=0;
	// for(int np=1; np<nparts+1; ++np){
		// xadjFaces[np] = xadjFaces[np-1] + nadjFaces[np-1];
	// }
	// //
	
	// // int nadjCells[nparts];
	// // std::fill_n(nadjCells, nparts, 0);
	// // vector<int> adjncyCells;
	// // for(int np=0; np<nparts; ++np){
		// // for(int i=0; i<mesh.cells.size(); ++i){
			// // if(np==npart[i]){
				// // adjncyCells.push_back(i);
				// // ++nadjCells[np];
			// // }
		// // }
	// // }
	// // int xadjCells[nparts+1];
	// // xadjCells[0]=0;
	// // for(int np=1; np<nparts+1; ++np){
		// // xadjCells[np] = xadjCells[np-1] + nadjCells[np-1];
	// // }
	
	
	
	// int nadjOwner[nparts];
	// int nadjNeighbour[nparts];
	// std::fill_n(nadjOwner, nparts, 0);
	// std::fill_n(nadjNeighbour, nparts, 0);
	// vector<int> adjncyOwner;
	// vector<int> adjncyNeighbour;
	// for(int np=0; np<nparts; ++np){
		// int temporder = 0;
		// int forOwner[mesh.cells.size()];
		// // int forNeighbour[mesh.cells.size()];
		// bool boolOwner[mesh.cells.size()];
		// // bool boolNeighbour[mesh.cells.size()];
		// std::fill_n(boolOwner, mesh.cells.size(), false);
		// // std::fill_n(boolNeighbour, mesh.cells.size(), false);
		// for(int i=0; i<mesh.faces.size(); ++i){
			// if(tempFaces[i][np]==true){
				// if(boolOwner[ mesh.faces[i].owner ] == false){
					// forOwner[ mesh.faces[i].owner ] = temporder;
					// boolOwner[ mesh.faces[i].owner ] = true;
					// ++temporder;
				// }
			// }
		// }
		
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// if(tempFaces[i][np]==true){
				// adjncyOwner.push_back( forOwner[ mesh.faces[i].owner ] );
				// ++nadjOwner[np];
				// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
					// adjncyNeighbour.push_back( forOwner[ mesh.faces[i].neighbour ] );
					// ++nadjNeighbour[np];
				// }
			// }
		// }
	// }
	// int xadjOwner[nparts+1];
	// int xadjNeighbour[nparts+1];
	// xadjOwner[0]=0;
	// xadjNeighbour[0]=0;
	// for(int np=1; np<nparts+1; ++np){
		// xadjOwner[np] = xadjOwner[np-1] + nadjOwner[np-1];
		// xadjNeighbour[np] = xadjNeighbour[np-1] + nadjNeighbour[np-1];
	// }
	
	
	
	// // int nadjBoundary[nparts];
	// // std::fill_n(nadjBoundary, nparts, 0);
	// // vector<vector<string>> adjncyBoundaryNames;
	// // vector<vector<int>> adjncyBoundaryNFaces;
	// // vector<vector<int>> adjncyBoundaryStartFace;
	// int myadjncyBoundaryNFaces[nparts][nparts];
	// int myadjncyBoundaryStartFace[nparts][nparts];
	// for(int np=0; np<nparts; ++np){
		// int temporder = 0;
		// int forBoundary[mesh.faces.size()];
		// bool boolBoundary[mesh.faces.size()];
		// std::fill_n(boolBoundary, mesh.faces.size(), false);
		// for(int i=0; i<mesh.faces.size(); ++i){
			// if(tempFaces[i][np] == true){
				// if(boolBoundary[i] == false){
					// forBoundary[i] = temporder;
					// ++temporder;
					// boolBoundary[i] = true;
				// }
			// }
		// }
		
		
		// bool bcstartface[nparts];
		// std::fill_n(bcstartface, nparts, false);
		// std::fill_n(myadjncyBoundaryStartFace[np], nparts, 0);
		// std::fill_n(myadjncyBoundaryNFaces[np], nparts, 0);
		// for(int i=0; i<mesh.faces.size(); ++i){
			// if(tempFaces[i][np]==true){
				// // vector<string> myadjncyBoundaryNames;
				// // vector<int> myadjncyBoundaryStartFace;
				// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
					// if( npart[mesh.faces[i].owner] != npart[mesh.faces[i].neighbour] ){
						
						
						// if( np == npart[mesh.faces[i].owner]){
							// // cout << mesh.faces[i].neighbour << endl;
							// // cout << ncells << endl;
							// // cout << mesh.cells.size() << endl;
							// // cout << npart[mesh.faces[i].neighbour] << endl;
							// ++myadjncyBoundaryNFaces[np][ npart[mesh.faces[i].neighbour] ];
							// if( bcstartface[ npart[mesh.faces[i].neighbour] ] == false ){
								// myadjncyBoundaryStartFace[np][ npart[mesh.faces[i].neighbour] ] = forBoundary[i];
								// bcstartface[ npart[mesh.faces[i].neighbour] ] = true;
							// }
						// }
						// if( np == npart[mesh.faces[i].neighbour]){
							// ++myadjncyBoundaryNFaces[np][ npart[mesh.faces[i].owner] ];
							// if( bcstartface[ npart[mesh.faces[i].owner] ] == false ){
								// // cout << np << forBoundary[i] << endl;
								// myadjncyBoundaryStartFace[np][ npart[mesh.faces[i].owner] ] = forBoundary[i];
								// bcstartface[ npart[mesh.faces[i].owner] ] = true;
							// }
						// }
					// }
				// }
				// // else if(mesh.faces[i].getType() == SEMO_Types::BOUNDARY_FACE){
					// // for(int j=0; j<mesh.boundary.size(); ++j){
						// // if( mesh.boundary[j].name
					// // }
				// // }
			// }
		// }
	// }
	
	
	
	// // moving cells each processes
	// // ofstream decomPointsName[nparts];
	// // ofstream decomFacesName[nparts];
	// // ofstream decomOwnerName[nparts];
	// // ofstream decomNeighbourName[nparts];
	// // ofstream decomBoundaryName[nparts];
	// // int filenametempnum = 0;
	// // for(int i=0; i<nparts; ++i){
		// // string pointsFilename = "./grid/decomposit/points." + to_string(i);
		// // string facesFilename = "./grid/decomposit/faces." + to_string(i);
		// // string ownerFilename = "./grid/decomposit/owner." + to_string(i);
		// // string neighbourFilename = "./grid/decomposit/neighbour." + to_string(i);
		// // string boundaryFilename = "./grid/decomposit/boundary." + to_string(i);
		// // decomPointsName[i].open(pointsFilename);
		// // decomFacesName[i].open(facesFilename);
		// // decomOwnerName[i].open(ownerFilename);
		// // decomNeighbourName[i].open(neighbourFilename);
		// // decomBoundaryName[i].open(boundaryFilename);
	// // }
	
	// for(int np=0; np<nparts; ++np){
		// ofstream decomPointsName;
		// string sFilename = "./grid/decomposit/points." + to_string(np);
		// decomPointsName.open(sFilename);
		// decomPointsName << nadjPoints[np] << endl;
		// decomPointsName << "(" << endl;
		// for(int i=xadjPoints[np]; i<xadjPoints[np+1]; ++i){
			// decomPointsName << "(";
			// decomPointsName << mesh.points[adjncyPoints[i]].x << " ";
			// decomPointsName << mesh.points[adjncyPoints[i]].y << " ";
			// decomPointsName << mesh.points[adjncyPoints[i]].z << ")" << endl;
		// }
		// decomPointsName << ")";
		// decomPointsName.close();
	// }
	
	// for(int np=0; np<nparts; ++np){
		// ofstream decomPointsName;
		// string sFilename = "./grid/decomposit/faces." + to_string(np);
		// decomPointsName.open(sFilename);
		// decomPointsName << nadjFaces[np] << endl;
		// decomPointsName << "(" << endl;
		// for(int i=xadjFaces[np]; i<xadjFaces[np+1]; ++i){
			// decomPointsName << adjncyFaces[i][0];
			// decomPointsName << "(";
			// for(int j=0; j<adjncyFaces[i][0]-1; ++j){
				// decomPointsName << adjncyFaces[i][j+1] << " ";
			// }
			// decomPointsName << adjncyFaces[i][adjncyFaces[i][0]] << ")" << endl;
		// }
		// decomPointsName << ")";
		// decomPointsName.close();
	// }
	
	// for(int np=0; np<nparts; ++np){
		// ofstream decomPointsName;
		// string sFilename = "./grid/decomposit/owner." + to_string(np);
		// decomPointsName.open(sFilename);
		// decomPointsName << nadjOwner[np] << endl;
		// decomPointsName << "(" << endl;
		// for(int i=xadjOwner[np]; i<xadjOwner[np+1]; ++i){
			// decomPointsName << adjncyOwner[i] << endl;
		// }
		// decomPointsName << ")";
		// decomPointsName.close();
	// }
	
	// for(int np=0; np<nparts; ++np){
		// ofstream decomPointsName;
		// string sFilename = "./grid/decomposit/neighbour." + to_string(np);
		// decomPointsName.open(sFilename);
		// decomPointsName << nadjNeighbour[np] << endl;
		// decomPointsName << "(" << endl;
		// for(int i=xadjNeighbour[np]; i<xadjNeighbour[np+1]; ++i){
			// decomPointsName << adjncyNeighbour[i] << endl;
		// }
		// decomPointsName << ")";
		// decomPointsName.close();
	// }
	
	
	// for(int np=0; np<nparts; ++np){
		// ofstream decomPointsName;
		// string sFilename = "./grid/decomposit/boundary." + to_string(np);
		// decomPointsName.open(sFilename);
		// // decomPointsName << nadjBoundary[np] << endl;
		// decomPointsName << "(" << endl;
		// for(int i=0; i<nparts; ++i){
			// if(np==i) continue;
			// if(myadjncyBoundaryNFaces[np][i]==0) continue;
			// string bcnames = "   procBoundary" + to_string(np) + "to" + to_string(i);
			// decomPointsName << bcnames << endl;
			// decomPointsName << "   {" << endl;
			// decomPointsName << "      type processor" << endl;
			// decomPointsName << "      nFaces " << myadjncyBoundaryNFaces[np][i] << endl;
			// decomPointsName << "      startFace " << myadjncyBoundaryStartFace[np][i] << endl;
			// decomPointsName << "      myProcNo " << np << endl;
			// decomPointsName << "      neighbProcNo " << i << endl;
			// decomPointsName << "   }" << endl;
		// }
		// decomPointsName << ")";
		// decomPointsName.close();
	// }
	
	
	
	
	







	// int METIS_STATUS;
    // idx_t objval;
	// idx_t nparts=size;
	// idx_t wgtflag=0;
	// idx_t numflag=0;
	// idx_t ncon=1;
	// real_t tpwgts[nparts*ncon];
	// for(int i=0;i<nparts*ncon;++i){
		// tpwgts[i]=1.0/nparts;
	// }
	// real_t ubvec[ncon];
	// // std::fill(ubvec, ubvec+ncon, 1.05);
	// // std::fill(ubvec[0], ubvec[ncon], 1.05);
	// std::fill_n(ubvec, ncon, 1.05);
	// real_t itr = 1000.0;
	// idx_t options[METIS_NOPTIONS];
	// METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	
	// MPI_Comm comm;
	// MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	// idx_t vtxdist[4];
	// idx_t adjncy[2];
	// idx_t xadj[2];
	// idx_t npart[1];
	
	// vtxdist[0] = 0;
	// vtxdist[1] = 1;
	// vtxdist[2] = 6;
	// vtxdist[3] = 9;
	
	// adjncy[0] = 1;
	// adjncy[1] = 3;
	
	// xadj[0] = 0;
	// xadj[1] = 2;
	
	
	// ParMETIS_V3_PartKway(
		// vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
		// &ncon,       &nparts, tpwgts, ubvec,
		// options, &objval, npart, &comm);


	// vtxdist

	// idx_t ncommon = 30;
	
	// idx_t eptr[ncells+1];
	// eptr[0] = 0;
	// int tempnumbering=0;
	// int tempnumbering2=0;
	// for(auto& cell : mesh.cells){
		// ++tempnumbering;
		// eptr[tempnumbering] = eptr[tempnumbering-1] + cell.points.size();
		// tempnumbering2 += cell.points.size();
	// }
	// idx_t eind[tempnumbering2];
	// tempnumbering=0;
	// for(auto& cell : mesh.cells){
		// for(auto& point : cell.points){
			// eind[tempnumbering] = point;
		// }
	// }
	
	// idx_t epart[ncells];











#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "build.h"
#include <cmath>
#include <array>
#include "mpi.h"

#include "parmetis.h" 
#include "scotch.h" 

#include "../mpi/build.h"
#include "../math/utility.h"

void SEMO_Mesh_Builder::distributeOneToAll(string type){

		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	SEMO_Mesh_Builder& mesh = *this;
	
	SEMO_Utility_Math utility;
	
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	int nBlocks = size;
		
		
		
		
		
		
	mesh.distributeEvenlyOneToAll();
		
		
		
		
		
		
		
	mesh.check();
	
	mesh.buildCells();
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// // check list
	// mesh.checkLists();
	
	// cell's faces connection
	mesh.connectCelltoFaces();
	
	// cell's points connection
	mesh.connectCelltoPoints();
	
	// set processor face counts
	mesh.setCountsProcFaces();
	
	// set processor face displacements
	mesh.setDisplsProcFaces(); 
	
	
	// mesh.checkQualities();
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	int idBlockCell[mesh.cells.size()];
	
	
	mesh.parMETIS(nBlocks, idBlockCell);
	
	//======================================================================
	

	vector<int> sendValues(0,0);
	vector<int> recvValues(0,0);
	vector<int> sendCounts(nBlocks,0);
	vector<int> recvCounts(nBlocks,0);
	vector<int> sendDisps(nBlocks,0);
	vector<int> recvDisps(nBlocks,0);
	sendDisps[0] = 0;
	recvDisps[0] = 0;
	int sendSize = 0;
	int recvSize = 0;
	
	vector<SEMO_Mesh_Builder> newMesh(nBlocks);
	

	vector<int> rank_Send;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			rank_Send.push_back(rank);
		}
	}
	vector<int> rank_Recv(rank_Send.size(),0);
	MPI_Alltoallv( rank_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   rank_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	// each points block id , local cell id , local N cells
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0)); // point block id (copies)
	vector<int> nCellsLocal(nBlocks,0); // local total cells
	vector<int> idCellLocal(mesh.cells.size(),0); // local cell id
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
	
	int nTotalLocalPointSize = 0; // total all # of local point
	for(auto& i : idBlockPoint) nTotalLocalPointSize += i.size();

	// point distribution (CSR format)
	//
	// / gP : global points / lP : local points / BL : block id /
	//
	//     gP[0]     gP[1]      gP[2]      gP[3]   gP[4] ...
	// - - - - - - - - - - - - - - - - - - - - - - - - - ....
	// |  lP[0,5]  |lP[6,9]| lP[10,15] |  ....   |     |
	// |  BL[0,5]  |BL[6,9]| BL[10,15] |  ....   |     |
	// |           |       |           |         |     |
	// strPoints[0]|  strPoints[2]  strPoints[3] |  strPoints[5]  ...
	//         strPoints[1]                  strPoints[4]
	//
	vector<int> nPointsLocal_Send(nBlocks,0); // local total points
	vector<int> idPointLocal(nTotalLocalPointSize,0); // local point id
	vector<int> idBlockPointLocal(nTotalLocalPointSize,0); // local point block id
	vector<int> strPoints(mesh.points.size()+1,0); // start of each global point
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(int j=0; j<idBlockPoint[i].size(); ++j){
			int idBlock = idBlockPoint[i][j];
			idPointLocal[nIndex] = nPointsLocal_Send[ idBlock ];
			++nPointsLocal_Send[ idBlock ];
			idBlockPointLocal[nIndex] = idBlock;
			++nIndex;
		}
	}
	strPoints[mesh.points.size()] = nIndex;
	
	
	// MPI local npoints
	//     local point's x y z
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> nPointsLocal_Recv(nBlocks,0);
	
	MPI_Alltoallv( nPointsLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nPointsLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	

	for(int i=0; i<nBlocks; ++i) sendCounts[i] = nPointsLocal_Send[i]*3;
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendSize = 0;
	for(auto& i : sendCounts) sendSize += i;
	double* xyz_Send = new double[sendSize]; 
	int nDisplPoint[nBlocks];
	std::fill_n(nDisplPoint, nBlocks, 0);
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			int displPoint = sendDisps[idBlock];
			displPoint += nDisplPoint[idBlock]*3;
			if(abs(mesh.points[i].x) < 1.e-300) mesh.points[i].x = 0.0;
			if(abs(mesh.points[i].y) < 1.e-300) mesh.points[i].y = 0.0;
			if(abs(mesh.points[i].z) < 1.e-300) mesh.points[i].z = 0.0;
			xyz_Send[displPoint] = mesh.points[i].x;
			xyz_Send[displPoint+1] = mesh.points[i].y;
			xyz_Send[displPoint+2] = mesh.points[i].z;
			++nDisplPoint[idBlock];
		}
	}
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nPointsLocal_Recv[i]*3;
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	vector<double> xyz_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0.0);
	
	
	MPI_Alltoallv( xyz_Send, sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   xyz_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
			
	for(int ip=0; ip<nBlocks; ++ip){
		int i=recvDisps[ip];
		int num = 0;
		while(i<recvDisps[ip] + recvCounts[ip]){
			// cout << recvDisps[ip] + recvCounts[ip] << " " << i << endl;
			newMesh[ip].addPoint();
			newMesh[ip].points[num].x = xyz_Recv[i]; ++i;
			newMesh[ip].points[num].y = xyz_Recv[i]; ++i;
			newMesh[ip].points[num].z = xyz_Recv[i]; ++i;
			++num;
		}
	}
	delete[] xyz_Send; 
	xyz_Recv.clear();
	
	
	
	// local processor face id order
	
	vector<int> idBlockCell_Send;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idBlockCell_Send.push_back(idBlockCell[face.owner]);
		}
	}
	vector<int> idBlockCell_Recv(idBlockCell_Send.size(),0);
	MPI_Alltoallv( idBlockCell_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   idBlockCell_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.neighbour = proc_num;
			++proc_num;
		}
	}

	// ====================   Internal face 
	// face setting
	proc_num = 0;
	vector<int> nFacesLocal_Send(nBlocks,0); // total local faces
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal_Send[idBlockOwner];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				// face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				++nFacesLocal_Send[idBlockNeighbour];
			}
		}
	}
	
	
	
	
	
	vector<int> nFacesLocal_Recv(nBlocks,0);

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	MPI_Alltoallv( nFacesLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFacesLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	vector<int> idFaceLocal_Send(nBlocks,0);  // temporary local face id
	vector<vector<int>> idFacePoint_Send(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE || 
		   // face.getType() == SEMO_Types::TO_BE_INTERNAL_FACE){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlock = idBlockCell[face.owner];
			++idFaceLocal_Send[idBlock];
			
			vector<int> idFacePoint; // face's points id
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			
			// save face's points size , local points id
			idFacePoint_Send[idBlock].push_back(face.points.size());
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
		}
	}
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		sendValues[i] = idFacePoint_Send[i].size();
	}
	recvValues.resize(nBlocks,0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	// sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);

	MPI_Alltoallv( sendValues.data(),       sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idFacePoint_Recv(nBlocks+1,0);
	str_idFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// ====================   boundary face 
	
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
	
	// boundary faces
	vector<int> nFacesBoundaryLocal(nBlocks,0); // local boundary face size
	vector<vector<int>> nFacesEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0));
	vector<vector<int>> nStartFaceEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0));
	vector<vector<bool>> boolStartFaceEachBoundaryLocal(nBlocks,vector<bool>(mesh.boundary.size(),false));
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int idBlock = idBlockCell[face.owner];
			++nFacesBoundaryLocal[ idBlock ];
			
			// mesh.getBoundaryNumber();
			++nFacesEachBoundaryLocal[ idBlock ][ face.getTypeBC() ];
			
			// if(idBlock==1) cout << face.getTypeBC() << endl;
			
			if(boolStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] == false){
				nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = idFaceLocal_Send[ idBlock ];
				boolStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = true;
			}
			++idFaceLocal_Send[ idBlock ];
			
			// cout << face.getTypeBC() << endl;
			
			// wirte bc
			vector<int> idFacePoint;
			idFacePoint.clear();
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			
			// save face's points size , local points id
			idFacePoint_Send[idBlock].push_back(face.points.size());
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
		
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		sendValues[i] = idFacePoint_Send[i].size();
	}
	recvValues.resize(nBlocks,0);

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];

	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idBoundaryFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idBoundaryFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idBoundaryFacePoint_Recv(nBlocks+1,0);
	str_idBoundaryFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idBoundaryFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];


	// //===============
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
	
	// PROCESSOR_FACE
	vector<int> nFacesOriginProcessorLocal(nBlocks,0); 
	int sendCounts_temp[nBlocks][nBlocks]; // sending size from each processor to each processor , i<->j
	for(int i=0; i<nBlocks; ++i){ for(int j=0; j<nBlocks; ++j) { sendCounts_temp[i][j]=0; } }
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::PROCESSOR_FACE ){
			int idBlock = idBlockCell[face.owner];
			++nFacesOriginProcessorLocal[idBlock];
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell_Recv[face.neighbour];
			
			++sendCounts_temp[idBlockOwner][idBlockNeighbour];
			// ++sendCounts_temp[idBlockNeighbour][idBlockOwner];
			
			vector<int> idFacePoint;
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
		
			// save face's points size , local points id
			idFacePoint_Send[idBlock].push_back(face.points.size());
			
			//######################################
			if( 
			idBlock == idBlockCell_Recv[face.neighbour] &&
			rank > rank_Send[face.neighbour]
			){
				std::reverse(idFacePoint.begin(),idFacePoint.end());
			}
			
			//######################################
			
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
	}
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
	recvValues.resize(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idOriginProcessorFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),                sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idOriginProcessorFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idOriginProcessorFacePoint_Recv(nBlocks+1,0);
	str_idOriginProcessorFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idOriginProcessorFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];


	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	// //====================
	
	
	
	// ====================   processor face 
	// TO_BE_PROCESSOR_FACE
	vector<int> nFacesProcessorLocal(nBlocks,0); // local processor face size
	int sendCounts_temp2[nBlocks][nBlocks]; // sending size from each processor to each processor , i<->j
	for(int i=0; i<nBlocks; ++i){ for(int j=0; j<nBlocks; ++j) { sendCounts_temp2[i][j]=0; } }
	vector<int> idFacesProcessor(0,0); // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::TO_BE_PROCESSOR_FACE ){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			++sendCounts_temp2[idBlockOwner][idBlockNeighbour];
			++sendCounts_temp2[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
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
			
				idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				for(int i=0; i<idFacePoint.size(); ++i){
					// cout << idFacePoint[i] << endl;
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
				}
				
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
			
				idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				std::reverse(idFacePoint.begin(),idFacePoint.end());
				for(int i=0; i<idFacePoint.size(); ++i){
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
				}
					
			}
		}
	}
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
	recvValues.resize(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idProcessorFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),                sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idProcessorFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idProcessorFacePoint_Recv(nBlocks+1,0);
	str_idProcessorFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idProcessorFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// //====================
	
	
	
	
	// write owner of internal face
	vector<vector<int>> idOwnerLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::INTERNAL_FACE ){
			int j = face.owner;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(int i=0; i<nBlocks; ++i){
		for(int j=0; j<idFacesProcessor.size(); ++j){
			int k = idFacesProcessor[j];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			if(n==i) {
				idOwnerLocal[m].push_back( idCellLocal[ mesh.faces[k].owner ] );
			}
			else if(m==i) {
				idOwnerLocal[n].push_back( idCellLocal[ mesh.faces[k].neighbour ] );
			}
		}
	}
	
	
	 
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idOwnerLocal[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	sendValues.clear();
	for(auto& i : idOwnerLocal){ for(auto& j : i){ sendValues.push_back(j); } }
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nFacesLocal_Recv[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idOwnerLocal_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idOwnerLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	vector<int> str_idOwnerLocal_Recv(nBlocks+1,0);
	str_idOwnerLocal_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idOwnerLocal_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// //====================

	// write # of neighbour
	vector<vector<int>> idNeighbourLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int j = face.neighbour;
			int i = idBlockCell[j];
			idNeighbourLocal[i].push_back(idCellLocal[j]);
		}
	}
	
	vector<int> nNeighbourLocal(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		nNeighbourLocal[i] = nFacesLocal_Send[i] 
		- nFacesBoundaryLocal[i] - nFacesOriginProcessorLocal[i] - nFacesProcessorLocal[i];
	}
	// cout<< idNeighbourLocal[0].size() << " " << nNeighbourLocal[0] << endl;
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : nNeighbourLocal){ sendValues.push_back(i); }
	
	vector<int> nNeighbourLocal_Recv(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(),           sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nNeighbourLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idNeighbourLocal[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	sendValues.clear();
	for(auto& i : idNeighbourLocal){ for(auto& j : i){ sendValues.push_back(j); } }
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nNeighbourLocal_Recv[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idNeighbourLocal_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),            sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idNeighbourLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	vector<int> str_idNeighbourLocal_Recv(nBlocks+1,0);
	str_idNeighbourLocal_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idNeighbourLocal_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// ==================
	
	
	
	// write of boundary faces
	int bound_num=0;
	for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
		if(mesh.boundary[ibcs].myProcNo == -1){
			++bound_num;
		}
	}
	
	int nFaceBoundaryFaceLocal_Send[nBlocks*bound_num];
	int nStartBoundaryFaceLocal_Send[nBlocks*bound_num];
	int temp_bound=0;
	for(int i=0; i<nBlocks; ++i){
		for(int ibcs=0; ibcs<bound_num; ++ibcs){
			if( nFacesEachBoundaryLocal[i][ibcs] > 0) {
				nFaceBoundaryFaceLocal_Send[temp_bound] = nFacesEachBoundaryLocal[i][ibcs];
				nStartBoundaryFaceLocal_Send[temp_bound] = nStartFaceEachBoundaryLocal[i][ibcs];
			}
			else{
				nFaceBoundaryFaceLocal_Send[temp_bound] = 0;
				nStartBoundaryFaceLocal_Send[temp_bound] = 0;
			}
			++temp_bound;
		}
	}
		
	
	std::fill(sendCounts.begin(),sendCounts.end(),bound_num);
	std::fill(recvCounts.begin(),recvCounts.end(),bound_num);
	
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : nFaceBoundaryFaceLocal_Send){ 
		sendValues.push_back(i); 
		// cout << i << endl;
	}
	
	vector<int> nFaceBoundaryFaceLocal_Recv(nBlocks*bound_num,0);

	MPI_Alltoallv( sendValues.data(),                  sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceBoundaryFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	sendValues.clear();
	for(auto& i : nStartBoundaryFaceLocal_Send){ 
		sendValues.push_back(i);
		
		// cout << i << endl;
	}
	
	vector<int> startFaceBoundaryFaceLocal_Recv(nBlocks*bound_num,0);
	

	MPI_Alltoallv( sendValues.data(),                      sendCounts.data(), sendDisps.data(), MPI_INT, 
				   startFaceBoundaryFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	
	// boundary name save
	vector<string> boundaryNames; 
	string test; 
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1) {
			

			// boundary.name.erase(std::find_if(boundary.name.rbegin(), boundary.name.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), boundary.name.end());
			// boundary.name.erase(boundary.name.begin(), std::find_if(boundary.name.begin(), boundary.name.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			// boundary.name.erase(remove(boundary.name.begin(), boundary.name.end(), ' '), boundary.name.end());
			boundaryNames.push_back(boundary.name);

			// strrrr.erase(std::find_if(strrrr.rbegin(), strrrr.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), strrrr.end());
			// strrrr.erase(strrrr.begin(), std::find_if(strrrr.begin(), strrrr.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// cout << boundary.name.length() << endl;
		}
	}
	// for(int i=0; i<boundaryNames.size(); ++i){
		// cout << boundaryNames[i] << boundaryNames[i].length() << endl;
	// }
	// for(auto& i : boundaryNames){
		// cout << i.data() << endl;
		
	// } 
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int ibcs=0; ibcs<bound_num; ++ibcs){
			int numm = ip*bound_num + ibcs;
			// if(rank==2) cout << rank << " " << ip << " " << ibcs << " " << nBlocks << endl;
			newMesh[ip].addBoundary();
			newMesh[ip].boundary.back().name = mesh.boundary[ibcs].name;
			newMesh[ip].boundary.back().nFaces = nFaceBoundaryFaceLocal_Recv[numm];
			newMesh[ip].boundary.back().startFace = startFaceBoundaryFaceLocal_Recv[numm];
			newMesh[ip].boundary.back().myProcNo = -1;
			newMesh[ip].boundary.back().neighbProcNo = -1;
			// // cout << nFaceBoundaryFaceLocal_Recv[numm] << endl;
			// if(rank==2) cout << rank << " " << numm << " " << newMesh[ip].points.size() << " " << bound_num << " " << startFaceBoundaryFaceLocal_Recv[numm] << " " << nFaceBoundaryFaceLocal_Recv[numm] << endl;
		}
	}
	
	
	//=======================
	
	
	// write of processor faces, PROCESSOR_FACE
	vector<int> nFaceProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> nStartProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> myProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> neighbProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	int temp_proc=0;
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal_Send[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp[i][j];
		}
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp2[i][j];
		}
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts_temp[i][j] > 0) {
				nFaceProcessorFaceLocal_Send[temp_proc] = sendCounts_temp[i][j];
				nStartProcessorFaceLocal_Send[temp_proc] = n;
				myProcNoProcessorFaceLocal_Send[temp_proc] = i;
				neighbProcNoProcessorFaceLocal_Send[temp_proc] = j;
				
				// cout << "========= " << sendCounts_temp[i][j] << " " << n << endl;
				
				n += sendCounts_temp[i][j];
			}
			else {
				nFaceProcessorFaceLocal_Send[temp_proc]=0;
				nStartProcessorFaceLocal_Send[temp_proc]=0;
				myProcNoProcessorFaceLocal_Send[temp_proc]=i;
				neighbProcNoProcessorFaceLocal_Send[temp_proc]=j;
			}
			++temp_proc;
		}
	}
	
	

	// write of processor faces, TO_BE_PROCESSOR_FACE
	vector<int> nFaceToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> nStartToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> myProcNoToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> neighbProcNoToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	temp_proc=0;
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal_Send[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp2[i][j];
		}
		// cout << n << endl;
		
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts_temp2[i][j] > 0) {
				nFaceToBeProcessorFaceLocal_Send[temp_proc] = sendCounts_temp2[i][j];
				nStartToBeProcessorFaceLocal_Send[temp_proc] = n;
				myProcNoToBeProcessorFaceLocal_Send[temp_proc] = i;
				neighbProcNoToBeProcessorFaceLocal_Send[temp_proc] = j;
				
				
				n += sendCounts_temp2[i][j];
			}
			else {
				nFaceToBeProcessorFaceLocal_Send[temp_proc]=0;
				nStartToBeProcessorFaceLocal_Send[temp_proc]=0;
				myProcNoToBeProcessorFaceLocal_Send[temp_proc]=i;
				neighbProcNoToBeProcessorFaceLocal_Send[temp_proc]=j;
			}
			++temp_proc;
		}
	}
	
	
	
	std::fill(sendCounts.begin(),sendCounts.end(),nBlocks);
	std::fill(recvCounts.begin(),recvCounts.end(),nBlocks);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];

	vector<int> nFaceProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nFaceProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> nStartProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nStartProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nStartProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> myProcNoProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( myProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   myProcNoProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> neighbProcNoProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( neighbProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   neighbProcNoProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);


	vector<vector<int>> ibcProcessorFace(nBlocks,vector<int>(0,0));
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=0; i<nBlocks; ++i){
			int numm = ip*nBlocks + i;
			if(nFaceProcessorFaceLocal_Recv[numm] > 0){
				
				ibcProcessorFace[ip].push_back( newMesh[ip].boundary.size() );
				
				newMesh[ip].addBoundary();
				string bcnames = "procBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[numm]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[numm]);
				newMesh[ip].boundary.back().name = bcnames;
				newMesh[ip].boundary.back().nFaces = nFaceProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().startFace = nStartProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().myProcNo = myProcNoProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().neighbProcNo = neighbProcNoProcessorFaceLocal_Recv[numm];
			}
		}
	}
	

	
	vector<int> nFaceToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nFaceToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> nStartToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nStartToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nStartToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> myProcNoToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( myProcNoToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   myProcNoToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> neighbProcNoToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( neighbProcNoToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   neighbProcNoToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);


	vector<vector<int>> ibcToBeProcessorFace(nBlocks,vector<int>(0,0));

	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=0; i<nBlocks; ++i){
			int numm = ip*nBlocks + i;
			if(nFaceToBeProcessorFaceLocal_Recv[numm] > 0){
				
				ibcToBeProcessorFace[ip].push_back( newMesh[ip].boundary.size() );
				
				newMesh[ip].addBoundary();
				string bcnames = "toBeProcBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[numm]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[numm]);
				newMesh[ip].boundary.back().name = bcnames;
				
				// cout << newMesh[ip].boundary.back().name << endl;
				
				newMesh[ip].boundary.back().nFaces = nFaceToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().startFace = nStartToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().myProcNo = myProcNoToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().neighbProcNo = neighbProcNoToBeProcessorFaceLocal_Recv[numm];
			}
		}
	}
	
	
	// startFace re calc.
	
	for(int ip=0; ip<nBlocks; ++ip){
		int nbc=0;
		int startProcsFace=0;
		for(auto& ibc : newMesh[ip].boundary){
			if(ibc.neighbProcNo == -1) {
				if(ibc.nFaces>0) startProcsFace = ibc.startFace + ibc.nFaces;
				++nbc;
			}
		}
		
		
		for(int ibc=nbc; ibc<newMesh[ip].boundary.size(); ++ibc){
			int ostart = newMesh[ip].boundary[ibc-1].startFace;
			int onface = newMesh[ip].boundary[ibc-1].nFaces;
			int nstart = ostart + onface;
			if(ibc==nbc) nstart = startProcsFace;
			newMesh[ip].boundary[ibc].startFace = nstart;
		}
	}
	
	
	
	
	
	
	
	
	//==========================================
	
	
	// push face points
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idFacePoint_Recv[ip]; i<str_idFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idFacePoint_Recv[i] );
			}
		}
	}
	

	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idBoundaryFacePoint_Recv[ip]; i<str_idBoundaryFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idBoundaryFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idBoundaryFacePoint_Recv[i] );
				// cout << num1 << " " << idBoundaryFacePoint_Recv[i] << endl;
			}
		}
	}


	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idOriginProcessorFacePoint_Recv[ip]; i<str_idOriginProcessorFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idOriginProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idOriginProcessorFacePoint_Recv[i] );
			}
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idProcessorFacePoint_Recv[ip]; i<str_idProcessorFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idProcessorFacePoint_Recv[i] );
				// cout << num1 << " " << idProcessorFacePoint_Recv[i] << endl;
			}
		}
	}

// if(rank==1){
	// cout << rank << " " << newMesh[0].faces.size() << endl;
	// cout << rank << " " << newMesh[1].faces.size() << endl;
	// cout << rank << " " << newMesh[0].boundary[0].nFaces << " " <<  newMesh[1].boundary[0].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[0].startFace << " " <<  newMesh[1].boundary[0].startFace << endl;
	
	// cout << rank << " " << newMesh[0].boundary[1].nFaces << " " <<  newMesh[1].boundary[1].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[1].startFace << " " <<  newMesh[1].boundary[1].startFace << endl;
	
	// cout << rank << " " << newMesh[0].boundary[2].nFaces << " " <<  newMesh[1].boundary[2].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[2].startFace << " " <<  newMesh[1].boundary[2].startFace << endl;

// }
	//==========================================
	
	
	// push owner 
	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idOwnerLocal_Recv[ip]; i<str_idOwnerLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].owner = idOwnerLocal_Recv[i];
			
			// if(idOwnerLocal_Recv[i]>100000000)  cout<< idOwnerLocal_Recv[i] << endl;
			// cout << ip << " " << i << " " << idOwnerLocal_Recv[i] << endl;
			++numm;
		}
	}
	
	
	// push neighbour 
	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idNeighbourLocal_Recv[ip]; i<str_idNeighbourLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].neighbour = idNeighbourLocal_Recv[i];
			// cout << ip << " " << i << " " << idNeighbourLocal_Recv[i] << endl;
			++numm;
		}
	}
	
	
	// cout << newMesh[1].points.size() << 
	// " " << newMesh[1].faces.size() <<
	// " " <<	newMesh[1].boundary.size() <<
	// " " <<	newMesh[1].boundary[0].myProcNo <<
	// " " <<	newMesh[1].boundary[1].myProcNo <<
	// " " <<	newMesh[1].boundary[2].myProcNo  << endl;
	// // cout << newMesh[1].boundary[2].name  << endl;
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idOwnerLocal_Recv[ip]; i<str_idOwnerLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].owner = idOwnerLocal_Recv[i];
			if( newMesh[ip].faces[numm].neighbour == -1 ){
				newMesh[ip].faces[numm].setType(SEMO_Types::BOUNDARY_FACE);
			}
			else{
				newMesh[ip].faces[numm].setType(SEMO_Types::INTERNAL_FACE);
			}
			++numm;
		}
		
		// newMesh[ip].check();
		
		// newMesh[ip].buildCells();
		
		// newMesh[ip].setFaceTypes();
	
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// if(rank==1){
	for(int ip=0; ip<size; ++ip){
		// int ip=1;
		
		// for(auto& i : newMesh[ip].faces){
			// cout << endl;
			// for(auto& j : i.points){
				// cout << " face's points : " << j << endl;
			// }
		// }
		// for(auto& i : newMesh[ip].faces){
			// cout << " face's owner : " << i.owner << endl;
		// }
		// for(auto& i : newMesh[ip].faces){
			// cout << " face's neighbour : " << i.neighbour << endl;
		// }
		if(rank==2){
			int ibc = 0;
			for(auto& i : newMesh[ip].boundary){
				cout << ip << " " << ibc << endl;
				cout << " boundary name : " << i.name << endl;
				cout << " boundary startFace : " << i.startFace << endl;
				cout << " boundary nFaces : " << i.nFaces << endl;
				cout << " boundary neighbProcNo : " << i.neighbProcNo << endl;
				++ibc;
			}
		}
		
		newMesh[ip].check();
		
		newMesh[ip].buildCells();
		
		newMesh[ip].setFaceTypes();
		
		// create list
		newMesh[ip].buildLists();
		
		// check list
		// newMesh[ip].checkLists();
		
		// cell's faces connection
		newMesh[ip].connectCelltoFaces();
		
		// cell's points connection
		newMesh[ip].connectCelltoPoints();
		
		// set processor face counts
		newMesh[ip].setCountsProcFaces();
		
		// set processor face displacements
		newMesh[ip].setDisplsProcFaces(); 

		// newMesh[ip].saveFile("vtu");
		
		// newMesh[ip].checkQualities(); 
	
	}
	
	// if(rank==2) newMesh[3].saveFile("vtu");
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	//==========================================
	
	// if(rank==0) newMesh[0].checkQualities();
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	
	
	

	//==========================================
	// connection new meshes
	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute connection partitioned meshes ... ";
	
	vector<int> startToBeInternalFace_Send(nBlocks,0);
	vector<vector<int>> idToBeInternalFace(nBlocks,vector<int>(0,0));
	vector<vector<int>> neighbMeshNoToBeInternalFace(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::INTERNAL_FACE ){
			int j = face.owner;
			int i = idBlockCell[j];
			++startToBeInternalFace_Send[i];
		}
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			++startToBeInternalFace_Send[i];
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			int k = idBlockCell_Recv[face.neighbour];
			if(i == k){
				idToBeInternalFace[i].push_back( startToBeInternalFace_Send[i] );
				neighbMeshNoToBeInternalFace[i].push_back( rank_Recv[face.neighbour] );
			}
			++startToBeInternalFace_Send[i];
		}
	}
	
	
	
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	sendValues.clear();
	sendValues.resize(nBlocks,0);
	recvValues.clear();
	recvValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idToBeInternalFace[i].size();
	
	
	mpi.exchangeDatas(sendCounts, recvCounts, sendValues, recvValues);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = sendValues[i];
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	
	sendValues.clear();
	recvValues.clear();
	
	vector<int> idToBeInternalFace_Send;
	vector<int> idToBeInternalFace_Recv;
	for(int ip=0; ip<nBlocks; ++ip) {for(auto& i : idToBeInternalFace[ip]) { idToBeInternalFace_Send.push_back(i); } } 

	mpi.exchangeDatas(sendCounts, recvCounts, idToBeInternalFace_Send, idToBeInternalFace_Recv);
			
			
	vector<int> neighbMeshNoToBeInternalFace_Send;
	vector<int> neighbMeshNoToBeInternalFace_Recv;	 
	for(int ip=0; ip<nBlocks; ++ip) {for(auto& i : neighbMeshNoToBeInternalFace[ip]) { neighbMeshNoToBeInternalFace_Send.push_back(i); } }  
	
	mpi.exchangeDatas(sendCounts, recvCounts, neighbMeshNoToBeInternalFace_Send, neighbMeshNoToBeInternalFace_Recv);
	
	
	vector<int> str_idToBeInternalFace_Recv(nBlocks+1,0);
	for(int i=1; i<nBlocks+1; ++i) 
		str_idToBeInternalFace_Recv[i] = str_idToBeInternalFace_Recv[i-1] + recvCounts[i-1];
	
	
	
	vector<int> nCellsEach(nBlocks,0);
	vector<int> startCells(nBlocks+1,0);
	
	for(int ip=0; ip<nBlocks; ++ip){
		int maxCells=-1;
		for(auto& face : newMesh[ip].faces){
			maxCells = max(maxCells, face.owner);
		}
		if(maxCells==-1){
			nCellsEach[ip] = 0; 
		}
		else{
			nCellsEach[ip] = maxCells+1; 
		}
	}
	
	for(int ip=1; ip<nBlocks+1; ++ip) startCells[ip] = startCells[ip-1] + nCellsEach[ip-1];
	
	
	// # of physical boundary  
	int nnbcs = 0;
	for(int ibcs=0; ibcs<newMesh[0].boundary.size(); ++ibcs){
		if(newMesh[0].boundary[ibcs].myProcNo == -1){
			++nnbcs;
		}
	}
	vector<vector<int>> temp_bcfacepoint(nnbcs,vector<int>(0,0));
	vector<vector<int>> temp_procfacepoint(nBlocks,vector<int>(0,0));
	vector<vector<int>> temp_Toprocfacepoint(nBlocks,vector<int>(0,0));
	
	vector<vector<int>> temp_bcfacepointOwner(nnbcs,vector<int>(0,0));
	vector<vector<int>> temp_procfacepointOwner(nBlocks,vector<int>(0,0));
	vector<vector<int>> temp_ToprocfacepointOwner(nBlocks,vector<int>(0,0));
	
	
	// delete original mesh
	mesh.points.clear();
	mesh.faces.clear();
	mesh.cells.clear();
	mesh.boundary.clear();
	
	
	
	// point id reset
	
	
	if(rank==2){
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=0; i<nBlocks; ++i){
			int numm = ip*nBlocks + i;
			if(nFaceToBeProcessorFaceLocal_Recv[numm] > 0){
				cout << nFaceToBeProcessorFaceLocal_Recv[numm] << " " << ip << " " << i << endl;
			}
		}
	}
	}



	// vector<int> faceId_Send;
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// faceId_Send.push_back(rank);
		// }
	// }
	// vector<int> faceId_Recv(faceId_Send.size(),0);
	// MPI_Alltoallv( faceId_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // faceId_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // MPI_COMM_WORLD);



	// vector<vector<int>> localIdFaces; 
	// vector<vector<int>> localProcFaces; 
	// for(int ip=0; ip<nBlocks; ++ip){
		// for(auto& boundary : newMesh[ip].boundary){
			// if(boundary.neighbProcNo != -1){
				
			// }
		// }
	// }



	vector<int> localIdPoints;
	vector<int> localProcPoints;
	vector<int> strGlobalPoints(nBlocks,0);
	vector<vector<int>> globalIdPoints(nBlocks,vector<int>(0,0));
	vector<vector<bool>> globalErasePoints(nBlocks,vector<bool>(0,false));
	int idGlobalPoints=0;
	for(int ip=0; ip<nBlocks; ++ip){
	// if(rank==2){
		// cout << "AAAAAAAAAAA" << endl;
		// cout << newMesh[ip].boundary[1].startFace << " " << newMesh[ip].boundary[1].nFaces << endl;
		// cout << newMesh[ip].boundary[2].startFace << " " << newMesh[ip].boundary[2].nFaces << endl;
		// cout << newMesh[ip].boundary[3].startFace << " " << newMesh[ip].boundary[3].nFaces << endl;
		// cout << newMesh[ip].boundary[4].startFace << " " << newMesh[ip].boundary[4].nFaces << endl;
	// }
		vector<int> replaceIdPoints(newMesh[ip].points.size(),-1);
		vector<int> replaceProcPoints(newMesh[ip].points.size(),-1);
		
		vector<bool> executePoints(newMesh[ip].points.size(),true);
		for(int i=str_idToBeInternalFace_Recv[ip]; i<str_idToBeInternalFace_Recv[ip+1]; ++i){
			int j = idToBeInternalFace_Recv[i]; // id
			int ipNgb = neighbMeshNoToBeInternalFace_Recv[i]; // neighbour face's Processor
			if(ip>ipNgb){
				for(auto& point : newMesh[ip].faces[j].points){
					executePoints[point] = false;
				}
			}
		}
		
		for(int jp=0; jp<ip; ++jp){
			int n=0;
			for(auto& ipoint : newMesh[ip].points){
				int nn=0;
				for(auto& jpoint : newMesh[jp].points){
					// int cx = utility.CompareDoubleAbsoulteAndUlps(ipoint.x, jpoint.x, 1.0e-8, 4);
					// int cy = utility.CompareDoubleAbsoulteAndUlps(ipoint.y, jpoint.y, 1.0e-8, 4);
					// int cz = utility.CompareDoubleAbsoulteAndUlps(ipoint.z, jpoint.z, 1.0e-8, 4); 
					// if(cx==0 && cy==0 && cz==0) executePoints[n] = false;
					// double aaa;
					bool cx = utility.approximatelyEqualAbsRel(ipoint.x, jpoint.x, 1.e-12, 1.e-8);
					bool cy = utility.approximatelyEqualAbsRel(ipoint.y, jpoint.y, 1.e-12, 1.e-8);
					bool cz = utility.approximatelyEqualAbsRel(ipoint.z, jpoint.z, 1.e-12, 1.e-8);
					if(cx && cy && cz) {
						replaceIdPoints[n] = nn;
						replaceProcPoints[n] = jp;
						executePoints[n] = false; 
						break;
					}
					++nn;
					// cout << cx << " " << cy << " " << cz << endl;
				}
				++n;
			}
		}
		
		
		
		for(int point=0; point<newMesh[ip].points.size(); ++point){
			if(executePoints[point]==true){
				localIdPoints.push_back(point);
				localProcPoints.push_back(ip);
				
				globalErasePoints[ip].push_back(false);
				globalIdPoints[ip].push_back(localIdPoints.size()-1);
			}
			else{
				int id = replaceIdPoints[point];
				int proc = replaceProcPoints[point];
				
				globalErasePoints[ip].push_back(true);
				globalIdPoints[ip].push_back(strGlobalPoints[proc] + id);
			}
		}
		
		strGlobalPoints[ip] = idGlobalPoints;
		idGlobalPoints = localIdPoints.size();
	}
	
	
	// if(rank==2){
		
	// for(int ip=0; ip<nBlocks; ++ip){
		// for(auto& i : globalIdPoints[ip]){
			// cout << i << endl;
		// }
	// }
		
		
		// for(int i=0; i<localIdPoints.size(); ++i){
			// // cout << localIdPoints[i] << " " << localProcPoints[i] << endl;
			// int proc = localProcPoints[i];
			// int id = localIdPoints[i];
			// cout << newMesh[proc].points[id].x << " " << newMesh[proc].points[id].y << " " << newMesh[proc].points[id].z << endl;
		// }
		

		// for(int ip=0; ip<nBlocks; ++ip){
			// for(int i=0; i<newMesh[ip].points.size(); ++i){
				// if(globalErasePoints[ip][i] == false){
						
					// // mesh.addPoint();
					// // mesh.points.back().x = newMesh[ip].points[i].x;
					// // mesh.points.back().y = newMesh[ip].points[i].y;
					// // mesh.points.back().z = newMesh[ip].points[i].z;
					// // replacePointsId[i] = globalIdPoints[ip][i];
					// cout << globalIdPoints[ip][i] << endl;
					
				// }
				// else{
					// // replacePointsId[i] = globalIdPoints[ip][i];
				// }
			// }
		// }
		
	// }
	


	vector<vector<int>> erasePoints(nBlocks,vector<int>(0,0));
	// vector<vector<int>> eraseFaces(nBlocks,vector<int>(0,0));
	vector<vector<int>> replacePoints(nBlocks,vector<int>(0,0));
	vector<vector<int>> replacePointsRank(nBlocks,vector<int>(0,0));
	
	vector<int> iii_temp(nBlocks,0);
	for(int ip=0; ip<nBlocks; ++ip){
		
		// points
		vector<int> replacePointsId(newMesh[ip].points.size(),0);
		// int iipoint = 0;
		for(int i=0; i<newMesh[ip].points.size(); ++i){
			
			// cout << globalErasePoints[ip][i] << " " << (std::find(erasePoints[ip].begin(), erasePoints[ip].end(), i) == erasePoints[ip].end()) << endl;
			if(std::find(erasePoints[ip].begin(), erasePoints[ip].end(), i) 
				== erasePoints[ip].end()){
			// if(globalErasePoints[ip][i] == false){
					
				mesh.addPoint();
				mesh.points.back().x = newMesh[ip].points[i].x;
				mesh.points.back().y = newMesh[ip].points[i].y;
				mesh.points.back().z = newMesh[ip].points[i].z;
				replacePointsId[i] = mesh.points.size()-1;
				
			}
			else{
				int l = std::find(erasePoints[ip].begin(), erasePoints[ip].end(), i) 
				        - erasePoints[ip].begin();
				if(replacePointsRank[ip][l] >= ip){
					cout << " ERROR " << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				replacePointsId[i] = replacePoints[ip][l];
			}
			// if(globalErasePoints[ip][i] == false){
					
				// mesh.addPoint();
				// mesh.points.back().x = newMesh[ip].points[i].x;
				// mesh.points.back().y = newMesh[ip].points[i].y;
				// mesh.points.back().z = newMesh[ip].points[i].z;
				// replacePointsId[i] = globalIdPoints[ip][i];
				// // if(rank==2) cout << replacePointsId[i] << endl;
				
			// }
			// else{
				// replacePointsId[i] = globalIdPoints[ip][i];
			// }
		}
		
		for(int i=str_idToBeInternalFace_Recv[ip]; i<str_idToBeInternalFace_Recv[ip+1]; ++i){
			int j = idToBeInternalFace_Recv[i]; // id
			int ipNgb = neighbMeshNoToBeInternalFace_Recv[i]; // neighbour face's Processor
			if( ip < ipNgb ) {
				int ii = str_idToBeInternalFace_Recv[ipNgb] + iii_temp[ipNgb];
				int jj = idToBeInternalFace_Recv[ii];
				
				// int owner = newMesh[ip].faces[j].owner;
				int ownerNgb = newMesh[ipNgb].faces[jj].owner;
				
				// if(owner < 0 || ownerNgb < 0 ){
					// cout << "#ERROR : owner < 0 || ownerNgb < 0" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				newMesh[ip].faces[j].neighbour = startCells[ipNgb] - startCells[ip] + ownerNgb;
				
				
				int opoint = 0;
				for(auto& point : newMesh[ipNgb].faces[jj].points){
					// there is not 'point' at erasePoints
					if(std::find(erasePoints[ipNgb].begin(), erasePoints[ipNgb].end(), point) 
						== erasePoints[ipNgb].end()){
					// if(globalErasePoints[ipNgb][point] == false){
							
						// reverse connection for each other
						int reverse_opoint = newMesh[ip].faces[j].points.size()-1 - opoint;
						int addIPoint = newMesh[ip].faces[j].points[ reverse_opoint ];
						int addIPointReplace = replacePointsId[ addIPoint ];
						// int addIPointReplace = globalIdPoints[ip][ addIPoint ];
						
						erasePoints[ipNgb].push_back(point);
						replacePoints[ipNgb].push_back( addIPointReplace );
						replacePointsRank[ipNgb].push_back(ip);
					}
					++opoint;
				}
				// eraseFaces[ipNgb].push_back(jj);
				newMesh[ip].faces[j].setType(SEMO_Types::INTERNAL_FACE);
				newMesh[ipNgb].faces[jj].setType(SEMO_Types::TO_BE_DELETE_FACE);
				++iii_temp[ipNgb];
				
				
			}
		}
		
			
			
		// faces , faces points
		for(auto& face : newMesh[ip].faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				mesh.addFace();
				for(auto& point : face.points){
					mesh.faces.back().points.push_back( replacePointsId[point] );
					// mesh.faces.back().points.push_back( globalIdPoints[ip][replacePointsId[point]] );
				}
				
				// owner
				mesh.faces.back().owner = startCells[ip] + face.owner;
				
				// neighbour
				mesh.faces.back().neighbour = startCells[ip] + face.neighbour;
			}
		}
		
		
		//boundary
		int nbcs = 0;
		for(int ibcs=0; ibcs<newMesh[ip].boundary.size(); ++ibcs){
			if(newMesh[ip].boundary[ibcs].myProcNo == -1){
				// mesh.addBoundary();
				++nbcs;
			}
		}
		// cout << "aaaa" << nbcs << endl;
		for(int ibcs=0; ibcs<nbcs; ++ibcs){
			int strF = newMesh[ip].boundary[ibcs].startFace;
			int endF = strF + newMesh[ip].boundary[ibcs].nFaces;
			for(int i=strF; i<endF; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
					temp_bcfacepointOwner[ibcs].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				temp_bcfacepoint[ibcs].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_bcfacepoint[ibcs].push_back( replacePointsId[point] );
				}
				
			}
		}
		
	

		//procs
		for(auto& ibcs : ibcProcessorFace[ip]){
			
			int strf = newMesh[ip].boundary[ibcs].startFace;
			int endf = strf + newMesh[ip].boundary[ibcs].nFaces;
			int negh = newMesh[ip].boundary[ibcs].neighbProcNo;
			for(int i=strf; i<endf; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
				temp_procfacepointOwner[negh].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				temp_procfacepoint[negh].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_procfacepoint[negh].push_back( replacePointsId[point] );
					// if(rank==2) cout << replacePointsId[point] << endl;
					// if(rank==1) cout << replacePointsId[point] << endl;
				}
			}
		}
		
		for(auto& ibcs : ibcToBeProcessorFace[ip]){
			int strf = newMesh[ip].boundary[ibcs].startFace;
			int endf = strf + newMesh[ip].boundary[ibcs].nFaces;
			int negh = newMesh[ip].boundary[ibcs].neighbProcNo;
	// if(rank==0) cout << "aaa " << strf << " " << endf <<  endl;
			for(int i=strf; i<endf; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
				temp_ToprocfacepointOwner[negh].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				// if(rank==0) cout << startCells[ip] << endl;
				
				temp_Toprocfacepoint[negh].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_Toprocfacepoint[negh].push_back( replacePointsId[point] );
					// if(rank==2) cout << replacePointsId[point] << endl;
				}
			}
		}
	// if(rank==2) {
		// for(auto& i : replacePointsId){
			// cout << ip << " " << i << endl;
		// }
	// }
		
	}

// if(rank==2){
	// for(int ip=0; ip<nBlocks; ++ip){
		
		// for(int i=0; i<newMesh[ip].points.size(); ++i){
			
			// cout << globalErasePoints[ip][i] << " " << (std::find(erasePoints[ip].begin(), erasePoints[ip].end(), i) == erasePoints[ip].end()) << endl;
		// }
	// }
// }
	
	// for(auto& i : boundaryNames){
	for(int i=0; i<boundaryNames.size(); ++i){
		mesh.addBoundary();
		
		mesh.boundary.back().name = boundaryNames[i];
	}
	int nbcs = boundaryNames.size();
	
	
	
		
	for(int ibcs=0; ibcs<nbcs; ++ibcs){
		
		mesh.boundary[ibcs].startFace = mesh.faces.size();
		int nFaces = 0;
		for(int i=0; i<temp_bcfacepoint[ibcs].size(); ++i){
			mesh.addFace();
			
			
			int nnn = temp_bcfacepoint[ibcs][i];
			for(int j=0; j<nnn; ++j){
				++i;
				mesh.faces.back().points.push_back( temp_bcfacepoint[ibcs][i] );
			}
			++nFaces;
			
			
		}
		mesh.boundary[ibcs].nFaces = nFaces;
		mesh.boundary[ibcs].myProcNo = -1;
		mesh.boundary[ibcs].neighbProcNo = -1;
		
		
		// owner
		int iii=0;
		for(int i=mesh.faces.size()-temp_bcfacepointOwner[ibcs].size(); 
		i<mesh.faces.size(); ++i){
			mesh.faces[i].owner = temp_bcfacepointOwner[ibcs][iii];
			++iii;
		}
		
	}
	

	// processor face
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(auto& i : temp_Toprocfacepoint[ip]){
			temp_procfacepoint[ip].push_back( i );
		}
		temp_Toprocfacepoint[ip].clear();
		
		for(auto& i : temp_ToprocfacepointOwner[ip]){
			temp_procfacepointOwner[ip].push_back( i );
		}
		temp_ToprocfacepointOwner[ip].clear();
	}
	
	
	// if(rank==1){
		// cout << endl;
		// for(auto& i : temp_procfacepoint[0]){
			// cout << i << endl;
		// }
	// }
	
	
	
	int iii_proc=0;
	for(int ip=0; ip<nBlocks; ++ip){
		
		// if(ip==rank) continue;
		
		if(temp_procfacepoint[ip].size()>0){
			
			mesh.addBoundary();
		
			mesh.boundary[nbcs+iii_proc].startFace = mesh.faces.size();
			int nFaces = 0;
			for(int i=0; i<temp_procfacepoint[ip].size(); ++i){
				mesh.addFace();
				int nnn = temp_procfacepoint[ip][i];
				for(int j=0; j<nnn; ++j){
					++i;
					// if(rank==2) cout << temp_procfacepoint[ip][i] << endl;
					mesh.faces.back().points.push_back( temp_procfacepoint[ip][i] );
				}
				// if(rank==2) cout << endl;
				
				// std::reverse(mesh.faces.back().points.begin(), mesh.faces.back().points.end());
				
				++nFaces;
			}
			
			string bcnames = "procBoundary" + to_string(rank) + "to" + to_string(ip);
			
			mesh.boundary[nbcs+iii_proc].name = bcnames;
			mesh.boundary[nbcs+iii_proc].nFaces = nFaces;
			mesh.boundary[nbcs+iii_proc].myProcNo = rank;
			mesh.boundary[nbcs+iii_proc].neighbProcNo = ip;
			
			// // owner
			// mesh.faces.back().owner = 
			
			++iii_proc;
			
			// owner
			int iii=0;
			for(int i=mesh.faces.size()-temp_procfacepointOwner[ip].size(); 
			i<mesh.faces.size(); ++i){
				mesh.faces[i].owner = temp_procfacepointOwner[ip][iii];
				++iii;
			}
			
		}
	}
	
	//==========================================
	

	
	// startFace re calc.
	
	for(int ip=0; ip<nBlocks; ++ip){
		int nbc=0;
		for(auto& ibc : mesh.boundary){
			if(ibc.neighbProcNo == -1) ++nbc;
		}
		
		for(int ibc=nbc; ibc<mesh.boundary.size(); ++ibc){
			int ostart = mesh.boundary[ibc-1].startFace;
			int onface = mesh.boundary[ibc-1].nFaces;
			int nstart = ostart + onface;
			mesh.boundary[ibc].startFace = nstart;
		}
	}
	
	
	

	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	
	
	

	// if(rank==1){
		// int ip=0;
		
		// for(auto& i : temp_procfacepoint[ip]){
			// cout << i << endl;
		// }
		
		
		// int face = 0;
		// for(auto& i : mesh.faces){
			// cout << " face : " << face << endl;
			// for(auto& j : i.points){
				// cout << " face's points : " << j << endl;
			// }
			// ++face;
		// }
		// for(auto& i : mesh.faces){
			// cout << " face's owner : " << i.owner << endl;
		// }
		// for(auto& i : mesh.faces){
			// cout << " face's neighbour : " << i.neighbour << endl;
		// }
		// for(auto& i : mesh.boundary){
			// cout << " boundary name : " << i.name << endl;
			// cout << " boundary startFace : " << i.startFace << endl;
			// cout << " boundary nFaces : " << i.nFaces << endl;
			// cout << " boundary neighbProcNo : " << i.neighbProcNo << endl;
		// }
	// }
	
		newMesh.clear();

		mesh.check();
		
		mesh.buildCells();
		
		mesh.setFaceTypes();
		
		// create list
		mesh.buildLists();
		
		// check list
		// mesh.checkLists();
		
		// cell's faces connection
		mesh.connectCelltoFaces();
		
		// cell's points connection
		mesh.connectCelltoPoints();
		
		// set processor face counts
		mesh.setCountsProcFaces();
		
		// set processor face displacements
		mesh.setDisplsProcFaces(); 

		mesh.saveFile("vtu");
	
		mesh.informations();
	
	// for(auto& i : mesh.faces){
		// cout << i.points.size() << endl;
	// }
	

	
	
		
	

	mesh.checkQualities();
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	// if(rank==3){
		// for(auto& i : newMesh[3].points){
			// cout << i.x << " " << i.y << " " << i.z << endl;
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// if(rank==0){
	// // for(int ip=3; ip<4; ++ip){
		// int ip = 1;
		
		// newMesh[ip].check();
		
		// newMesh[ip].buildCells();
		
		// newMesh[ip].setFaceTypes();
		
		// // create list
		// newMesh[ip].buildLists();
		
		// // check list
		// newMesh[ip].checkLists();
		
		// // cell's faces connection
		// newMesh[ip].connectCelltoFaces();
		
		// // cell's points connection
		// newMesh[ip].connectCelltoPoints();
		
		// // set processor face counts
		// newMesh[ip].setCountsProcFaces();
		
		// // set processor face displacements
		// newMesh[ip].setDisplsProcFaces(); 

		// newMesh[ip].saveFile("vtu");
	
	// // }
	// }
		

		
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// //======================================================================
	
	
	// // mesh.saveFile("vtu");
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}


// void SEMO_Mesh_Builder::partitionOneToAll(){
	
	
	
	
	
	
	
	
	
	
	
	
	
// }



void SEMO_Mesh_Builder::partitionInit(string type){
	
	int ncells = (*this).cells.size();
	int idBlockCell[ncells];
	
	
	int nBlocks = (*this).mpi.getSize();
	
	// if(type=="METIS"){
		SEMO_Mesh_Builder::parMETIS(nBlocks, idBlockCell);
	// }


	SEMO_Mesh_Builder::distribute(nBlocks, idBlockCell);
	
}




void SEMO_Mesh_Builder::partition(string type){
	
	int ncells = (*this).cells.size();
	int idBlockCell[ncells];
	
	
	int nBlocks = 2;
	
	// if(type=="METIS"){
		SEMO_Mesh_Builder::parMETIS(nBlocks, idBlockCell);
	// }


	SEMO_Mesh_Builder::distribute(nBlocks, idBlockCell);
	
}




void SEMO_Mesh_Builder::parMETIS(int nBlocks, int idBlockCell[]){
		
		
	SEMO_Mesh_Builder& mesh = *this;
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_DBGLVL]=0;
	options[0] = 0;
	
	// int idBlockCell[ncells];
	// int dummy1[npoints];

	// cout << ncells << endl;
	

	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	// cout << "111" << endl;
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	// tpwgts[0]=0.4;
	// tpwgts[1]=0.9;
	
	real_t ubvec[ncon];
	std::fill_n(ubvec, ncon, 1.05);
	
	
	
	int temp_vtxdist[nBlocks];
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	int vtxdist[nBlocks+1];
	vtxdist[0] = 0;
	for(int i=1; i<nBlocks+1; ++i){
		vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	}
	
	// int str_Cell = 0;
	// for(int i=0; i<rank; ++i){
		// str_Cell += 
	// }
	
	int xadj[ncells+1];
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		for(auto& face : cell.faces){
			if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE){
				++numt2;
			}
			
			
			else if(mesh.faces[face].getType() == SEMO_Types::PROCESSOR_FACE) {
				++numt2;
			}
		}
		++numt;
		xadj[numt] = xadj[numt-1] + numt2;
	}
	
	
	
	
	
	
	int nSize = mesh.displsProcFaces[nBlocks-1] + mesh.countsProcFaces[nBlocks-1];
	vector<int> idCell_Send;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			idCell_Send.push_back(face.owner);
		}
	}
	// cout << nSize << endl;
	// int nnn = 0;
	// for(auto& i : mesh.boundary){
		// ++nnn;
		// cout << nnn << " " << i.neighbProcNo << endl;
	// }
	
	// cout << nSize << " " << idCell_Send.size() << 
	// " " << mesh.displsProcFaces[0] << 
	// " " << mesh.displsProcFaces[1] << 
	// " " << mesh.displsProcFaces[2] << 
	// " " << mesh.displsProcFaces[3] << 
	// " " << mesh.countsProcFaces[0] << 
	// " " << mesh.countsProcFaces[1] << 
	// " " << mesh.countsProcFaces[2] << 
	// " " << mesh.countsProcFaces[3] << 
	// endl;
	
	
	// if(nSize==0) nSize = 1;
	vector<int> idCell_Recv(nSize,0);
	MPI_Alltoallv( idCell_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   idCell_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
				   
				   
	int rank = mesh.mpi.getRank();
	
	vector<int> iRank_Send;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			iRank_Send.push_back(rank);
		}
	}
	// if(nSize==0) nSize = 1;
	vector<int> iRank_Recv(nSize,0);;
	MPI_Alltoallv( iRank_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   iRank_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	int adjncy[ xadj[ncells] ];
	
	
	numt = 0;
	
	
	// cout << xadj[ncells] << endl;
	// cout << idCell_Recv.size() << endl;
	// cout << iRank_Recv.size() << endl;
	
	int num_ncell[ncells];
	std::fill_n(num_ncell, ncells, 0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			int tnum = xadj[face.owner] + num_ncell[face.owner];
			adjncy[tnum] = vtxdist[rank] + face.neighbour;
			++num_ncell[face.owner];
			
			tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			adjncy[tnum] = vtxdist[rank] + face.owner;
			++num_ncell[face.neighbour];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			int tnum = xadj[face.owner] + num_ncell[face.owner];
			int rank_neighbour = iRank_Recv[proc_num];
			// cout << proc_num << " " << iRank_Recv.size() << endl;
			adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			++num_ncell[face.owner];
			
			++proc_num;
		}
	}
	
	
	// for(auto& cell : mesh.cells){
		// for(auto& face : cell.faces){
			// if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE) {
				// if( &cell == &mesh.cells[ mesh.faces[face].owner ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].neighbour;
				// }
				// else if( &cell == &mesh.cells[ mesh.faces[face].neighbour ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].owner;
				// }
				// else {
					// cerr << "#error, not matching cell != face owner or neighbour" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				// ++numt;
			// }
			
			
			// else if(mesh.faces[face].getType() == SEMO_Types::PROCESSOR_FACE) {
				
				// if( &cell == &mesh.cells[ mesh.faces[face].owner ] ){
					// // int idNeighbour = mesh.faces[face].neighbour;
					// int rank_neighbour = idCell_Recv[i];
					
					// adjncy[numt] = vtxdist[rank_neighbour] + idNeighbour;
					
					
				// }
				// else if( &cell == &mesh.cells[ mesh.faces[face].neighbour ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].owner;
				// }
				// else {
					// cerr << "#error, not matching cell != face owner or neighbour" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				
				// ++numt;
			// }
			
			
		// }
	// }
	


	ParMETIS_V3_PartKway(
		vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
		&ncon, &nBlocks, tpwgts, ubvec,
		options, &objval, idBlockCell, &comm);

	
}







void SEMO_Mesh_Builder::distribute(int nBlocks, int idBlockCell[]){
	
	SEMO_MPI_Builder mpi;
	mpi.setRank();
	mpi.setSize();
	
	SEMO_Mesh_Builder& mesh = *this;
	
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0)); // point block id (copies)
	int nCellsLocal[nBlocks]; // local total cells
	int idCellLocal[mesh.cells.size()]; // local cell id
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
	
	int nTotalLocalPointSize = 0; // total all # of local point
	for(auto& i : idBlockPoint) 
		nTotalLocalPointSize += i.size();

	// point distribution (CSR format)
	//
	// / gP : global points / lP : local points / BL : block id /
	//
	//     gP[0]     gP[1]      gP[2]      gP[3]   gP[4] ...
	// - - - - - - - - - - - - - - - - - - - - - - - - - ....
	// |  lP[0,5]  |lP[6,9]| lP[10,15] |  ....   |     |
	// |  BL[0,5]  |BL[6,9]| BL[10,15] |  ....   |     |
	// |           |       |           |         |     |
	// strPoints[0]|  strPoints[2]  strPoints[3] |  strPoints[5]  ...
	//         strPoints[1]                  strPoints[4]
	//
	int nPointsLocal[nBlocks]; // local total points
	std::fill_n(nPointsLocal, nBlocks, 0);
	int idPointLocal[nTotalLocalPointSize]; // local point id
	int idBlockPointLocal[nTotalLocalPointSize]; // local point block id
	int strPoints[mesh.points.size()+1]; // start of each global point
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(int j=0; j<idBlockPoint[i].size(); ++j){
			int idBlock = idBlockPoint[i][j];
			idPointLocal[nIndex] = nPointsLocal[ idBlock ];
			++nPointsLocal[ idBlock ];
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
	
	
	
	// //==========================================
	// // @@@ at MPI
	// int* idBlockCell_Send;
	// int* idBlockCell_Recv;
	// if(mpi.getSize() > 0){
		// // MPI_Allgather(&ncells, 1, MPI_INT, myvtxdist, 1, MPI_INT, MPI_COMM_WORLD);
		// // int numb_RecvNode[size];
		// // MPI_Alltoall(&numb_SendNode, 1, MPI_INT, numb_RecvNode, 1, MPI_INT, MPI_COMM_WORLD);
		// idBlockCell_Send = new int[mesh.mpi.nSend];
		// idBlockCell_Recv = new int[mesh.mpi.nRecv];
		// MPI_Alltoallv( idBlockCell_Send, mesh.mpi.sendCounts, mesh.mpi.sendDispls, MPI_INT, 
					   // idBlockCell_Recv, mesh.mpi.recvCounts, mesh.mpi.recvDispls, MPI_INT, 
					   // MPI_COMM_WORLD);
	// }
	// //==========================================
				   
	
	
	// face setting
	int nFacesLocal[nBlocks]; // total local faces
	std::fill_n(nFacesLocal, nBlocks, 0);
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal[idBlockOwner];
		
		int idBlockNeighbour;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				face.setType(SEMO_Types::PROCESSOR_FACE); // set processor face
				++nFacesLocal[idBlockNeighbour];
			}
		}
		// //==========================================
		// // @@@ at MPI 
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// idBlockNeighbour = idBlockCell_Recv[face.neighbour]; // MPI neighbour
			// if(idBlockOwner != idBlockNeighbour){
				// face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				// ++nFacesLocal[idBlockNeighbour];
			// }
			// else {
				// face.setType(SEMO_Types::TO_BE_INTERNAL_FACE); // set unknown face
				// ++nFacesLocal[idBlockNeighbour];
			// }
		// }
		// //==========================================
		
	}
	
	// write # of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/faces." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	
	int idFaceLocal[nBlocks];  // temporary local face id
	std::fill_n(idFaceLocal, nBlocks, 0);
	
	// internal faces
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlock = idBlockCell[face.owner];
			++idFaceLocal[idBlock];
			
			vector<int> idFacePoint; // face's points id
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
		
			// save face's points size , local points id
			decomPointsName[idBlock] << face.points.size() << "(";
			decomPointsName[idBlock] << idFacePoint[0];
			for(int i=1; i<idFacePoint.size(); ++i){
				decomPointsName[idBlock] << " " << idFacePoint[i];
			}
			decomPointsName[idBlock] << ")" << endl;
		}
	}
	
	
	// //==========================================
	// // @@@ at MPI
	// // TO_BE_INTERNAL_FACE 
	
	// vector<int> idBlockCell_Recv;
	// if(mpi.getSize() > 0){
		
		// int nBlockCell_Send[nBlocks];
		// int nBlockCell_Recv[nBlocks];
		// std::fill_n(nBlockCell_Send, nBlocks, 0);
		// std::fill_n(nBlockCell_Recv, nBlocks, 0);
		
		// vector<int> idBlockCell_vec;
		// for(auto& face : mesh.faces){
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// int idBlock = idBlockCell[face.owner];
				// idBlockCell_vec.push_back(idBlock);
				// ++nBlockCell_Send[idBlock];
				// ++nBlockCell_Recv[idBlock];
			// }
		// }
		
		// int disBlockCell_Send[nBlocks];
		// int disBlockCell_Recv[nBlocks];
		// disBlockCell_Send[0]=0;
		// disBlockCell_Recv[0]=0;
		// for(int i=1; i<nBlocks; ++i){
			// disBlockCell_Send[i] = disBlockCell_Send[i-1] + nBlockCell_Send[i-1];
			// disBlockCell_Recv[i] = disBlockCell_Recv[i-1] + nBlockCell_Recv[i-1];
		// }
		// int idBlockCell_Send[idBlockCell_vec.size()];
		// // int idBlockCell_Recv[idBlockCell_vec.size()];
		
		// for(int i=0; i<idBlockCell_vec.size(); ++i){
			// idBlockCell_Send[i] = idBlockCell_vec[i];
		// }
		
		// MPI_Alltoallv( idBlockCell_Send, nBlockCell_Send, disBlockCell_Send, MPI_INT, 
					   // idBlockCell_Recv.data(), nBlockCell_Recv, disBlockCell_Recv, MPI_INT, 
					   // MPI_COMM_WORLD);
		
		
	// }
	
	
	
	
	// vector<int> idFaceToBeInt_Send;
	// int nFaceToBeInt_Send[nBlocks];
	// int nFaceToBeInt_Recv[nBlocks];
	// std::fill_n(nFaceToBeInt_Send, nBlocks, 0);
	// std::fill_n(nFaceToBeInt_Recv, nBlocks, 0);
	// int temp_num_to_be_internal = 0;
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::TO_BE_INTERNAL_FACE){
			
			// int idBlockOwner = idBlockCell[face.owner];
			// int idBlockNeighbour = idBlockCell_Recv[face.neighbour];
			
			// if( idBlockOwner != mpi.getRank() ){
				// idFaceToBeInt_Send.push_back(temp_num_to_be_internal);
				// // local send face size [from block]
				// ++nFaceToBeInt_Send[idBlockNeighbour]; 
			// }
			// else {
				// // local recv face size [from block]
				// ++nFaceToBeInt_Recv[idBlockNeighbour]; 
			// }
			
		// }
		// ++temp_num_to_be_internal;
	// }
	
	// vector<int> idCellLocal_Recv;
	// if(mpi.getSize() > 0){
		// // int sendCounts[nBlocks];
		// // int recvCounts[nBlocks];
		// // std::fill_n(sendCounts, nBlocks, 0);
		// // std::fill_n(recvCounts, nBlocks, 0);
		// // MPI_Allgather(&ncells, 1, MPI_INT, myvtxdist, 1, MPI_INT, MPI_COMM_WORLD);
		// // MPI_Alltoall(&numb_SendNode, 1, MPI_INT, numb_RecvNode, 1, MPI_INT, MPI_COMM_WORLD);
		// // for(int i=0; i<mpi.getSize(); ++i){
			// // sendCounts[i] = nFaceNeighbourLocalToBeInternal[i];
		// // }
		// // int nSend = nFaceNeighbourLocalToBeInternal;
		// // int* idBlockCell_Send = new int[mesh.mpi.nSend];
		// // int* idBlockCell_Recv = new int[mesh.mpi.nRecv];
		// int sendDispls[nBlocks];
		// int recvDispls[nBlocks];
		// sendDispls[0]=0;
		// recvDispls[0]=0;
		// for(int i=1; i<nBlocks; ++i){
			// sendDispls[i] = sendDispls[i-1] + nFaceToBeInt_Send[i-1];
			// recvDispls[i] = recvDispls[i-1] + nFaceToBeInt_Recv[i-1];
		// }
		// int nSizeSend = sendDispls[nBlocks-1]+nFaceToBeInt_Send[nBlocks-1];
		// int nSizeRecv = recvDispls[nBlocks-1]+nFaceToBeInt_Recv[nBlocks-1];
		// int idCellLocal_Send[nSizeSend];
		
		// for(int i=0; i<nBlocks; ++i){
			// int strB=sendDispls[i];
			// int endB=strB+nFaceToBeInt_Send[i];
			// for(int j=strB; j<endB; ++j){
				// int idFace = idFaceToBeInt_Send[j];
				// int idOwner = mesh.faces[idFace].owner;
				// int idOwnerLocal = idCellLocal[idOwner];
				// idCellLocal_Send[j] = idOwnerLocal;
			// }
		// }
		
		// MPI_Alltoallv( idCellLocal_Send, nFaceToBeInt_Send, sendDispls, MPI_INT, 
					   // idCellLocal_Recv.data(), nFaceToBeInt_Recv, recvDispls, MPI_INT, 
					   // MPI_COMM_WORLD);
					   
		// // for(int i=1; i<nBlocks; ++i){
			// // int strB=recvDispls[i];
			// // int endB=strB+nFaceToBeInt_Recv[i];
			// // for(int j=strB; j<endB; ++j){
				// // int idOwnerLocal = idCellLocal_Recv[j];
			// // }
		// // }
		
	// }
	
	
	
	// // idCellLocal_Recv : TO_BE_INTERNAL_FACE connection neighbour local cell
	
	
	// //==========================================
	
	
	
	// boundary faces
	int nFacesBoundaryLocal[nBlocks]; // local boundary face size
	int nFacesEachBoundaryLocal[nBlocks][mesh.boundary.size()]; // local each boundary face size
	int nStartFaceEachBoundaryLocal[nBlocks][mesh.boundary.size()]; // local each boundary face start
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
	int nFacesProcessorLocal[nBlocks]; // local processor face size
	int sendCounts[nBlocks][nBlocks]; // sending size from each processor to each processor 
	for(int i=0; i<nBlocks; ++i){
		nFacesProcessorLocal[i] = 0;
		for(int j=0; j<nBlocks; ++j)
			sendCounts[i][j]=0;
	}
	vector<int> idFacesProcessor; // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
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
				
                // inverse face vector direction, because ngbr cell index convert to owner cell index
				std::reverse(idFacePoint.begin(),idFacePoint.end());
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
	
	
	
	
	//==========================================
	// @@@ at MPI
	// TO_BE_PROCESSOR_FACE
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::TO_BE_PROCESSOR_FACE){
			
		}
	}
	
	
	
	//==========================================
	
	
	// write # of owner, same as number of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/owner." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// write owner of internal face
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// write owner of boundary face
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// write owner of processor face
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
	
	// write neighbour of internal face
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
	
	
	
	// write # of boundary and processor faces
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
	
	// write of boundary faces
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
	
	
	// write of processor faces
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts[i][j];
		}
		
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
	
	
	// close boundary & processor face files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	
}