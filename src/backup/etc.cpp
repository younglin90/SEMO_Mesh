
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

