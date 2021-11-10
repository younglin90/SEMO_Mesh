#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string>
using namespace std;

#include "../mesh/build.h" 

void load_vtu(SEMO_Mesh_Builder &mesh, string folder);
SEMO_Mesh_Builder mesh;

int main(int argc, char* argv[]) {
	

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();


	load_vtu(mesh, "./save/0/");
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
}








void load_vtu(SEMO_Mesh_Builder &mesh, string folder){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load vtu data files ... ";
	}
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rank);
	
	ifstream inputFile;
	string openFileName;
	
	// points
	openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	inputFile.open(openFileName);
	if(inputFile.fail()){
		cerr << "Unable to open file for reading : " << openFileName << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	vector<int> connectivity;
	vector<int> offsets;
	vector<int> faces;
	vector<int> faceoffsets;
	
	vector<int> owner;
	vector<int> neighbour;
	vector<string> bcName;
	vector<int> bcStartFace;
	vector<int> bcNFaces;
	vector<int> bcNeighbProcNo;
	
	string nextToken;
	bool startPoints=false;
	bool startConnectivity=false;
	bool startOffsets=false;
	bool startTypes=false;
	bool startFaces=false;
	bool startFaceoffsets=false;
	bool startowner=false;
	bool startneighbour=false;
	bool startbcName=false;
	bool startbcStartFace=false;
	bool startbcNFaces=false;
	bool startbcNeighbProcNo=false;

	while(getline(inputFile, nextToken)){
		string asignToken;

		if(startPoints){
			if(nextToken.find("</DataArray>") != string::npos){
				startPoints=false;
			}
			else{
				double xyz[3];
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
		else if(startConnectivity){
			if(nextToken.find("</DataArray>") != string::npos){
				startConnectivity=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					connectivity.push_back(tempint);
				}
			}
		}
		else if(startOffsets){
			if(nextToken.find("</DataArray>") != string::npos){
				startOffsets=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					offsets.push_back(tempint);
				}
			}
		}
		else if(startTypes){
			if(nextToken.find("</DataArray>") != string::npos){
				startTypes=false;
			}
		}
		else if(startFaces){
			if(nextToken.find("</DataArray>") != string::npos){
				startFaces=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					faces.push_back(tempint);
				}
			}
		}
		else if(startFaceoffsets){
			if(nextToken.find("</DataArray>") != string::npos){
				startFaceoffsets=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					faceoffsets.push_back(tempint);
				}
			}
		}
		// additional
		else if(startowner){
			if(nextToken.find("</owner>") != string::npos){
				startowner=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					owner.push_back(tempint);
				}
			}
		}
		else if(startneighbour){
			if(nextToken.find("</neighbour>") != string::npos){
				startneighbour=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					neighbour.push_back(tempint);
				}
			}
		}
		else if(startbcName){
			if(nextToken.find("</bcName>") != string::npos){
				startbcName=false;
			}
			else{
				istringstream iss(nextToken);
				string tempint;
				while(iss >> tempint){
					bcName.push_back(tempint);
				}
			}
		}
		else if(startbcStartFace){
			if(nextToken.find("</bcStartFace>") != string::npos){
				startbcStartFace=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					bcStartFace.push_back(tempint);
				}
			}
		}
		else if(startbcNFaces){
			if(nextToken.find("</bcNFaces>") != string::npos){
				startbcNFaces=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					bcNFaces.push_back(tempint);
				}
			}
		}
		else if(startbcNeighbProcNo){
			if(nextToken.find("</bcNeighbProcNo>") != string::npos){
				startbcNeighbProcNo=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					bcNeighbProcNo.push_back(tempint);
				}
			}
		}
		else{
			if( nextToken.find("\"NodeCoordinates\"") != string::npos ){
				startPoints=true;
			}
			else if( nextToken.find("\"connectivity\"") != string::npos ){
				startConnectivity=true;
			}
			else if( nextToken.find("\"offsets\"") != string::npos ){
				startOffsets=true;
			}
			else if( nextToken.find("\"types\"") != string::npos ){
				startTypes=true;
			}
			else if( nextToken.find("\"faces\"") != string::npos ){
				startFaces=true;
			}
			else if( nextToken.find("\"faceoffsets\"") != string::npos ){
				startFaceoffsets=true;
			}
			// additional
			else if( nextToken.find("owner") != string::npos ){
				startowner=true;
			}
			else if( nextToken.find("neighbour") != string::npos ){
				startneighbour=true;
			}
			else if( nextToken.find("bcName") != string::npos ){
				startbcName=true;
			}
			else if( nextToken.find("bcStartFace") != string::npos ){
				startbcStartFace=true;
			}
			else if( nextToken.find("bcNFaces") != string::npos ){
				startbcNFaces=true;
			}
			else if( nextToken.find("bcNeighbProcNo") != string::npos ){
				startbcNeighbProcNo=true;
			}
		}
		
	}
	// cout << "------------------------------------" << endl;
	// cout << "rank : " << rank << " / point x,y,z size : " << mesh.points.size() << endl;
	inputFile.close();
	
	
	
	int n=0;
	int ncells=-1;
	mesh.faces.clear();
	for(auto& i : owner){
		mesh.addFace();
		mesh.faces.back().owner = i;
		ncells = max(ncells, mesh.faces.back().owner);
	}
	owner.clear();
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();
	}
	
	n=0;
	for(auto& i : neighbour){
		mesh.faces[n].neighbour = i;
		++n;
	}
	neighbour.clear();
	
	int m=0;
	n=0;
	for(auto& i : offsets){
		for(int j=n; j<i; ++j){
			int point = connectivity[j];
			mesh.cells[m].points.push_back( point );
			// if(rank==1) cout << point << endl;
		}
		n=i;
		++m;
	}
	
	
	
	n=0;
	int nFacesInt=0;
	for(auto& face : mesh.faces){
		if(face.neighbour != -1){
			mesh.cells[ face.owner ].faces.push_back( n );
			mesh.cells[ face.neighbour ].faces.push_back( n );
			++nFacesInt;
		}
		else{
			mesh.cells[ face.owner ].faces.push_back( n );
		}
		++n;
	}
	
	
	m=0;
	n=0;
	for(auto& i : faceoffsets){
		// if(faces[n]>5) cout << faces[n] << endl;
		int N=0;
		int face_size = faces[m+N];
		for(int j=0; j<face_size; ++j){
			int face = mesh.cells[n].faces[j];
			++N;
			int point_size = faces[m+N];
			for(int k=0; k<point_size; ++k){
				++N;
				int point = faces[m+N];
				if(mesh.faces[ face ].points.size() == point_size) continue;
				mesh.faces[ face ].points.push_back( point );
				// if(rank==1) cout << point << endl;
			}
		}
		m=i;
		++n;
	}
	faces.clear();
	faceoffsets.clear();
	
	n=0;
	for(auto& startFace : bcStartFace){
		
		mesh.addBoundary();
		
		// trim;
		bcName[n].erase(std::find_if(bcName[n].rbegin(), bcName[n].rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName[n].end());
		bcName[n].erase(bcName[n].begin(), std::find_if(bcName[n].begin(), bcName[n].end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		
		mesh.boundary.back().name = bcName[n];
		mesh.boundary.back().startFace = bcStartFace[n];
		mesh.boundary.back().nFaces = bcNFaces[n];
		mesh.boundary.back().neighbProcNo = bcNeighbProcNo[n];
		mesh.boundary.back().myProcNo = rank;
		if(bcNeighbProcNo[n] < 0){
			mesh.boundary.back().myProcNo = -1;
		}
		
		++n;
		
	}
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
		
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	mesh.informations();
	
	
	// mesh.buildCells();
	
		
	mesh.check();
	
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// mesh.checkLists();
	
	// // cell's points, faces connection
	// mesh.connectCelltoFaces();
	
	// mesh.connectCelltoPoints();


	// mesh.saveFile("vtu");
	
	
	// mesh.cellsGlobal();
		
	
	
	
}






