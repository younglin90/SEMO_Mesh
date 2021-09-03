#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <mpi.h>
using namespace std;

#include "load.h" 
#include "build.h"  
#include "../controls/build.h" 
#include "../utility/read.h" 
#include "../solvers/build.h" 

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



void SEMO_Mesh_Load::OpenFoam(SEMO_Mesh_Builder &mesh, string folder){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute load OpenFoam files ... ";
	
	if(rank == 0){
		
		string gridFolderName = folder;
		string pointsName = "points";
		string facesName = "faces";
		string ownerName = "owner";
		string neighbourName = "neighbour";
		string boundaryName = "boundary";
		
		ifstream inputFile;
		string openFileName;
		
		// points
		openFileName = gridFolderName + "/" + pointsName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		string nextToken;
		bool startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					double xyz[3];
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
		// cout << "------------------------------------" << endl;
		// cout << "point x,y,z size : " << mesh.points.size() << endl;
		inputFile.close();
		
		
		

		// faces
		openFileName = gridFolderName + "/" + facesName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
		
		
		
		
		// cout << "face size : " << mesh.faces.size() << endl;
		inputFile.close();
		
		
		
		// owner
		openFileName = gridFolderName + "/" + ownerName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
		// cout << "owner size : " << temp_num << endl;
		inputFile.close();
		
		
		
		
		// neighbour
		openFileName = gridFolderName + "/" + neighbourName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
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
						// cout << tempint << endl;
						
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
		// cout << "neighbour size : " << temp_num << endl;
		inputFile.close();
		
		
		
		// boundary
		openFileName = gridFolderName + "/" + boundaryName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		vector<string> boundary_name;
		vector<char> boundary_type;
		vector<int> boundary_nFaces;
		vector<int> boundary_startFace;
		vector<int> boundary_myProcNo;
		vector<int> boundary_neighbProcNo;
		
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
					trim(names);
					
			// cout << names.length() << endl;
					
					boundary_name.push_back(names);
					int l=0;
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
						if( item2.find("myProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_myProcNo.push_back(temptempint);
							++l;
						}
						if( item2.find("neighbProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_neighbProcNo.push_back(temptempint);
							++l;
						}
					}
					if(l==0) {
						boundary_myProcNo.push_back(-1);
						boundary_neighbProcNo.push_back(-1);
					}
					
					startInput=false;
					setToken2.clear();
				}
				setToken2.push_back(item.c_str());
			}
			else{ 
				if( item.find("{") != string::npos ){
					names.clear();
					names = backToken;
					startInput=true;
				}
			}
			backToken = item;
		}
		// cout << mesh.faces.size() << " " << mesh.cells.size() << " " << mesh.points.size() << " " << mesh.boundary.size() << endl;
		
		inputFile.close();
		mesh.boundary.clear();
		int nbcs=boundary_name.size();
		for (int i=0; i<boundary_name.size(); ++i) {
			// cout << endl;
			// cout << endl;
			// cout << trim(boundary_name[i]) << endl;
			// cout << boundary_nFaces[i] << endl;
			// cout << boundary_startFace[i] << endl;
			// cout << boundary_myProcNo[i] << endl;
			// cout << boundary_neighbProcNo[i] << endl;
			
			mesh.addBoundary();
		}
		for (int i=0; i<mesh.boundary.size(); ++i) {
			mesh.boundary[i].name = trim(boundary_name[i]);
			mesh.boundary[i].nFaces = boundary_nFaces[i];
			mesh.boundary[i].startFace = boundary_startFace[i];
			mesh.boundary[i].myProcNo = boundary_myProcNo[i];
			mesh.boundary[i].neighbProcNo = boundary_neighbProcNo[i];
			
			
			// cout << mesh.boundary.back().startFace << endl;
		}
	}
	
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	

	mesh.check();
	
	mesh.buildCells();
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// mesh.checkLists();
	
	// cell's points, faces connection
	mesh.connectCelltoFaces();
	
	mesh.connectCelltoPoints();
		
	mesh.informations();
		
		
		
	// mesh.checkQualities();
	
	
	
}



void SEMO_Mesh_Load::vtu(
	SEMO_Mesh_Builder &mesh, 
	string folder){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load vtu data files ... ";
	// }
		
	// string saveFolderName = folder;
	// string saveFileName = "plot";
	// string saveRankName = to_string(rank);
	
	// ifstream inputFile;
	// string openFileName;
	
	// // points
	// openFileName = saveFolderName + "/" + saveFileName + "." + saveRankName + ".vtu";
	// inputFile.open(openFileName);
	// if(inputFile.fail()){
		// cerr << "Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	// vector<int> connectivity;
	// vector<int> offsets;
	// vector<int> faces;
	// vector<int> faceoffsets;
	
	// vector<int> owner;
	// vector<int> neighbour;
	// vector<string> bcName;
	// vector<int> bcStartFace;
	// vector<int> bcNFaces;
	// vector<int> bcNeighbProcNo;
	
	// string nextToken;
	// bool startPoints=false;
	// bool startConnectivity=false;
	// bool startOffsets=false;
	// bool startTypes=false;
	// bool startFaces=false;
	// bool startFaceoffsets=false;
	// bool startowner=false;
	// bool startneighbour=false;
	// bool startbcName=false;
	// bool startbcStartFace=false;
	// bool startbcNFaces=false;
	// bool startbcNeighbProcNo=false;

	// while(getline(inputFile, nextToken)){
		// string asignToken;

		// if(startPoints){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startPoints=false;
			// }
			// else{
				// double xyz[3];
				// stringstream sstream(nextToken);
				// string word;
				// char del = ' ';
				// int num=0;
				// while (getline(sstream, word, del)){
					// xyz[num] = stod(word);
					// ++num;
				// }
				// mesh.addPoint();
				// mesh.points.back().x = xyz[0];
				// mesh.points.back().y = xyz[1];
				// mesh.points.back().z = xyz[2];
				
			// }
		// }
		// else if(startConnectivity){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startConnectivity=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// connectivity.push_back(tempint);
				// }
			// }
		// }
		// else if(startOffsets){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startOffsets=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// offsets.push_back(tempint);
				// }
			// }
		// }
		// else if(startTypes){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startTypes=false;
			// }
		// }
		// else if(startFaces){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startFaces=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// faces.push_back(tempint);
				// }
			// }
		// }
		// else if(startFaceoffsets){
			// if(nextToken.find("</DataArray>") != string::npos){
				// startFaceoffsets=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// faceoffsets.push_back(tempint);
				// }
			// }
		// }
		// // additional
		// else if(startowner){
			// if(nextToken.find("</owner>") != string::npos){
				// startowner=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// owner.push_back(tempint);
				// }
			// }
		// }
		// else if(startneighbour){
			// if(nextToken.find("</neighbour>") != string::npos){
				// startneighbour=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// neighbour.push_back(tempint);
				// }
			// }
		// }
		// else if(startbcName){
			// if(nextToken.find("</bcName>") != string::npos){
				// startbcName=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// string tempint;
				// while(iss >> tempint){
					// bcName.push_back(tempint);
				// }
			// }
		// }
		// else if(startbcStartFace){
			// if(nextToken.find("</bcStartFace>") != string::npos){
				// startbcStartFace=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// bcStartFace.push_back(tempint);
				// }
			// }
		// }
		// else if(startbcNFaces){
			// if(nextToken.find("</bcNFaces>") != string::npos){
				// startbcNFaces=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// bcNFaces.push_back(tempint);
				// }
			// }
		// }
		// else if(startbcNeighbProcNo){
			// if(nextToken.find("</bcNeighbProcNo>") != string::npos){
				// startbcNeighbProcNo=false;
			// }
			// else{
				// istringstream iss(nextToken);
				// int tempint;
				// while(iss >> tempint){
					// bcNeighbProcNo.push_back(tempint);
				// }
			// }
		// }
		// else{
			// if( nextToken.find("\"NodeCoordinates\"") != string::npos ){
				// startPoints=true;
			// }
			// else if( nextToken.find("\"connectivity\"") != string::npos ){
				// startConnectivity=true;
			// }
			// else if( nextToken.find("\"offsets\"") != string::npos ){
				// startOffsets=true;
			// }
			// else if( nextToken.find("\"types\"") != string::npos ){
				// startTypes=true;
			// }
			// else if( nextToken.find("\"faces\"") != string::npos ){
				// startFaces=true;
			// }
			// else if( nextToken.find("\"faceoffsets\"") != string::npos ){
				// startFaceoffsets=true;
			// }
			// // additional
			// else if( nextToken.find("owner") != string::npos ){
				// startowner=true;
			// }
			// else if( nextToken.find("neighbour") != string::npos ){
				// startneighbour=true;
			// }
			// else if( nextToken.find("bcName") != string::npos ){
				// startbcName=true;
			// }
			// else if( nextToken.find("bcStartFace") != string::npos ){
				// startbcStartFace=true;
			// }
			// else if( nextToken.find("bcNFaces") != string::npos ){
				// startbcNFaces=true;
			// }
			// else if( nextToken.find("bcNeighbProcNo") != string::npos ){
				// startbcNeighbProcNo=true;
			// }
		// }
		
	// }
	// // cout << "------------------------------------" << endl;
	// // cout << "rank : " << rank << " / point x,y,z size : " << mesh.points.size() << endl;
	// inputFile.close();
	
	
	
	// int n=0;
	// int ncells=-1;
	// mesh.faces.clear();
	// for(auto& i : owner){
		// mesh.addFace();
		// mesh.faces.back().owner = i;
		// ncells = max(ncells, mesh.faces.back().owner);
	// }
	// owner.clear();
	
	// mesh.cells.clear();
	// for(int i=0; i<ncells+1; ++i){
		// mesh.addCell();
	// }
	
	// n=0;
	// for(auto& i : neighbour){
		// mesh.faces[n].neighbour = i;
		// ++n;
	// }
	// neighbour.clear();
	
	// int m=0;
	// n=0;
	// for(auto& i : offsets){
		// for(int j=n; j<i; ++j){
			// int point = connectivity[j];
			// mesh.cells[m].points.push_back( point );
			// // if(rank==1) cout << point << endl;
		// }
		// n=i;
		// ++m;
	// }
	
	
	
	// n=0;
	// int nFacesInt=0;
	// for(auto& face : mesh.faces){
		// if(face.neighbour != -1){
			// mesh.cells[ face.owner ].faces.push_back( n );
			// mesh.cells[ face.neighbour ].faces.push_back( n );
			// ++nFacesInt;
		// }
		// else{
			// mesh.cells[ face.owner ].faces.push_back( n );
		// }
		// ++n;
	// }
	
	
	// m=0;
	// n=0;
	// for(auto& i : faceoffsets){
		// // if(faces[n]>5) cout << faces[n] << endl;
		// int N=0;
		// int face_size = faces[m+N];
		// for(int j=0; j<face_size; ++j){
			// int face = mesh.cells[n].faces[j];
			// ++N;
			// int point_size = faces[m+N];
			// for(int k=0; k<point_size; ++k){
				// ++N;
				// int point = faces[m+N];
				// if(mesh.faces[ face ].points.size() == point_size) continue;
				// mesh.faces[ face ].points.push_back( point );
				// // if(rank==1) cout << point << endl;
			// }
		// }
		// m=i;
		// ++n;
	// }
	// faces.clear();
	// faceoffsets.clear();
	
	// n=0;
	// for(auto& startFace : bcStartFace){
		
		// mesh.addBoundary();
		
		// // trim;
		// bcName[n].erase(std::find_if(bcName[n].rbegin(), bcName[n].rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName[n].end());
		// bcName[n].erase(bcName[n].begin(), std::find_if(bcName[n].begin(), bcName[n].end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		
		// mesh.boundary.back().name = bcName[n];
		// mesh.boundary.back().startFace = bcStartFace[n];
		// mesh.boundary.back().nFaces = bcNFaces[n];
		// mesh.boundary.back().neighbProcNo = bcNeighbProcNo[n];
		// mesh.boundary.back().myProcNo = rank;
		// if(bcNeighbProcNo[n] < 0){
			// mesh.boundary.back().myProcNo = -1;
		// }
		
		// ++n;
		
	// }
	// bcName.clear();
	// bcStartFace.clear();
	// bcNFaces.clear();
	// bcNeighbProcNo.clear();
	
		
		
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	

	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load boundary property files ... ";
	// }
	
	// vector<string> saveNameP;
	// vector<string> saveTypeP;
	// vector<double> saveValueP;
	// this->boundaryProperties("./constant/pressure", mesh, saveNameP, saveTypeP, saveValueP);
	
	// vector<string> saveNameT;
	// vector<string> saveTypeT;
	// vector<double> saveValueT;
	// this->boundaryProperties("./constant/temperature", mesh, saveNameT, saveTypeT, saveValueT);
	

	// vector<vector<string>> saveNameVF;
	// vector<vector<string>> saveTypeVF;
	// vector<vector<double>> saveValueVF;
	// for(int ns=0; ns<controls.nSp-1; ++ns){
		// vector<string> saveNameVF0;
		// vector<string> saveTypeVF0;
		// vector<double> saveValueVF0;
		// string tmpname = "./constant/" + species[ns].name;
		// this->boundaryProperties(tmpname, mesh, saveNameVF0, saveTypeVF0, saveValueVF0);
		// saveNameVF.push_back(saveNameVF0);
		// saveTypeVF.push_back(saveTypeVF0);
		// saveValueVF.push_back(saveValueVF0);
	// }
	
	
	// vector<string> saveNameVel;
	// vector<string> saveTypeVel;
	// vector<vector<double>> saveValueVel;
	// this->boundaryVelocities("./constant/velocities", mesh, saveNameVel, saveTypeVel, saveValueVel);
	
	// int neq = 6;
	// for(int i=0; i<mesh.boundary.size(); ++i){
		
		// if(mesh.boundary[i].neighbProcNo == -1){
			
			// mesh.boundary[i].type.resize(neq,"");
			// mesh.boundary[i].var.resize(neq,0.0);
			
			// bool nameMatching = false;
			// for(int j=0; j<saveNameP.size(); ++j){
				// if(mesh.boundary[i].name == saveNameP[j]){
					// mesh.boundary[i].type[0] = trim(saveTypeP[j]);
					// mesh.boundary[i].var[0] = saveValueP[j];
					// nameMatching = true;
					// break;
				// }
			// }
			// if(nameMatching==false){
				// cerr << "| #Error : boundary name matching failure" << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
			
			// nameMatching = false;
			// for(int j=0; j<saveNameVel.size(); ++j){
				// if(mesh.boundary[i].name == saveNameVel[j]){
					// mesh.boundary[i].type[1] = trim(saveTypeVel[j]);
					// mesh.boundary[i].type[2] = trim(saveTypeVel[j]);
					// mesh.boundary[i].type[3] = trim(saveTypeVel[j]);
					// mesh.boundary[i].var[1] = saveValueVel[j][0];
					// mesh.boundary[i].var[2] = saveValueVel[j][1];
					// mesh.boundary[i].var[3] = saveValueVel[j][2];
					// nameMatching = true;
					// break;
				// }
			// }
			// if(nameMatching==false){
				// cerr << "| #Error : boundary name matching failure" << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
			
			// nameMatching = false;
			// for(int j=0; j<saveNameT.size(); ++j){
				// if(mesh.boundary[i].name == saveNameT[j]){
					// mesh.boundary[i].type[4] = trim(saveTypeT[j]);
					// mesh.boundary[i].var[4] = saveValueT[j];
					// nameMatching = true;
					// break;
				// }
			// }
			// if(nameMatching==false){
				// cerr << "| #Error : boundary name matching failure" << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
			

			// for(int ns=0; ns<controls.nSp-1; ++ns){
				// nameMatching = false;
				// for(int j=0; j<saveNameVF[ns].size(); ++j){
					// if(mesh.boundary[i].name == saveNameVF[ns][j]){
						// mesh.boundary[i].type[5+ns] = trim(saveTypeVF[ns][j]);
						// mesh.boundary[i].var[5+ns] = saveValueVF[ns][j];
						// nameMatching = true;
						// break;
					// }
				// }
				// if(nameMatching==false){
					// cerr << "| #Error : boundary name matching failure, VF" << endl;
					// for(auto& j : saveNameVF[ns]){
						// cerr << j << endl;
					// }
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
			// }
			
			
			
			// // if(rank==0) cout << mesh.boundary[i].type[1];
			
		// }
	// }
	
	
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	// // mesh.buildCells();
	
		
	// mesh.check();
	
	
	// mesh.setFaceTypes();

	// // create list
	// mesh.buildLists();
	
	// // mesh.checkLists();
	
	// // cell's points, faces connection
	// // mesh.connectCelltoFaces();
	
	// // mesh.connectCelltoPoints();
	
		
	// // set processor face counts
	// mesh.setCountsProcFaces();
	
	// // set processor face displacements
	// mesh.setDisplsProcFaces(); 
	

	// mesh.informations();

	// // mesh.saveFile("vtu");
		
	
	
	
}














void SEMO_Mesh_Load::vtu(
	string folder, 
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species){
	
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
	bool startTimeValue=false;
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
	
	bool startPressure=false;
	vector<double> pressure;
	bool startVelocity=false;
	vector<double> xVelocity;
	vector<double> yVelocity;
	vector<double> zVelocity;
	bool startTemperature=false;
	vector<double> temperature;
	vector<bool> startVolFrac(1,false);
	vector<vector<double>> volFrac(1,vector<double>(0,0.0));
	
	bool startPointLevels=false;
	vector<int> pointLevels;
	bool startCellLevels=false;
	vector<int> cellLevels;
	bool startCellGroups=false;
	vector<int> cellGroups;
	bool startFaceLevels=false;
	vector<int> faceLevels;
	bool startFaceGroups=false;
	vector<int> faceGroups;
	
	vector<string> volFracName;
	string tmpVolFracName;
	tmpVolFracName = "volumeFraction-" + species[0].name;
	volFracName.push_back(tmpVolFracName);
	
	// cout << volFracName[0] << endl;
	
	double timevalue = 0.0;

	while(getline(inputFile, nextToken)){
		string asignToken;

		if(startTimeValue){
			// cout << " start " << endl;
			if(nextToken.find("</DataArray>") != string::npos){
				startTimeValue=false;
			}
			else{
				istringstream iss(nextToken);
				double tempint;
				while(iss >> tempint){
					timevalue = tempint;
				}
			}
		}
		else if(startPointLevels){
			if(nextToken.find("</DataArray>") != string::npos){
				startPointLevels=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					pointLevels.push_back(tempint);
				}
			}
		}
		else if(startPressure){
			// cout << " start " << endl;
			if(nextToken.find("</DataArray>") != string::npos){
				startPressure=false;
			}
			else{
				istringstream iss(nextToken);
				double tempint;
				while(iss >> tempint){
					pressure.push_back(tempint);
				}
			}
		}
		else if(startVelocity){
			if(nextToken.find("</DataArray>") != string::npos){
				startVelocity=false;
			}
			else{
				istringstream iss(nextToken);
				double tempint;
				int order = 0;
				while(iss >> tempint){
					if(order==0){
						xVelocity.push_back(tempint);
						++order;
					}
					else if(order==1){
						yVelocity.push_back(tempint);
						++order;
					}
					else if(order==2){
						zVelocity.push_back(tempint);
						order=0;
					}
				}
			}
		}
		else if(startTemperature){
			if(nextToken.find("</DataArray>") != string::npos){
				startTemperature=false;
			}
			else{
				istringstream iss(nextToken);
				double tempint;
				while(iss >> tempint){
					temperature.push_back(tempint);
				}
			}
		}
		else if(startVolFrac[0]){
			if(nextToken.find("</DataArray>") != string::npos){
				startVolFrac[0]=false;
			}
			else{
				istringstream iss(nextToken);
				double tempint;
				while(iss >> tempint){
					volFrac[0].push_back(tempint);
				}
			}
		}
		else if(startCellLevels){
			if(nextToken.find("</DataArray>") != string::npos){
				startCellLevels=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					cellLevels.push_back(tempint);
				}
			}
		}
		else if(startCellGroups){
			if(nextToken.find("</DataArray>") != string::npos){
				startCellGroups=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					cellGroups.push_back(tempint);
				}
			}
		}
		else if(startPoints){
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
		else if(startFaceLevels){
			if(nextToken.find("</faceLevels>") != string::npos){
				startFaceLevels=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					faceLevels.push_back(tempint);
				}
			}
		}
		else if(startFaceGroups){
			if(nextToken.find("</faceGroups>") != string::npos){
				startFaceGroups=false;
			}
			else{
				istringstream iss(nextToken);
				int tempint;
				while(iss >> tempint){
					faceGroups.push_back(tempint);
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
					// cout << tempint << endl;
					bcName.push_back(tempint);
					// cout << bcName.back() << endl;
					
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
			if( nextToken.find("\"TimeValue\"") != string::npos ){
				startTimeValue=true;
			}
			else if( nextToken.find("\"pointLevels\"") != string::npos ){
				startPointLevels=true;
			}
			else if( nextToken.find("\"pressure\"") != string::npos ){
				startPressure=true;
			}
			else if( nextToken.find("\"velocity\"") != string::npos ){
				startVelocity=true;
			}
			else if( nextToken.find("\"temperature\"") != string::npos ){
				startTemperature=true;
			}
			else if( nextToken.find(volFracName[0]) != string::npos ){
				startVolFrac[0]=true;
			}
			else if( nextToken.find("\"cellLevels\"") != string::npos ){
				startCellLevels=true;
			}
			else if( nextToken.find("\"cellGroups\"") != string::npos ){
				startCellGroups=true;
			}
			else if( nextToken.find("\"NodeCoordinates\"") != string::npos ){
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
			// else if( nextToken.find("faceLevels") != string::npos ){
				// startFaceLevels=true;
			// }
			// else if( nextToken.find("faceGroups") != string::npos ){
				// startFaceGroups=true;
			// }
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
		// cout << pressure.size() << endl;
	
	
	
	controls.time = timevalue;
	
	
	
	// SEMO_Solvers_Builder solver;
	
	for(int i=0; i<mesh.points.size(); ++i){
		mesh.points[i].level = pointLevels[i];
		
		// cout << pointLevels[i] << endl;
	}
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();

		mesh.cells.back().var.resize(controls.nTotalCellVar,0.0);

		mesh.cells.back().var[controls.P] = pressure[i];
		mesh.cells.back().var[controls.U] = xVelocity[i];
		mesh.cells.back().var[controls.V] = yVelocity[i];
		mesh.cells.back().var[controls.W] = zVelocity[i];
		mesh.cells.back().var[controls.T] = temperature[i];
		
		if(controls.nSp>1){
			mesh.cells.back().var[controls.VF[0]] = volFrac[0][i];
			for(int is=1; is<controls.nSp; ++is){
				mesh.cells.back().var[controls.VF[is]] = volFrac[0][is];
			}
		}
		else{
			mesh.cells.back().var[controls.VF[0]] = 1.0;
		}
		
		
		
		
		mesh.cells.back().level = cellLevels[i];
		mesh.cells.back().group = cellGroups[i];
		
		// mesh.cells.back().level = 0;
		// mesh.cells.back().group = i;
		
		// // EOS
		// vector<double> volumeFractions;
		// double VFnSp = 0.0;
		// for(int ns=0; ns<controls.nSp-1; ++ns){
			// volumeFractions.push_back(mesh.cells.back().var[controls.VF[ns]]);
			// VFnSp += mesh.cells.back().var[controls.VF[ns]];
		// }
		// volumeFractions.push_back(1.0 - VFnSp);
		
		// solver.getValuesFromEOSVF(
			// species,
			// mesh.cells.back().var[controls.P], 
			// mesh.cells.back().var[controls.U], 
			// mesh.cells.back().var[controls.V], 
			// mesh.cells.back().var[controls.W], 
			// mesh.cells.back().var[controls.T], 
			// volumeFractions, 
			// mesh.cells.back().var[controls.Rho], 
			// mesh.cells.back().var[controls.C], 
			// mesh.cells.back().var[controls.Ht]);



		
	}
	
	
	int tmpNum = 0;
	for(auto& face : mesh.faces){
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		face.var.resize(controls.nTotalFaceVar,0.0);
		
		// face.level = faceLevels[tmpNum];
		// face.group = faceGroups[tmpNum];
		// face.level = 0;
		// face.group = tmpNum;
		
		++tmpNum;
	}
	
	
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
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
	
	// n=0;
	// for(auto& startFace : bcStartFace){
		// cout << bcName[n] << endl;
		// ++n;
	// }
	
	// SEMO_Utility_Read read;
	
	n=0;
	for(auto& startFace : bcStartFace){
		
		mesh.addBoundary();
		
		// trim;
		// read.trim(mesh.boundary.back().name);
		
		mesh.boundary.back().name = trim(bcName[n]);
		mesh.boundary.back().startFace = bcStartFace[n];
		mesh.boundary.back().nFaces = bcNFaces[n];
		mesh.boundary.back().neighbProcNo = bcNeighbProcNo[n];
		mesh.boundary.back().myProcNo = rank;
		if(bcNeighbProcNo[n] < 0){
			mesh.boundary.back().myProcNo = -1;
		}
		
		++n;
		
		// cout << mesh.boundary.back().name << endl;
		
	}
	
	
	int maxBCNum = mesh.boundary.size()-1;
	mesh.boundary[maxBCNum].startFace = mesh.faces.size()-mesh.boundary[maxBCNum].nFaces;
	for(int i=maxBCNum-1; i>=0; --i){
		mesh.boundary[i].startFace = mesh.boundary[i+1].startFace-mesh.boundary[i].nFaces;
	}
	
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
		
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load boundary property files ... ";
	}
	
	

	// map<string,string> mapBCP;
	// read.file("./constant/pressure", "fluxScheme", mapBCP);
	
	
	vector<string> saveNameP;
	vector<string> saveTypeP;
	vector<double> saveValueP;
	this->boundaryProperties("./constant/pressure", mesh, saveNameP, saveTypeP, saveValueP);
	
	vector<string> saveNameT;
	vector<string> saveTypeT;
	vector<double> saveValueT;
	this->boundaryProperties("./constant/temperature", mesh, saveNameT, saveTypeT, saveValueT);
	
	vector<vector<string>> saveNameVF;
	vector<vector<string>> saveTypeVF;
	vector<vector<double>> saveValueVF;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		vector<string> saveNameVF0;
		vector<string> saveTypeVF0;
		vector<double> saveValueVF0;
		string tmpname = "./constant/" + species[ns].name;
		this->boundaryProperties(tmpname, mesh, saveNameVF0, saveTypeVF0, saveValueVF0);
		saveNameVF.push_back(saveNameVF0);
		saveTypeVF.push_back(saveTypeVF0);
		saveValueVF.push_back(saveValueVF0);
	}
	
	vector<string> saveNameVel;
	vector<string> saveTypeVel;
	vector<vector<double>> saveValueVel;
	this->boundaryVelocities("./constant/velocities", mesh, saveNameVel, saveTypeVel, saveValueVel);
	
	int neq = 6;
	for(int i=0; i<mesh.boundary.size(); ++i){
		
		if(mesh.boundary[i].neighbProcNo == -1){
			
			mesh.boundary[i].type.resize(neq,"");
			mesh.boundary[i].var.resize(neq,0.0);
			
			
			bool nameMatching = false;
			for(int j=0; j<saveNameP.size(); ++j){
				if(mesh.boundary[i].name == saveNameP[j]){
					mesh.boundary[i].type[0] = trim(saveTypeP[j]);
					mesh.boundary[i].var[0] = saveValueP[j];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure, P" << endl;
				for(auto& j : saveNameP){
					cerr << j << endl;
				}
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			nameMatching = false;
			for(int j=0; j<saveNameVel.size(); ++j){
				if(mesh.boundary[i].name == saveNameVel[j]){
					mesh.boundary[i].type[1] = trim(saveTypeVel[j]);
					mesh.boundary[i].type[2] = trim(saveTypeVel[j]);
					mesh.boundary[i].type[3] = trim(saveTypeVel[j]);
					mesh.boundary[i].var[1] = saveValueVel[j][0];
					mesh.boundary[i].var[2] = saveValueVel[j][1];
					mesh.boundary[i].var[3] = saveValueVel[j][2];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure, Vel" << endl;
				for(auto& j : saveNameVel){
					cerr << j << endl;
				}
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			nameMatching = false;
			for(int j=0; j<saveNameT.size(); ++j){
				if(mesh.boundary[i].name == saveNameT[j]){
					mesh.boundary[i].type[4] = trim(saveTypeT[j]);
					mesh.boundary[i].var[4] = saveValueT[j];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure, T" << endl;
				for(auto& j : saveNameT){
					cerr << j << endl;
				}
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			for(int ns=0; ns<controls.nSp-1; ++ns){
				nameMatching = false;
				for(int j=0; j<saveNameVF[ns].size(); ++j){
					if(mesh.boundary[i].name == saveNameVF[ns][j]){
						mesh.boundary[i].type[5+ns] = trim(saveTypeVF[ns][j]);
						mesh.boundary[i].var[5+ns] = saveValueVF[ns][j];
						nameMatching = true;
						break;
					}
				}
				if(nameMatching==false){
					cerr << "| #Error : boundary name matching failure, VF" << endl;
					for(auto& j : saveNameVF[ns]){
						cerr << j << endl;
					}
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
			}
			
			
			
			// if(rank==0) cout << mesh.boundary[i].type[1];
			
		}
	}
	
	
	// if(rank==0){
	// for(int i=0; i<mesh.boundary.size(); ++i){
		// if(mesh.boundary[i].neighbProcNo == -1){
		// cout << i << endl;
		// cout << mesh.boundary[i].name << endl;
		// cout << mesh.boundary[i].type[1] << endl;
		// }
	// }
	// }
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	// mesh.buildCells();
	
		
	mesh.check();
	
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// mesh.checkLists();
	
	// cell's points, faces connection
	// mesh.connectCelltoFaces();
	
	// mesh.connectCelltoPoints();
	
		
	// set processor face counts
	mesh.setCountsProcFaces();
	
	// set processor face displacements
	mesh.setDisplsProcFaces(); 
	
	
	
	// face level
	vector<int> cLevel_recv;
	if(size>1){
		vector<int> cLevel_send;
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cLevel_send.push_back(mesh.cells[face.owner].level);
			}
		}
		
		cLevel_recv.clear();
		cLevel_recv.resize(cLevel_send.size(),0);

		MPI_Alltoallv( cLevel_send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   cLevel_recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
					   MPI_COMM_WORLD);
	}
	
	int proc_num=0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			face.level = max(
				mesh.cells[face.owner].level, mesh.cells[face.neighbour].level);
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.level = max(
				mesh.cells[face.owner].level, cLevel_recv[proc_num]);
			++proc_num;
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			face.level = mesh.cells[face.owner].level;
		}
	}
	
	
	
	

	mesh.informations();

	// mesh.saveFile("vtu");
		
	
	
	
}
















void SEMO_Mesh_Load::boundaryProperties(string file, 
SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<double>& value){
	
	name.clear();
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;
	bool startReading2 = false;
	

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading || startReading2){
				if(nextToken.find("{") != string::npos) startReading=true;
				if(!startReading && nextToken.find("}") != string::npos) startReading2=false;
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("boundaryField") != string::npos ){
					startReading=true;
					startReading2 = true;
					// cout << nextToken << endl;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			for(int i=0; i<mesh.boundary.size(); ++i){
				if(read[line].find(mesh.boundary[i].name) != string::npos){
					
					name.push_back(mesh.boundary[i].name);
					
					for(int line2=line; line2<read.size(); ++line2){
						if(read[line2].find("type") != string::npos){
							istringstream iss2(read[line2]);
							string tempstring2;
							string tempstring3;
							iss2 >> tempstring2 >> tempstring3;
							tempstring3.erase(tempstring3.find(";"),1); 
							type.push_back(tempstring3);
							if(
							type.back() == "fixedValue" ||
							type.back() == "inletOutlet" ||
							type.back() == "switch"){
								istringstream iss2(read[line2+1]);
								string tempstring2;
								string tempstring3;
								iss2 >> tempstring2 >> tempstring3;
								tempstring3.erase(tempstring3.find(";"),1); 
								value.push_back(stod(tempstring3));
							}
							else{
								value.push_back(0.0);
							}
							break;
						}
					}
					
				}
			}
		}
	}
	inputFile.close();
	
}


void SEMO_Mesh_Load::boundaryVelocities(string file, 
SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<vector<double>>& value){
	
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;
	bool startReading2 = false;
	

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading || startReading2){
				if(nextToken.find("{") != string::npos) startReading=true;
				if(!startReading && nextToken.find("}") != string::npos) startReading2=false;
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("boundaryField") != string::npos ){
					startReading=true;
					startReading2 = true;
					// cout << nextToken << endl;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			for(int i=0; i<mesh.boundary.size(); ++i){
				if(read[line].find(mesh.boundary[i].name) != string::npos){
					
					name.push_back(mesh.boundary[i].name);
					
					for(int line2=line; line2<read.size(); ++line2){
						if(read[line2].find("type") != string::npos){
							istringstream iss2(read[line2]);
							string tempstring2;
							string tempstring3;
							iss2 >> tempstring2 >> tempstring3;
							tempstring3.erase(tempstring3.find(";"),1); 
							type.push_back(tempstring3);
							if(type.back() == "fixedValue"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								lineee.erase(lineee.find("("),1); 
								lineee.erase(lineee.find(")"),1); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0] >> tmpD[1] >> tmpD[2];
								value.push_back(tmpD);
								break;
							}
							else if(type.back() == "surfaceNormalFixedValue"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0];
								value.push_back(tmpD);
								break;
							}
							else if(type.back() == "inletOutlet"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0];
								value.push_back(tmpD);
								break;
							}
							else{
								vector<double> tmpD(3,0.0);
								value.push_back(tmpD);
								break;
							}
						}
					}
					
				}
			}
		}
	}
	else{
		cerr << "| #Warning : Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	inputFile.close();
	
}
