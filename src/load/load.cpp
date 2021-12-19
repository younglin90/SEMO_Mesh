
using namespace std;

#include "load.h" 
#include "../mesh/build.h"  
#include "../controls/build.h" 
#include "../utility/read.h" 
#include "../solvers/build.h" 


//앞에 있는 개행 문자 제거 
string &SEMO_Mesh_Load::ltrim(std::string &s) { 
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
	return s; 
}

//뒤에 있는 개행 문자 제거 
string &SEMO_Mesh_Load::rtrim(std::string &s) { 
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end()); 
	return s; 
}

//양쪽 끝의 개행 문자 제거 
string &SEMO_Mesh_Load::trim(std::string &s) { 
	return ltrim(rtrim(s)); 
}
	

template void SEMO_Mesh_Load::loadDatasAtVTU<int>(
	ifstream& inputFile, string dataName, vector<int>& outData);
template void SEMO_Mesh_Load::loadDatasAtVTU<double>(
	ifstream& inputFile, string dataName, vector<double>& outData);
template void SEMO_Mesh_Load::loadDatasAtVTU<string>(
	ifstream& inputFile, string dataName, vector<string>& outData);
template<typename T>
void SEMO_Mesh_Load::loadDatasAtVTU(
	ifstream& inputFile, string dataName, vector<T>& outData){
	
	string nextToken;
	string combDataName = "\"" + dataName + "\"";
	trim(combDataName);
	
	inputFile.clear();
	// inputFile.seekg( 0, std::ios_base::beg);
	outData.clear();
	
	bool startValue = false;
	bool boolBinary = false;
	while(getline(inputFile, nextToken)){
		string asignToken;

		if(startValue){
			if(nextToken.find("</DataArray>") != string::npos){
				startValue=false;
				boolBinary=false;
				break;
			}
			else{
				// istringstream iss(nextToken);
				// T tempint;
				// while(iss >> tempint){
					// outData.push_back(tempint);
				// }
				stringstream sstream(nextToken);
				string word;
				
				
				char del = ' ';
				if(boolBinary==false){
					while (getline(sstream, word, del)){
						istringstream iss(word);
						T tempint;
						iss >> tempint;
						outData.push_back(tempint);
					}
				}
				else{
					if(boolCompress==false){
						while (getline(sstream, word, del)){
							readBinary(word, outData);
						}
					}
					else{
						while (getline(sstream, word, del)){
							readCompress(word, outData);
						}
					}
				}
				
			}
		}
		else{
			if( nextToken.find(combDataName) != string::npos ){
				startValue=true;
				if( nextToken.find("format=\"binary\"") != string::npos ){
					boolBinary=true;
				}
			}
		}
		
	}
	
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
	
	
	string nextToken;
	// boolBinary = false;
	boolCompress = false;
	while(getline(inputFile, nextToken)){
		if( nextToken.find("VTKFile") != string::npos ){
			if( nextToken.find("vtkZLibDataCompressor") != string::npos ){
				boolCompress = true;
			}
			break;
		}
	}
	inputFile.clear();
	
	
	
	vector<string> volFracName;
	// string tmpVolFracName;
	// tmpVolFracName = "volumeFraction-" + species[0].name;
	// volFracName.push_back(tmpVolFracName);
	volFracName.push_back(controls.name[controls.VF[0]]);
	
	vector<string> massFracName;
	// string tmpMassFracName;
	// tmpMassFracName = "massFraction-" + species[0].name;
	// massFracName.push_back(tmpMassFracName);
	massFracName.push_back(controls.name[controls.MF[0]]);
	
	// vector<double> NodeCoordinates;
	// loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	
	// cout << NodeCoordinates.size() << endl;
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	vector<double> timevalue;
	loadDatasAtVTU(inputFile, "TimeValue", timevalue);
	if(timevalue.size()==0) timevalue.push_back(0.0);
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	vector<int> pointLevels;
	loadDatasAtVTU(inputFile, "pointLevels", pointLevels);
	
	
	vector<double> pressure;
	loadDatasAtVTU(inputFile, "pressure", pressure);
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	vector<double> velocity;
	loadDatasAtVTU(inputFile, "velocity", velocity);
	
	vector<double> temperature;
	loadDatasAtVTU(inputFile, "temperature", temperature);
	
	vector<vector<double>> volFrac(1,vector<double>(0,0.0));
	loadDatasAtVTU(inputFile, volFracName[0], volFrac[0]);
	
	vector<vector<double>> massFrac(1,vector<double>(0,0.0));
	loadDatasAtVTU(inputFile, massFracName[0], massFrac[0]);
	
	vector<int> cellLevels;
	loadDatasAtVTU(inputFile, "cellLevels", cellLevels);
	
	vector<int> cellGroups;
	loadDatasAtVTU(inputFile, "cellGroups", cellGroups);
	
	vector<double> NodeCoordinates;
	loadDatasAtVTU(inputFile, "NodeCoordinates", NodeCoordinates);
	// cout << NodeCoordinates.size() << endl;
	
	vector<int> connectivity;
	loadDatasAtVTU(inputFile, "connectivity", connectivity);
	
	vector<int> offsets;
	loadDatasAtVTU(inputFile, "offsets", offsets);
	
	vector<int> faces;
	loadDatasAtVTU(inputFile, "faces", faces);
	
	vector<int> faceoffsets;
	loadDatasAtVTU(inputFile, "faceoffsets", faceoffsets);
	
	vector<int> owner;
	loadDatasAtVTU(inputFile, "owner", owner);
	
	vector<int> neighbour;
	loadDatasAtVTU(inputFile, "neighbour", neighbour);
	
	vector<string> bcName;
	loadDatasAtVTU(inputFile, "bcName", bcName);
	
	vector<int> bcStartFace;
	loadDatasAtVTU(inputFile, "bcStartFace", bcStartFace);
	
	vector<int> bcNFaces;
	loadDatasAtVTU(inputFile, "bcNFaces", bcNFaces);
	
	vector<int> bcNeighbProcNo;
	loadDatasAtVTU(inputFile, "bcNeighbProcNo", bcNeighbProcNo);


	inputFile.close();
	
	
	
	
	
	for(int i=0; i<NodeCoordinates.size()/3; ++i){
		mesh.addPoint();
		mesh.points.back().x = NodeCoordinates[i*3+0];
		mesh.points.back().y = NodeCoordinates[i*3+1];
		mesh.points.back().z = NodeCoordinates[i*3+2];
	}
	NodeCoordinates.clear();
	
	// cout << "AAAAAAAA" << NodeCoordinates.size()/3 << endl;
	
	
	int n=0;
	int ncells=-1;
	mesh.faces.clear();
	for(auto& i : owner){
		mesh.addFace();
		mesh.faces.back().owner = i;
		ncells = max(ncells, mesh.faces.back().owner);
	}
	owner.clear();
	
	
	controls.time = timevalue[0];
	
	
	
	for(int i=0; i<mesh.points.size(); ++i){
		mesh.points[i].level = pointLevels[i];
	}
	
	mesh.cells.clear();
	for(int i=0; i<ncells+1; ++i){
		mesh.addCell();

		mesh.cells.back().var.resize(controls.nTotalCellVar,0.0);

		mesh.cells.back().var[controls.P] = pressure[i];
		
		mesh.cells.back().var[controls.U] = velocity[i*3+0];
		mesh.cells.back().var[controls.V] = velocity[i*3+1];
		mesh.cells.back().var[controls.W] = velocity[i*3+2];
		// mesh.cells.back().var[controls.U] = xVelocity[i];
		// mesh.cells.back().var[controls.V] = yVelocity[i];
		// mesh.cells.back().var[controls.W] = zVelocity[i];
		
		mesh.cells.back().var[controls.T] = temperature[i];
	
		for(int is=0; is<controls.nSp; ++is){
			mesh.cells.back().var[controls.VF[is]] = 1.0;
			mesh.cells.back().var[controls.MF[is]] = 1.0;
		}
		
		if(controls.nSp>1){
			if(i<volFrac[0].size()){
				double tmp_sum = 0.0;
				for(int is=0; is<controls.nSp-1; ++is){
					mesh.cells.back().var[controls.VF[is]] = volFrac[is][i];
					tmp_sum += volFrac[is][i];
				}
				mesh.cells.back().var[controls.VF[controls.nSp-1]] = max(0.0,min(1.0,1.0-tmp_sum));
			}
			
			if(i<massFrac[0].size()){
				double tmp_sum = 0.0;
				for(int is=0; is<controls.nSp-1; ++is){
					mesh.cells.back().var[controls.MF[is]] = massFrac[is][i];
					tmp_sum += massFrac[is][i];
				}
				mesh.cells.back().var[controls.MF[controls.nSp-1]] = max(0.0,min(1.0,1.0-tmp_sum));
			}
		}
		
		mesh.cells.back().level = cellLevels[i];
		mesh.cells.back().group = cellGroups[i];
		
		
	}
	
	
	// 추가적인 셀 값들
	mesh.cellsProcVar.resize(controls.nTotalCellVar,vector<double>());
	mesh.cellsProcGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	mesh.cellsGradientVar.resize(controls.nTotalCellVar,vector<vector<double>>());
	
	
	
	
	int tmpNum = 0;
	for(auto& face : mesh.faces){
		face.var.resize(controls.nTotalFaceVar,0.0);
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		++tmpNum;
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
	
	
	// cout << "BBBBBBB" << endl;
	// cout << bcName.size() << endl;
	// cout << bcStartFace.size() << endl;
	// cout << bcNFaces.size() << endl;
	// cout << bcNeighbProcNo.size() << endl;
	
	
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
	
	mesh.cellsGlobal();
	
	
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

	
	
}




