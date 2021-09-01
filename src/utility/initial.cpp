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
#include <stdexcept>
#include <iomanip>
#include <sys/stat.h>
#include <cstring>
#include <vector>
// #include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>

using namespace std;
#include "../mesh/build.h" 
#include "../mesh/load.h" 
#include "../mesh/geometric.h" 

#include "../controls/build.h" 
 
void initialConditionSetting(string file,
	vector<string>& type, vector<double>& value);
void initialConditionSettingVel(string file,
	vector<string>& type, vector<double>& value);
void saveInitialField(string folder);

void loadOnlyMeshVtu(
	SEMO_Mesh_Builder &mesh, 
	string folder);
	
SEMO_Mesh_Builder mesh;
SEMO_Mesh_Load load;
SEMO_Utility_Read read;
SEMO_Controls_Builder controls;
vector<SEMO_Species> species;
SEMO_Mesh_Save save;

int main(int argc, char* argv[]) {
	
	// cout << "AAA" << endl;

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	

	// cout << "BBB" << endl;
	
		
	controls.readSpecies(species);
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	SEMO_MPI_Builder mpi;
	

	// cout << "BBB" << endl; 
	
	// load.vtu("./grid/0/", mesh, controls, species);

	loadOnlyMeshVtu(mesh, "./grid/0/");
	
	
	for(auto& point : mesh.points){
		point.level = 0;
	}
	

	vector<vector<double>> cellXYZ(mesh.cells.size(),vector<double>(3,0.0));
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.level = 0;
		cell.group = i;
		
		for(auto& j : mesh.cells[i].points){
			cellXYZ[i][0] += mesh.points[j].x;
			cellXYZ[i][1] += mesh.points[j].y;
			cellXYZ[i][2] += mesh.points[j].z;
		}
		cellXYZ[i][0] /= (double)mesh.cells[i].points.size();
		cellXYZ[i][1] /= (double)mesh.cells[i].points.size();
		cellXYZ[i][2] /= (double)mesh.cells[i].points.size();
		
		cell.var.resize(6,0.0);
	}
	
	
	
	// 중복 포인트, 면 > 3 찾기
	int noAMR_Cells = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		for(auto& p0 : cell.points){
			int p0Num = 0;
			for(auto& j : cell.faces){
				auto& face = mesh.faces[j];
				auto& facePoints = face.points;
				if( std::find(facePoints.begin(),facePoints.end(),p0)!=facePoints.end() ){
					++p0Num;
				}
			}
			if(p0Num>3){
				cell.level = -1;
				// cout << "| #WARNING : 1 point share " << p0Num << " faces" << endl;
			}
		}
		if(cell.level==-1) ++noAMR_Cells;
		
	}
	cout << "| #WARNING : " << noAMR_Cells << " cells can NOT AMR" << endl;
	
	
	
	
	vector<string> type;
	vector<double> values;
	
	
	//====== initial conditions ====
	
	
	initialConditionSetting("./constant/pressure", type, values);
	
	
	if(type.back() == "fixedValue"){
		for(auto& cell : mesh.cells){
			cell.var[0] = values[0];
		}
	}
	else if(type.back() == "boxToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			if(
			(cellXYZ[i][0] >= values[1] && cellXYZ[i][0] <= values[4]) &&
			(cellXYZ[i][1] >= values[2] && cellXYZ[i][1] <= values[5]) &&
			(cellXYZ[i][2] >= values[3] && cellXYZ[i][2] <= values[6]) ){
				cell.var[0] = values[7];
			}
			else{
				cell.var[0] = values[0];
			}
		}
	}
	
	// cout << "AAA" << endl;
	
	
	initialConditionSettingVel("./constant/velocities", type, values);
	// cout << type.back() << endl;
	if(type.back() == "fixedValue"){
		// cout << values[0] << endl;
		for(auto& cell : mesh.cells){
			cell.var[1] = values[0];
			cell.var[2] = values[1];
			cell.var[3] = values[2];
		}
	}
	else if(type.back() == "boxToCell"){
		// cout << "AAAAAAAAA" << endl;
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			if(
			(cellXYZ[i][0] >= values[3] && cellXYZ[i][0] <= values[6]) &&
			(cellXYZ[i][1] >= values[4] && cellXYZ[i][1] <= values[7]) &&
			(cellXYZ[i][2] >= values[5] && cellXYZ[i][2] <= values[8]) ){
				cell.var[1] = values[9];
				cell.var[2] = values[10];
				cell.var[3] = values[11];
			}
			else{
				cell.var[1] = values[0];
				cell.var[2] = values[1];
				cell.var[3] = values[2];
			}
		}
		// cout << "BBBBBBBB" << endl;
	}
		
	initialConditionSetting("./constant/temperature", type, values);
	
	if(type.back() == "fixedValue"){
		for(auto& cell : mesh.cells){
			cell.var[4] = values[0];
		}
	}
	else if(type.back() == "boxToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			if(
			(cellXYZ[i][0] >= values[1] && cellXYZ[i][0] <= values[4]) &&
			(cellXYZ[i][1] >= values[2] && cellXYZ[i][1] <= values[5]) &&
			(cellXYZ[i][2] >= values[3] && cellXYZ[i][2] <= values[6]) ){
				cell.var[4] = values[7];
			}
			else{
				cell.var[4] = values[0];
			}
		}
	}
	
	// cout << "AAA" << endl;
	
	string speciesFileName = "./constant/" + species[0].name;
	initialConditionSetting(speciesFileName, type, values);
	
	// cout << "BBB" << endl;
	if(type.back() == "fixedValue"){
		for(auto& cell : mesh.cells){
			cell.var[5] = values[0];
		}
	}
	else if(type.back() == "boxToCell"){
		// cout << "boxToCell" << endl;
		// for(auto& i : values){
			// if(rank==0) cout << i << endl;
		// }
		// cout << values[4] << values[5] << values[6] << endl;
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			if(
			(cellXYZ[i][0] >= values[1] && cellXYZ[i][0] <= values[4]) &&
			(cellXYZ[i][1] >= values[2] && cellXYZ[i][1] <= values[5]) &&
			(cellXYZ[i][2] >= values[3] && cellXYZ[i][2] <= values[6]) ){
				// cout << values[7] << endl;
				cell.var[5] = values[7];
					// cout.precision(20);
				// cout << "YES " << cellXYZ[i][1] << endl;
			}
			else{
				// if( cellXYZ[i][1] < 0.31 && cellXYZ[i][1] > 0.29 ){
					// cout.precision(20);
					// cout << "NO " << cellXYZ[i][1] << endl;
				// }
				cell.var[5] = values[0];
			}
		}
	}
	else if(type.back() == "cylinderToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double radius = values[7];
			double deValue = values[0];
			double value = values[8];
			
			double pVec[3];
			double vec0[3];
			double vec1[3];
			
			pVec[0] = values[4]-values[1];
			pVec[1] = values[5]-values[2];
			pVec[2] = values[6]-values[3];
			
			vec0[0] = cellXYZ[i][0]-values[1];
			vec0[1] = cellXYZ[i][1]-values[2];
			vec0[2] = cellXYZ[i][2]-values[3];
			
			vec1[0] = cellXYZ[i][0]-values[4];
			vec1[1] = cellXYZ[i][1]-values[5];
			vec1[2] = cellXYZ[i][2]-values[6];
			
			double pVecMag = sqrt(pow(pVec[0],2.0)+pow(pVec[1],2.0)+pow(pVec[2],2.0));
			double vec0Mag = sqrt(pow(vec0[0],2.0)+pow(vec0[1],2.0)+pow(vec0[2],2.0));
			double vec1Mag = sqrt(pow(vec1[0],2.0)+pow(vec1[1],2.0)+pow(vec1[2],2.0));
			
			double cosVec0 = (vec0[0]*pVec[0]+vec0[1]*pVec[1]+vec0[2]*pVec[2])/pVecMag/(vec0Mag+1.e-200);
			double cosVec1 = -(vec1[0]*pVec[0]+vec1[1]*pVec[1]+vec1[2]*pVec[2])/pVecMag/(vec1Mag+1.e-200);
			
			double unitPVec[3];
			unitPVec[0] = pVec[0]/(pVecMag);
			unitPVec[1] = pVec[1]/(pVecMag);
			unitPVec[2] = pVec[2]/(pVecMag);
			
			bool boolInternal = true;
			if(cosVec0 <= 0.0 || cosVec1 <= 0.0) boolInternal = false;
			
			// cout << cosVec0 << " " << cosVec1 << endl;
			
			if(boolInternal){
				double distance = 
				sqrt( pow(vec0Mag,2.0)-pow(unitPVec[0]*vec0[0]+unitPVec[1]*vec0[1]+unitPVec[2]*vec0[2] , 2.0) );
				
				if(distance >= radius) boolInternal = false;
			}
			
			if(boolInternal){
				cell.var[5] = value;
			}
			else{
				cell.var[5] = deValue;
			}
		}
	}
	
	

	// SEMO_Utility_Math math;
	// SEMO_Mesh_Geometric geometric;
	// geometric.init(mesh);
	
	saveInitialField("./save/0/");
	// save.vtu("./save/0/", mesh, controls, species);
	

	if(rank==0) cout << "| completed save initial files : ./save/0/" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
}







void initialConditionSetting(string file,
	vector<string>& type, vector<double>& value){
	
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading){
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("initialField") != string::npos ){
					startReading=true;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			if(read[line].find("type") != string::npos){
				istringstream iss2(read[line]);
				string dummy;
				string stmp;
				iss2 >> dummy >> stmp;
				stmp.erase(stmp.find(";"),1); 
				// cout << stmp << endl;
				type.push_back(stmp);
				if(type.back() == "fixedValue"){
					istringstream iss2(read[line+1]);
					double dtmp;
					iss2 >> dummy >> dtmp;
					// dtmp.erase(dtmp.find(";"),1); 
					value.push_back(dtmp);
					// value.push_back(stod(dtmp));
				}
				else if(type.back() == "boxToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream box(read[line+2]);
					istringstream boxValue(read[line+3]);
					double dDefaultVal;
					double xyzMin[3];
					double xyzMax[3];
					double dBoxValue;
					defaultVal >> dummy >> dDefaultVal;
					// box >> dummy >> xyzMin[0] >> xyzMin[1] >> xyzMin[2] >> xyzMax[0] >> xyzMax[1] >> xyzMax[2];
					
					
					string dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
					box >> dummy >> dummy1 >> dummy2 >> dummy3 >> dummy4 >> dummy5 >> dummy6;
					
					
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy4.erase(dummy4.find("("),1);
					dummy6.erase(dummy6.find(")"),1);
					dummy6.erase(dummy6.find(";"),1);
					// cout << dummy4 << dummy5 << dummy6 << endl;
					
					xyzMin[0] = stod(dummy1);
					xyzMin[1] = stod(dummy2);
					xyzMin[2] = stod(dummy3);
					xyzMax[0] = stod(dummy4);
					xyzMax[1] = stod(dummy5);
					xyzMax[2] = stod(dummy6);
					
					// tmpdummy >> xyzMin[0] >> xyzMin[1] >> xyzMin[2] >> xyzMax[0] >> xyzMax[1] >> xyzMax[2];
					
					boxValue >> dummy >> dBoxValue;
					
					value.push_back(dDefaultVal);
					value.push_back(xyzMin[0]);
					value.push_back(xyzMin[1]);
					value.push_back(xyzMin[2]);
					value.push_back(xyzMax[0]);
					value.push_back(xyzMax[1]);
					value.push_back(xyzMax[2]);
					value.push_back(dBoxValue);
				}
				else if(type.back() == "cylinderToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream points(read[line+2]);
					istringstream radius(read[line+3]);
					istringstream cylinderValue(read[line+4]);
					double dDefaultVal;
					double xyzMin[3];
					double xyzMax[3];
					double radiusVal;
					double cylinderVal;
					defaultVal >> dummy >> dDefaultVal;
					// box >> dummy >> xyzMin[0] >> xyzMin[1] >> xyzMin[2] >> xyzMax[0] >> xyzMax[1] >> xyzMax[2];
					
					
					string dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
					points >> dummy >> dummy1 >> dummy2 >> dummy3 >> dummy4 >> dummy5 >> dummy6;
					
					
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy4.erase(dummy4.find("("),1);
					dummy6.erase(dummy6.find(")"),1);
					dummy6.erase(dummy6.find(";"),1);
					// cout << dummy4 << dummy5 << dummy6 << endl;
					
					xyzMin[0] = stod(dummy1);
					xyzMin[1] = stod(dummy2);
					xyzMin[2] = stod(dummy3);
					xyzMax[0] = stod(dummy4);
					xyzMax[1] = stod(dummy5);
					xyzMax[2] = stod(dummy6);
					
					// tmpdummy >> xyzMin[0] >> xyzMin[1] >> xyzMin[2] >> xyzMax[0] >> xyzMax[1] >> xyzMax[2];
					
					radius >> dummy >> radiusVal;
					cylinderValue >> dummy >> cylinderVal;
					
					value.push_back(dDefaultVal);
					value.push_back(xyzMin[0]);
					value.push_back(xyzMin[1]);
					value.push_back(xyzMin[2]);
					value.push_back(xyzMax[0]);
					value.push_back(xyzMax[1]);
					value.push_back(xyzMax[2]);
					value.push_back(radiusVal);
					value.push_back(cylinderVal);
				}
				else {
					cout << "| #Warning : not exist, type" << endl;
				}
				break;
			}
		}
	}
	else{
		cout << "| #Warning : File does not exist, file name = " << file << endl;
	}
	inputFile.close();
	
}






void initialConditionSettingVel(string file,
	vector<string>& type, vector<double>& value){
	
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading){
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("initialField") != string::npos ){
					startReading=true;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			if(read[line].find("type") != string::npos){
				istringstream iss2(read[line]);
				string dummy;
				string stmp;
				iss2 >> dummy >> stmp;
				stmp.erase(stmp.find(";"),1); 
				type.push_back(stmp);
				if(type.back() == "fixedValue"){
					istringstream iss2(read[line+1]);
					string dummy1;
					string dummy2;
					string dummy3;
					iss2 >> dummy >> dummy1 >> dummy2 >> dummy3;
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy3.erase(dummy3.find(";"),1);
					// cout << dtmp0 << dtmp1 << dtmp2 << endl;
					// dtmp.erase(dtmp.find(";"),1); 
					// value.push_back(dtmp0);
					value.push_back(stod(dummy1));
					value.push_back(stod(dummy2));
					value.push_back(stod(dummy3));
				}
				else if(type.back() == "boxToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream box(read[line+2]);
					istringstream boxValue(read[line+3]);
					string dDefaultVal1;
					string dDefaultVal2;
					string dDefaultVal3;
					double xyzMin[3];
					double xyzMax[3];
					string dBoxValue1;
					string dBoxValue2;
					string dBoxValue3;
					
					defaultVal >> dummy >> dDefaultVal1 >> dDefaultVal2 >> dDefaultVal3 ;
					// cout << dDefaultVal1 << endl;
					dDefaultVal1.erase(dDefaultVal1.find("("),1); 
					dDefaultVal3.erase(dDefaultVal3.find(")"),1); 
					dDefaultVal3.erase(dDefaultVal3.find(";"),1);
					
					string dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
					box >> dummy >> dummy1 >> dummy2 >> dummy3 >> dummy4 >> dummy5 >> dummy6;
					
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy4.erase(dummy4.find("("),1);
					dummy6.erase(dummy6.find(")"),1);
					dummy6.erase(dummy6.find(";"),1);
					// cout << dummy1 << dummy2 << dummy3 << endl;
					
					xyzMin[0] = stod(dummy1);
					xyzMin[1] = stod(dummy2);
					xyzMin[2] = stod(dummy3);
					xyzMax[0] = stod(dummy4);
					xyzMax[1] = stod(dummy5);
					xyzMax[2] = stod(dummy6);
					
					
					boxValue >> dummy >> dBoxValue1 >> dBoxValue2 >> dBoxValue3;
					// cout << "AAAAAAAAAAAAAA" << endl;
					dBoxValue1.erase(dBoxValue1.find("("),1); 
					dBoxValue3.erase(dBoxValue3.find(")"),1); 
					dBoxValue3.erase(dBoxValue3.find(";"),1);
					
					// cout << dBoxValue1 << endl;
					
					value.push_back(stod(dDefaultVal1));
					value.push_back(stod(dDefaultVal2));
					value.push_back(stod(dDefaultVal3));
					value.push_back(xyzMin[0]);
					value.push_back(xyzMin[1]);
					value.push_back(xyzMin[2]);
					value.push_back(xyzMax[0]);
					value.push_back(xyzMax[1]);
					value.push_back(xyzMax[2]);
					value.push_back(stod(dBoxValue1));
					value.push_back(stod(dBoxValue2));
					value.push_back(stod(dBoxValue3));
				}
				else {
					cout << "| #Warning : not exist, type" << endl;
				}
				break;
			}
		}
	}
	else{
		cout << "| #Warning : File does not exist, file name = " << file << endl;
	}
	inputFile.close();
	
}






void saveInitialmkdirs(char *dir_path){
	char buff[1024];
	char *p_dir = buff;
	
	strcpy(buff, dir_path);
	buff[1024-1] = '\0';
	
	while(*p_dir){
		if('/'==*p_dir){
			*p_dir = '\0';
			mkdir(buff,0777);
			*p_dir = '/';
		}
		p_dir ++;
	}
	
}


void saveInitialField(string folder){


	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	saveInitialmkdirs(folder_name);

	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << folder_name << "plot...) ... ";
	}
	
	outputFile.open(filenamePlot);
	
	outputFile.precision(20);
	
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"pointLevels\" format=\"ascii\">" << endl;
	for(auto& point : mesh.points) outputFile << point.level << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "    </PointData>" << endl;
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[0] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[1] << " " 
		<< cell.var[2] << " " << cell.var[3] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[4] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[5] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellLevels\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.level << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellGroups\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.group << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	
	
	outputFile << "    </CellData>" << endl;
	
	
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.points){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto& cell : mesh.cells){
		outputFile << cell.faces.size() << endl;
		for(auto& i : cell.faces){
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
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.faces.size();
		for(auto& i : cell.faces){
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
	

	// additional informations
	outputFile << " <owner>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.owner << " ";
	}
	outputFile << endl;
	outputFile << " </owner>" << endl;
	
	outputFile << " <neighbour>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.neighbour << " ";
	}
	outputFile << endl;
	outputFile << " </neighbour>" << endl;
	
	// outputFile << " <faceLevels>" << endl;
	// for(auto& face : mesh.faces){
		// outputFile << face.level << " ";
	// }
	// outputFile << endl;
	// outputFile << " </faceLevels>" << endl;
	
	outputFile << " <faceGroups>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.group << " ";
	}
	outputFile << endl;
	outputFile << " </faceGroups>" << endl;
	
	
	outputFile << " <bcName>" << endl;
	for(auto& boundary : mesh.boundary){
		// cout << boundary.name << endl;
		// read.trim;
		string bcName = boundary.name;
		
		bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
		bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		
		
		// outputFile << boundary.name << " ";
		outputFile << bcName << " ";
	}
	outputFile << endl;
	outputFile << " </bcName>" << endl;
	
	outputFile << " <bcStartFace>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.startFace << " ";
	}
	outputFile << endl;
	outputFile << " </bcStartFace>" << endl;
	
	outputFile << " <bcNFaces>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.nFaces << " ";
	}
	outputFile << endl;
	outputFile << " </bcNFaces>" << endl;
	
	outputFile << " <bcNeighbProcNo>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.neighbProcNo << " ";
	}
	outputFile << endl;
	outputFile << " </bcNeighbProcNo>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	

	if(rank==0){
		string filenamePvtu = folder + "../plot.0.pvtu";
		
		outputFile.open(filenamePvtu);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		// string out_line;
		outputFile << "<?xml version=\"1.0\"?>" << endl;
		outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		outputFile << "  <PUnstructuredGrid>" << endl;

		outputFile << "   <PFieldData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"TimeValue\"/>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0; ip<size; ++ip){
			string filenamevtus = "./0/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"pointLevels\"/>" << endl;
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellLevels\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellGroups\"/>" << endl;
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}








void loadOnlyMeshVtu(
	SEMO_Mesh_Builder &mesh, 
	string folder){
	
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
	
	// cout << "AAAAAAAA" << endl;
	
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
	
	// cout << "BBBBBBBB" << endl;
	n=0;
	for(auto& startFace : bcStartFace){
		
		mesh.addBoundary();
		
		// read.trim;
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
	
	
	
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute load boundary property files ... ";
	}
	
	vector<string> saveNameP;
	vector<string> saveTypeP;
	vector<double> saveValueP;
	load.boundaryProperties("./constant/pressure", mesh, saveNameP, saveTypeP, saveValueP);
	
	vector<string> saveNameT;
	vector<string> saveTypeT;
	vector<double> saveValueT;
	load.boundaryProperties("./constant/temperature", mesh, saveNameT, saveTypeT, saveValueT);
	

	vector<vector<string>> saveNameVF;
	vector<vector<string>> saveTypeVF;
	vector<vector<double>> saveValueVF;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		vector<string> saveNameVF0;
		vector<string> saveTypeVF0;
		vector<double> saveValueVF0;
		string tmpname = "./constant/" + species[ns].name;
		load.boundaryProperties(tmpname, mesh, saveNameVF0, saveTypeVF0, saveValueVF0);
		saveNameVF.push_back(saveNameVF0);
		saveTypeVF.push_back(saveTypeVF0);
		saveValueVF.push_back(saveValueVF0);
	}
	
	
	vector<string> saveNameVel;
	vector<string> saveTypeVel;
	vector<vector<double>> saveValueVel;
	load.boundaryVelocities("./constant/velocities", mesh, saveNameVel, saveTypeVel, saveValueVel);
	
	int neq = 6;
	for(int i=0; i<mesh.boundary.size(); ++i){
		
		if(mesh.boundary[i].neighbProcNo == -1){
			
			mesh.boundary[i].type.resize(neq,"");
			mesh.boundary[i].var.resize(neq,0.0);
			
			bool nameMatching = false;
			for(int j=0; j<saveNameP.size(); ++j){
				if(mesh.boundary[i].name == saveNameP[j]){
					mesh.boundary[i].type[0] = read.trim(saveTypeP[j]);
					mesh.boundary[i].var[0] = saveValueP[j];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			nameMatching = false;
			for(int j=0; j<saveNameVel.size(); ++j){
				if(mesh.boundary[i].name == saveNameVel[j]){
					mesh.boundary[i].type[1] = read.trim(saveTypeVel[j]);
					mesh.boundary[i].type[2] = read.trim(saveTypeVel[j]);
					mesh.boundary[i].type[3] = read.trim(saveTypeVel[j]);
					mesh.boundary[i].var[1] = saveValueVel[j][0];
					mesh.boundary[i].var[2] = saveValueVel[j][1];
					mesh.boundary[i].var[3] = saveValueVel[j][2];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			nameMatching = false;
			for(int j=0; j<saveNameT.size(); ++j){
				if(mesh.boundary[i].name == saveNameT[j]){
					mesh.boundary[i].type[4] = read.trim(saveTypeT[j]);
					mesh.boundary[i].var[4] = saveValueT[j];
					nameMatching = true;
					break;
				}
			}
			if(nameMatching==false){
				cerr << "| #Error : boundary name matching failure" << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			

			for(int ns=0; ns<controls.nSp-1; ++ns){
				nameMatching = false;
				for(int j=0; j<saveNameVF[ns].size(); ++j){
					if(mesh.boundary[i].name == saveNameVF[ns][j]){
						mesh.boundary[i].type[5+ns] = read.trim(saveTypeVF[ns][j]);
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
	

	mesh.informations();

	// mesh.saveFile("vtu");
		
	
	
	
}





