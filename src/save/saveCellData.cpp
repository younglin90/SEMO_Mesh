// #include <iostream>
#include <cstring>
#include <fstream>
// #include <sstream>
#include <string>
#include <vector>
// #include <algorithm>
// #include <stdexcept>
// #include <iomanip>
#include <mpi.h>
using namespace std;


#include "save.h" 
#include "../mesh/build.h" 
#include "../controls/build.h" 

void SEMO_Mesh_Save::cellDataToMemoryOverTime(
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls){
	

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int extSize = controls.extractNames.size();
	
	// 셀 중심 포인트 x,y,z 추출
	// 구 내에 있는 점들 중, 중심하고 가장 가까운 점의 셀 추출
	vector<double> extCenterX(extSize,0.0);
	vector<double> extCenterY(extSize,0.0);
	vector<double> extCenterZ(extSize,0.0);
	vector<int> cellID(extSize,0);
	vector<double> distMin(extSize,1.e15);
	for(int j=0; j<extSize; ++j){
		double extX = controls.extractCenterPoints[j][0];
		double extY = controls.extractCenterPoints[j][1];
		double extZ = controls.extractCenterPoints[j][2];
		double extR = controls.extractRadii[j];
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			double X = mesh.cells[i].x;
			double Y = mesh.cells[i].y;
			double Z = mesh.cells[i].z;
			double dist = (X-extX)*(X-extX);
			dist += (Y-extY)*(Y-extY);
			dist += (Z-extZ)*(Z-extZ);
			if(dist <= extR*extR){
				if(dist<distMin[j]){
					extCenterX[j] = cell.x;
					extCenterY[j] = cell.y;
					extCenterZ[j] = cell.z;
					cellID[j] = i;
				}
			}
		}
	}
	
	
	// 셀 값 추출 및 저장 (시간, 데이터들...)
	int fieldSize = controls.extractFieldDatas.size();
	int averageSize = controls.extractAverageDatas.size();
	int valueSize = controls.extractCellValueTargets.size();
	
	controls.extractDatas.push_back(vector<double>());
	controls.extractDatas.back().push_back(controls.time);
	
	// 필드 데이터
	for(auto& data : controls.extractFieldDatas){
		if(data == "residual") 
			controls.extractDatas.back().push_back(controls.residual);
	}
	
	
	// 평균 데이터
	double cellVolume = 0.0;
	vector<double> averageData(averageSize,0.0);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		int num = 0;
		cellVolume += cell.volume;
		for(auto& data : controls.extractAverageDatas){
			if(data == "pressure") averageData[num] += cell.var[controls.P]*cell.volume;
			if(data == "x-velocity") averageData[num] += cell.var[controls.U]*cell.volume;
			if(data == "y-velocity") averageData[num] += cell.var[controls.V]*cell.volume;
			if(data == "z-velocity") averageData[num] += cell.var[controls.W]*cell.volume;
			if(data == "temperature") averageData[num] += cell.var[controls.T]*cell.volume;
			if(data == controls.name[controls.VF[0]]) averageData[num] += cell.var[controls.VF[0]]*cell.volume;
			if(data == controls.name[controls.MF[0]]) averageData[num] += cell.var[controls.MF[0]]*cell.volume;
			++num;
		}
	}
	vector<double> averageData_global(averageSize,0.0);
	MPI_Allreduce(averageData.data(), averageData_global.data(), averageSize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double totalCellVolume;
	MPI_Allreduce(&cellVolume, &totalCellVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	for(auto& data : averageData_global) 
		controls.extractDatas.back().push_back(data/totalCellVolume);
	
	// 셀 데이터
	for(auto& i : cellID){
		auto& cell = mesh.cells[i];
		for(auto& data : controls.extractCellValueTargets){
			if(data == "pressure") 
				controls.extractDatas.back().push_back(cell.var[controls.P]);
			if(data == "x-velocity") 
				controls.extractDatas.back().push_back(cell.var[controls.U]);
			if(data == "y-velocity") 
				controls.extractDatas.back().push_back(cell.var[controls.V]);
			if(data == "z-velocity") 
				controls.extractDatas.back().push_back(cell.var[controls.W]);
			if(data == "temperature") 
				controls.extractDatas.back().push_back(cell.var[controls.T]);
		}
	}

	// cout << valueSize << endl;
	
	
}


void SEMO_Mesh_Save::cellDataOverTime(
	string folder, 
	SEMO_Controls_Builder &controls){
		
		
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	int fieldSize = controls.extractFieldDatas.size();
	int averageSize = controls.extractAverageDatas.size();
	int extSize = controls.extractNames.size();
	int valueSize = controls.extractCellValueTargets.size();
	int totalSize = 1 + fieldSize + averageSize + valueSize;
	int extSizeTime = controls.extractDatas.size();
	
	int formatPrecision = 22;

// cout << controls.extractNames[0] << endl;

	if(rank==0){
		string fileName = "field_datas.csv";
		ofstream outputFile;
		string filenamePlot = folder + fileName;
		if(rank==0){
			cout << "| execute save file (" << folder_name << "extract field data...) ... " << endl;
		}
		outputFile.open(filenamePlot);
		outputFile.precision( formatPrecision );
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		outputFile << "time";
		for(auto& data : controls.extractFieldDatas){
			if(data == "residual") outputFile << ", " << "residual";
		}
		outputFile << endl;
		for(int j=0; j<extSizeTime; ++j) {
			outputFile << scientific << controls.extractDatas[j][0];
			for(int k=0; k<fieldSize; ++k){
				int tmp_num = k+1;
				outputFile << scientific << ", " << controls.extractDatas[j][tmp_num];
			}
			outputFile << endl;
		}
		outputFile.close();
	}

	if(rank==0){
		string fileName = "average_datas.csv";
		ofstream outputFile;
		string filenamePlot = folder + fileName;
		if(rank==0){
			cout << "| execute save file (" << folder_name << "extract average data...) ... " << endl;
		}
		outputFile.open(filenamePlot);
		outputFile.precision( formatPrecision );
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		outputFile << "time";
		for(auto& data : controls.extractAverageDatas){
			if(data == "pressure") 
				outputFile << ", " << "pressure";
			if(data == "x-velocity") 
				outputFile << ", " << "x-velocity";
			if(data == "y-velocity") 
				outputFile << ", " << "y-velocity";
			if(data == "z-velocity") 
				outputFile << ", " << "z-velocity";
			if(data == "temperature") 
				outputFile << ", " << "temperature";
			if(data == controls.name[controls.VF[0]]) 
				outputFile << ", " << controls.name[controls.VF[0]];
			if(data == controls.name[controls.MF[0]]) 
				outputFile << ", " << controls.name[controls.MF[0]];
		}
		outputFile << endl;
		for(int j=0; j<extSizeTime; ++j) {
			outputFile << scientific << controls.extractDatas[j][0];
			for(int k=0; k<averageSize; ++k){
				int tmp_num = k+1 + fieldSize;
				outputFile << scientific << ", " << controls.extractDatas[j][tmp_num];
			}
			outputFile << endl;
		}
		outputFile.close();
		
	}



	for(int i=0, num=0; i<extSize; ++i){
		string fileName;
		std::ostringstream streamXYZ;
		// streamXYZ << "point_";
		// streamXYZ << i;
		streamXYZ << controls.extractNames[i];
		fileName = streamXYZ.str() + ".csv";
		{
			ofstream outputFile;
			string filenamePlot = folder + fileName;
			if(rank==0){
				cout << "| execute save file (" << folder_name << "extract cell data...) ... " << endl;
			}
			outputFile.open(filenamePlot);
			outputFile.precision( formatPrecision );
			if(outputFile.fail()){
				cerr << "Unable to write file for writing." << endl;
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			outputFile << "time";
			for(auto& data : controls.extractCellValueTargets){
				if(data == "pressure") 
					outputFile << ", " << "pressure";
				if(data == "x-velocity") 
					outputFile << ", " << "x-velocity";
				if(data == "y-velocity") 
					outputFile << ", " << "y-velocity";
				if(data == "z-velocity") 
					outputFile << ", " << "z-velocity";
				if(data == "temperature") 
					outputFile << ", " << "temperature";
			}
			outputFile << endl;
			for(int j=0; j<extSizeTime; ++j) {
				outputFile << scientific << controls.extractDatas[j][0];
				for(int k=0; k<valueSize; ++k){
					int tmp_num = k+1 + fieldSize + averageSize + num*valueSize;
					outputFile << scientific << ", " << controls.extractDatas[j][tmp_num];
				}
				outputFile << endl;
			}
			outputFile.close();
		}

		++num;
	}
	
	controls.extractDatas.clear();
	
	MPI_Barrier(MPI_COMM_WORLD); 
	
}