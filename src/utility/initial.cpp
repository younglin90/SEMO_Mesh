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
#include <dlfcn.h>

using namespace std;
#include "../mesh/build.h" 
#include "../load/load.h" 
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
// SEMO_Utility_Read read;
SEMO_Controls_Builder controls;
vector<SEMO_Species> species;
SEMO_Mesh_Save save;
SEMO_Mesh_Geometric geometric;
SEMO_Utility_Math math;

int main(int argc, char* argv[]) {
	
	// cout << "AAA" << endl;

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	

	// cout << "BBB" << endl;
	
		
	// cout << "AAA" << endl;
	controls.readSpecies(species);
	
	// cout << "AAA" << endl;
	controls.readConfigures();
	
	controls.setValues(species);
	
	SEMO_MPI_Builder mpi;
	

	// cout << "BBB" << endl; 
	
	// load.vtu("./grid/0/", mesh, controls, species);

	loadOnlyMeshVtu(mesh, "./grid/0/");

	// geometric.init(mesh);
	
	// math.initLeastSquare(mesh);
	// math.initLeastSquare2nd(mesh); 
	
	
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
		
		cell.var.resize(controls.nTotalCellVar,0.0);
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
	
	
	// cout << "AAA" << endl;
	
	
	vector<string> type;
	vector<double> values;
	
	
	//====== initial conditions ====
	
	
	initialConditionSetting("./constant/pressure", type, values);
	
	// cout << "BBB" << endl;
	
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
	// UDF 펑션
	else if(type.back() == "function")
	{
		void *handle = dlopen("./constant/initialFunctions.so", RTLD_NOW);
		char *error = nullptr;
		if (handle) {
			typedef void (*setFunc_t)(double, double, double, double, double&);
			setFunc_t setFunction;
			*(void **) (&setFunction) = dlsym(handle, "setFunctionPressure");
			for(int i=0; i<mesh.cells.size(); ++i){
				SEMO_Cell& cell = mesh.cells[i];
				double phi;
				setFunction(0.0, cell.x, cell.y, cell.z, phi);
				cell.var[0] = phi;
			}
		}
		dlclose(handle);
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
	else if(type.back() == "cylinderToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double radius = values[9];
			double deValue1 = values[0];
			double deValue2 = values[1];
			double deValue3 = values[2];
			double value1 = values[10];
			double value2 = values[11];
			double value3 = values[12];
			
			double pVec[3];
			double vec0[3];
			double vec1[3];
			
			pVec[0] = values[6]-values[3];
			pVec[1] = values[7]-values[4];
			pVec[2] = values[8]-values[5];
			
			vec0[0] = cellXYZ[i][0]-values[3];
			vec0[1] = cellXYZ[i][1]-values[4];
			vec0[2] = cellXYZ[i][2]-values[5];
			
			vec1[0] = cellXYZ[i][0]-values[6];
			vec1[1] = cellXYZ[i][1]-values[7];
			vec1[2] = cellXYZ[i][2]-values[8];
			
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
			
				cell.var[1] = value1;
				cell.var[2] = value2;
				cell.var[3] = value3;
			}
			else{
				cell.var[1] = deValue1;
				cell.var[2] = deValue2;
				cell.var[3] = deValue3;
			}
		}
	}
	else if(type.back() == "sphereToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double radius = values[6];
			double deValue1 = values[0];
			double deValue2 = values[1];
			double deValue3 = values[2];
			double value1 = values[7];
			double value2 = values[8];
			double value3 = values[9];
			
			double origin0[3];
			origin0[0] = cellXYZ[i][0]-values[3];
			origin0[1] = cellXYZ[i][1]-values[4];
			origin0[2] = cellXYZ[i][2]-values[5];
			
			bool boolInternal = true;
			if(pow(origin0[0],2.0)+pow(origin0[1],2.0)+pow(origin0[2],2.0) > pow(radius,2.0)) boolInternal = false;
			
			if(boolInternal){
				cell.var[1] = value1;
				cell.var[2] = value2;
				cell.var[3] = value3;
			}
			else{
				cell.var[1] = deValue1;
				cell.var[2] = deValue2;
				cell.var[3] = deValue3;
			}
		}
	}
	// UDF 펑션
	else if(type.back() == "function")
	{
		void *handle = dlopen("./constant/initialFunctions.so", RTLD_NOW);
		char *error = nullptr;
		if (handle) {
			typedef void (*setFunc_t)(double, double, double, double, double&);
			setFunc_t setFunctionU, setFunctionV, setFunctionW;
			*(void **) (&setFunctionU) = dlsym(handle, "setFunctionUvelocity");
			*(void **) (&setFunctionV) = dlsym(handle, "setFunctionVvelocity");
			*(void **) (&setFunctionW) = dlsym(handle, "setFunctionWvelocity");
			for(int i=0; i<mesh.cells.size(); ++i){
				SEMO_Cell& cell = mesh.cells[i];
				double phiU, phiV, phiW;
				setFunctionU(0.0, cell.x, cell.y, cell.z, phiU);
				setFunctionV(0.0, cell.x, cell.y, cell.z, phiV);
				setFunctionW(0.0, cell.x, cell.y, cell.z, phiW);
				cell.var[1] = phiU;
				cell.var[2] = phiV;
				cell.var[3] = phiW;
			}
		}
		dlclose(handle);
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
	// UDF 펑션
	else if(type.back() == "function")
	{
		void *handle = dlopen("./constant/initialFunctions.so", RTLD_NOW);
		char *error = nullptr;
		if (handle) {
			typedef void (*setFunc_t)(double, double, double, double, double&);
			setFunc_t setFunction;
			*(void **) (&setFunction) = dlsym(handle, "setFunctionTemperature");
			for(int i=0; i<mesh.cells.size(); ++i){
				SEMO_Cell& cell = mesh.cells[i];
				double phi;
				setFunction(0.0, cell.x, cell.y, cell.z, phi);
				cell.var[4] = phi;
			}
		}
		dlclose(handle);
	}
	
	
	
	// cout << "AAA" << endl;
	
	string speciesFileName = "./constant/" + species[0].name;
	initialConditionSetting(speciesFileName, type, values);
	
	
	// cout << "BBB" << endl;
	if(type.back() == "fixedValue"){
		for(auto& cell : mesh.cells){
			cell.var[controls.MF[0]] = values[0];
			cell.var[controls.VF[0]] = values[0];
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
				cell.var[controls.MF[0]] = values[7];
				cell.var[controls.VF[0]] = values[7];
					// cout.precision(20);
				// cout << "YES " << cellXYZ[i][1] << endl;
			}
			else{
				// if( cellXYZ[i][1] < 0.31 && cellXYZ[i][1] > 0.29 ){
					// cout.precision(20);
					// cout << "NO " << cellXYZ[i][1] << endl;
				// }
				cell.var[controls.MF[0]] = values[0];
				cell.var[controls.VF[0]] = values[0];
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
				cell.var[controls.MF[0]] = value;
				cell.var[controls.VF[0]] = value;
			}
			else{
				cell.var[controls.MF[0]] = deValue;
				cell.var[controls.VF[0]] = deValue;
			}
		}
	}
	else if(type.back() == "sphereToCell"){
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double radius = values[4];
			double deValue = values[0];
			double value = values[5];
			
			double origin0[3];
			origin0[0] = cellXYZ[i][0]-values[1];
			origin0[1] = cellXYZ[i][1]-values[2];
			origin0[2] = cellXYZ[i][2]-values[3];
			
			bool boolInternal = true;
			if(pow(origin0[0],2.0)+pow(origin0[1],2.0)+pow(origin0[2],2.0) > pow(radius,2.0)) boolInternal = false;
			
			if(boolInternal){
				cell.var[controls.MF[0]] = value;
				cell.var[controls.VF[0]] = value;
			}
			else{
				cell.var[controls.MF[0]] = deValue;
				cell.var[controls.VF[0]] = deValue;
			}
		}
	}
	// UDF 펑션
	else if(type.back() == "function")
	{
		void *handle = dlopen("./constant/initialFunctions.so", RTLD_NOW);
		char *error = nullptr;
		if (handle) {
			typedef void (*setFunc_t)(double, double, double, double, double&);
			setFunc_t setFunction;
			*(void **) (&setFunction) = dlsym(handle, "setFunctionMassFraction");
			for(int i=0; i<mesh.cells.size(); ++i){
				SEMO_Cell& cell = mesh.cells[i];
				double phi;
				setFunction(0.0, cell.x, cell.y, cell.z, phi);
				cell.var[controls.MF[0]] = phi;
				cell.var[controls.VF[0]] = phi;
			}
		}
		dlclose(handle);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		// // 2차원 오실레이션 드랍
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			
			// bool boolLiquid = false;
			
			
			// // double X = cellXYZ[i][0]-0.500300000000000000000001; double Y = cellXYZ[i][1]-0.500000000000000000001;
			// // double X = cellXYZ[i][0]-0.025; double Y = cellXYZ[i][1]-0.025;
			// // double X = cellXYZ[i][0]-0.5*0.00125; double Y = cellXYZ[i][1]-0.5*0.00125;
			// double X = cellXYZ[i][0]-0.000625375; double Y = cellXYZ[i][1]-0.000625;
			// double theta = atan2(abs(Y),abs(X));
			// if(X>0.0 && Y>0.0){
			// }
			// else if(X<=0.0 && Y>0.0){
				// theta = 3.141592-theta;
			// }
			// else if(X<=0.0 && Y<=0.0){
				// theta = 3.141592+theta;
			// }
			// else if(X>0.0 && Y<=0.0){
				// theta = 2.0*3.141592-theta;
			// }
			
			// // double R0 = 0.2;
			// // double R0 = 0.1;
			// // double R0 = 0.005;
			// double R0 = 1.25e-4;
			// // double Eta = 0.2;
			// double Eta = 0.05;
			// double N = 2.0;
			
			// // double dagak = sqrt(pow(cell.x-0.5,2.0) + pow(cell.y-0.5,2.0)+1.e-200);
			// // double theta = acos(cell.x/dagak);
			// double R = R0*(1.0+Eta*cos(N*theta)-0.25*Eta*Eta);
			// // cout << X << " " << Y << " " << theta << " " << R << endl;
			// if( pow(X,2.0) + pow(Y,2.0) <= R*R ){
				// boolLiquid = true;
			// }
			
			
			// if(boolLiquid){
				// cell.var[controls.MF[0]] = 1.0;
				// cell.var[controls.VF[0]] = 1.0;
			// }
			// else{
				// cell.var[controls.MF[0]] = 0.0;
				// cell.var[controls.VF[0]] = 0.0;
			// }
		// }
		
		
		
		// // 3차원 오실레이션 드랍
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			
			// bool boolLiquid = false;
			
			
			// double X = cellXYZ[i][0]-0.5; double Y = cellXYZ[i][1]-0.5;
			// double Z = cellXYZ[i][2]-0.5;
			
			// double R0 = 0.2 * 1.192;
			// double R1 = 0.2 * 1.784;
			
			// if( (pow(Z,2.0) + pow(X,2.0))/R0/R0 + pow(Y,2.0)/R1/R1 < 1.0 ){
				// boolLiquid = true;
			// }
			
			// if(boolLiquid){
				// cell.var[controls.MF[0]] = 1.0;
				// cell.var[controls.VF[0]] = 1.0;
			// }
			// else{
				// cell.var[controls.MF[0]] = 0.0;
				// cell.var[controls.VF[0]] = 0.0;
			// }
		// }
	
	
	
	
		// for(int i=0; i<mesh.cells.size(); ++i){
			// SEMO_Cell& cell = mesh.cells[i];
			
			// double radius = values[4];
			// double deValue = values[0];
			// double value = values[5];
			
			// radius = 0.005;
			
			// double origin0[3];
			// origin0[0] = cellXYZ[i][0]-0.0;
			// origin0[1] = cellXYZ[i][1]-0.025;
			// origin0[2] = cellXYZ[i][2]-0.0;
			
			// bool boolInternal = false;
			// if(pow(origin0[0],2.0)+pow(origin0[1],2.0)+pow(origin0[2],2.0) <= pow(radius,2.0)) boolInternal = true;
			
			// origin0[0] = cellXYZ[i][0]-0.0085;
			// origin0[1] = cellXYZ[i][1]-0.01;
			// origin0[2] = cellXYZ[i][2]-0.0;
			// if(pow(origin0[0],2.0)+pow(origin0[1],2.0)+pow(origin0[2],2.0) <= pow(radius,2.0)) boolInternal = true;
			
			// if(boolInternal){
				// cell.var[5] = value;
			// }
			// else{
				// cell.var[5] = deValue;
			// }
			
		// }

	// SEMO_Utility_Math math;
	// SEMO_Mesh_Geometric geometric;
	// geometric.init(mesh);
	
	// cout << "AAAAAA" << endl;
	
	saveInitialField("./save/0/");
	// save.vtu("./save/0/", mesh, controls, species);
	
	// cout << "BBBBBB" << endl;

	if(rank==0) cout << "| completed save initial files : ./save/0/" << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	MPI_Finalize();
	// exit(EXIT_FAILURE);
	exit(EXIT_SUCCESS);
	// return 0;
	return EXIT_SUCCESS;
	
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
				else if(type.back() == "sphereToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream points(read[line+2]);
					istringstream radius(read[line+3]);
					istringstream sphereValue(read[line+4]);
					double dDefaultVal;
					double xyzOrigin[3];
					double radiusVal;
					double sphereVal;
					defaultVal >> dummy >> dDefaultVal;
					
					string dummy1, dummy2, dummy3;
					points >> dummy >> dummy1 >> dummy2 >> dummy3;
					
					
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy3.erase(dummy3.find(";"),1);
					// cout << dummy4 << dummy5 << dummy6 << endl;
					
					xyzOrigin[0] = stod(dummy1);
					xyzOrigin[1] = stod(dummy2);
					xyzOrigin[2] = stod(dummy3);
					
					radius >> dummy >> radiusVal;
					sphereValue >> dummy >> sphereVal;
					
					value.push_back(dDefaultVal);
					value.push_back(xyzOrigin[0]);
					value.push_back(xyzOrigin[1]);
					value.push_back(xyzOrigin[2]);
					value.push_back(radiusVal);
					value.push_back(sphereVal);
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
				else if(type.back() == "cylinderToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream points(read[line+2]);
					istringstream radius(read[line+3]);
					istringstream cylinderValue(read[line+4]);
					double dDefaultVal;
					string dDefaultVal1;
					string dDefaultVal2;
					string dDefaultVal3;
					double xyzMin[3];
					double xyzMax[3];
					double radiusVal;
					double cylinderVal;
					string cylinderVal1;
					string cylinderVal2;
					string cylinderVal3;
					// defaultVal >> dummy >> dDefaultVal;
					defaultVal >> dummy >> dDefaultVal1 >> dDefaultVal2 >> dDefaultVal3 ;
					// cout << dDefaultVal1 << endl;
					dDefaultVal1.erase(dDefaultVal1.find("("),1); 
					dDefaultVal3.erase(dDefaultVal3.find(")"),1); 
					dDefaultVal3.erase(dDefaultVal3.find(";"),1);
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
					// cylinderValue >> dummy >> cylinderVal;
					cylinderValue >> dummy >> cylinderVal1 >> cylinderVal2 >> cylinderVal3 ;
					cylinderVal1.erase(cylinderVal1.find("("),1); 
					cylinderVal3.erase(cylinderVal3.find(")"),1); 
					cylinderVal3.erase(cylinderVal3.find(";"),1);
					
					value.push_back(stod(dDefaultVal1));
					value.push_back(stod(dDefaultVal2));
					value.push_back(stod(dDefaultVal3));
					value.push_back(xyzMin[0]);
					value.push_back(xyzMin[1]);
					value.push_back(xyzMin[2]);
					value.push_back(xyzMax[0]);
					value.push_back(xyzMax[1]);
					value.push_back(xyzMax[2]);
					value.push_back(radiusVal);
					value.push_back(stod(cylinderVal1));
					value.push_back(stod(cylinderVal2));
					value.push_back(stod(cylinderVal3));
				}
				else if(type.back() == "sphereToCell"){
					istringstream defaultVal(read[line+1]);
					istringstream points(read[line+2]);
					istringstream radius(read[line+3]);
					istringstream sphereValue(read[line+4]);
					string dDefaultVal1,dDefaultVal2,dDefaultVal3;
					double xyzOrigin[3];
					double radiusVal;
					string sphereVal1,sphereVal2,sphereVal3;
					defaultVal >> dummy >> dDefaultVal1 >> dDefaultVal2 >> dDefaultVal3;
					dDefaultVal1.erase(dDefaultVal1.find("("),1); 
					dDefaultVal3.erase(dDefaultVal3.find(")"),1); 
					dDefaultVal3.erase(dDefaultVal3.find(";"),1);
					
					string dummy1, dummy2, dummy3, dummy4, dummy5, dummy6;
					points >> dummy >> dummy1 >> dummy2 >> dummy3;
					
					
					dummy1.erase(dummy1.find("("),1); 
					dummy3.erase(dummy3.find(")"),1); 
					dummy6.erase(dummy3.find(";"),1);
					// cout << dummy4 << dummy5 << dummy6 << endl;
					
					xyzOrigin[0] = stod(dummy1);
					xyzOrigin[1] = stod(dummy2);
					xyzOrigin[2] = stod(dummy3);
					
					radius >> dummy >> radiusVal;
					sphereValue >> dummy >> sphereVal1 >> sphereVal2 >> sphereVal3 ;
					sphereVal1.erase(sphereVal1.find("("),1); 
					sphereVal3.erase(sphereVal3.find(")"),1); 
					sphereVal3.erase(sphereVal3.find(";"),1);
					
					value.push_back(stod(dDefaultVal1));
					value.push_back(stod(dDefaultVal2));
					value.push_back(stod(dDefaultVal3));
					value.push_back(xyzOrigin[0]);
					value.push_back(xyzOrigin[1]);
					value.push_back(xyzOrigin[2]);
					value.push_back(radiusVal);
					value.push_back(stod(sphereVal1));
					value.push_back(stod(sphereVal2));
					value.push_back(stod(sphereVal3));
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





void saveInitialField(string folder){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_Utility_Math math;
	
	int compressSize = controls.saveCompression;
	// int compressSize = 6;
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	save.mkdirs(folder_name);

	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(rank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << folder_name << "plot...) ... ";
	}
	
	outputFile.open(filenamePlot);
	
	outputFile.precision( controls.writePrecision );
	
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;

	string saveFormat = controls.saveFormat;
	if(controls.saveFormat == "ascii"){
		outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	}
	else if(controls.saveFormat == "binary"){
		if(controls.saveCompression==0){
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		}
		else{
			outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
		}
	}
	else{
		cout << endl;
		cout << endl;
		cout << "| warning : not defined saveFormat at controlDic file" << endl;
		cout << endl;
		cout << endl;
	}
	
	outputFile << "  <UnstructuredGrid>" << endl;
	

	// Field data
	outputFile << "    <FieldData>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		values.push_back(controls.time);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"pointLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& point : mesh.points) values.push_back(point.level);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "    </PointData>" << endl;
	
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& cell : mesh.cells) values.push_back(cell.var[controls.P]);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& cell : mesh.cells){
			values.push_back(cell.var[controls.U]); 
			values.push_back(cell.var[controls.V]); 
			values.push_back(cell.var[controls.W]); 
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& cell : mesh.cells) values.push_back(cell.var[controls.T]);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" << controls.name[controls.VF[0]] << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& cell : mesh.cells) values.push_back(cell.var[controls.VF[0]]);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"" << controls.name[controls.MF[0]] << "\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& cell : mesh.cells) values.push_back(cell.var[controls.MF[0]]);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// cell levels
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellLevels\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells) values.push_back(cell.level);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// cell groups
	{
		outputFile << "     <DataArray type=\"Int32\" Name=\"cellGroups\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells) values.push_back(cell.group);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	
	
	outputFile << "    </CellData>" << endl;
	
	
	// Points
	outputFile << "    <Points>" << endl;
	{
		outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"" << saveFormat << "\">" << endl;
		vector<double> values;
		for(auto& point : mesh.points) {
			values.push_back(point.x);
			values.push_back(point.y);
			values.push_back(point.z);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			for(auto i : cell.points){
				values.push_back(i);
			}
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// offsets (cell's points offset)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFaceOffset = 0;
		for(auto& cell : mesh.cells){
			cellFaceOffset += cell.points.size();
			values.push_back(cellFaceOffset);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// types (cell's type, 42 = polyhedron)
	{
		outputFile << "    <DataArray type=\"Int32\" Name=\"types\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values(mesh.cells.size(),42);
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faces\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& cell : mesh.cells){
			values.push_back(cell.faces.size());
			for(auto& i : cell.faces){
				values.push_back(mesh.faces[i].points.size());
				for(auto& j : mesh.faces[i].points){
					values.push_back(j);
				}
			}
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	// faceoffsets (cell's face offset)
	{
		outputFile << "    <DataArray type=\"Int32\" IdType=\"1\" Name=\"faceoffsets\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		int cellFacePointOffset = 0;
		for(auto& cell : mesh.cells){
			int numbering = 1 + cell.faces.size();
			for(auto& i : cell.faces){
				numbering += mesh.faces[i].points.size();
			}
			cellFacePointOffset += numbering;
			values.push_back(cellFacePointOffset);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << "     </DataArray>" << endl;
	}
	
	outputFile << "   </Cells>" << endl;
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	

	// additional informations
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			values.push_back(face.owner);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << " </DataArray>" << endl;
	}
	{
		outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		vector<int> values;
		for(auto& face : mesh.faces){
			values.push_back(face.neighbour);
		}
		save.writeDatasAtVTU(controls, outputFile, values);
		outputFile << " </DataArray>" << endl;
	}
	
	

	// boundary informations
	{
		outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundary){
			// cout << boundary.name << endl;
			// trim;
			string bcName = boundary.name;
			
			bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// outputFile << boundary.name << " ";
			outputFile << bcName << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundary){
			outputFile << boundary.startFace << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundary){
			outputFile << boundary.nFaces << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundary){
			outputFile << boundary.neighbProcNo << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
	}
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	
	
	// pvtu file
	if(rank==0){
		string filenamePvtu = "./save/plot.";
		string stime = folder;
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		filenamePvtu += stime;
		filenamePvtu += ".pvtu";
		
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
			string filenamevtus = "./" + stime;
			filenamevtus += "/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;

		outputFile << "    <PDataArray type=\"Int32\" Name=\"pointLevels\"/>" << endl;
		
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"density\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"" << controls.name[controls.VF[0]] << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"" << controls.name[controls.MF[0]] << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int32\" Name=\"cellLevels\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int32\" Name=\"cellGroups\"/>" << endl;
		// }
		outputFile << "   </PCellData>" << endl;
		outputFile << "  </PUnstructuredGrid>" << endl;
		outputFile << "</VTKFile>" << endl;
		
		
		outputFile.close();
		
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
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
	
	// cout << "AAA" << endl;
	int neq = 6;
	for(int i=0; i<mesh.boundary.size(); ++i){
		
		if(mesh.boundary[i].neighbProcNo == -1){
			
			mesh.boundary[i].type.resize(neq,"");
			mesh.boundary[i].var.resize(neq,0.0);
			
			bool nameMatching = false;
			for(int j=0; j<saveNameP.size(); ++j){
				if(mesh.boundary[i].name == saveNameP[j]){
					mesh.boundary[i].type[0] = load.trim(saveTypeP[j]);
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
					mesh.boundary[i].type[1] = load.trim(saveTypeVel[j]);
					mesh.boundary[i].type[2] = load.trim(saveTypeVel[j]);
					mesh.boundary[i].type[3] = load.trim(saveTypeVel[j]);
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
					mesh.boundary[i].type[4] = load.trim(saveTypeT[j]);
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
						mesh.boundary[i].type[5+ns] = load.trim(saveTypeVF[ns][j]);
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





