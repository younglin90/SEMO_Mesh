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
void saveField(
	string folder, SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species,
	int rankTarget);
void loadVTU(
	string folder,
	SEMO_Mesh_Builder &mesh,
	int rankTarget);
void loadSaveVTU(
	string folder, 
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species,
	int rankTarget);

// SEMO_Mesh_Builder mesh;
// SEMO_Mesh_Builder meshObj;
SEMO_Mesh_Load load;
SEMO_Utility_Read read;
SEMO_Controls_Builder controls;
vector<SEMO_Species> species;
SEMO_Mesh_Save save;

int main(int argc, char* argv[]) {
	
	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	
	
	map<string,string> mapArgv;
	
	for(int i=1; i<argc; i+=2){
		string first = argv[i];
		string second;
		if(i+1==argc){
			second = "nan";
		}
		else{
			second = argv[i+1];
		}
		mapArgv.insert(make_pair(first, second));
		// cout <<  first<< " " << second << endl;
	}
	
	if( 
	mapArgv.find("-help") != mapArgv.end() ||
	mapArgv.find("-h") != mapArgv.end() ||
	mapArgv.size() == 0
	){
		if(rank==0){
			cout << endl;
			cout << endl;
			cout << "┌─────── MapField helper ─────────────────────────────── " << endl;
			cout << "| -s \"source grid par num.\"  : # of par" << endl;
			cout << "| -sd \"source grid dir.\"     : directory of par" << endl;
			cout << "| -t \"target grid par num.\"  : # of par" << endl;
			cout << "| -td \"target grid dir.\"     : directory of par" << endl;
			cout << "| -r \"distance coeff.\"       : cell to cell distance coeff." << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	
	
	int sizeSource = stoi(mapArgv["-s"]);
	int sizeTarget = stoi(mapArgv["-t"]);
	
	double coeffDist = stod(mapArgv["-r"]);
	
	
	if( sizeTarget != size ){
		if(rank==0){
			cout << endl;
			cout << endl;
			cout << "| #Error : not matching \"target grid size\" /= \"mpi size\" " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute MapField " << endl;
	}
	
	
	controls.readSpecies(species);
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	// SEMO_Solvers_Builder solvers;
	
	SEMO_MPI_Builder mpi;
	
	
	double starttime = stod(controls.startFrom);
	string foldername;
	std::ostringstream streamObj;
	streamObj << starttime;
	// foldername = "./save/" + streamObj.str() + "/";
	// if(starttime == 0.0){
		// foldername = "./save/0/";
	// }
	
	foldername = mapArgv["-sd"];
	
	
	
	SEMO_Mesh_Builder meshTarget;
	loadVTU(mapArgv["-td"], meshTarget, rank);
	
	SEMO_Mesh_Geometric geometric;
	
	geometric.init(meshTarget);
	
	int meshTargetSize = meshTarget.cells.size();
	vector<double> distance(meshTargetSize,1.e+10);
	vector<bool> boolBreakDist(meshTargetSize,false);
	vector<double> distanceMax(meshTargetSize,0.0);
	for(int i=0; i<meshTargetSize; ++i){
		auto& cellTarget = meshTarget.cells[i];
		cellTarget.var.resize(controls.nEq);
		distanceMax[i] = coeffDist*pow(cellTarget.volume,0.333333);
	}
	
	for(int size1=0; size1<sizeSource; ++size1){
		
		
		SEMO_Mesh_Builder meshSource;
		loadSaveVTU(foldername, meshSource, controls, species, size1);
		
		// cout << rank << " " << size1 << endl;
		int meshSourceSize = meshSource.cells.size();
		
		for(int j=0; j<meshSourceSize; ++j){
			auto& cellSource = meshSource.cells[j];
			
			cellSource.x = 0.0;
			cellSource.y = 0.0;
			cellSource.z = 0.0;
			for(auto& j : cellSource.points){
				cellSource.x += meshSource.points[j].x;
				cellSource.y += meshSource.points[j].y;
				cellSource.z += meshSource.points[j].z;
			}
			double pointSize = (double)cellSource.points.size();
			cellSource.x /= pointSize;
			cellSource.y /= pointSize;
			cellSource.z /= pointSize;
		}
		
		
		for(int i=0; i<meshTargetSize; ++i){
			
			if(i % 100 == 0){
				if(rank==0) {
					cout << size1 << " / " << sizeSource << " | " << i << " / " << meshTargetSize << endl;
				}
			}
			
			if(boolBreakDist[i] == true) continue;
		
			auto& cellTarget = meshTarget.cells[i];
			double X = cellTarget.x;
			double Y = cellTarget.y;
			double Z = cellTarget.z;
		
			for(int j=0; j<meshSourceSize; ++j){
				auto& cellSource = meshSource.cells[j];
				double nX = cellSource.x;
				double nY = cellSource.y;
				double nZ = cellSource.z;
				
				double tmpDist = sqrt( (X-nX)*(X-nX) + (Y-nY)*(Y-nY) + (Z-nZ)*(Z-nZ) );
				// double tmpDist = sqrt( pow(X-nX,2.0) + pow(Y-nY,2.0) + pow(Z-nZ,2.0) );
					
				if(tmpDist<distance[i]){
					for(int num=0; num<controls.nEq; ++num){
						cellTarget.var[num] = cellSource.var[num];
					}
					distance[i] = tmpDist;
					
					if(tmpDist < distanceMax[i]) boolBreakDist[i] = true;
				}
			}
		}
	}
	
	
	// // cout << "| completed calc. distance, rank = " << rank << endl;



	saveField("./save/target/", meshTarget, controls, species, rank);
	
	
	// if(rank==0){
		// string folder = "./save/target/";
		// string filenamePvtu = folder + "../plot.target.pvtu";
		
		// ofstream outputFile;
		// outputFile.open(filenamePvtu);
		// if(outputFile.fail()){
			// cerr << "Unable to write file for writing." << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		// // string out_line;
		// outputFile << "<?xml version=\"1.0\"?>" << endl;
		// outputFile << " <VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
		// outputFile << "  <PUnstructuredGrid>" << endl;

		// outputFile << "   <PFieldData>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"TimeValue\"/>" << endl;
		// outputFile << "   </PFieldData>" << endl;
		
		// outputFile << "   <PPoints>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		// outputFile << "   </PPoints>" << endl;
		// for(int ip=0; ip<sizeTarget; ++ip){
			// string filenamevtus = "./target/plot.";
			// filenamevtus += to_string(ip);
			// filenamevtus += ".vtu";
			// outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		// }
		// outputFile << "   <PPointData>" << endl;
		// outputFile << "   </PPointData>" << endl;
		// outputFile << "   <PCellData>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"density\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\"/>" << endl;
		// outputFile << "   </PCellData>" << endl;
		// outputFile << "  </PUnstructuredGrid>" << endl;
		// outputFile << "</VTKFile>" << endl;
		
		// outputFile.close();
	// }
	

	if(rank==0){
		cout << "| completed save vtu files : ./save/target/" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	// if(rank==0) cout << endl << "| completed save vtu files : ./save/target/" << endl << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
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


void saveField(
	string folder, SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species,
	int rankTarget){



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
	
	outputFile.precision( controls.writePrecision );
	
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	

	// Field data
	outputFile << "    <FieldData>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">" << endl;
	outputFile << controls.time << endl;
	outputFile << "     </DataArray>" << endl;
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	
	
	

	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.P] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[controls.U] << " " 
		<< cell.var[controls.V] << " " << cell.var[controls.W] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.T] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.Rho] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	

	vector<string> volFracName;
	string tmpVolFracName;
	tmpVolFracName = "volumeFraction-" + species[0].name;
	volFracName.push_back(tmpVolFracName);
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"" << volFracName[0] << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.VF[0]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	// cell volume
	outputFile << "     <DataArray type=\"Float64\" Name=\"volume\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.volume << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	// UDV
	outputFile << "     <DataArray type=\"Float64\" Name=\"UDV0\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[0]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"UDV1\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[1]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"UDV2\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[2]] << " ";
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
	
	outputFile << " <bcName>" << endl;
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
	
	









	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// char folder_name[1000];
	// strcpy(folder_name, folder.c_str());
	// saveInitialmkdirs(folder_name);

	// ofstream outputFile;
	// string filenamePlot = folder + "plot." + to_string(rankTarget) + ".vtu";
	
	// // if(rank==0){
		// // cout << "┌────────────────────────────────────────────────────" << endl;
		// // cout << "| execute save file (" << folder_name << "plot...) ... ";
	// // }
	
	// outputFile.open(filenamePlot);
	
	// outputFile.precision(20);
	
	// if(outputFile.fail()){
		// cerr << "Unable to write file for writing." << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
	// // string out_line;
	// outputFile << "<?xml version=\"1.0\"?>" << endl;
	// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	// outputFile << "  <UnstructuredGrid>" << endl;
	// outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// // Points data
	// outputFile << "    <PointData>" << endl;
	// outputFile << "    </PointData>" << endl;
	
	
	// // Cells data
	// outputFile << "    <CellData>" << endl;
	
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[0] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells){
		// outputFile << scientific << cell.var[1] << " " 
		// << cell.var[2] << " " << cell.var[3] << " ";
	// }
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[4] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[5] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	
	
	
	
	// outputFile << "    </CellData>" << endl;
	
	
	// // Points
	// outputFile << "    <Points>" << endl;
	// // }
	// outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	// stringstream streamXYZ;
	// for(auto& point : mesh.points){
		// outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	// }
	
	// outputFile << "    </DataArray>" << endl;
	// outputFile << "   </Points>" << endl;
	
	// // cells
	// outputFile << "   <Cells>" << endl; 
	// // connectivity (cell's points)
	// outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	// for(auto& cell : mesh.cells){
		// for(auto i : cell.points){
			// outputFile << i << " ";
		// }
		// outputFile << endl;
	// }
	
	// outputFile << "    </DataArray>" << endl;
	
	// // offsets (cell's points offset)
	// int cellFaceOffset = 0;
	// outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	// cellFaceOffset = 0;
	// for(auto& cell : mesh.cells){
		// cellFaceOffset += cell.points.size();
		// outputFile << cellFaceOffset << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // types (cell's type, 42 = polyhedron)
	// outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	// for(auto& cell : mesh.cells){
		// outputFile << "42" << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // faces (cell's faces number, each face's point number, cell's faces's points)
	// outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// // outputFile << mesh.faces.size() << endl;
	// for(auto& cell : mesh.cells){
		// outputFile << cell.faces.size() << endl;
		// for(auto& i : cell.faces){
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
	// for(auto& cell : mesh.cells){
		// int numbering = 1 + cell.faces.size();
		// for(auto& i : cell.faces){
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
	

	// // additional informations
	// outputFile << " <owner>" << endl;
	// for(auto& face : mesh.faces){
		// outputFile << face.owner << " ";
	// }
	// outputFile << endl;
	// outputFile << " </owner>" << endl;
	
	// outputFile << " <neighbour>" << endl;
	// for(auto& face : mesh.faces){
		// outputFile << face.neighbour << " ";
	// }
	// outputFile << endl;
	// outputFile << " </neighbour>" << endl;
	
	// outputFile << " <bcName>" << endl;
	// for(auto& boundary : mesh.boundary){
		// // cout << boundary.name << endl;
		// // read.trim;
		// string bcName = boundary.name;
		
		// bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
		// bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
		
		
		// // outputFile << boundary.name << " ";
		// outputFile << bcName << " ";
	// }
	// outputFile << endl;
	// outputFile << " </bcName>" << endl;
	
	// outputFile << " <bcStartFace>" << endl;
	// for(auto& boundary : mesh.boundary){
		// outputFile << boundary.startFace << " ";
	// }
	// outputFile << endl;
	// outputFile << " </bcStartFace>" << endl;
	
	// outputFile << " <bcNFaces>" << endl;
	// for(auto& boundary : mesh.boundary){
		// outputFile << boundary.nFaces << " ";
	// }
	// outputFile << endl;
	// outputFile << " </bcNFaces>" << endl;
	
	// outputFile << " <bcNeighbProcNo>" << endl;
	// for(auto& boundary : mesh.boundary){
		// outputFile << boundary.neighbProcNo << " ";
	// }
	// outputFile << endl;
	// outputFile << " </bcNeighbProcNo>" << endl;
	
	
	// outputFile << "</VTKFile>" << endl;
	
	// outputFile.close();
	
	

	
	// // if(rank==0){
		// // cout << "-> completed" << endl;
		// // cout << "└────────────────────────────────────────────────────" << endl;
	// // }
	
	
}








void loadVTU(
	string folder,
	SEMO_Mesh_Builder &mesh,
	int rankTarget){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load vtu data files ... ";
	// }
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rankTarget);
	
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
		mesh.boundary.back().myProcNo = rankTarget;
		if(bcNeighbProcNo[n] < 0){
			mesh.boundary.back().myProcNo = -1;
		}
		
		++n;
		
	}
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
		
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	

	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load boundary property files ... ";
	// }
	
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
	
	
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
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
	

	// mesh.informations();

	// mesh.saveFile("vtu");
		
	
	
	
}







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




void loadSaveVTU(
	string folder, 
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species,
	int rankTarget){
		
	SEMO_Mesh_Load load;
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load vtu data files ... ";
	// }
		
	string saveFolderName = folder;
	string saveFileName = "plot";
	string saveRankName = to_string(rankTarget);
	
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
	
	
	for(auto& face : mesh.faces){
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		face.var.resize(controls.nTotalFaceVar,0.0);
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
		mesh.boundary.back().myProcNo = rankTarget;
		if(bcNeighbProcNo[n] < 0){
			mesh.boundary.back().myProcNo = -1;
		}
		
		++n;
		
		// cout << mesh.boundary.back().name << endl;
		
	}
	bcName.clear();
	bcStartFace.clear();
	bcNFaces.clear();
	bcNeighbProcNo.clear();
	
		
		
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
	

	// if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute load boundary property files ... ";
	// }
	
	

	// map<string,string> mapBCP;
	// read.file("./constant/pressure", "fluxScheme", mapBCP);
	
	
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
	
	
	// if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	// }
	
	
	
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
	

	// mesh.informations();

	// mesh.saveFile("vtu");
		
	
	
	
}



