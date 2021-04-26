#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// #include <algorithm>
#include <mpi.h>
using namespace std;

#include "save.h" 
#include "build.h" 

void SEMO_Mesh_Save::vtu(SEMO_Mesh_Builder &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	if(rank==0) cout << "| execute save file (.vtu format) ... ";
	ofstream outputFile;
	string filenamePlot = "./save/plot." + to_string(rank) + ".vtu";
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.listPoints.size() << "\" NumberOfCells=\"" << mesh.listCells.size() << "\">" << endl;
	
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
	for(auto iter=mesh.listPoints.begin(); iter!=mesh.listPoints.end(); iter++){
		outputFile << scientific << (*iter)->x << " " << (*iter)->y << " " << (*iter)->z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
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
	for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		cellFaceOffset += (*iter)->points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
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
	for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
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
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "------------------------------------" << endl;
	
	
	// }
	
	
	
}




