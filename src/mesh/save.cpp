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
		outputFile << boundary.name << " ";
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
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "------------------------------------" << endl;
	
	
	// }
	
	
	
}




