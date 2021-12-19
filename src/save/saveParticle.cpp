#include <fstream>
#include <cstring>
#include <mpi.h>
using namespace std;

#include "save.h" 
#include "../mesh/build.h" 
#include "../controls/build.h" 

void SEMO_Mesh_Save::particles(
	string folder, 
	SEMO_Mesh_Builder &mesh, 
	SEMO_Controls_Builder &controls,
	vector<SEMO_Species>& species){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);

	ofstream outputFile;
	string filenamePlot = folder + "particles_cell." + to_string(rank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save particles file (" << folder_name << "plot...) ... ";
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
	
	outputFile << "    <FieldData>" << endl;
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.cells.size() << "\" NumberOfCells=\"" << 0 << "\">" << endl;
	
	
	outputFile << "    <PointData>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"Diameter\" format=\"ascii\">" << endl;
	
	double diameter = 1.0;
	for(auto& cell : mesh.cells) {
		outputFile << scientific << diameter << " ";
	}
	outputFile << endl;
	
	outputFile << "     </DataArray>" << endl;
	outputFile << "    </PointData>" << endl;
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	
	
	// Points
	outputFile << "    <Points>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.x << " " << cell.y << " " << cell.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	outputFile << "   <Cells>" << endl; 
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
	
	
	
	
	
	
	// faces
	string filenamePlot2 = folder + "particles_face." + to_string(rank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save particles file (" << folder_name << "plot...) ... ";
	}
	
	outputFile.open(filenamePlot2);
	
	outputFile.precision( controls.writePrecision );
	
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	
	outputFile << "    <FieldData>" << endl;
	outputFile << "    </FieldData>" << endl;
	
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.faces.size() << "\" NumberOfCells=\"" << 0 << "\">" << endl;
	
	
	outputFile << "    <PointData>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"Diameter\" format=\"ascii\">" << endl;
	
	// double diameter2 = 1.0;
	for(auto& face : mesh.faces) {
		outputFile << scientific << diameter << " ";
	}
	outputFile << endl;
	
	outputFile << "     </DataArray>" << endl;
	outputFile << "    </PointData>" << endl;
	
	
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	
	
	// Points
	outputFile << "    <Points>" << endl;
	outputFile << "     <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	for(auto& face : mesh.faces){
		outputFile << scientific << face.x << " " << face.y << " " << face.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	outputFile << "   <Cells>" << endl; 
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	// }
	
	
	
	
	
}


