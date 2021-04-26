#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>

#include "parmetis.h" 
#include "scotch.h" 

#include "./mesh/build.h" 
#include "./mesh/load.h" 
#include "./mesh/geometric.h" 
// #include "./mpi/build.h"

int main(int argc, char* argv[]) {

	// printf("%d\n",argc);
	// printf("%s\n",argv[1]);
	
	// MPI initializing
	// int rank, size;
	// MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    // MPI_Comm_size(MPI_COMM_WORLD, &size); 
	
	// MPI::Status status;
	// MPI::Init(); 
    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();
	
	// SEMO_MPI_Builder mpi;
	// mpi.init();
	// mpi.setRank();
	// mpi.setSize();
	
	SEMO_Mesh_Builder mesh;
	
	mesh.mpi.init();
	mesh.mpi.setSize();
	mesh.mpi.setRank();
	
	
	if(mesh.mpi.getRank() == 0){
		
		mesh.loadFile("OpenFoam");
		
		mesh.check();
		
		mesh.buildCells();
		
		mesh.setFaceTypes();
		

		// create list
		mesh.buildLists();
		
		mesh.checkLists();
		
		// cell's points, faces connection
		mesh.connectCelltoFaces();
		
		mesh.connectCelltoPoints();
		
		
		// for(auto& i : mesh.faces){
			// cout << i.owner << endl;
		// }
		
		
	}
	
	bool boolPartitioning = true;
	bool boolPlotting = false;
	bool boolGeometric = false;
	
	// partitioning
	if(boolPartitioning){
		
		mesh.distributeOneToAll("EVENLY");
		
		// mesh.partitionInit("METIS");
		
		// mesh.partition("METIS");
		
	}
	
	
	
	
	// // geometric
	// if(boolGeometric){
	// cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	// SEMO_Mesh_Geometric geometric;

	
	// // polygon face normal vectors & polygon face area
	// // polygon face center x,y,z
	// // 3D Polygon Areas : https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
	// for(auto iter=mesh.faces.begin(); iter!=mesh.faces.end(); iter++){
		// vector<double> Vx, Vy, Vz;
		// for(auto iPoint : iter->points){
			// Vx.push_back(mesh.points[iPoint].x);
			// Vy.push_back(mesh.points[iPoint].y);
			// Vz.push_back(mesh.points[iPoint].z);
		// }
		
		// geometric.calcUnitNormals_Area3dPolygon(
		// iter->points.size(), Vx,Vy,Vz,
		// iter->unitNormals, iter->area );
		
		// iter->x = accumulate(Vx.begin(), Vx.end(), 0.0) / iter->points.size();
		// iter->y = accumulate(Vy.begin(), Vy.end(), 0.0) / iter->points.size();
		// iter->z = accumulate(Vz.begin(), Vz.end(), 0.0) / iter->points.size();
	// }
	
	
	// // polyhedron cell volume (Green-Gauss Theorem.)
	// for(auto& cell : mesh.cells)
		// cell.volume = 0.0;
	
	// for(auto& face : mesh.faces){
		
		// // cout << iter->owner << endl;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// mesh.cells[face.owner].volume += 
				// (face.x*face.unitNormals[0]+
				 // face.y*face.unitNormals[1]+
				 // face.z*face.unitNormals[2])*face.area;
			
			// mesh.cells[face.neighbour].volume -= 
				// (face.x*face.unitNormals[0]+
				 // face.y*face.unitNormals[1]+
				 // face.z*face.unitNormals[2])*face.area;
		// }
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// mesh.cells[face.owner].volume += 
				// (face.x*face.unitNormals[0]+
				 // face.y*face.unitNormals[1]+
				 // face.z*face.unitNormals[2])*face.area;
		// }
		
	// }
	// for(auto iter=mesh.cells.begin(); iter!=mesh.cells.end(); iter++){
		// iter->volume /= 3.0;
		
		// if(iter->volume < std::numeric_limits<double>::min()) {
			// cerr << endl;
			// cerr << "#error, from calc cell volume, cell volume = " << iter->volume << " < cpu_min_val " << endl;
			// cerr << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		// // cout << iter->volume << endl;
	// }
	// cout << "-> completed" << endl;
	
	// }
	
	
	
	// // plotting
	
	
	// if(boolPlotting){
	
	// cout << "| execute save file (.vtu format) ... ";
	// ofstream outputFile;
	// outputFile.open("./save/plot.vtu");
	// // outputFile.open("./save/plot.vtu",ios::binary);
	// // outputFile.open("./save/plot.vtu",ios::out | ios::binary);
	// if(outputFile.fail()){
		// cerr << "Unable to write file for writing." << endl;
		// return 1;
	// }
	
	// // string out_line;
	// outputFile << "<?xml version=\"1.0\"?>" << endl;
	// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	// outputFile << "  <UnstructuredGrid>" << endl;
	// outputFile << "   <Piece NumberOfPoints=\"" << mesh.listPoints.size() << "\" NumberOfCells=\"" << mesh.listCells.size() << "\">" << endl;
	
	// // Points data
	// outputFile << "    <PointData>" << endl;
	// outputFile << "    </PointData>" << endl;
	// // Cells data
	// outputFile << "    <CellData>" << endl;
	// outputFile << "    </CellData>" << endl;
	// // Points
	// outputFile << "    <Points>" << endl;
	// // }
	// outputFile << "     <DataArray type=\"Float32\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	// stringstream streamXYZ;
	// // for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	// for(auto iter=mesh.listPoints.begin(); iter!=mesh.listPoints.end(); iter++){
		// outputFile << scientific << (*iter)->x << " " << (*iter)->y << " " << (*iter)->z << endl;

	// }
	
	// outputFile << "    </DataArray>" << endl;
	// outputFile << "   </Points>" << endl;
	
	// // cells
	// outputFile << "   <Cells>" << endl; 
	// // connectivity (cell's points)
	// outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	// for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		// for(auto i : (*iter)->points){
			// outputFile << i << " ";
		// }
		// outputFile << endl;
	// }
	
	// outputFile << "    </DataArray>" << endl;
	
	// // offsets (cell's points offset)
	// int cellFaceOffset = 0;
	// outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	// cellFaceOffset = 0;
	// for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		// cellFaceOffset += (*iter)->points.size();
		// outputFile << cellFaceOffset << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // types (cell's type, 42 = polyhedron)
	// outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	// for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		// outputFile << "42" << " ";
	// }
	// outputFile << endl;
	
	// outputFile << "    </DataArray>" << endl;
	
	// // faces (cell's faces number, each face's point number, cell's faces's points)
	// outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// // outputFile << mesh.faces.size() << endl;
	// for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		// outputFile << (*iter)->faces.size() << endl;
		// for(auto& i : (*iter)->faces){
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
	// for(auto iter=mesh.listCells.begin(); iter!=mesh.listCells.end(); iter++){
		// int numbering = 1 + (*iter)->faces.size();
		// for(auto& i : (*iter)->faces){
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
	// outputFile << "</VTKFile>" << endl;
	
	// outputFile.close();
	// cout << "-> completed" << endl;
	// cout << "------------------------------------" << endl;
	

	// }
	

	return EXIT_SUCCESS;
}














