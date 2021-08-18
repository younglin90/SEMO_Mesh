#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// #include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>
using namespace std;

#include <sys/stat.h>
// #include <filesystem> // C++17 standard header file name
// #include <experimental/filesystem> // Header file for pre-standard implementation
// namespace fs = std::experimental::filesystem;

#include <zlib.h>

#include "save.h" 
#include "build.h" 
#include "../controls/build.h" 

void mkdirs(char *dir_path){
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

/** Compress a STL string using zlib with given compression level and return
  * the binary data. */
std::string compress_string(const std::string& str,
                            int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();           // set the z_stream's input

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // retrieve the compressed bytes blockwise
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = deflate(&zs, Z_FINISH);

        if (outstring.size() < zs.total_out) {
            // append the block to the output string
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }
    } while (ret == Z_OK);

    deflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib compression: (" << ret << ") " << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}

/** Decompress an STL string using zlib and return the original data. */
std::string decompress_string(const std::string& str)
{
    z_stream zs;                        // z_stream is zlib's control structure
    memset(&zs, 0, sizeof(zs));

    if (inflateInit(&zs) != Z_OK)
        throw(std::runtime_error("inflateInit failed while decompressing."));

    zs.next_in = (Bytef*)str.data();
    zs.avail_in = str.size();

    int ret;
    char outbuffer[32768];
    std::string outstring;

    // get the decompressed bytes blockwise using repeated calls to inflate
    do {
        zs.next_out = reinterpret_cast<Bytef*>(outbuffer);
        zs.avail_out = sizeof(outbuffer);

        ret = inflate(&zs, 0);

        if (outstring.size() < zs.total_out) {
            outstring.append(outbuffer,
                             zs.total_out - outstring.size());
        }

    } while (ret == Z_OK);

    inflateEnd(&zs);

    if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
        std::ostringstream oss;
        oss << "Exception during zlib decompression: (" << ret << ") "
            << zs.msg;
        throw(std::runtime_error(oss.str()));
    }

    return outstring;
}


void SEMO_Mesh_Save::vtu(SEMO_Mesh_Builder &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (.vtu format) ... ";
	}
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
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

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
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}




void SEMO_Mesh_Save::vtu(string folder, int myRank, SEMO_Mesh_Builder &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	// fs::path Paths(folder);
	// if( ! exists(Paths) ){
		// create_directories(Paths);
	// }
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(myRank) + ".vtu";
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << filenamePlot << ") ... ";
	}
	
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
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;  
	outputFile.precision(20);  
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
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}







void SEMO_Mesh_Save::vtu(
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
	
	
	
	// Points data
	
	// SEMO_Solvers_Builder solvers;
	// solvers.setCompValuesLeftRightFace(mesh, controls, species);
	
	outputFile << "    <PointData>" << endl;
	// for(auto& point : mesh.points){
		// point.var.resize(7,0.0);
	// }
	// for(auto& face : mesh.faces){
		// if(
		// face.getType() == SEMO_Types::INTERNAL_FACE ||
		// face.getType() == SEMO_Types::PROCESSOR_FACE
		// ){
			// for(auto& point : face.points){
				
				// mesh.points[point].var[0] += face.varL[controls.fP];
				// mesh.points[point].var[0] += face.varR[controls.fP];
				
				// mesh.points[point].var[1] += face.varL[controls.fU];
				// mesh.points[point].var[1] += face.varR[controls.fU];
				
				// mesh.points[point].var[2] += face.varL[controls.fV];
				// mesh.points[point].var[2] += face.varR[controls.fV];
				
				// mesh.points[point].var[3] += face.varL[controls.fW];
				// mesh.points[point].var[3] += face.varR[controls.fW];
				
				// mesh.points[point].var[4] += face.varL[controls.fT];
				// mesh.points[point].var[4] += face.varR[controls.fT];
				
				// mesh.points[point].var[5] += face.varL[controls.fRho];
				// mesh.points[point].var[5] += face.varR[controls.fRho];
				
				// // mesh.points[point].var[5] += face.varL[controls.fVF[0]];
				// // mesh.points[point].var[5] += face.varR[controls.fVF[0]];
				
				// mesh.points[point].var[6] += 2.0;
			// }
		// }
	// }
	// for(auto& point : mesh.points){
		// point.var[0] /= point.var[6];
		// point.var[1] /= point.var[6];
		// point.var[2] /= point.var[6];
		// point.var[3] /= point.var[6];
		// point.var[4] /= point.var[6];
		// point.var[5] /= point.var[6];
	// }
	// for(auto& face : mesh.faces){
		// if(
		// face.getType() == SEMO_Types::BOUNDARY_FACE
		// ){
			// for(auto& point : face.points){
				
				// mesh.points[point].var[0] = face.varL[controls.fP];
				// mesh.points[point].var[1] = face.varL[controls.fU];
				// mesh.points[point].var[2] = face.varL[controls.fV];
				// mesh.points[point].var[3] = face.varL[controls.fW];
				// mesh.points[point].var[4] = face.varL[controls.fT];
				// mesh.points[point].var[5] = face.varL[controls.fRho];
				// // mesh.points[point].var[5] = face.varL[controls.fVF[0]];
			// }
		// }
	// }
	// outputFile << "     <DataArray type=\"Float64\" Name=\"pointP\" format=\"ascii\">" << endl;
	// for(auto& point : mesh.points) outputFile << scientific << point.var[0] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"pointVel\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	// for(auto& point : mesh.points) outputFile << scientific << point.var[1] << " "<< point.var[2] << " "<< point.var[3] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;

	// outputFile << "     <DataArray type=\"Float64\" Name=\"pointT\" format=\"ascii\">" << endl;
	// for(auto& point : mesh.points) outputFile << scientific << point.var[4] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"pointRho\" format=\"ascii\">" << endl;
	// for(auto& point : mesh.points) outputFile << scientific << point.var[5] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// for(auto& point : mesh.points){
		// point.var.clear();
	// }
	outputFile << "    </PointData>" << endl;
	
	
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
	
	
	// cell levels
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellLevels\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.level << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	// cell groups
	outputFile << "     <DataArray type=\"Int64\" Name=\"cellGroups\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << cell.group << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	// // cell volume
	// outputFile << "     <DataArray type=\"Float64\" Name=\"volume\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.volume << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// // UDV
	// outputFile << "     <DataArray type=\"Float64\" Name=\"UDV0\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[0]] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"UDV1\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[1]] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	// outputFile << "     <DataArray type=\"Float64\" Name=\"UDV2\" format=\"ascii\">" << endl;
	// for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.UDV[2]] << " ";
	// outputFile << endl;
	// outputFile << "     </DataArray>" << endl;
	
	
	
	
	
	
	
	
	
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
	
	
	outputFile << " <faceLevels>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.level << " ";
	}
	outputFile << endl;
	outputFile << " </faceLevels>" << endl;
	
	
	outputFile << " <faceGroups>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.group << " ";
	}
	outputFile << endl;
	outputFile << " </faceGroups>" << endl;
	
	
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
	
	
	
	
	

	if(rank==0){
		// string filenamePvtu = folder + "../plot.";
		string filenamePvtu = "./save/plot.";
		string stime = folder;
		stime.erase(stime.find("./save/"),7);
		stime.erase(stime.find("/"),1);
		// filenamePvtu += to_string(controls.time);
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

		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pointP\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pointVel\" NumberOfComponents=\"3\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pointT\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"pointRho\"/>" << endl;
		
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"density\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"" << volFracName[0] << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellLevels\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellGroups\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"volume\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"UDV0\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"UDV1\"/>" << endl;
		// outputFile << "    <PDataArray type=\"Float64\" Name=\"UDV2\"/>" << endl;
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
	
	
	// }
	
	
	
	
	
}





void SEMO_Mesh_Save::vtuZlib(SEMO_Mesh_Builder &mesh, SEMO_Controls_Builder &controls){
	
	int rank = MPI::COMM_WORLD.Get_rank();



	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (.vtu format) ... ";
	}
	ofstream outputFile;
	string filenamePlot = "./save/plot." + to_string(rank) + ".vtu";
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	// outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\" compressor=\"vtkZLibDataCompressor\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	for(auto& point : mesh.points){
		point.var.resize(2,0.0);
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(auto& point : face.points){
				mesh.points[point].var[0] += mesh.cells[face.owner].var[0];
				mesh.points[point].var[0] += mesh.cells[face.neighbour].var[0];
				mesh.points[point].var[1] += 2.0;
			}
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			for(auto& point : face.points){
				mesh.points[point].var[0] += 101325.0;
				mesh.points[point].var[1] += 1.0;
			}
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			for(auto& point : face.points){
				mesh.points[point].var[0] += face.varR[0];
				mesh.points[point].var[1] += 1.0;
			}
		}
	}
	outputFile << "     <DataArray type=\"Float64\" Name=\"pressure\" format=\"appended\" offset=\"0\">" << endl;
	
    // // original string len = 36
    // char a[50] = "Hello Hello Hello Hello Hello Hello!"; 

    // char b[50];
    // z_stream defstream;
    // defstream.zalloc = Z_NULL;
    // defstream.zfree = Z_NULL;
    // defstream.opaque = Z_NULL;
    // // setup "a" as the input and "b" as the compressed output
    // defstream.avail_in = (uInt)strlen(a)+1; // size of input, string + terminator
    // defstream.next_in = (Bytef *)a; // input char array
    // defstream.avail_out = (uInt)sizeof(b); // size of output
    // defstream.next_out = (Bytef *)b; // output char array
    
    // // the actual compression work.
    // deflateInit(&defstream, Z_DEFAULT_COMPRESSION);
    // deflate(&defstream, Z_FINISH);
    // deflateEnd(&defstream);
	// string allinput;
	for(auto& point : mesh.points) {
		// allinput.append(inbuffer, std::cin.gcount());
		// string allinput = to_string(point.var[0]/point.var[1]) + " ";
		// std::string cstr = compress_string( allinput );
		// outputFile << cstr;
		// outputFile << scientific << point.var[0]/point.var[1] << " ";
	}
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	// for(auto& point : mesh.points){
		// point.var.clear();
	// }
	outputFile << "    </PointData>" << endl;
	
	
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
	
	

	outputFile << " <AppendedData encoding=\"raw\">" << endl;
	string allinput = "_";
	for(auto& point : mesh.points) {
		allinput.append(to_string(point.var[0]/point.var[1]));
		allinput.append(" ");
		// string allinput = "_" + to_string(point.var[0]/point.var[1]);
		// string allinput = to_string(point.var[0]/point.var[1]) + " ";
		// std::string cstr = compress_string( allinput );
		// outputFile << cstr << endl;
		// outputFile << scientific << point.var[0]/point.var[1] << " ";
		// std::string cstr = compress_string( allinput );
		// outputFile << cstr << endl;
	}
	cout << allinput << endl;
	std::string cstr = compress_string( allinput );
	outputFile << cstr << endl;
	outputFile << "  </AppendedData>" << endl;
	
	

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
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
}













