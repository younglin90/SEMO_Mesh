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
#include <stdexcept>
#include <iomanip>
#include <mpi.h>

using namespace std;

#include "../mesh/build.h" 
#include "../mesh/load.h" 
#include "../mesh/geometric.h" 

#include "../solvers/build.h" 
#include "../controls/build.h" 
#include "../variables/build.h" 


void saveInitialField(string folder);
	
SEMO_Mesh_Builder mesh;
vector<SEMO_Species> species;
SEMO_Controls_Builder controls;
SEMO_MPI_Builder mpi;
SEMO_Mesh_Load load;
SEMO_Utility_Math math;
SEMO_Mesh_Geometric geometric;
	
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
	mapArgv.find("-h") != mapArgv.end()
	){
		if(rank==0){
			cout << endl;
			cout << endl;
			cout << "┌─────── Potential Solver's helper ─────────────────────────────── " << endl;
			cout << "| -n \"int num.\" : iteration of solver" << endl;
			cout << "| -r \"real num.\"    : under relaxation factor of P, U, V, W" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	
	
	int iPotentialMax = stoi(mapArgv["-n"]);
	double relaxationFactor = stod(mapArgv["-r"]);
	
	
	

		
	controls.readSpecies(species);
	
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	
	SEMO_Solvers_Builder solvers;
	
	
	
	load.vtu("./save/0/", mesh, controls, species);
	
	solvers.calcCellEOSVF(mesh, controls, species);
	

		
	
	geometric.init(mesh);
	
	math.initLeastSquare(mesh);
	
	
	
	// calc velocity potential
	vector<vector<double>> gradP;
	math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	
	
	
	
	for(int iPotential = 0; iPotential<iPotentialMax; ++iPotential){
	
	
		solvers.setIncomValuesLeftRightFace(mesh, controls, species);
		
		
		
		vector<double> linA(mesh.cells.size(),0.0);
		vector<double> linAL(mesh.faces.size(),0.0);
		vector<double> linAR(mesh.faces.size(),0.0);
		vector<double> linB(mesh.cells.size(),0.0);

		
		vector<double> linA_P(mesh.cells.size(),0.0);
		vector<double> linAL_P(mesh.faces.size(),0.0);
		vector<double> linAR_P(mesh.faces.size(),0.0);
		vector<double> linB_P(mesh.cells.size(),0.0);
		
		
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double area = face.area;
			vector<double> nvec(3,0.0);
			nvec[0] = face.unitNormals[0];
			nvec[1] = face.unitNormals[1];
			nvec[2] = face.unitNormals[2];
			
			vector<double> distanceCells(3,0.0);
			distanceCells[0] = face.distCells[0];
			distanceCells[1] = face.distCells[1];
			distanceCells[2] = face.distCells[2];
		
			// if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				// distanceCells[0] *= 2.0;
				// distanceCells[1] *= 2.0;
				// distanceCells[2] *= 2.0;
			// }
			double UL = face.varL[controls.fU];
			double VL = face.varL[controls.fV];
			double WL = face.varL[controls.fW];
			double UR = face.varR[controls.fU];
			double VR = face.varR[controls.fV];
			double WR = face.varR[controls.fW];
			
			double PL = face.varL[controls.fP];
			double PR = face.varR[controls.fP];
			
			double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
			double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
			
			double RhoL = face.varL[controls.fRho];
			double RhoR = face.varR[controls.fRho];
			
		
			double UnF = 0.5*(UnL+UnR);
			
					
			double dPN = sqrt(pow(distanceCells[0],2.0)+pow(distanceCells[1],2.0)+pow(distanceCells[2],2.0));
			double dPNEff = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
			dPNEff = area/dPN *(dPNEff/dPN);
			
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				linAL[i] = -1.0*dPNEff;
				linAR[i] = -1.0*dPNEff;
				
				linA[face.owner] += (-linAL[i]);
				linA[face.neighbour] += (-linAR[i]);
				
				// convective term
				linB[face.owner] -= UnF*area;
				linB[face.neighbour] += UnF*area;
				
				
				
				// // PRESSURE
				// linAL[i] = -1.0*dPNEff;
				// linAR[i] = -1.0*dPNEff;
				
				// linA[face.owner] += (-linAL[i]);
				// linA[face.neighbour] += (-linAR[i]);
				
				// // convective term
				// linB[face.owner] -= UnF*area;
				// linB[face.neighbour] += UnF*area;
				
				
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				linAL[i] = -1.0*dPNEff;
				linAR[i] = -1.0*dPNEff;
				
				linA[face.owner] += (-linAL[i]);
				
				// convective term
				linB[face.owner] -= UnF*area;
				
				
				
				// // PRESSURE
				// linAL[i] = -1.0*dPNEff;
				// linAR[i] = -1.0*dPNEff;
				
				// linA[face.owner] += (-linAL[i]);
				
				// // // convective term
				// // linB[face.owner] -= UnF*area;
				
				
			}
			else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
				
				// convective term
				linB[face.owner] -= UnF*area;
				
				
				
				// // PRESSURE
				// // convective term
				// linB[face.owner] -= UnF*area;
				
			}
			
		}
		
		// boundary
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				if(boundary.type[0] == "fixedValue"){
					for(int i=str; i<end; ++i){
						auto& face = mesh.faces[i];
						
						double area = face.area;
						vector<double> nvec(3,0.0);
						nvec[0] = face.unitNormals[0];
						nvec[1] = face.unitNormals[1];
						nvec[2] = face.unitNormals[2];
						
						vector<double> distanceCells(3,0.0);
						distanceCells[0] = 2.0*face.distCells[0];
						distanceCells[1] = 2.0*face.distCells[1];
						distanceCells[2] = 2.0*face.distCells[2];
					
						// double RhoL = face.varL[controls.fRho];
						// double RhoR = face.varR[controls.fRho];
					
						// double tmp1 = 1.0/RhoL*controls.timeStep;
						// double tmp1 = mesh.cells[face.owner].volume/linAD[face.owner]/RhoL;
					
						// double dPN = 
							// sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
													// std::begin(distanceCells), 0.0));
						// double dPNEff = area/dPN/dPN*
							// (std::inner_product(std::begin(nvec), std::end(nvec), 
												// std::begin(distanceCells), 0.0));
														
						double dPN = sqrt(pow(distanceCells[0],2.0)+pow(distanceCells[1],2.0)+pow(distanceCells[2],2.0));
						double dPNEff = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
						dPNEff = area/dPN *(dPNEff/dPN);
						
						
						linA[face.owner] += 1.0*dPNEff;
					}
				}
				
			}
			
		}
		
		

		vector<double> resiVar(mesh.cells.size(),0.0);
		solvers.solvePETSc(mesh, resiVar, linA, linAL, linAR, linB,
			"tfqmr", 1.e-7, 1.e-7, "mg", 50);
		
		
		
		
		// // velocities correction
		// vector<vector<double>> gradResiP(mesh.cells.size(),vector<double>(3,0.0));
		
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// double PF = 0.5*(resiVar[face.owner]+resiVar[face.neighbour]);
				
				// for(int j=0; j<3; ++j){
					// gradResiP[face.owner][j] += 
						// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
					// gradResiP[face.neighbour][j] -= 
						// PF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
				// }
			// }
		// }

		
		// // boundary
		// for(auto& boundary : mesh.boundary){
			
			// if(boundary.neighbProcNo == -1){
				
				// int str = boundary.startFace;
				// int end = str + boundary.nFaces;
				
				// if(boundary.type[0] == "fixedValue"){
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						
						// double PF = 0.5*(resiVar[face.owner]);
						// for(int j=0; j<3; ++j){
							// gradResiP[face.owner][j] += 
								// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
						// }
						
					// }
				// }
				// else{
					// for(int i=str; i<end; ++i){
						// auto& face = mesh.faces[i];
						
						// double PF = resiVar[face.owner];
						// for(int j=0; j<3; ++j){
							// gradResiP[face.owner][j] += 
								// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
						// }
						
					// }
				// }
			// }
		// }
		
		
		// // processor faces
		// vector<double> sendValues;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// sendValues.push_back(resiVar[face.owner]);
			// }
		// }
		
		
		// vector<double> recvValues;
		// mpi.setProcsFaceDatasDouble(
					// sendValues, recvValues,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);

		// int num=0;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// double PF = 0.5*(resiVar[face.owner]+recvValues[num]);
				// for(int j=0; j<3; ++j){
					// gradResiP[face.owner][j] += 
						// PF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				// }
				// ++num;
			// }
		// }
		

		// calc gradient
		
		vector<vector<double>> gradResiP;
		
		gradResiP.clear();
		gradResiP.resize(mesh.cells.size(),vector<double>(3,0.0));
		
		
		for(auto& face : mesh.faces){
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				double DVar = 
					resiVar[face.neighbour] - resiVar[face.owner];
				
				gradResiP[face.owner][0] += face.distCells[0] * DVar;
				gradResiP[face.owner][1] += face.distCells[1] * DVar;
				gradResiP[face.owner][2] += face.distCells[2] * DVar;
				
				gradResiP[face.neighbour][0] += face.distCells[0] * DVar;
				gradResiP[face.neighbour][1] += face.distCells[1] * DVar;
				gradResiP[face.neighbour][2] += face.distCells[2] * DVar;
				
			}
			
		}
		
		
		// boundary
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				if(boundary.type[0] == "fixedValue"){
					for(int i=str; i<end; ++i){
						auto& face = mesh.faces[i];
						
						double DVar = 
							0.0 - resiVar[face.owner];
						
						gradResiP[face.owner][0] += face.distCells[0] * DVar;
						gradResiP[face.owner][1] += face.distCells[1] * DVar;
						gradResiP[face.owner][2] += face.distCells[2] * DVar;
						
						
					}
				}
			}
		}
		
		
		// processor faces
		vector<double> sendValues;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues.push_back(resiVar[face.owner]);
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		vector<double> recvValues;
		mpi.setProcsFaceDatas(
					sendValues, recvValues,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);

		int num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				double DVar = 
					recvValues[num] - resiVar[face.owner];
				
				gradResiP[face.owner][0] += face.distCells[0] * DVar;
				gradResiP[face.owner][1] += face.distCells[1] * DVar;
				gradResiP[face.owner][2] += face.distCells[2] * DVar;
						
				++num;
			}
		}
		
		
		
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			double tmp0 = 
				cell.coeffLeastSquare[0] * gradResiP[i][0] +
				cell.coeffLeastSquare[1] * gradResiP[i][1] +
				cell.coeffLeastSquare[2] * gradResiP[i][2];
				
			double tmp1 = 
				cell.coeffLeastSquare[1] * gradResiP[i][0] +
				cell.coeffLeastSquare[3] * gradResiP[i][1] +
				cell.coeffLeastSquare[4] * gradResiP[i][2];
				
			double tmp2 = 
				cell.coeffLeastSquare[2] * gradResiP[i][0] +
				cell.coeffLeastSquare[4] * gradResiP[i][1] +
				cell.coeffLeastSquare[5] * gradResiP[i][2];
				
			gradResiP[i][0] = tmp0;
			gradResiP[i][1] = tmp1;
			gradResiP[i][2] = tmp2;
			
		}
		
		
		double residual = 0.0;

		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			// cell.var[controls.P] += 0.1 * cell.var[controls.Rho] * resiVar[i];
			cell.var[controls.P] += relaxationFactor * resiVar[i];
			
			cell.var[controls.U] -= relaxationFactor * gradResiP[i][0];
			cell.var[controls.V] -= relaxationFactor * gradResiP[i][1];
			cell.var[controls.W] -= relaxationFactor * gradResiP[i][2];
			
			residual += abs(resiVar[i]);
		}
		

		double residualReduced;
		MPI_Allreduce(&residual, &residualReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		if(rank==0) cout << "| iteration = " << iPotential << " | residual = " << residualReduced << endl;
	}
	
	saveInitialField("./save/0_Potential/");
	

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
	
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.VF[0]] << " ";
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
		string filenamePvtu = folder + "../plot.0_Potential.pvtu";
		
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
			string filenamevtus = "./0_Potential/plot.";
			filenamevtus += to_string(ip);
			filenamevtus += ".vtu";
			outputFile << "    <Piece Source=\"" << filenamevtus << "\"/>" << endl;
		}
		outputFile << "   <PPointData>" << endl;
		outputFile << "   </PPointData>" << endl;
		outputFile << "   <PCellData>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"pressure\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"temperature\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"volumeFraction-" << species[0].name << "\"/>" << endl;
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




