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
#include "../load/load.h" 
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
	
	
	
	int iLaplaceMax = stoi(mapArgv["-n"]);
	double relaxationFactor = stod(mapArgv["-r"]);
	
	
	

		
	controls.readSpecies(species);
	
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	
	SEMO_Solvers_Builder solvers;
	
	
	
	load.vtu("./save/0/", mesh, controls, species);
	
	solvers.calcCellEOSVF(mesh, controls, species);
	

	
	
	geometric.init(mesh);
	
	math.initLeastSquare(mesh);
	
	
	// analytic solution
	double pi = 3.141592;
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		double x = cell.x;
		double y = cell.y;
		
		double value = 0.0;
		for(int j=1; j<200; ++j){
			double jj = (double)j;
			double temp_value = 0.0;
			temp_value = -( 4.0 * pow(-1.0,jj) * cos(0.5*pi*x*(-1.0+2.0*jj)) * cosh(0.5*pi*y*(-1.0+2.0*jj)) /
				cosh(pi*(-0.5+jj)) )/(-pi+2.0*pi*jj);
			if(abs(temp_value)<1.e-10) temp_value = 0.0;
			value += temp_value;
			if(abs(value)<1.e-10){
				value = 0.0;
				break;
			}
		}
		
		
		cell.var[controls.T] = value;
	}
	
	
	std::clock_t c_start = std::clock();


	for(int iLaplace = 0; iLaplace<iLaplaceMax; ++iLaplace){


		solvers.setIncomValuesLeftRightFace(mesh, controls, species);
	
		
		// diagonal terms
		int B_n = 1;
		int A_n = B_n * B_n;
		
		
		int proc_num=0;
		
		// gradient P
		vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
		vector<double> gradPx_recv;
		vector<double> gradPy_recv;
		vector<double> gradPz_recv;
		
		// linear solver
		vector<double> resiVar(B_n*mesh.cells.size(),0.0);
		
	
		// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
		{
			vector<double> dummyVec;
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.P, controls.fP, dummyVec, gradP);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.P, controls.fP, dummyVec, gradP);
		}
		
		
		if(size>1){
			// processor faces
			// gradP , 
			vector<double> gradPx_send;
			vector<double> gradPy_send;
			vector<double> gradPz_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					gradPx_send.push_back(gradP[face.owner][0]);
					gradPy_send.push_back(gradP[face.owner][1]);
					gradPz_send.push_back(gradP[face.owner][2]);
				}
			}
			// SEMO_MPI_Builder mpi;
			
			mpi.setProcsFaceDatas(
						gradPx_send, gradPx_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradPy_send, gradPy_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradPz_send, gradPz_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			gradPx_send.clear();
			gradPy_send.clear();
			gradPz_send.clear();
		}
		
		
		vector<int> A_rows(mesh.cells.size()*A_n, 0);
		vector<int> A_cols(mesh.cells.size()*A_n, 0);
		vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
		vector<double> B(mesh.cells.size()*B_n, 0.0);
		

		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];

			int id = 0;
			int i_glo = 0;
			int j_glo = 0;
			int i_loc = i;
			int str_glo = mesh.startCellGlobal*B_n;
			int step_loc = mesh.cells.size();
			int step_glo = mesh.ncellsTotal;
			
			id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
			A_rows[id] = i_glo; A_cols[id] = j_glo;
			
		}
		
		
		proc_num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
			
			int ijStartL_local = B_n*(face.owner) - 1;
			int ijStartR_local = B_n*(face.neighbour) - 1;
			
			int ijStartL = B_n*mesh.startCellGlobal + ijStartL_local;
			int ijStartR = B_n*mesh.startCellGlobal + ijStartR_local;
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				ijStartR = B_n*mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] + 
						   B_n*(mesh.procNeighbCellNo[proc_num]) - 1;
			}

			
			double area = face.area;
			vector<double> nvec(3,0.0);
			nvec[0] = face.unitNormals[0];
			nvec[1] = face.unitNormals[1];
			nvec[2] = face.unitNormals[2];
			
			// vector<double> distanceCells(3,0.0);
			// distanceCells[0] = face.distCells[0];
			// distanceCells[1] = face.distCells[1];
			// distanceCells[2] = face.distCells[2];
			
			double wCL = face.wC;
			// double wCL = face.wVC;
			// double wCL = 0.5;
			double wCR = 1.0-wCL;
			
			double PL = face.varL[controls.fP];
			double PR = face.varR[controls.fP];
				
			double dPN = face.magPN;
			double alpha = face.alphaF;
		
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = 0.0;
			vector<double> gradPL(3,0.0);
			gradPL[0] = gradP[face.owner][0];
			gradPL[1] = gradP[face.owner][1];
			gradPL[2] = gradP[face.owner][2];
			
			vector<double> gradPR(3,0.0);
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				orgPR = mesh.cells[face.neighbour].var[controls.P];
				gradPR[0] = gradP[face.neighbour][0];
				gradPR[1] = gradP[face.neighbour][1];
				gradPR[2] = gradP[face.neighbour][2];
				
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				orgPR = PR;
				gradPR[0] = gradPx_recv[proc_num];
				gradPR[1] = gradPy_recv[proc_num];
				gradPR[2] = gradPz_recv[proc_num];
			}
			
			
		
			double nonOrtholimiter = 1.0;
			// double cosAlpha = 1.0/alpha;
			// if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
				// nonOrtholimiter = 0.5;
			// }
			// else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
				// nonOrtholimiter = 0.333;
			// }
			// else if( cosAlpha < 0.342 ){
				// nonOrtholimiter = 0.0;
			// }
			
			
			double magGradPF = 0.0;
			magGradPF += alpha * (orgPR-orgPL)/dPN * area;
			// non-orthogonal
			for(int ii=0; ii<3; ++ii){
				magGradPF += nonOrtholimiter * 
					(wCL*gradPL[ii] + wCR*gradPR[ii])*
					(nvec[ii] - alpha*face.unitNomalsPN[ii]) * area;			
			}
			
			// // pressure correction (skewness)
			// for(int ii=0; ii<3; ++ii){
				// magGradPF += alpha * (gradPR[ii]-gradPL[ii])*face.vecSkewness[ii]/dPN * area;
			// }
			
			double diff_magGradPF = alpha / dPN * area;
			// double diff_magGradPF = alpha_grad * alpha / dPN * area;




			int str_glo_L = mesh.startCellGlobal*B_n;
			int str_glo_R = mesh.startCellGlobal*B_n;
			int step_loc_L = mesh.cells.size();
			int step_loc_R = mesh.cells.size();
			int i_loc_L = face.owner;
			int i_loc_R = face.neighbour;
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				str_glo_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]]*B_n;
				step_loc_R = mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]+1] - 
							 mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]];
				i_loc_R = mesh.procNeighbCellNo[proc_num];
			}
			
			int id = -1;
			int i_glo = -1;
			int j_glo = -1;
			double conv_flux = 0.0;
			double diff_flux = 0.0;
			
			
			//------------------------
			id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
			A_vals[id] += ( - diff_magGradPF );
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			A_vals.push_back(( + diff_magGradPF ));

			if(face.getType() == SEMO_Types::INTERNAL_FACE) {
				id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
				A_vals[id] -= ( + diff_magGradPF );
				A_rows.push_back(i_glo); A_cols.push_back(j_glo);
				A_vals.push_back(-( - diff_magGradPF ));
			}
			
	 
			// ----------------------------
			B[step_loc_L*0 + i_loc_L] -= magGradPF;
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				B[step_loc_R*0 + i_loc_R] += magGradPF;
			}
			
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
			
			
		}
		
		
		
		// boundary
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec(3,0.0);
					nvec[0] = face.unitNormals[0];
					nvec[1] = face.unitNormals[1];
					nvec[2] = face.unitNormals[2];
			
					double PF = 0.5*face.varL[controls.fP] + 0.5*face.varR[controls.fP];
				
					// pressure coefficient
					double coeffP = 0.0;
					double coeffP_diff = 0.0;
					if( boundary.type[controls.P] == "fixedValue" ){
						coeffP = 0.0;
						coeffP_diff = 1.0;
					}
					else if( boundary.type[controls.P] == "zeroGradient" ){
						coeffP = 1.0;
						coeffP_diff = 0.0;
						
					}
					else if( boundary.type[controls.P] == "switch" ){
						double machNum = 
							sqrt(pow(mesh.cells[face.owner].var[controls.U],2.0)+
								 pow(mesh.cells[face.owner].var[controls.V],2.0)+
								 pow(mesh.cells[face.owner].var[controls.W],2.0))/
								  mesh.cells[face.owner].var[controls.C];
						coeffP = (machNum > 1.0) ? 0.0 : 1.0;
						coeffP_diff = (machNum > 1.0) ? 1.0 : 0.0;
					}
					else {
						cerr << "| #Error : not defined B.C., var = " << controls.P << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					
					
					double dPN = 0.5 * face.magPN;
					double alpha = face.alphaF;
					
					double orgPL = mesh.cells[face.owner].var[controls.P];
					double orgPR = face.varR[controls.fP];
					

					int str_glo = mesh.startCellGlobal*B_n;
					int step_loc = mesh.cells.size();
					int step_glo = mesh.ncellsTotal;
					int i_loc_L = face.owner;
					int id=-1;
					
					
					double nonOrtholimiter = 1.0;

					double alpha_grad = 4.0/3.0;
					double magGradPF = 0.0;
					magGradPF += coeffP_diff * alpha * (orgPR-orgPL)/dPN * area;
					// magGradPF += coeffP_diff * alpha_grad * alpha * (orgPR-orgPL)/dPN * area;
					// non-orthogonal
					for(int ii=0; ii<3; ++ii){
						magGradPF += coeffP_diff * nonOrtholimiter * 
							gradP[face.owner][ii]*(nvec[ii] - alpha*face.unitNomalsPN[ii]) * area; 
						// magGradPF += coeffP_diff * gradP[face.owner][ii]*nvec[ii] * area; 
					}
			
						
					double diff_magGradPF = alpha / dPN * area;
					// double diff_magGradPF = alpha_grad * alpha / dPN * area;
			
			
			
					
					//
					id = step_loc*(B_n*0+0) + i_loc_L;
					A_vals[id] -= coeffP_diff * diff_magGradPF;
					
					// 
					B[step_loc*0 + i_loc_L] -= ( magGradPF );
					
				}
				
			}
			
		}
		
		
		
		solvers.solveAMGCL("pressure", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);


		
		// update
		double residual = 0.0;
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			cell.var[controls.P] += relaxationFactor * resiVar[i];
			
			cell.var[controls.UDV[0]] = gradP[i][0];
			cell.var[controls.UDV[1]] = gradP[i][1];
			cell.var[controls.UDV[2]] = gradP[i][2];
			
			residual += abs(resiVar[i]);
		}
	


		double residualReduced;
		MPI_Allreduce(&residual, &residualReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		if(rank==0) cout << "| iteration = " << iLaplace << " | residual = " << residualReduced << endl;
	}
	
	
	std::clock_t c_end = std::clock();

	long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << endl;
	std::cout << "| CPU time used: " << time_elapsed_ms << " ms\n";
	cout << endl;
	
	saveInitialField("./save/0_Laplace/");
	

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
	
	outputFile.precision(20);
	
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
	outputFile << "     <DataArray type=\"Float64\" Name=\"totalError\" NumberOfTuples=\"1\" format=\"ascii\">" << endl;
	double totalError = 0.0;
	for(auto& cell : mesh.cells){
		totalError += abs( cell.var[controls.P] - cell.var[controls.T] )*cell.volume;
	}
	outputFile << totalError << endl;
	outputFile << "     </DataArray>" << endl;
	outputFile << "    </FieldData>" << endl;
	
	
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
		outputFile << scientific << cell.var[controls.U] << " " << cell.var[controls.V] << " " << cell.var[controls.W] << " ";
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
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"massFraction-" << species[0].name << "\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells) outputFile << scientific << cell.var[controls.MF[0]] << " ";
	outputFile << endl;
	outputFile << "     </DataArray>" << endl;
	
	outputFile << "     <DataArray type=\"Float64\" Name=\"grad-pressure\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << cell.var[controls.UDV[0]] << " " << cell.var[controls.UDV[1]] << " " << cell.var[controls.UDV[2]] << " ";
	}
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



	outputFile << "     <DataArray type=\"Float64\" Name=\"error\" format=\"ascii\">" << endl;
	for(auto& cell : mesh.cells){
		outputFile << scientific << abs(cell.var[controls.P]-cell.var[controls.T]) << " ";
	}
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
		string filenamePvtu = folder + "../plot.0_Laplace.pvtu";
		
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
		outputFile << "    <PDataArray type=\"Float64\" Name=\"totalError\"/>" << endl;
		outputFile << "   </PFieldData>" << endl;
		
		outputFile << "   <PPoints>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Points\"/>" << endl;
		outputFile << "   </PPoints>" << endl;
		for(int ip=0; ip<size; ++ip){
			string filenamevtus = "./0_Laplace/plot.";
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
		outputFile << "    <PDataArray type=\"Float64\" Name=\"massFraction-" << species[0].name << "\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"grad-pressure\" NumberOfComponents=\"3\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellLevels\"/>" << endl;
		outputFile << "    <PDataArray type=\"Int64\" Name=\"cellGroups\"/>" << endl;
		outputFile << "    <PDataArray type=\"Float64\" Name=\"error\"/>" << endl;
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
