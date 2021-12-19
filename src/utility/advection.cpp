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

#include <zlib.h>

using namespace std;

#include "../mesh/build.h" 
#include "../load/load.h" 
#include "../mesh/geometric.h" 

#include "../solvers/build.h" 
#include "../controls/build.h" 
#include "../variables/build.h" 
#include "../mesh/polyAMR.h"


void saveInitialField(string folder);
	
SEMO_Mesh_Builder mesh;
vector<SEMO_Species> species;
SEMO_Controls_Builder controls;
SEMO_MPI_Builder mpi;
SEMO_Mesh_Load load;
SEMO_Utility_Math math;
SEMO_Mesh_Geometric geometric;
SEMO_Mesh_Save save;



void initialCD_VF(){

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		// // 1 test
		// if( pow(cell.x-0.25, 2.0) + pow(cell.y-0.25, 2.0) < 0.15*0.15 ){
			// cell.var[controls.VF[0]] = 1.0;
		// }
		// else{
			// cell.var[controls.VF[0]] = 0.0;
		// }
		
		// if( cell.x-0.25 > -0.075*0.5 && cell.x-0.25 < 0.075*0.5 &&
			// cell.y-0.25 < 0.1 ){
			// cell.var[controls.VF[0]] = 0.0;
		// }
		
		
		// // advection of Zalesak's disk in a vortex
		// if( pow(cell.x-0.5, 2.0) + pow(cell.y-0.75, 2.0) < 0.15*0.15 ){
			// cell.var[controls.VF[0]] = 1.0;
		// }
		// else{
			// cell.var[controls.VF[0]] = 0.0;
		// }
		
		// if( cell.x-0.5 > -0.075*0.5 && cell.x-0.5 < 0.075*0.5 &&
			// cell.y-0.75 < 0.1 ){
			// cell.var[controls.VF[0]] = 0.0;
		// }
		
		
		// Advection of a circular interface in shear flow
		if( pow(cell.x-0.5, 2.0) + pow(cell.y-0.25, 2.0) < 0.2*0.2 ){
			cell.var[controls.VF[0]] = 1.0;
		}
		else{
			cell.var[controls.VF[0]] = 0.0;
		}
		
		
	}
	
}

void initialCD_VEL(){

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		// // 1 test
		// cell.var[controls.U] = 1.0;
		// cell.var[controls.V] = 1.0;
		// cell.var[controls.W] = 0.0;
		
		
		// // advection of Zalesak's disk in a vortex
		// cell.var[controls.U] = -2.0*3.141592*(cell.y-0.5);
		// cell.var[controls.V] = 2.0*3.141592*(cell.x-0.5);
		// cell.var[controls.W] = 0.0;
		
		
		// Advection of a circular interface in shear flow
		if(controls.time < 2.5){
			cell.var[controls.U] = sin(3.141592*cell.x)*cos(3.141592*cell.y);
			cell.var[controls.V] = -cos(3.141592*cell.x)*sin(3.141592*cell.y);
			cell.var[controls.W] = 0.0;
		}
		else{
			cell.var[controls.U] = -sin(3.141592*cell.x)*cos(3.141592*cell.y);
			cell.var[controls.V] = cos(3.141592*cell.x)*sin(3.141592*cell.y);
			cell.var[controls.W] = 0.0;
		}
		
		
	}
	
}
		
	
	
	
	
	
	
	
	
	
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
			cout << "| -m \"real num.\" : max time" << endl;
			cout << "| -t \"real num.\" : save time interval" << endl;
			cout << "| -c \"real num.\" : cfl" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	
	
	double timeMax = stod(mapArgv["-m"]);
	double saveInterval = stod(mapArgv["-t"]);
	double cflInput = stod(mapArgv["-c"]);
	// int isave = stoi(mapArgv["-n"]);
	// inpTimeStep
	
	

		
	controls.readSpecies(species);
	
	
	controls.readConfigures();
	
	controls.setValues(species);
	
	
	SEMO_Solvers_Builder solvers;
	
	
	
	load.vtu("./save/0/", mesh, controls, species);
	
	solvers.calcCellEOSVF(mesh, controls, species);
	

		
	
	geometric.init(mesh);
	
	math.initLeastSquare(mesh);
	
	
	
	
	// initial
	for(int jjj=0; jjj<3; ++jjj){
		
		
		initialCD_VF();
		initialCD_VEL();

		controls.iterReal = controls.intervalRefine-1;
		if
		( (controls.iterReal+1) % controls.intervalRefine == 0 ||
		  (controls.iterReal+1) % controls.intervalUnrefine == 0)
		{ 
			int numnum=0; 
			SEMO_Poly_AMR_Builder AMR;
			AMR.polyAMR(mesh, controls, species, 0);
			mesh.cellsGlobal();
			solvers.calcIncomCellEOSVF(mesh, controls, species);
		} 
	}
	

	initialCD_VF();
	initialCD_VEL();
	
	
	
	
	
	// initial save
	{
		string foldername;
		std::ostringstream streamObj;
		streamObj << 0.0;
		foldername = "./save/" + streamObj.str() + "_Advection/";
		
		// solvers.setCompValuesLeftRightFace(mesh, controls, species);
		
		save.vtu(foldername, mesh, controls, species);
		
	}
	
	// cout << "AAAAAA" << endl;

	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		cell.var[controls.Qm[0]] = cell.var[controls.VF[0]];
		cell.var[controls.oldVF[0]] = cell.var[controls.VF[0]];
	}
	
	controls.time = 0.0;
	int iTime = 0;
	controls.iterReal = -1;
	
	while(controls.time<timeMax){
		
		++iTime;
		++controls.iterReal;
		
		
		// timestep 계산
		double CFL = cflInput;
		double dtMin = 10000.0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			double normalMagVel = mesh.cells[face.owner].var[controls.U]*face.unitNormals[0];
			normalMagVel += mesh.cells[face.owner].var[controls.V]*face.unitNormals[1];
			normalMagVel += mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
			normalMagVel = abs(normalMagVel);
			
			double normalFlux = normalMagVel*face.area;
			
			double dt = mesh.cells[face.owner].volume / (normalFlux + 1.e-100);
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				dt = min(dt, mesh.cells[face.neighbour].volume / (normalFlux + 1.e-100));
				
				normalMagVel = mesh.cells[face.neighbour].var[controls.U]*face.unitNormals[0];
				normalMagVel += mesh.cells[face.neighbour].var[controls.V]*face.unitNormals[1];
				normalMagVel += mesh.cells[face.neighbour].var[controls.W]*face.unitNormals[2];
				normalMagVel = abs(normalMagVel);
				normalFlux = normalMagVel*face.area;
				
				dt = min(dt, mesh.cells[face.owner].volume / (normalFlux + 1.e-100));
				dt = min(dt, mesh.cells[face.neighbour].volume / (normalFlux + 1.e-100));
			}
			dtMin = min(dtMin, dt);
			
		}
		double dtMinReduced;
		MPI_Allreduce(&dtMin, &dtMinReduced, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		controls.old2TimeStep = controls.oldTimeStep;
		controls.oldTimeStep = controls.timeStep;
		controls.timeStep = CFL * dtMinReduced;
		if(iTime==1){
			controls.old2TimeStep = controls.timeStep;
			controls.oldTimeStep = controls.timeStep;
		}
		
		
		
		
		controls.time += controls.timeStep;
		// controls.timeStep = inpTimeStep;



		if
		( (controls.iterReal+1) % controls.intervalRefine == 0 ||
		  (controls.iterReal+1) % controls.intervalUnrefine == 0)
		{ 
			int numnum=0; 
			SEMO_Poly_AMR_Builder AMR;
			AMR.polyAMR(mesh, controls, species, 0);
			mesh.cellsGlobal();
			solvers.calcIncomCellEOSVF(mesh, controls, species);
		} 
		initialCD_VEL();




		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			cell.var[controls.Qm[0]] = cell.var[controls.oldVF[0]];
			cell.var[controls.oldVF[0]] = cell.var[controls.VF[0]];
		}



	



		// gradient VF
		vector<vector<double>> gradVF(mesh.cells.size(),vector<double>(3,0.0));
		
		{
			vector<double> dummyVec;
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.VF[0], controls.fVF[0], dummyVec, gradVF);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.VF[0], controls.fVF[0], dummyVec, gradVF);
		}
		
		vector<double> gradVFx_recv;
		vector<double> gradVFy_recv;
		vector<double> gradVFz_recv;
		if(size>1){
			// processor faces
			vector<double> gradVFx_send;
			vector<double> gradVFy_send;
			vector<double> gradVFz_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					gradVFx_send.push_back(gradVF[face.owner][0]);
					gradVFy_send.push_back(gradVF[face.owner][1]);
					gradVFz_send.push_back(gradVF[face.owner][2]);
				}
			}
			// SEMO_MPI_Builder mpi;
			
			mpi.setProcsFaceDatas(
						gradVFx_send, gradVFx_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradVFy_send, gradVFy_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradVFz_send, gradVFz_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			gradVFx_send.clear();
			gradVFy_send.clear();
			gradVFz_send.clear();	
			
			int proc_num=0;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					vector<double> tmp;
					tmp.push_back(gradVFx_recv[proc_num]);
					tmp.push_back(gradVFy_recv[proc_num]);
					tmp.push_back(gradVFz_recv[proc_num]);
					gradVF.push_back(tmp);
					++proc_num;
				}
			}
		}
		
		
		
		
		
		
		
		
		
		

		// calc recon. zero order
		solvers.reconIncomZeroOrder(mesh, controls, species);





		// gradient U, V, W
		vector<vector<double>> gradU(mesh.cells.size(),vector<double>(3,0.0));
		vector<vector<double>> gradV(mesh.cells.size(),vector<double>(3,0.0));
		vector<vector<double>> gradW(mesh.cells.size(),vector<double>(3,0.0));
		// math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
		// math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
		// math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
		// math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
		// math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
		// math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
		// math.calcMGG(mesh, controls.U, controls.fU, 10, 1.e-8, gradU);
		// math.calcMGG(mesh, controls.V, controls.fV, 10, 1.e-8, gradV);
		// math.calcMGG(mesh, controls.W, controls.fW, 10, 1.e-8, gradW);
		
		{
			vector<double> dummyVec;
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.U, controls.fU, dummyVec, gradU);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.V, controls.fV, dummyVec, gradV);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.W, controls.fW, dummyVec, gradW);
				
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.U, controls.fU, dummyVec, gradU);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.V, controls.fV, dummyVec, gradV);
			math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
				controls.W, controls.fW, dummyVec, gradW);
		}
		
		vector<double> gradUx_recv, gradUy_recv, gradUz_recv;
		vector<double> gradVx_recv, gradVy_recv, gradVz_recv;
		vector<double> gradWx_recv, gradWy_recv, gradWz_recv;
		if(size>1){
			// processor faces
			vector<double> gradUx_send, gradUy_send, gradUz_send;
			vector<double> gradVx_send, gradVy_send, gradVz_send;
			vector<double> gradWx_send, gradWy_send, gradWz_send;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					gradUx_send.push_back(gradU[face.owner][0]);
					gradUy_send.push_back(gradU[face.owner][1]);
					gradUz_send.push_back(gradU[face.owner][2]);
					
					gradVx_send.push_back(gradV[face.owner][0]);
					gradVy_send.push_back(gradV[face.owner][1]);
					gradVz_send.push_back(gradV[face.owner][2]);
					
					gradWx_send.push_back(gradW[face.owner][0]);
					gradWy_send.push_back(gradW[face.owner][1]);
					gradWz_send.push_back(gradW[face.owner][2]);
				}
			}
			// SEMO_MPI_Builder mpi;
			
			mpi.setProcsFaceDatas(
						gradUx_send, gradUx_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradUy_send, gradUy_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradUz_send, gradUz_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			gradUx_send.clear();
			gradUy_send.clear();
			gradUz_send.clear();
						
			mpi.setProcsFaceDatas(
						gradVx_send, gradVx_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradVy_send, gradVy_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradVz_send, gradVz_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			gradVx_send.clear();
			gradVy_send.clear();
			gradVz_send.clear();
						
			mpi.setProcsFaceDatas(
						gradWx_send, gradWx_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradWy_send, gradWy_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			mpi.setProcsFaceDatas(
						gradWz_send, gradWz_recv,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			gradWx_send.clear();
			gradWy_send.clear();
			gradWz_send.clear();
		}
		



		// cout << "AAAAAAAAA" << endl;
		
		
		// NVD, MSTACS
		solvers.calcMSTACS(mesh, controls, controls.VF[0], controls.fVF[0], gradVF, controls.fVF[0]);
		// solvers.calcQUICK(mesh, controls, controls.VF[0], controls.fVF[0], gradVF);
		
			
		
		// cout << "BBBBBBBB" << endl;
	
	
	

	vector<double> maxU(mesh.cells.size(),-1.e15);
	vector<double> minU(mesh.cells.size(),+1.e15);
	vector<double> maxV(mesh.cells.size(),-1.e15);
	vector<double> minV(mesh.cells.size(),+1.e15);
	vector<double> maxW(mesh.cells.size(),-1.e15);
	vector<double> minW(mesh.cells.size(),+1.e15);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		maxU[i] = cell.var[controls.U];
		minU[i] = cell.var[controls.U];
		maxV[i] = cell.var[controls.V];
		minV[i] = cell.var[controls.V];
		maxW[i] = cell.var[controls.W];
		minW[i] = cell.var[controls.W];
		for(auto j : cell.stencil){
			auto& cellSten = mesh.cells[j];
			maxU[i] = max(maxU[i],cellSten.var[controls.U]);
			minU[i] = min(minU[i],cellSten.var[controls.U]);
			maxV[i] = max(maxV[i],cellSten.var[controls.V]);
			minV[i] = min(minV[i],cellSten.var[controls.V]);
			maxW[i] = max(maxW[i],cellSten.var[controls.W]);
			minW[i] = min(minW[i],cellSten.var[controls.W]);
		}
	}
	if(size>1){
		// processor faces
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				maxU[face.owner] = max(maxU[face.owner],face.varR[controls.fU]);
				minU[face.owner] = min(minU[face.owner],face.varR[controls.fU]);
				maxV[face.owner] = max(maxV[face.owner],face.varR[controls.fV]);
				minV[face.owner] = min(minV[face.owner],face.varR[controls.fV]);
				maxW[face.owner] = max(maxW[face.owner],face.varR[controls.fW]);
				minW[face.owner] = min(minW[face.owner],face.varR[controls.fW]);
			}
		}
	}
	vector<double> maxU_recv;
	vector<double> minU_recv;
	vector<double> maxV_recv;
	vector<double> minV_recv;
	vector<double> maxW_recv;
	vector<double> minW_recv;
	if(size>1){
		// processor faces
		vector<double> maxU_send;
		vector<double> minU_send;
		vector<double> maxV_send;
		vector<double> minV_send;
		vector<double> maxW_send;
		vector<double> minW_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				maxU_send.push_back(maxU[face.owner]);
				minU_send.push_back(minU[face.owner]);
				maxV_send.push_back(maxV[face.owner]);
				minV_send.push_back(minV[face.owner]);
				maxW_send.push_back(maxW[face.owner]);
				minW_send.push_back(minW[face.owner]);
			}
		}
		
		mpi.setProcsFaceDatas(
					maxU_send, maxU_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					minU_send, minU_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		maxU_send.clear();
		minU_send.clear();
		
		mpi.setProcsFaceDatas(
					maxV_send, maxV_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					minV_send, minV_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		maxV_send.clear();
		minV_send.clear();
		
		mpi.setProcsFaceDatas(
					maxW_send, maxW_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					minW_send, minW_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		maxW_send.clear();
		minW_send.clear();
	}
	
	
	
		
		// diagonal terms
		int B_n = 1;
		int A_n = B_n * B_n;
		
			
		int proc_num=0;
		


		vector<int> A_rows(mesh.cells.size()*A_n, 0);
		vector<int> A_cols(mesh.cells.size()*A_n, 0);
		vector<double> A_vals(mesh.cells.size()*A_n, 0.0);
		vector<double> B(mesh.cells.size()*B_n, 0.0);
		

		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];

			int ijStart_local = B_n*(i) - 1;
			int ijStart_global = B_n*mesh.startCellGlobal + ijStart_local;
			int Astart = A_n*(i) - 1;
			// int id = Astart;
			
			int id = 0;
			int i_glo = 0;
			int j_glo = 0;
			int i_loc = i;
			int str_glo = mesh.startCellGlobal*B_n;
			int step_loc = mesh.cells.size();
			int step_glo = mesh.ncellsTotal;
			
			double volume = cell.volume;
			double timeStep = controls.timeStep;
			
			double VF = cell.var[controls.VF[0]];
			double oldVF = cell.var[controls.oldVF[0]];
			double old2VF = cell.var[controls.Qm[0]];
			
			
			
			

			double beta_phi = 1.0;
			double phi_SOU = ( oldVF - old2VF ) /
						(abs(VF-old2VF)+1.e-200)*( VF>old2VF ? 1.0 : -1.0 );
			if(0.0 < phi_SOU && phi_SOU > 1.0){
				beta_phi = 0.0;
			}
			else{
				beta_phi = 1.0;
			}
			double coeff_old1 = controls.oldTimeStep/(controls.timeStep+controls.oldTimeStep);
			double coeff_old2 = controls.old2TimeStep/(controls.oldTimeStep+controls.old2TimeStep);
			double cons_np12 = VF + beta_phi * coeff_old1 * (VF-oldVF);
			double cons_nm12 = oldVF + beta_phi * coeff_old2 * (oldVF-old2VF);
			
			// cons_np12 = max(0.0,min(1.0,cons_np12));
			// cons_nm12 = max(0.0,min(1.0,cons_nm12));
			
			// volume fraction
			id = step_loc*(B_n*0+0) + i_loc; i_glo = str_glo + step_loc*0 + i_loc; j_glo = str_glo + step_loc*0 + i_loc;
			A_rows[id] = i_glo; A_cols[id] = j_glo;
			// A_vals[id] = volume/timeStep;
			A_vals[id] = 1.5*volume/timeStep;
			// A_vals[id] = (1.0+beta_phi*coeff_old1)*volume/timeStep;
			
			
			// B[step_loc*0 + i_loc] = -(VF - oldVF)*volume/timeStep;
			B[step_loc*0 + i_loc] = -(1.5*VF - 2.0*oldVF + 0.5*old2VF)*volume/timeStep;
			// B[step_loc*0 + i_loc] = -(cons_np12 - cons_nm12)*volume/timeStep;
			
			
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
			
			double wCL = face.wC;
			// double wCL = 0.5;
			double wCR = 1.0-wCL;
		
			double UL = face.varL[controls.fU];
			double VL = face.varL[controls.fV];
			double WL = face.varL[controls.fW];
			double VFL = face.varL[controls.fVF[0]];
			
			double UR = face.varR[controls.fU];
			double VR = face.varR[controls.fV];
			double WR = face.varR[controls.fW];
			double VFR = face.varR[controls.fVF[0]];
			
			double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
			double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
			
			double UnF = wCL*UnL+wCR*UnR;
			
			
			double dPN = face.magPN;
			double alpha = face.alphaF;
			
			// properties
			double orgUL = mesh.cells[face.owner].var[controls.U];
			double orgVL = mesh.cells[face.owner].var[controls.V];
			double orgWL = mesh.cells[face.owner].var[controls.W];
			double orgUR = 0.0;
			double orgVR = 0.0;
			double orgWR = 0.0;
			vector<double> gradUL(3,0.0);
			vector<double> gradVL(3,0.0);
			vector<double> gradWL(3,0.0);
			vector<double> gradUR(3,0.0);
			vector<double> gradVR(3,0.0);
			vector<double> gradWR(3,0.0);
			
			vector<double> gradVFL(3,0.0);
			gradVFL[0] = gradVF[face.owner][0];
			gradVFL[1] = gradVF[face.owner][1];
			gradVFL[2] = gradVF[face.owner][2];
			vector<double> gradVFR(3,0.0);
			
			double minUf = minU[face.owner];
			double maxUf = maxU[face.owner];
			double minVf = minV[face.owner];
			double maxVf = maxV[face.owner];
			double minWf = minW[face.owner];
			double maxWf = maxW[face.owner];
			
			for(int ii=0; ii<3; ++ii){
				gradUL[ii] = gradU[face.owner][ii];
				gradVL[ii] = gradV[face.owner][ii];
				gradWL[ii] = gradW[face.owner][ii];
			}
			
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				gradVFR[0] = gradVF[face.neighbour][0];
				gradVFR[1] = gradVF[face.neighbour][1];
				gradVFR[2] = gradVF[face.neighbour][2];
				
				orgUR = mesh.cells[face.neighbour].var[controls.U];
				orgVR = mesh.cells[face.neighbour].var[controls.V];
				orgWR = mesh.cells[face.neighbour].var[controls.W];
				
				for(int ii=0; ii<3; ++ii){
					gradUR[ii] = gradU[face.neighbour][ii];
					gradVR[ii] = gradV[face.neighbour][ii];
					gradWR[ii] = gradW[face.neighbour][ii];
				}
				
				minUf = min(minUf,minU[face.neighbour]);
				maxUf = max(maxUf,maxU[face.neighbour]);
				minVf = min(minVf,minV[face.neighbour]);
				maxVf = max(maxVf,maxV[face.neighbour]);
				minWf = min(minWf,minW[face.neighbour]);
				maxWf = max(maxWf,maxW[face.neighbour]);
				
			}
			else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				gradVFR[0] = gradVFx_recv[proc_num];
				gradVFR[1] = gradVFy_recv[proc_num];
				gradVFR[2] = gradVFz_recv[proc_num];
				
				orgUR = UR;
				orgVR = VR;
				orgWR = WR;
				
				gradUR[0] = gradUx_recv[proc_num]; gradUR[1] = gradUy_recv[proc_num]; gradUR[2] = gradUz_recv[proc_num];
				gradVR[0] = gradVx_recv[proc_num]; gradVR[1] = gradVy_recv[proc_num]; gradVR[2] = gradVz_recv[proc_num];
				gradWR[0] = gradWx_recv[proc_num]; gradWR[1] = gradWy_recv[proc_num]; gradWR[2] = gradWz_recv[proc_num];
				
				minUf = min(minUf,minU_recv[proc_num]);
				maxUf = max(maxUf,maxU_recv[proc_num]);
				minVf = min(minVf,minV_recv[proc_num]);
				maxVf = max(maxVf,maxV_recv[proc_num]);
				minWf = min(minWf,minW_recv[proc_num]);
				maxWf = max(maxWf,maxW_recv[proc_num]);
				
			}

			// recon.
			double reconUL = UL;
			double reconVL = VL;
			double reconWL = WL;
			double reconUR = UR;
			double reconVR = VR;
			double reconWR = WR;
			// for(int ii=0; ii<3; ++ii){
				// reconUL += gradUL[ii]*face.vecPF[ii];
				// reconVL += gradVL[ii]*face.vecPF[ii];
				// reconWL += gradWL[ii]*face.vecPF[ii];
				// reconUR += gradUR[ii]*face.vecNF[ii];
				// reconVR += gradVR[ii]*face.vecNF[ii];
				// reconWR += gradWR[ii]*face.vecNF[ii];
			// }
			for(int ii=0; ii<3; ++ii){
				reconUL += gradUL[ii]*face.vecSkewness[ii];
				reconVL += gradVL[ii]*face.vecSkewness[ii];
				reconWL += gradWL[ii]*face.vecSkewness[ii];
				reconUR += gradUR[ii]*face.vecSkewness[ii];
				reconVR += gradVR[ii]*face.vecSkewness[ii];
				reconWR += gradWR[ii]*face.vecSkewness[ii];
			}
			// reconUL = max(minUf,min(maxUf,reconUL));
			// reconUR = max(minUf,min(maxUf,reconUR));
			// reconVL = max(minVf,min(maxVf,reconVL));
			// reconVR = max(minVf,min(maxVf,reconVR));
			// reconWL = max(minWf,min(maxWf,reconWL));
			// reconWR = max(minWf,min(maxWf,reconWR));
			

			
			UnF = wCL*(reconUL*nvec[0]+reconVL*nvec[1]+reconWL*nvec[2]);
			UnF += wCR*(reconUR*nvec[0]+reconVR*nvec[1]+reconWR*nvec[2]);
			// UnF = 0.5*(reconUL*nvec[0]+reconVL*nvec[1]+reconWL*nvec[2]);
			// UnF += 0.5*(reconUR*nvec[0]+reconVR*nvec[1]+reconWR*nvec[2]);
			

			// // skewness
			// for(int ii=0; ii<3; ++ii){
				// UnF += wCL*gradUL[ii]*face.vecSkewness[ii] * nvec[0];
				// UnF += wCL*gradVL[ii]*face.vecSkewness[ii] * nvec[1];
				// UnF += wCL*gradWL[ii]*face.vecSkewness[ii] * nvec[2];
				// UnF += wCR*gradUR[ii]*face.vecSkewness[ii] * nvec[0];
				// UnF += wCR*gradVR[ii]*face.vecSkewness[ii] * nvec[1];
				// UnF += wCR*gradWR[ii]*face.vecSkewness[ii] * nvec[2];
			// }

			double weightL = (UnF > 0.0) ? 1.0 : 0.0;
			double weightR = 1.0 - weightL;
			
			double VFF = weightL*VFL + weightR*VFR;
			
			// // skewness
			// VFF += weightL*gradVFL[0]*face.vecSkewness[0];
			// VFF += weightL*gradVFL[1]*face.vecSkewness[1];
			// VFF += weightL*gradVFL[2]*face.vecSkewness[2];
			// VFF += weightR*gradVFR[0]*face.vecSkewness[0];
			// VFF += weightR*gradVFR[1]*face.vecSkewness[1];
			// VFF += weightR*gradVFR[2]*face.vecSkewness[2];
			// VFF = max(0.0,min(1.0,VFF));
			// // if(VFF>max(VFL,VFR)) VFF = max(VFL,VFR);
			// // if(VFF<min(VFL,VFR)) VFF = min(VFL,VFR);
			
			
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
			// continuity
			// VF'
			id = step_loc_L*(B_n*0+0) + i_loc_L; i_glo = str_glo_L + step_loc_L*0 + i_loc_L; j_glo = str_glo_R + step_loc_R*0 + i_loc_R;
			// A_vals[id] += ( weightL * UnF * area );
			// A_vals[id] += ( 0.5 * weightL * UnF * area );
			// A_vals[id] += ( ( weightL * face.varL[controls.fVF_NVD] + weightR * face.varR[controls.fVF_NVD] ) * UnF * area );
			
			A_rows.push_back(i_glo); A_cols.push_back(j_glo);
			// A_vals.push_back(( weightR * UnF * area ));
			// A_vals.push_back(( 0.5 * weightR * UnF * area ));
			// A_vals.push_back(( ( weightL * (1.0-face.varL[controls.fVF_NVD]) + weightR * (1.0-face.varR[controls.fVF_NVD]) ) * UnF * area ));
			
			if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
				A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
			}

			if(face.getType() == SEMO_Types::INTERNAL_FACE) {
				id = step_loc_R*(B_n*0+0) + i_loc_R; i_glo = str_glo_R + step_loc_R*0 + i_loc_R; j_glo = str_glo_L + step_loc_L*0 + i_loc_L;
				// A_vals[id] -= ( weightR * UnF * area );
				// A_vals[id] -= ( 0.5 * weightR * UnF * area );
				// A_vals[id] -= ( ( weightL * (1.0-face.varL[controls.fVF_NVD]) + weightR * (1.0-face.varR[controls.fVF_NVD]) ) * UnF * area );
				
				A_rows.push_back(i_glo); A_cols.push_back(j_glo);
				// A_vals.push_back(-( weightL * UnF * area ));
				// A_vals.push_back(-( 0.5 * weightL * UnF * area ));
				// A_vals.push_back(-( ( weightL * face.varL[controls.fVF_NVD] + weightR * face.varR[controls.fVF_NVD] ) * UnF * area ));
				
				if( abs(A_vals.back()) < std::numeric_limits<double>::min() ){
					A_rows.pop_back(); A_cols.pop_back(); A_vals.pop_back();
				}
			}
			// ----------------------------
			B[step_loc_L*0 + i_loc_L] -= ( VFF * UnF * area );

			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				B[step_loc_R*0 + i_loc_R] += ( VFF * UnF * area );
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
					
					// int ijStartL_local = B_n*(face.owner) - 1;
					// vector<int> id(B_n,0);
					// id[0] = A_n*(face.owner) - 1;
					// for(int j=1; j<B_n; ++j){
						// id[j] = id[j-1] + 4;
					// }

					
					double area = face.area;
			
					double UnL = face.varL[controls.fU]*face.unitNormals[0] + 
								 face.varL[controls.fV]*face.unitNormals[1] + 
								 face.varL[controls.fW]*face.unitNormals[2];
					double UnR = face.varR[controls.fU]*face.unitNormals[0] + 
								 face.varR[controls.fV]*face.unitNormals[1] + 
								 face.varR[controls.fW]*face.unitNormals[2];
								
					double UnF = 0.5*UnL+0.5*UnR;
					
					double RhoL = face.varL[controls.fRho];
					double tmp1 = 1.0/RhoL*controls.timeStep;
					
					// for(int ii=0; ii<3; ++ii){
						// UnF += tmp1 * gradP[face.owner][ii]*face.unitNormals[ii];
					// }
					
					double VFF = 0.5*face.varL[controls.fVF[0]] + 0.5*face.varR[controls.fVF[0]];
				
					// volume-fraction coefficient
					double coeffVF = 0.0;
					if( boundary.type[controls.VF[0]] == "fixedValue" ){
						coeffVF = 0.0;
					}
					else if( boundary.type[controls.VF[0]] == "zeroGradient" ){
						coeffVF = 1.0;
						
					}
					else if( boundary.type[controls.VF[0]] == "inletOutlet" ){
						double ownNorVel =  
							mesh.cells[face.owner].var[controls.U]*face.unitNormals[0] +
							mesh.cells[face.owner].var[controls.V]*face.unitNormals[1] +
							mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
						coeffVF = (ownNorVel > 1.0) ? 1.0 : 0.0;
					}
					else {
						cerr << "| #Error : not defined B.C., var = " << controls.VF[0] << endl;
						MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					}
					
					
					int str_glo = mesh.startCellGlobal*B_n;
					int step_loc = mesh.cells.size();
					int step_glo = mesh.ncellsTotal;
					int i_loc_L = face.owner;
					int id=-1;
					
					// continuity
					id = step_loc*(B_n*0+0) + i_loc_L;
					// A_vals[id] += coeffVF * ( UnF * area );
					
					// convective term
					B[step_loc*0 + i_loc_L] -= ( VFF * UnF * area );
				}
				
				
			}
			
		}
		
		
		
		

		vector<double> resiVar(B_n*mesh.cells.size(),0.0);
		
		
		// solvers.solveAMGCL("volume_fraction", mesh, B_n, A_rows, A_cols, A_vals, B, resiVar);
		
		
		// update
		double residual = 0.0;
		for(int i=0; i<mesh.cells.size(); ++i){
			SEMO_Cell& cell = mesh.cells[i];
			
			resiVar[i] = B[i]/A_vals[i];
			
			cell.var[controls.VF[0]] += resiVar[i];
			// cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
			
			residual += abs(resiVar[i]);
		}
	


		double residualReduced;
		MPI_Allreduce(&residual, &residualReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		controls.residual = residualReduced;
		
		if(rank==0) cout << "| time = " << controls.time << " | iteration = " << iTime << " | residual = " << residualReduced << endl;
		
		
		

		// if(iTime % isave == 0 || time>=timeMax){
			// string foldername;
			// std::ostringstream streamObj;
			// streamObj << time;
			// foldername = "./save/" + streamObj.str() + "_Advection/";
			
			// // solvers.setCompValuesLeftRightFace(mesh, controls, species);
			
			// save.vtu(foldername, mesh, controls, species);
			
		// }
		
	
		{
			int jung = controls.time / saveInterval;
			double namuji = controls.time - (double)jung * saveInterval;
			
			
			save.cellDataToMemoryOverTime(mesh, controls);
			
			// cout << namuji << " " << controls.timeStep << endl;
			
			if(
			namuji < controls.timeStep &&
			namuji >= 0.0
			){
				string foldername;
				std::ostringstream streamObj;
				streamObj << controls.time;
				foldername = "./save/" + streamObj.str() + "_Advection/";
				
				// solvers.setCompValuesLeftRightFace(mesh, controls, species);
				
				save.vtu(foldername, mesh, controls, species);
				
				save.cellDataOverTime(foldername, controls);
				
				// save.particles(foldername, mesh, controls, species);
				
			}
		}
		
		
	}
	
	
	// cout << "GGGGGGGGGGG" << endl;
	// // saveInitialField("./save/0_Advection/");
	// {
		// string foldername;
		// std::ostringstream streamObj;
		// streamObj << time;
		// foldername = "./save/" + streamObj.str() + "_Advection/";
		
		// // solvers.setCompValuesLeftRightFace(mesh, controls, species);
		
		// save.vtu(foldername, mesh, controls, species);
		
	// }
	// cout << "GGGGGGGGGGG2" << endl;
	






	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
}














