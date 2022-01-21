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
SEMO_Solvers_Builder solvers;


	
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
			cout << "| -c \"real num.\" : CFL" << endl;
			cout << "| -i \"int. num.\" : max iteration number" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		MPI_Finalize();
		// exit(EXIT_FAILURE);
		exit(EXIT_SUCCESS);
		return EXIT_SUCCESS;
	}
	
	
	
	double CFL = stod(mapArgv["-c"]);
	int iterMax = stoi(mapArgv["-i"]);
	// double cflInput = stod(mapArgv["-c"]);
	// int isave = stoi(mapArgv["-n"]);
	// inpTimeStep
	
	
	// 초기화
	controls.readSpecies(species);
	controls.readConfigures();
	controls.setValues(species);
	
	
	// 메쉬 및 데이터 읽기
	{
		
		SEMO_Mesh_Load load;
		double starttime = stod(controls.startFrom);
		string foldername;
		std::ostringstream streamObj;
		streamObj << starttime;
		foldername = "./save/" + streamObj.str() + "/";
		if(starttime == 0.0) foldername = "./save/0/";
		
		load.vtu(foldername, mesh, controls, species);
		
		solvers.calcCellEOSVF(mesh, controls, species);
		
		// for(auto& cell : mesh.cells){
			// cell.var[controls.Qm[0]] = cell.var[controls.U];
			// cell.var[controls.Qm[1]] = cell.var[controls.V];
			// cell.var[controls.Qm[2]] = cell.var[controls.W];
			// for(int i=0; i<controls.nSp-1; ++i){
				// cell.var[controls.Qm[3+i]] = cell.var[controls.VF[i]];
			// }
			
			// cell.var[controls.Qn[0]] = cell.var[controls.U];
			// cell.var[controls.Qn[1]] = cell.var[controls.V];
			// cell.var[controls.Qn[2]] = cell.var[controls.W];
			// for(int i=0; i<controls.nSp-1; ++i){
				// cell.var[controls.Qn[3+i]] = cell.var[controls.VF[i]];
			// }
		// }
	}
	
	// geometric 계산
	{
		geometric.init(mesh);
		math.initLeastSquare(mesh);
	}

	
	// 체적분율 압축 방정식 솔버, OpenFoam compression term 참고
	{
		
		
		int gradIterMax_LG = 1;
		int gradIterMax_GG = 1;
		

		// 로컬 시간스템 구하기
		double tau_ls = 1.e15;
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			double maxA=0.0;
			double minA=1000000.0;
			for(auto& j : cell.faces){
				maxA = max(maxA, mesh.faces[j].area);
				minA = min(minA, mesh.faces[j].area);
			}
			double minX = cell.volume / maxA;
			tau_ls = min(tau_ls,min(pow(cell.volume,0.3), minX));
		}
		double tau_ls_glob;
		MPI_Allreduce(&tau_ls, &tau_ls_glob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		tau_ls = tau_ls_glob;
		
	
		
		
		for(int outerIter=0; outerIter<iterMax; ++outerIter){
			
			vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
			for(int iter=0; iter<gradIterMax_LG; ++iter){
				vector<double> dummyVec;
				math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
					controls.VF[0], controls.fVF[0], dummyVec, gradAi);
			}
			// for(int iter=0; iter<gradIterMax_GG; ++iter){
				// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
			// }
			
			
			vector<double> UN(mesh.cells.size());
			vector<double> VN(mesh.cells.size());
			vector<double> WN(mesh.cells.size());
			vector<double> alphaN(mesh.cells.size());
			for(int i=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				
				vector<double> sufNorVec(3,0.0);
				double magGrad = sqrt(
					pow(gradAi[i][0],2.0)+
					pow(gradAi[i][1],2.0)+
					pow(gradAi[i][2],2.0));
				if(magGrad != 0.0){
					for(int ii=0; ii<3; ++ii){
						sufNorVec[ii] = gradAi[i][ii]/magGrad;
					}
				}
				
				// double magVel = sqrt(
					// pow(cell.var[controls.U],2.0)+
					// pow(cell.var[controls.V],2.0)+
					// pow(cell.var[controls.W],2.0));
				double magVel = 1.0;
				UN[i] = magVel*sufNorVec[0];
				VN[i] = magVel*sufNorVec[1];
				WN[i] = magVel*sufNorVec[2];
				
				double alpha = cell.var[controls.VF[0]];
				alphaN[i] = alpha*(1.0-alpha);
				
			}
			// processor faces
			vector<double> UN_recv, VN_recv, WN_recv, alphaN_recv;
			if(size>1){
				vector<double> UN_send, VN_send, WN_send, alphaN_send;
				for(int i=0; i<mesh.faces.size(); ++i){
					auto& face = mesh.faces[i];
					if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						UN_send.push_back(UN[face.owner]);
						VN_send.push_back(VN[face.owner]);
						WN_send.push_back(WN[face.owner]);
						alphaN_send.push_back(alphaN[face.owner]);
					}
				}
				mpi.setProcsFaceDatas(
							UN_send, UN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							VN_send, VN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							WN_send, WN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
				mpi.setProcsFaceDatas(
							alphaN_send, alphaN_recv,
							mesh.countsProcFaces, mesh.countsProcFaces, 
							mesh.displsProcFaces, mesh.displsProcFaces);
			}
			
			
			
			vector<double> resi(mesh.cells.size(),0.0);
			for(int i=0, ip=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				vector<double> nvec(3,0.0);
				nvec[0] = face.unitNormals[0];
				nvec[1] = face.unitNormals[1];
				nvec[2] = face.unitNormals[2];
			
				double wCL = face.wC;
				double wCR = 1.0-wCL;
				
				double UnF = wCL*(UN[face.owner]*nvec[0]+VN[face.owner]*nvec[1]+WN[face.owner]*nvec[2]);
				double alphaL = alphaN[face.owner];
				double alphaR = 0.0;
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					UnF += wCR*(UN[face.neighbour]*nvec[0]+VN[face.neighbour]*nvec[1]+WN[face.neighbour]*nvec[2]);
					alphaR = alphaN[face.neighbour];
				}
				else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					UnF += wCR*(UN_recv[ip]*nvec[0]+VN_recv[ip]*nvec[1]+WN_recv[ip]*nvec[2]);
					alphaR = alphaN_recv[ip];
				}
				else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
					UnF = 0.0;
				}
				
				double weightL = (UnF > 0.0) ? 1.0 : 0.0;
				double weightR = 1.0 - weightL;
			
				double alphaF = weightL*alphaL + weightR*alphaR;
				
				if(alphaL*(1.0-alphaL) < 1.e-15) alphaF = 0.0;
				if(alphaR*(1.0-alphaR) < 1.e-15) alphaF = 0.0;
				
				double flux = alphaF*UnF*face.area;
				
				resi[face.owner] -= flux;
				if(face.getType() == SEMO_Types::INTERNAL_FACE){
					resi[face.neighbour] += flux;
				}
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
				
			}
			
			double residual =0.0;
			for(int i=0; i<mesh.cells.size(); ++i){
				auto& cell = mesh.cells[i];
				
				double updateVar = CFL * tau_ls * resi[i] / cell.volume;
				
				cell.var[controls.VF[0]] += updateVar;
				
				residual += (updateVar*updateVar);
				
				cell.var[controls.VF[0]] = max(0.0,min(1.0,cell.var[controls.VF[0]]));
				
			}
			residual = sqrt(residual);
			double residualReduced;
			MPI_Allreduce(&residual, &residualReduced, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			controls.residual = residualReduced;
			
			if(rank==0) cout << 
			" | iteration = " << outerIter << 
			" | residual = " << residualReduced << endl;
		
		}
		
		
	}
		
		
	// 저장
	{
		string foldername;
		std::ostringstream streamObj;
		streamObj << controls.time;
		foldername = "./save/" + streamObj.str() + "_CompressVF/";
		
		// solvers.setCompValuesLeftRightFace(mesh, controls, species);
		
		save.vtu(foldername, mesh, controls, species);
		
		save.cellDataOverTime(foldername, controls);
		
		// save.particles(foldername, mesh, controls, species);
	}



	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	return 0;
}














