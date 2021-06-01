#include "build.h"
#include <cmath>
#include <array>
#include <numeric>




void SEMO_Solvers_Builder::calcCurvature(
	SEMO_Mesh_Builder& mesh,
	int cn,
	vector<double>& kappa){
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	
	
	vector<double> smoothAi;
	vector<double> AiUp;
	vector<double> AiDown;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		smoothAi.push_back(mesh.cells[i].var[cn]);
		
	}
	
	//================================================
	// smoothing Ai
	for(int iter=0; iter<3; ++iter){
		
		AiUp.clear();
		AiDown.clear();
		for(int i=0; i<mesh.cells.size(); ++i){
			AiUp.push_back(0.0);
			AiDown.push_back(0.0);
		}
		
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
				double AiF = 
					face.wC*smoothAi[face.owner] +
					(1.0-face.wC)*smoothAi[face.neighbour];
				
				AiUp[face.owner] += AiF*face.area;
				AiUp[face.neighbour] += AiF*face.area;
				
				AiDown[face.owner] += face.area;
				AiDown[face.neighbour] += face.area;
			
			}
			
		}

		// boundary
		for(auto& boundary : mesh.boundary){
			
			if(boundary.neighbProcNo == -1){
				
				int str = boundary.startFace;
				int end = str + boundary.nFaces;
				
				for(int i=str; i<end; ++i){
					auto& face = mesh.faces[i];
					
					double AiF = smoothAi[face.owner];
					
					AiUp[face.owner] += AiF*face.area;
				
					AiDown[face.owner] += face.area;
					
				}
			}
		}
		
		if(size>1){
			// processor faces
			vector<double> sendValues;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){
					sendValues.push_back(smoothAi[face.owner]);
				}
			}
			vector<double> recvValues;
			mpi.setProcsFaceDatasDouble(
						sendValues, recvValues,
						mesh.countsProcFaces, mesh.countsProcFaces, 
						mesh.displsProcFaces, mesh.displsProcFaces);
			int num=0;
			for(int i=0; i<mesh.faces.size(); ++i){
				auto& face = mesh.faces[i];
				
				if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
					double AiF = face.wC*smoothAi[face.owner]+(1.0-face.wC)*recvValues[num];
					
					AiUp[face.owner] += AiF*face.area;
				
					AiDown[face.owner] += face.area;
					
					++num;
				}
			}
			
		}
	
		
		for(int i=0; i<mesh.cells.size(); ++i){
			smoothAi[i] = AiUp[i]/AiDown[i];
		}
		
	}
	
	
	//================================================
	// calc gauss-green gradient
	vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
	
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double varF = face.wC*smoothAi[face.owner]+(1.0-face.wC)*smoothAi[face.neighbour];
			
			for(int j=0; j<3; ++j){
				gradAi[face.owner][j] += 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				gradAi[face.neighbour][j] -= 
					varF*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			}
		}
	}

	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				
				double varF = smoothAi[face.owner];
				for(int j=0; j<3; ++j){
					gradAi[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				
			}
		}
	}
	
	if(size>1){
		// processor faces
		vector<double> sendValues;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues.push_back(smoothAi[face.owner]);
			}
		}
		vector<double> recvValues;
		mpi.setProcsFaceDatasDouble(
					sendValues, recvValues,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		int num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
				double varF = face.wC*smoothAi[face.owner]+(1.0-face.wC)*recvValues[num];
				
				for(int j=0; j<3; ++j){
					gradAi[face.owner][j] += 
						varF*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				++num;
			}
		}
		
	}
	


	//================================================
	// surfNormalVec = dAidx/|dAidx|
	vector<vector<double>> surfNormalVec;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		double magGrad = sqrt(
			pow(gradAi[i][0],2.0)+pow(gradAi[i][1],2.0)+pow(gradAi[i][2],2.0));
			
		vector<double> tmpVec;
		tmpVec.push_back(gradAi[i][0]/(magGrad+1.e-200));
		tmpVec.push_back(gradAi[i][1]/(magGrad+1.e-200));
		tmpVec.push_back(gradAi[i][2]/(magGrad+1.e-200));
		
		surfNormalVec.push_back(tmpVec);
		
	}
	
	//================================================
	// calc Curvature
	kappa.clear();
	kappa.resize(mesh.cells.size(),0.0);
	
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double varF[3];
			varF[0] = face.wC*surfNormalVec[face.owner][0]+(1.0-face.wC)*surfNormalVec[face.neighbour][0];
			varF[1] = face.wC*surfNormalVec[face.owner][1]+(1.0-face.wC)*surfNormalVec[face.neighbour][1];
			varF[2] = face.wC*surfNormalVec[face.owner][2]+(1.0-face.wC)*surfNormalVec[face.neighbour][2];
			
			for(int j=0; j<3; ++j){
				kappa[face.owner] += 
					varF[j]*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				kappa[face.neighbour] -= 
					varF[j]*face.unitNormals[j]*face.area/mesh.cells[face.neighbour].volume;
			}
		}
	}
	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				
				double varF[3];
				varF[0] = surfNormalVec[face.owner][0];
				varF[1] = surfNormalVec[face.owner][1];
				varF[2] = surfNormalVec[face.owner][2];
				for(int j=0; j<3; ++j){
					kappa[face.owner] += 
						varF[j]*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				
			}
		}
	}
	if(size>1){
		// processor faces
		vector<double> sendValues0, sendValues1, sendValues2;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues0.push_back(surfNormalVec[face.owner][0]);
				sendValues1.push_back(surfNormalVec[face.owner][1]);
				sendValues2.push_back(surfNormalVec[face.owner][2]);
			}
		}
		vector<double> recvValues0, recvValues1, recvValues2;
		mpi.setProcsFaceDatasDouble(
					sendValues0, recvValues0,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					sendValues1, recvValues1,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					sendValues2, recvValues2,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		int num=0;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){ 
				double varF[3];
				varF[0] = face.wC*surfNormalVec[face.owner][0]+(1.0-face.wC)*recvValues0[num];
				varF[1] = face.wC*surfNormalVec[face.owner][1]+(1.0-face.wC)*recvValues1[num];
				varF[2] = face.wC*surfNormalVec[face.owner][2]+(1.0-face.wC)*recvValues2[num];
				
				for(int j=0; j<3; ++j){
					kappa[face.owner] += 
						varF[j]*face.unitNormals[j]*face.area/mesh.cells[face.owner].volume;
				}
				++num;
			}
		}
		
	}
	
	
	
}