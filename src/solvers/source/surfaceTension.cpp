
#include "../build.h"
#include <cmath>
#include <array>
#include <numeric>

double SEMO_Solvers_Builder::calcSourceSurfaceTension(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	SEMO_Utility_Math math;
	
	
	
	// bool boolSkewnessCorrection = true;
	// bool boolSkewnessCorrection = false;
	int gradIterMax_LG = 1;
	int gradIterMax_GG = 1;
	
	
	// ===============================
	vector<double> kappa;
	this->calcCurvature(mesh, controls.VF[0], kappa);
	
	
	vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(9,0.0));
	for(int iter=0; iter<gradIterMax_LG; ++iter){
		// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
		vector<double> dummyVec;
		math.calcLeastSquare(mesh, "cellVertex", "2nd", "cell", 
			controls.VF[0], controls.fVF[0], dummyVec, gradAi);
	}
	for(int iter=0; iter<gradIterMax_GG; ++iter){
		math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	}
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		// cell.var[controls.sourceSurfaceTension[0]] = 
			// (species[0].sigma * 100.0 * gradAi[i][0]);
		// cell.var[controls.sourceSurfaceTension[1]] = 
			// (species[0].sigma * 100.0 * gradAi[i][1]);
		// cell.var[controls.sourceSurfaceTension[2]] = 
			// (species[0].sigma * 100.0 * gradAi[i][2]);
		
		cell.var[controls.sourceSurfaceTension[0]] = 
			(species[0].sigma * kappa[i] * gradAi[i][0]);
		cell.var[controls.sourceSurfaceTension[1]] = 
			(species[0].sigma * kappa[i] * gradAi[i][1]);
		cell.var[controls.sourceSurfaceTension[2]] = 
			(species[0].sigma * kappa[i] * gradAi[i][2]);
		cell.var[controls.sourceSurfaceTension[3]] = 
			species[0].sigma * kappa[i] * 
			(gradAi[i][0] * cell.var[controls.U] +
			 gradAi[i][1] * cell.var[controls.V] +
			 gradAi[i][2] * cell.var[controls.W]);
		
	}
	
	
	
	
	// mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int iter=0; iter<gradIterMax; ++iter){
		// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], mesh.cellsGradientVar[controls.VF[0]]);
	// }
	// if(size>1){
		// mpi.sendRecvTemporaryCellData(mesh, controls.VF[0], mesh.cellsProcVar[controls.VF[0]]);
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.VF[0]], mesh.cellsProcGradientVar[controls.VF[0]]);
		
	// }
	
	
	// vector<vector<double>> gradAi(mesh.cells.size(),vector<double>(3,0.0));
	

	// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		// double wCL = face.wC;
		// double wCR = 1.0-wCL;
	
		// double VFL = mesh.cells[face.owner].var[controls.VF[0]];
		// vector<double> gradVFL(3,0.0);
		
		// double VFR = 0.0;
		// vector<double> gradVFR(3,0.0);
		
		// for(int ii=0; ii<3; ++ii){
			// gradVFL[ii] = mesh.cellsGradientVar[controls.VF[0]][face.owner][ii];
		// }
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// VFR = mesh.cells[face.neighbour].var[controls.VF[0]];
			// for(int ii=0; ii<3; ++ii){
				// gradVFR[ii] = mesh.cellsGradientVar[controls.VF[0]][face.neighbour][ii];
			// }
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// VFR = mesh.cellsProcVar[controls.VF[0]][ip];
			// for(int ii=0; ii<3; ++ii){
				// gradVFR[ii] = mesh.cellsProcGradientVar[controls.VF[0]][ip][ii];
			// }
			
		// }
		
		// double VFF = wCL*VFL + wCR*VFR; 
		// if(boolSkewnessCorrection){
			// for(int ii=0; ii<3; ++ii){
				// VFF += wCL*gradVFL[ii]*face.vecSkewness[ii];
				// VFF += wCR*gradVFR[ii]*face.vecSkewness[ii];
			// }
			// VFF = max(0.0,min(1.0,VFF));
		// }
		
		// gradAi[face.owner][0] += (VFF * face.area * face.unitNormals[0]);
		// gradAi[face.owner][1] += (VFF * face.area * face.unitNormals[1]);
		// gradAi[face.owner][2] += (VFF * face.area * face.unitNormals[2]);
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// gradAi[face.neighbour][0] -= (VFF * face.area * face.unitNormals[0]);
			// gradAi[face.neighbour][1] -= (VFF * face.area * face.unitNormals[1]);
			// gradAi[face.neighbour][2] -= (VFF * face.area * face.unitNormals[2]);
		// }
		
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++ip;
		
		
	// }
	

	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// for(int i=str; i<end; ++i){
				// auto& face = mesh.faces[i];
				
				// double area = face.area;
				
				// double VFF = 0.5*face.varL[controls.fVF[0]] + 0.5*face.varR[controls.fVF[0]];
			
				// // convective term
				// gradAi[face.owner][0] += ( VFF * face.area * face.unitNormals[0]);
				// gradAi[face.owner][1] += ( VFF * face.area * face.unitNormals[1]);
				// gradAi[face.owner][2] += ( VFF * face.area * face.unitNormals[2]);
			// }
			
			
		// }
		
	// }
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// // gradAi[i][0] /= cell.volume;
		// // gradAi[i][1] /= cell.volume;
		// // gradAi[i][2] /= cell.volume;
		
		// cell.var[controls.sourceSurfaceTension[0]] = 
			// (species[0].sigma * kappa[i] * gradAi[i][0]);
		// cell.var[controls.sourceSurfaceTension[1]] = 
			// (species[0].sigma * kappa[i] * gradAi[i][1]);
		// cell.var[controls.sourceSurfaceTension[2]] = 
			// (species[0].sigma * kappa[i] * gradAi[i][2]);
		// cell.var[controls.sourceSurfaceTension[3]] = 
			// species[0].sigma * kappa[i] * 
			// (gradAi[i][0] * cell.var[controls.U] +
			 // gradAi[i][1] * cell.var[controls.V] +
			 // gradAi[i][2] * cell.var[controls.W]);
			 
			 
		// // cell.var[controls.U] = (gradAi[i][0]);
		// // cell.var[controls.V] = (gradAi[i][1]);
		// // cell.var[controls.W] = (gradAi[i][2]);
		// // cell.var[controls.T] = (kappa[i]);
		
	// }
	
	
	// ===============================
	
	
	return 0.0;
	
	
}