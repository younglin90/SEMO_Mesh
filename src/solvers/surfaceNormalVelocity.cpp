
#include "build.h"
#include <cmath>
#include <array>
#include <numeric>


// void SEMO_Solvers_Builder::calcSurfaceNormalVelocity(
	// double& UnF, 
	// vector<double>& gradUL, vector<double>& gradUR, 
	// vector<double>& gradVL, vector<double>& gradVR, 
	// vector<double>& gradWL, vector<double>& gradWR,
	// vector<double>& gradPL, vector<double>& gradPR, 
	// vector<double>& srcGRV_L, vector<double>& srcGRV_R, 
	// vector<double>& srcSFT_L, vector<double>& srcSFT_R,
	// double& orgPL, double& orgPR, 
	// double& RhoL, double& RhoR, 
	// double& wCL, double& wCR, 
	// double& timeStep,
	// vector<double>& gravityAcceleration,
	// double& sigma, 
	// double& kappaL, double& kappaR,
	// double& orgVFL, double& orgVFR,
	// double& alpha, double& dPN,
	// vector<double>& nvec,
	// vector<double>& vecSkewness,
	// vector<double>& unitNomalsPN
	// ){

	// double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
	// double tmp1 = timeStep / Rho_star;
	
	// // skewness correction
	// {
		// for(int ii=0; ii<3; ++ii){
			// UnF += wCL * gradUL[ii]*vecSkewness[ii] * nvec[0];
			// UnF += wCR * gradUR[ii]*vecSkewness[ii] * nvec[0];
			
			// UnF += wCL * gradVL[ii]*vecSkewness[ii] * nvec[1];
			// UnF += wCR * gradVR[ii]*vecSkewness[ii] * nvec[1];
			
			// UnF += wCL * gradWL[ii]*vecSkewness[ii] * nvec[2];
			// UnF += wCR * gradWR[ii]*vecSkewness[ii] * nvec[2];
		// }
	// }
	
	// // pressure correction (cell->face) 
	// {
		// UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
		// // non-orthogonal
		// for(int ii=0; ii<3; ++ii){
			// UnF -= tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			// UnF -= tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		// }
		// for(int ii=0; ii<3; ++ii){
			// UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			// UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
		// }
	// }
	
	
	// // gravity correction (cell->face)
	// {
		// for(int ii=0; ii<3; ++ii){
			// UnF += alpha * tmp1 * Rho_star * gravityAcceleration[ii] * alpha*unitNomalsPN[ii];		
		// }
		// // non-orthogonal
		// for(int ii=0; ii<3; ++ii){
			// UnF += tmp1 * Rho_star * wCL/RhoL*srcGRV_L[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			// UnF += tmp1 * Rho_star * wCR/RhoR*srcGRV_R[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		// }
		// for(int ii=0; ii<3; ++ii){
			// UnF -= alpha * tmp1 * Rho_star * wCL/RhoL* srcGRV_L[ii] * nvec[ii];
			// UnF -= alpha * tmp1 * Rho_star * wCR/RhoR* srcGRV_R[ii] * nvec[ii];
		// }
	// }

		
	// // surface tension correction (cell->face)
	// {
		// UnF += alpha * tmp1 * sigma * (wCL*kappaL + wCR*kappaR) * (orgVFR-orgVFL)/dPN;
		// // non-orthogonal
		// for(int ii=0; ii<3; ++ii){
			// UnF += tmp1 * Rho_star * wCL/RhoL*srcSFT_L[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			// UnF += tmp1 * Rho_star * wCR/RhoR*srcSFT_R[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		// }
		// for(int ii=0; ii<3; ++ii){
			// UnF -= alpha * tmp1 * Rho_star * wCL/RhoL* srcSFT_L[ii] * nvec[ii];
			// UnF -= alpha * tmp1 * Rho_star * wCR/RhoR* srcSFT_R[ii] * nvec[ii];
		// }
	// }
// }







void SEMO_Solvers_Builder::calcInterpolVelSkewness(
	double& UnF, 
	vector<double>& gradUL, vector<double>& gradUR, 
	vector<double>& gradVL, vector<double>& gradVR, 
	vector<double>& gradWL, vector<double>& gradWR,
	double& wCL, double& wCR, 
	vector<double>& nvec,
	vector<double>& vecSkewness
	){

	// skewness correction
	{
		for(int ii=0; ii<3; ++ii){
			UnF += wCL * gradUL[ii]*vecSkewness[ii] * nvec[0];
			UnF += wCR * gradUR[ii]*vecSkewness[ii] * nvec[0];
			
			UnF += wCL * gradVL[ii]*vecSkewness[ii] * nvec[1];
			UnF += wCR * gradVR[ii]*vecSkewness[ii] * nvec[1];
			
			UnF += wCL * gradWL[ii]*vecSkewness[ii] * nvec[2];
			UnF += wCR * gradWR[ii]*vecSkewness[ii] * nvec[2];
		}
	}
	
}
void SEMO_Solvers_Builder::calcInterpolVelPressure(
	double& UnF, 
	vector<double>& gradPL, vector<double>& gradPR, 
	double& orgPL, double& orgPR, 
	double& RhoL, double& RhoR, 
	double& wCL, double& wCR, 
	double& timeStep,
	double& alpha, double& dPN,
	vector<double>& nvec,
	vector<double>& unitNomalsPN
	){

	double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
	double tmp1 = timeStep / Rho_star;
	
	// pressure correction (cell->face) 
	{
		UnF -= alpha * tmp1*(orgPR-orgPL)/dPN;
		// non-orthogonal
		for(int ii=0; ii<3; ++ii){
			UnF -= tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			UnF -= tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		}
		for(int ii=0; ii<3; ++ii){
			UnF += tmp1 * Rho_star * wCL/RhoL*gradPL[ii]*nvec[ii];
			UnF += tmp1 * Rho_star * wCR/RhoR*gradPR[ii]*nvec[ii];
		}
	}
	
	
}
void SEMO_Solvers_Builder::calcInterpolVelGravity(
	double& UnF, 
	vector<double>& srcGRV_L, vector<double>& srcGRV_R, 
	double& RhoL, double& RhoR, 
	double& wCL, double& wCR, 
	double& timeStep,
	vector<double>& gravityAcceleration,
	double& alpha, 
	vector<double>& nvec,
	vector<double>& unitNomalsPN
	){

	double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
	double tmp1 = timeStep / Rho_star;
	
	
	// gravity correction (cell->face)
	{
		for(int ii=0; ii<3; ++ii){
			UnF += alpha * tmp1 * Rho_star * gravityAcceleration[ii] * alpha*unitNomalsPN[ii];		
		}
		// non-orthogonal
		for(int ii=0; ii<3; ++ii){
			UnF += tmp1 * Rho_star * wCL/RhoL*srcGRV_L[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			UnF += tmp1 * Rho_star * wCR/RhoR*srcGRV_R[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		}
		for(int ii=0; ii<3; ++ii){
			UnF -= alpha * tmp1 * Rho_star * wCL/RhoL* srcGRV_L[ii] * nvec[ii];
			UnF -= alpha * tmp1 * Rho_star * wCR/RhoR* srcGRV_R[ii] * nvec[ii];
		}
	}
}


void SEMO_Solvers_Builder::calcInterpolVelSurfTens(
	double& UnF, 
	vector<double>& srcSFT_L, vector<double>& srcSFT_R,
	double& RhoL, double& RhoR, 
	double& wCL, double& wCR, 
	double& timeStep,
	double& sigma, 
	double& kappaL, double& kappaR,
	double& orgVFL, double& orgVFR,
	double& alpha, double& dPN,
	vector<double>& nvec,
	vector<double>& unitNomalsPN
	){

	double Rho_star = 1.0 / (wCL/RhoL + wCR/RhoR);
	double tmp1 = timeStep / Rho_star;
	
		
	// surface tension correction (cell->face)
	{
		UnF += alpha * tmp1 * sigma * (wCL*kappaL + wCR*kappaR) * (orgVFR-orgVFL)/dPN;
		// non-orthogonal
		for(int ii=0; ii<3; ++ii){
			UnF += tmp1 * Rho_star * wCL/RhoL*srcSFT_L[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);
			UnF += tmp1 * Rho_star * wCR/RhoR*srcSFT_R[ii]*(nvec[ii] - alpha*unitNomalsPN[ii]);		
		}
		for(int ii=0; ii<3; ++ii){
			UnF -= alpha * tmp1 * Rho_star * wCL/RhoL* srcSFT_L[ii] * nvec[ii];
			UnF -= alpha * tmp1 * Rho_star * wCR/RhoR* srcSFT_R[ii] * nvec[ii];
		}
	}
}



void SEMO_Solvers_Builder::calcInterpolVelUnsteady(
	double& UnF, 
	double& wCL, double& wCR, 
	double& UnF_old, 
	double& UL_old, double& UR_old,
	double& VL_old, double& VR_old,
	double& WL_old, double& WR_old,
	vector<double>& nvec
	){

	// unsteady correction (cell->face)
	{
		UnF += UnF_old;
		UnF -= (wCL*UL_old + wCR*UR_old)*nvec[0];
		UnF -= (wCL*VL_old + wCR*VR_old)*nvec[1];
		UnF -= (wCL*WL_old + wCR*WR_old)*nvec[2];
	}
}




void SEMO_Solvers_Builder::saveSurfaceNormalVelocity(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){


    // int rank = MPI::COMM_WORLD.Get_rank(); 
    // int size = MPI::COMM_WORLD.Get_size();

	// SEMO_MPI_Builder mpi;
	
	// bool boolSkewnessCorrection = true;
	// // bool boolSkewnessCorrection = false;
	// int gradIterMax = 2;
	
	// // linear solver
	
	// vector<double> A_p_vals(mesh.cells.size(), 0.0);
	
	
	// SEMO_Utility_Math math;
	
	// setCellVarMinMax(mesh, controls.U, controls.maximumU, controls.minimumU);
	// setCellVarMinMax(mesh, controls.V, controls.maximumV, controls.minimumV);
	// setCellVarMinMax(mesh, controls.W, controls.maximumW, controls.minimumW);
	
	
	
	// mesh.cellsGradientVar[controls.P].resize(mesh.cells.size(),vector<double>(3,0.0));
	// mesh.cellsGradientVar[controls.U].resize(mesh.cells.size(),vector<double>(3,0.0));
	// mesh.cellsGradientVar[controls.V].resize(mesh.cells.size(),vector<double>(3,0.0));
	// mesh.cellsGradientVar[controls.W].resize(mesh.cells.size(),vector<double>(3,0.0));
	// mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int iter=0; iter<gradIterMax; ++iter){
		// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, mesh.cellsGradientVar[controls.P]);
		// math.calcLeastSquare2nd(mesh, controls.U, controls.fU, mesh.cellsGradientVar[controls.U]);
		// math.calcLeastSquare2nd(mesh, controls.V, controls.fV, mesh.cellsGradientVar[controls.V]);
		// math.calcLeastSquare2nd(mesh, controls.W, controls.fW, mesh.cellsGradientVar[controls.W]);
		// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], mesh.cellsGradientVar[controls.VF[0]]);
	// }
	
	// if(size>1){
		// mpi.sendRecvTemporaryCellData(mesh, controls.maximumU, mesh.cellsProcVar[controls.maximumU]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.maximumV, mesh.cellsProcVar[controls.maximumV]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.maximumW, mesh.cellsProcVar[controls.maximumW]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.minimumU, mesh.cellsProcVar[controls.minimumU]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.minimumV, mesh.cellsProcVar[controls.minimumV]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.minimumW, mesh.cellsProcVar[controls.minimumW]);
		
		// mpi.sendRecvTemporaryCellData(mesh, controls.P, mesh.cellsProcVar[controls.P]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.U, mesh.cellsProcVar[controls.U]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.V, mesh.cellsProcVar[controls.V]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.W, mesh.cellsProcVar[controls.W]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.Rho, mesh.cellsProcVar[controls.Rho]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.VF[0], mesh.cellsProcVar[controls.VF[0]]);
		
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[0], mesh.cellsProcVar[controls.sourceGravity[0]]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[1], mesh.cellsProcVar[controls.sourceGravity[1]]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceGravity[2], mesh.cellsProcVar[controls.sourceGravity[2]]);
		
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[0], mesh.cellsProcVar[controls.sourceSurfaceTension[0]]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[1], mesh.cellsProcVar[controls.sourceSurfaceTension[1]]);
		// mpi.sendRecvTemporaryCellData(mesh, controls.sourceSurfaceTension[2], mesh.cellsProcVar[controls.sourceSurfaceTension[2]]);
		
		// mpi.sendRecvTemporaryCellData(mesh, controls.kappa, mesh.cellsProcVar[controls.kappa]);
		
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.P], mesh.cellsProcGradientVar[controls.P]);
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.U], mesh.cellsProcGradientVar[controls.U]);
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.V], mesh.cellsProcGradientVar[controls.V]);
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.W], mesh.cellsProcGradientVar[controls.W]);
		// mpi.sendRecvTemporaryVectorData(mesh, mesh.cellsGradientVar[controls.VF[0]], mesh.cellsProcGradientVar[controls.VF[0]]);
		
	// }

	
	// for(int i=0, ip=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		// vector<double> nvec(3,0.0);
		// nvec[0] = face.unitNormals[0];
		// nvec[1] = face.unitNormals[1];
		// nvec[2] = face.unitNormals[2];
		
		// double wCL = face.wC;
		// double wCR = 1.0-wCL;
		
		// double RhoL = mesh.cells[face.owner].var[controls.Rho];
		// double PL = mesh.cells[face.owner].var[controls.P];
		// double UL = mesh.cells[face.owner].var[controls.U];
		// double VL = mesh.cells[face.owner].var[controls.V];
		// double WL = mesh.cells[face.owner].var[controls.W];
		// double VFL = mesh.cells[face.owner].var[controls.VF[0]];
		// double kappaL = mesh.cells[face.owner].var[controls.kappa];
		// double maximumUL = mesh.cells[face.owner].var[controls.maximumU];
		// double minimumUL = mesh.cells[face.owner].var[controls.minimumU];
		// double maximumVL = mesh.cells[face.owner].var[controls.maximumV];
		// double minimumVL = mesh.cells[face.owner].var[controls.minimumV];
		// double maximumWL = mesh.cells[face.owner].var[controls.maximumW];
		// double minimumWL = mesh.cells[face.owner].var[controls.minimumW];
		// vector<double> gradPL(3,0.0);
		// vector<double> gradUL(3,0.0);
		// vector<double> gradVL(3,0.0);
		// vector<double> gradWL(3,0.0);
		// vector<double> gradVFL(3,0.0);
		// vector<double> srcGRV_L(3,0.0);
		// vector<double> srcSFT_L(3,0.0);
		
		// double RhoR = 0.0;
		// double PR = 0.0;
		// double UR = 0.0;
		// double VR = 0.0;
		// double WR = 0.0;
		// double VFR = 0.0;
		// double kappaR = 0.0;
		// double maximumUR = 0.0;
		// double minimumUR = 0.0;
		// double maximumVR = 0.0;
		// double minimumVR = 0.0;
		// double maximumWR = 0.0;
		// double minimumWR = 0.0;
		// vector<double> gradPR(3,0.0);
		// vector<double> gradUR(3,0.0);
		// vector<double> gradVR(3,0.0);
		// vector<double> gradWR(3,0.0);
		// vector<double> gradVFR(3,0.0);
		// vector<double> srcGRV_R(3,0.0);
		// vector<double> srcSFT_R(3,0.0);
		
		// for(int ii=0; ii<3; ++ii){
			// gradPL[ii] = mesh.cellsGradientVar[controls.P][face.owner][ii];
			// gradUL[ii] = mesh.cellsGradientVar[controls.U][face.owner][ii];
			// gradVL[ii] = mesh.cellsGradientVar[controls.V][face.owner][ii];
			// gradWL[ii] = mesh.cellsGradientVar[controls.W][face.owner][ii];
			// gradVFL[ii] = mesh.cellsGradientVar[controls.VF[0]][face.owner][ii];
			// srcGRV_L[ii] = mesh.cells[face.owner].var[controls.sourceGravity[ii]];
			// srcSFT_L[ii] = mesh.cells[face.owner].var[controls.sourceSurfaceTension[ii]];
		// }
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// RhoR = mesh.cells[face.neighbour].var[controls.Rho];
			// PR = mesh.cells[face.neighbour].var[controls.P];
			// UR = mesh.cells[face.neighbour].var[controls.U];
			// VR = mesh.cells[face.neighbour].var[controls.V];
			// WR = mesh.cells[face.neighbour].var[controls.W];
			// VFR = mesh.cells[face.neighbour].var[controls.VF[0]];
			// kappaR = mesh.cells[face.neighbour].var[controls.kappa];
			// maximumUR = mesh.cells[face.neighbour].var[controls.maximumU];
			// minimumUR = mesh.cells[face.neighbour].var[controls.minimumU];
			// maximumVR = mesh.cells[face.neighbour].var[controls.maximumV];
			// minimumVR = mesh.cells[face.neighbour].var[controls.minimumV];
			// maximumWR = mesh.cells[face.neighbour].var[controls.maximumW];
			// minimumWR = mesh.cells[face.neighbour].var[controls.minimumW];
			// for(int ii=0; ii<3; ++ii){
				// gradPR[ii] = mesh.cellsGradientVar[controls.P][face.neighbour][ii];
				// gradUR[ii] = mesh.cellsGradientVar[controls.U][face.neighbour][ii];
				// gradVR[ii] = mesh.cellsGradientVar[controls.V][face.neighbour][ii];
				// gradWR[ii] = mesh.cellsGradientVar[controls.W][face.neighbour][ii];
				// gradVFR[ii] = mesh.cellsGradientVar[controls.VF[0]][face.neighbour][ii];
				// srcGRV_R[ii] = mesh.cells[face.neighbour].var[controls.sourceGravity[ii]];
				// srcSFT_R[ii] = mesh.cells[face.neighbour].var[controls.sourceSurfaceTension[ii]];
			// }
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// RhoR = mesh.cellsProcVar[controls.Rho][ip];
			// PR = mesh.cellsProcVar[controls.P][ip];
			// UR = mesh.cellsProcVar[controls.U][ip];
			// VR = mesh.cellsProcVar[controls.V][ip];
			// WR = mesh.cellsProcVar[controls.W][ip];
			// VFR = mesh.cellsProcVar[controls.VF[0]][ip];
			// maximumUR = mesh.cellsProcVar[controls.maximumU][ip];
			// minimumUR = mesh.cellsProcVar[controls.minimumU][ip];
			// maximumVR = mesh.cellsProcVar[controls.maximumV][ip];
			// minimumVR = mesh.cellsProcVar[controls.minimumV][ip];
			// maximumWR = mesh.cellsProcVar[controls.maximumW][ip];
			// minimumWR = mesh.cellsProcVar[controls.minimumW][ip];
			// kappaR = mesh.cellsProcVar[controls.kappa][ip];
			// for(int ii=0; ii<3; ++ii){
				// gradPR[ii] = mesh.cellsProcGradientVar[controls.P][ip][ii];
				// gradUR[ii] = mesh.cellsProcGradientVar[controls.U][ip][ii];
				// gradVR[ii] = mesh.cellsProcGradientVar[controls.V][ip][ii];
				// gradWR[ii] = mesh.cellsProcGradientVar[controls.W][ip][ii];
				// gradVFR[ii] = mesh.cellsProcGradientVar[controls.VF[0]][ip][ii];
				// srcGRV_R[ii] = mesh.cellsProcVar[controls.sourceGravity[ii]][ip];
				// srcSFT_R[ii] = mesh.cellsProcVar[controls.sourceSurfaceTension[ii]][ip];
			// }
			
		// }
		
		// double Rho_star = 1.0 / (wCR/RhoL + wCR/RhoR);
		// double tmp1 = controls.timeStep / Rho_star;
		
		// double dPN = face.magPN;
		// double alpha = face.alphaF;
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		// double UnF = wCL*UnL+wCR*UnR;
			
		// if(boolSkewnessCorrection){
			// // calcInterpolVelSkewness(UnF, gradUL, gradUR, gradVL, gradVR, gradWL, gradWR,
				// // wCL, wCR, 
				// // nvec,
				// // face.vecSkewness);
			// double UL_skew = UL;
			// double VL_skew = VL;
			// double WL_skew = WL;
			// double UR_skew = UR;
			// double VR_skew = VR;
			// double WR_skew = WR;
			// for(int ii=0; ii<3; ++ii){
				// UL_skew += gradUL[ii]*face.vecSkewness[ii];
				// UR_skew += gradUR[ii]*face.vecSkewness[ii];
				// VL_skew += gradVL[ii]*face.vecSkewness[ii];
				// VR_skew += gradVR[ii]*face.vecSkewness[ii];
				// WL_skew += gradWL[ii]*face.vecSkewness[ii];
				// WR_skew += gradWR[ii]*face.vecSkewness[ii];
			// }
			// UL_skew = max(minimumUL,min(maximumUL,UL_skew));
			// UR_skew = max(minimumUR,min(maximumUR,UR_skew));
			// VL_skew = max(minimumVL,min(maximumVL,VL_skew));
			// VR_skew = max(minimumVR,min(maximumVR,VR_skew));
			// WL_skew = max(minimumWL,min(maximumWL,WL_skew));
			// WR_skew = max(minimumWR,min(maximumWR,WR_skew));
			
			// UnL = UL_skew*nvec[0] + VL_skew*nvec[1] + WL_skew*nvec[2];
			// UnR = UR_skew*nvec[0] + VR_skew*nvec[1] + WR_skew*nvec[2];
			// UnF = wCL*UnL+wCR*UnR;
		// }
		// calcInterpolVelPressure(UnF, 
			// gradPL, gradPR, 
			// PL, PR, RhoL, RhoR, wCL, wCR, 
			// controls.timeStep,
			// alpha, dPN,
			// nvec,
			// face.unitNomalsPN);
		// calcInterpolVelGravity(UnF, srcGRV_L, srcGRV_R,
			// RhoL, RhoR, wCL, wCR, 
			// controls.timeStep,
			// controls.gravityAcceleration,
			// alpha,
			// nvec,
			// face.unitNomalsPN);
		// calcInterpolVelSurfTens(UnF, srcSFT_L, srcSFT_R,
			// RhoL, RhoR, wCL, wCR, 
			// controls.timeStep,
			// species[0].sigma, kappaL, kappaR,
			// VFL, VFR,
			// alpha, dPN,
			// nvec,
			// face.unitNomalsPN);
			
		// face.var[controls.Un] = UnF;
		
	// }
}