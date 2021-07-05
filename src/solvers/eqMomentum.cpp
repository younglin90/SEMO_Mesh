#include "build.h"
#include <cmath>
#include <array>
#include <numeric>




void SEMO_Solvers_Builder::calcMomentumEqs(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	vector<double>& linAD,
	vector<double>& linAOD,
	vector<double>& linAFL,
	vector<double>& linAFR,
	vector<double>& residuals){
		
		
	
	double coeffPf = 0.05;
	// double coeffPf = 0.0;
	double coeffDiss = 0.9;
	
		
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
		
	SEMO_MPI_Builder mpi;


	SEMO_Utility_Math math;
	
	// gradient P
	vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	// math.calcGaussGreen(mesh, controls.P, controls.fP, gradP);
	// math.calcGGLSQ(mesh, controls.P, controls.fP, gradP);
	// math.calcMGG(mesh, controls.P, controls.fP, 50, 1.e-8, gradP);
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	
	// vector<double> limGradP;
	// math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] *= limGradP[i]; 
		// }
	// }
	
	
	
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	// // math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// if(size>1){
		// // processor faces
		// // gradP , 
		// vector<double> gradPx_send;
		// vector<double> gradPy_send;
		// vector<double> gradPz_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gradPx_send.push_back(gradP[face.owner][0]);
				// gradPy_send.push_back(gradP[face.owner][1]);
				// gradPz_send.push_back(gradP[face.owner][2]);
			// }
		// }
		// // SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatasDouble(
					// gradPx_send, gradPx_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatasDouble(
					// gradPy_send, gradPy_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatasDouble(
					// gradPz_send, gradPz_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	// }
	

	
	// gradient U, V, W
	vector<vector<double>> gradU;
	vector<vector<double>> gradV;
	vector<vector<double>> gradW;
	// math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// math.calcGGLSQ(mesh, controls.U, controls.fU, gradU);
	// math.calcGGLSQ(mesh, controls.V, controls.fV, gradV);
	// math.calcGGLSQ(mesh, controls.W, controls.fW, gradW);
	// math.calcMGG(mesh, controls.U, controls.fU, 10, 1.e-8, gradU);
	// math.calcMGG(mesh, controls.V, controls.fV, 10, 1.e-8, gradV);
	// math.calcMGG(mesh, controls.W, controls.fW, 10, 1.e-8, gradW);
	math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	
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
		
		mpi.setProcsFaceDatasDouble(
					gradUx_send, gradUx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradUy_send, gradUy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradUz_send, gradUz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradUx_send.clear();
		gradUy_send.clear();
		gradUz_send.clear();
					
		mpi.setProcsFaceDatasDouble(
					gradVx_send, gradVx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradVy_send, gradVy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradVz_send, gradVz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradVx_send.clear();
		gradVy_send.clear();
		gradVz_send.clear();
					
		mpi.setProcsFaceDatasDouble(
					gradWx_send, gradWx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradWy_send, gradWy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradWz_send, gradWz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		gradWx_send.clear();
		gradWy_send.clear();
		gradWz_send.clear();
	}
	
	
		
	vector<double> linA;
	vector<double> linB0;
	vector<double> linB1;
	vector<double> linB2;
	
	
	
	vector<vector<double>> gradAi;
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// math.calcGGLSQ(mesh, controls.VF[0], controls.fVF[0], gradAi);
	// math.calcMGG(mesh, controls.VF[0], controls.fVF[0], 1, 1.e-8, gradAi);
	// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradAi);
	
	vector<double> kappa;
	this->calcCurvature(mesh, controls.VF[0], kappa);
	
	
	
	
	
	
	
	//=========================
	// Un initialization
	if( controls.iterPBs == 0 && controls.iterMom == 0){
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];

			vector<double> nvec;
			nvec.push_back(face.unitNormals[0]);
			nvec.push_back(face.unitNormals[1]);
			nvec.push_back(face.unitNormals[2]);
			
			double wCL = face.wC;
			double wCR = 1.0-face.wC;
		
			double UL = face.varL[controls.fU];
			double VL = face.varL[controls.fV];
			double WL = face.varL[controls.fW];
			
			double UR = face.varR[controls.fU];
			double VR = face.varR[controls.fV];
			double WR = face.varR[controls.fW];
			
			double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
			double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
			
			face.var[controls.Un] = wCL*UnL+wCR*UnR;
			
		}
	}
	//=========================
	
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
		// time marching -> first order euler
		// linA.push_back(cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		// linB0.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.U] - 1.0*cell.var[controls.oldU]
			// )*cell.volume/controls.timeStep);
			
		// linB1.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.V] - 1.0*cell.var[controls.oldV]
			// )*cell.volume/controls.timeStep);
			
		// linB2.push_back(-
			// cell.var[controls.Rho]*(
			// 1.0*cell.var[controls.W] - 1.0*cell.var[controls.oldW]
			// )*cell.volume/controls.timeStep);
		
		// time marching -> second order upwind euler
		linA.push_back(1.5*cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		linB0.push_back(-
			cell.var[controls.Rho]*(
			1.5*cell.var[controls.U] - 2.0*cell.var[controls.oldU] + 0.5*cell.var[controls.Qm[0]]
			)*cell.volume/controls.timeStep);
			
		linB1.push_back(-
			cell.var[controls.Rho]*(
			1.5*cell.var[controls.V] - 2.0*cell.var[controls.oldV] + 0.5*cell.var[controls.Qm[1]]
			)*cell.volume/controls.timeStep);
			
		linB2.push_back(-
			cell.var[controls.Rho]*(
			1.5*cell.var[controls.W] - 2.0*cell.var[controls.oldW] + 0.5*cell.var[controls.Qm[2]]
			)*cell.volume/controls.timeStep);
			
		// gravity force terms
		linB0.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[0];
		linB1.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[1];
		linB2.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[2];
			
		// surface tension force terms
		linB0.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][0]);
		linB1.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][1]);
		linB2.back() += cell.volume*(-species[0].sigma * kappa[i] * gradAi[i][2]);
		
		
	}
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	vector<double> linAL;
	vector<double> linAR;

	int proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		// face.area = 0.00025;
		
		double area = face.area;
		vector<double> nvec;
		nvec.push_back(face.unitNormals[0]);
		nvec.push_back(face.unitNormals[1]);
		nvec.push_back(face.unitNormals[2]);
		
		vector<double> distanceCells;
		distanceCells.push_back(face.distCells[0]);
		distanceCells.push_back(face.distCells[1]);
		distanceCells.push_back(face.distCells[2]);
		
		double wCL = face.wC;
		double wCR = 1.0-face.wC;
	
		double UL = face.varL[controls.fU];
		double VL = face.varL[controls.fV];
		double WL = face.varL[controls.fW];
		double PL = face.varL[controls.fP];
		double RhoL = face.varL[controls.fRho];
		
		double UR = face.varR[controls.fU];
		double VR = face.varR[controls.fV];
		double WR = face.varR[controls.fW];
		double PR = face.varR[controls.fP];
		double RhoR = face.varR[controls.fRho];
		
		double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		
		double UnF = wCL*UnL+wCR*UnR;
		// double UnF = 0.5*UnL+0.5*UnR;
		// double UnF = face.var[controls.Un];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			UnF = face.var[controls.Un];
		}
		

		double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  pow(distanceCells[1],2.0) + 
						  pow(distanceCells[2],2.0));
		// double dPNEff = dPN 
			// /( (nvec[0]*distanceCells[0] + 
			    // nvec[1]*distanceCells[1] + 
				// nvec[2]*distanceCells[2])/dPN );

		double nonOrtholimiter = 1.0;
		double cosAlpha = dPN_e/dPN;
		if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			nonOrtholimiter = 0.5;
		}
		else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			nonOrtholimiter = 0.333;
		}
		else if( cosAlpha < 0.342 ){
			nonOrtholimiter = 0.0;
		}
		// nonOrtholimiter = 0.0;
				
				
		double Ef[3];
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
				
		// over-relaxed approach
		// double overRelaxCoeff = 1.0;
		// double overRelaxCoeff = Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2];
		// double overRelaxCoeff = 1.0 / (Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2]);
		
		// non-orthogonal, over-relaxed approach
		// double Tf[3];
		// Tf[0] = nvec[0] - Ef[0]*overRelaxCoeff;
		// Tf[1] = nvec[1] - Ef[1]*overRelaxCoeff;
		// Tf[2] = nvec[2] - Ef[2]*overRelaxCoeff;
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		// double dPNEff = dPN / overRelaxCoeff;
		// double dPNEff = dPN;
		
		
		double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		double UF = UL*weightL + UR*weightR;
		double VF = VL*weightL + VR*weightR;
		double WF = WL*weightL + WR*weightR;
		
		double PF = wCL*PL + wCR*PR;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			double barA1 = RhoR*mesh.cells[face.neighbour].volume/controls.timeStep;
			
			// PF = (barA0*PL + barA1*PR)/(barA0 + barA1);
			
			PF += coeffPf * 0.5*barA0*barA1/(barA0+barA1)*(UnL-UnR)/area;
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			double barA1 = RhoR*mesh.cells[face.owner].volume/controls.timeStep;
			
			// PF = (barA0*PL + barA1*PR)/(barA0 + barA1);
			
			PF += coeffPf * 0.5*barA0*barA1/(barA0+barA1)*(UnL-UnR)/area;
		}
		// else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			
			// double barA0 = RhoL*mesh.cells[face.owner].volume/controls.timeStep;
			
			// PF += coeffPf * 0.5*barA0*(UnL-UnR)/area;
		// }
		
		
		
		
		
		linAL.push_back(-weightL*RhoL*UnF*area);
		linAR.push_back(+weightR*RhoR*UnF*area);
		
		// properties, viscous term
		double muL = face.varL[controls.fmu];
		double muR = face.varR[controls.fmu];
		double muF = wCL*muL + wCR*muR;
		
		linAL.back() -= muF/dPN_e*area;
		linAR.back() -= muF/dPN_e*area;
		
			
		// pressure term
		linB0.at(face.owner) -= PF*nvec[0]*area;
		linB1.at(face.owner) -= PF*nvec[1]*area;
		linB2.at(face.owner) -= PF*nvec[2]*area;
		
		
		gradP[face.owner][0] += PF*nvec[0]*area;
		gradP[face.owner][1] += PF*nvec[1]*area;
		gradP[face.owner][2] += PF*nvec[2]*area;
		
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			linA.at(face.owner)     += (-linAL.back());
			linA.at(face.neighbour) += (-linAR.back());
			
			// convective term
			linB0.at(face.owner) -= UF*RhoL*UnF*area;
			linB1.at(face.owner) -= VF*RhoL*UnF*area;
			linB2.at(face.owner) -= WF*RhoL*UnF*area;
		
			linB0.at(face.neighbour) += UF*RhoR*UnF*area;
			linB1.at(face.neighbour) += VF*RhoR*UnF*area;
			linB2.at(face.neighbour) += WF*RhoR*UnF*area;
			
			// pressure term
			linB0.at(face.neighbour) += PF*nvec[0]*area;
			linB1.at(face.neighbour) += PF*nvec[1]*area;
			linB2.at(face.neighbour) += PF*nvec[2]*area;
			
			
			
			gradP[face.neighbour][0] -= PF*nvec[0]*area;
			gradP[face.neighbour][1] -= PF*nvec[1]*area;
			gradP[face.neighbour][2] -= PF*nvec[2]*area;
			
			
			// viscous term
			double orgUL = mesh.cells[face.owner].var[controls.U];
			double orgUR = mesh.cells[face.neighbour].var[controls.U];
			double orgVL = mesh.cells[face.owner].var[controls.V];
			double orgVR = mesh.cells[face.neighbour].var[controls.V];
			double orgWL = mesh.cells[face.owner].var[controls.W];
			double orgWR = mesh.cells[face.neighbour].var[controls.W];
			
			linB0.at(face.owner) += muF*(orgUR-orgUL)/dPN_e*area;
			linB1.at(face.owner) += muF*(orgVR-orgVL)/dPN_e*area;
			linB2.at(face.owner) += muF*(orgWR-orgWL)/dPN_e*area;
		
			linB0.at(face.neighbour) -= muF*(orgUR-orgUL)/dPN_e*area;
			linB1.at(face.neighbour) -= muF*(orgVR-orgVL)/dPN_e*area;
			linB2.at(face.neighbour) -= muF*(orgWR-orgWL)/dPN_e*area;
			
			// viscous term, properties
			double delPhiECF;
			
			double gradUf[3];
			gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradU[face.neighbour][0]);
			gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradU[face.neighbour][1]);
			gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradU[face.neighbour][2]);
			// delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
			
			// gradUf[0] += nvec[0]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// gradUf[1] += nvec[1]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// gradUf[2] += nvec[2]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			
			double gradVf[3];
			gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradV[face.neighbour][0]);
			gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradV[face.neighbour][1]);
			gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradV[face.neighbour][2]);
			// delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
			
			// gradVf[0] += nvec[0]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// gradVf[1] += nvec[1]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// gradVf[2] += nvec[2]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			
			double gradWf[3];
			gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradW[face.neighbour][0]);
			gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradW[face.neighbour][1]);
			gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradW[face.neighbour][2]);
			// delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
			
			// gradWf[0] += nvec[0]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// gradWf[1] += nvec[1]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// gradWf[2] += nvec[2]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			
			
			// viscous term, non-orthogonal
			double vNonOrthFlux0 = 
				nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			double vNonOrthFlux1 = 
				nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			double vNonOrthFlux2 = 
				nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
			linB0[face.owner] += vNonOrthFlux0;
			linB1[face.owner] += vNonOrthFlux1;
			linB2[face.owner] += vNonOrthFlux2;
				
			linB0[face.neighbour] -= vNonOrthFlux0;
			linB1[face.neighbour] -= vNonOrthFlux1;
			linB2[face.neighbour] -= vNonOrthFlux2;
			
			
			// double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);;
			// linAL.back() -= muF/dPN*area*tmpCoeff2;
			// linAR.back() -= muF/dPN*area*tmpCoeff2;
			
			// linA.at(face.owner)     += muF/dPN*area*tmpCoeff2;
			// linA.at(face.neighbour) += muF/dPN*area*tmpCoeff2;
			
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			linA.at(face.owner) += (-linAL.back());

			// convective term
			linB0.at(face.owner) -= UF*RhoL*UnF*area;
			linB1.at(face.owner) -= VF*RhoL*UnF*area;
			linB2.at(face.owner) -= WF*RhoL*UnF*area;
			

			// viscous term
			double orgUL = mesh.cells[face.owner].var[controls.U];
			double orgUR = UR;
			double orgVL = mesh.cells[face.owner].var[controls.V];
			double orgVR = VR;
			double orgWL = mesh.cells[face.owner].var[controls.W];
			double orgWR = WR;
			
			linB0.at(face.owner) += muF*(orgUR-orgUL)/dPN_e*area;
			linB1.at(face.owner) += muF*(orgVR-orgVL)/dPN_e*area;
			linB2.at(face.owner) += muF*(orgWR-orgWL)/dPN_e*area;
		
			// viscous term, properties
			double delPhiECF;
			
			double gradUf[3];
			gradUf[0] = (wCL*gradU[face.owner][0]+wCR*gradUx_recv[proc_num]);
			gradUf[1] = (wCL*gradU[face.owner][1]+wCR*gradUy_recv[proc_num]);
			gradUf[2] = (wCL*gradU[face.owner][2]+wCR*gradUz_recv[proc_num]);
			// delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
			
			// gradUf[0] += nvec[0]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// gradUf[1] += nvec[1]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			// gradUf[2] += nvec[2]*overRelaxCoeff*((orgUR-orgUL)/dPN - delPhiECF);
			
			double gradVf[3];
			gradVf[0] = (wCL*gradV[face.owner][0]+wCR*gradVx_recv[proc_num]);
			gradVf[1] = (wCL*gradV[face.owner][1]+wCR*gradVy_recv[proc_num]);
			gradVf[2] = (wCL*gradV[face.owner][2]+wCR*gradVz_recv[proc_num]);
			// delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
			
			// gradVf[0] += nvec[0]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// gradVf[1] += nvec[1]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			// gradVf[2] += nvec[2]*overRelaxCoeff*((orgVR-orgVL)/dPN - delPhiECF);
			
			double gradWf[3];
			gradWf[0] = (wCL*gradW[face.owner][0]+wCR*gradWx_recv[proc_num]);
			gradWf[1] = (wCL*gradW[face.owner][1]+wCR*gradWy_recv[proc_num]);
			gradWf[2] = (wCL*gradW[face.owner][2]+wCR*gradWz_recv[proc_num]);
			// delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
			
			// gradWf[0] += nvec[0]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// gradWf[1] += nvec[1]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			// gradWf[2] += nvec[2]*overRelaxCoeff*((orgWR-orgWL)/dPN - delPhiECF);
			
			// viscous term, non-orthogonal
			double vNonOrthFlux0 = 
				nonOrtholimiter * muF*(gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
			double vNonOrthFlux1 = 
				nonOrtholimiter * muF*(gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
			double vNonOrthFlux2 = 
				nonOrtholimiter * muF*(gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
			linB0[face.owner] += vNonOrthFlux0;
			linB1[face.owner] += vNonOrthFlux1;
			linB2[face.owner] += vNonOrthFlux2;
			
			
			// double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);;
			// linAL.back() -= muF/dPN*area*tmpCoeff2;
			
			// linA.at(face.owner)     += muF/dPN*area*tmpCoeff2;
			
			
			
			
			++proc_num;
			
		}
	}
	
	
	// // pressure term
	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		// linB0[i] -= cell.volume*gradP[i][0];
		// linB1[i] -= cell.volume*gradP[i][1];
		// linB2[i] -= cell.volume*gradP[i][2];
		
		// // cout << cell.volume*gradP[i][1] << endl;
	// }
	
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// boundary
	for(auto& boundary : mesh.boundary){
		
		if(boundary.neighbProcNo == -1){
			
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			
			if(
			boundary.type[1] == "noSlip" ||
			boundary.type[1] == "slip" ){
				for(int i=str; i<end; ++i){
					SEMO_Face& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec;
					nvec.push_back(face.unitNormals[0]);
					nvec.push_back(face.unitNormals[1]);
					nvec.push_back(face.unitNormals[2]);
					
					vector<double> distanceCells;
					distanceCells.push_back(face.distCells[0]);
					distanceCells.push_back(face.distCells[1]);
					distanceCells.push_back(face.distCells[2]);
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
		
					double UnL = UL*face.unitNormals[0]+
					             VL*face.unitNormals[1]+
								 WL*face.unitNormals[2];
					double UnR = UR*face.unitNormals[0]+
					             VR*face.unitNormals[1]+
								 WR*face.unitNormals[2];
					
					double UnF = 0.5*(UnL+UnR);
					double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					double RhoL = face.varL[controls.fRho];
					double RhoR = face.varR[controls.fRho];
					
					// // convective term
					// linA.at(face.owner) += RhoL*UnL*area;
					
					linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					
				}
			}
			else if(
			boundary.type[1] == "fixedValue" ||
			boundary.type[1] == "surfaceNormalFixedValue"
			){
				for(int i=str; i<end; ++i){
					SEMO_Face& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec;
					nvec.push_back(face.unitNormals[0]);
					nvec.push_back(face.unitNormals[1]);
					nvec.push_back(face.unitNormals[2]);
					
					vector<double> distanceCells;
					distanceCells.push_back(face.distCells[0]);
					distanceCells.push_back(face.distCells[1]);
					distanceCells.push_back(face.distCells[2]);
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double orgUL = mesh.cells[face.owner].var[controls.U];
					double orgVL = mesh.cells[face.owner].var[controls.V];
					double orgWL = mesh.cells[face.owner].var[controls.W];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
					
					double orgUnL = orgUL*face.unitNormals[0]+
					             orgVL*face.unitNormals[1]+
								 orgWL*face.unitNormals[2];
		
					double UnL = UL*face.unitNormals[0]+
					             VL*face.unitNormals[1]+
								 WL*face.unitNormals[2];
					double UnR = UR*face.unitNormals[0]+
					             VR*face.unitNormals[1]+
								 WR*face.unitNormals[2];
					
					double UnF = 0.5*(UnL+UnR);
					double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					double RhoL = face.varL[controls.fRho];
					double RhoR = face.varR[controls.fRho];
					
					// convective term
					// linA.at(face.owner) += 0.5*RhoL*UnL*area;
					
					// linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					// linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					// linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					linB0[face.owner] -= UR*RhoR*UnF*area;
					linB1[face.owner] -= VR*RhoR*UnF*area;
					linB2[face.owner] -= WR*RhoR*UnF*area;
					// linB0[face.owner] -= 0.5*(orgUL+UR)*RhoL*0.5*(orgUnL+UnR)*area;
					// linB1[face.owner] -= 0.5*(orgVL+VR)*RhoL*0.5*(orgUnL+UnR)*area;
					// linB2[face.owner] -= 0.5*(orgWL+WR)*RhoL*0.5*(orgUnL+UnR)*area;
					
				}
			}
			else if(boundary.type[1] == "zeroGradient"){
				for(int i=str; i<end; ++i){
					SEMO_Face& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec;
					nvec.push_back(face.unitNormals[0]);
					nvec.push_back(face.unitNormals[1]);
					nvec.push_back(face.unitNormals[2]);
					
					vector<double> distanceCells;
					distanceCells.push_back(face.distCells[0]);
					distanceCells.push_back(face.distCells[1]);
					distanceCells.push_back(face.distCells[2]);
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
		
					double UnL = UL*face.unitNormals[0]+
					             VL*face.unitNormals[1]+
								 WL*face.unitNormals[2];
					double UnR = UR*face.unitNormals[0]+
					             VR*face.unitNormals[1]+
								 WR*face.unitNormals[2];
					
					double UnF = 0.5*(UnL+UnR);
					double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					double RhoL = face.varL[controls.fRho];
					double RhoR = face.varR[controls.fRho];
					
					double orgUL = mesh.cells[face.owner].var[controls.U];
					double orgVL = mesh.cells[face.owner].var[controls.V];
					double orgWL = mesh.cells[face.owner].var[controls.W];
					
					// convective term
					linA[face.owner] += RhoL*UnF*area;
						
					linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
					linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
					linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					
				}
			}
			else if(boundary.type[1] == "inletOutlet"){
				
				for(int i=str; i<end; ++i){
					SEMO_Face& face = mesh.faces[i];
					
					double area = face.area;
					vector<double> nvec;
					nvec.push_back(face.unitNormals[0]);
					nvec.push_back(face.unitNormals[1]);
					nvec.push_back(face.unitNormals[2]);
					
					vector<double> distanceCells;
					distanceCells.push_back(face.distCells[0]);
					distanceCells.push_back(face.distCells[1]);
					distanceCells.push_back(face.distCells[2]);
					
					double UL = face.varL[controls.fU];
					double VL = face.varL[controls.fV];
					double WL = face.varL[controls.fW];
					
					double orgUL = mesh.cells[face.owner].var[controls.U];
					double orgVL = mesh.cells[face.owner].var[controls.V];
					double orgWL = mesh.cells[face.owner].var[controls.W];
					
					double UR = face.varR[controls.fU];
					double VR = face.varR[controls.fV];
					double WR = face.varR[controls.fW];
					
					double orgUnL = orgUL*face.unitNormals[0]+
					             orgVL*face.unitNormals[1]+
								 orgWL*face.unitNormals[2];
		
					double UnL = UL*face.unitNormals[0]+
					             VL*face.unitNormals[1]+
								 WL*face.unitNormals[2];
					double UnR = UR*face.unitNormals[0]+
					             VR*face.unitNormals[1]+
								 WR*face.unitNormals[2];
					
					double UnF = 0.5*(UnL+UnR);
					double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					double RhoL = face.varL[controls.fRho];
					double RhoR = face.varR[controls.fRho];
					
					// linB0[face.owner] -= UR*RhoR*UnF*area;
					// linB1[face.owner] -= VR*RhoR*UnF*area;
					// linB2[face.owner] -= WR*RhoR*UnF*area;
					
					//=========
					if( UnF > 0.0 ){
						// outflow
						// zeroGradient
						linA[face.owner] += RhoL*UnF*area;
							
						linB0[face.owner] -= 0.5*(UL+UR)*RhoL*UnF*area;
						linB1[face.owner] -= 0.5*(VL+VR)*RhoL*UnF*area;
						linB2[face.owner] -= 0.5*(WL+WR)*RhoL*UnF*area;
					}
					else{
						// inflow
						// fixedValue
						
						linB0[face.owner] -= UR*RhoR*UnF*area;
						linB1[face.owner] -= VR*RhoR*UnF*area;
						linB2[face.owner] -= WR*RhoR*UnF*area;
					}
					
					
						
					
					
				}
					
			}
			else{
				
				if(rank==0) cerr << "| #Error : not defined velocity B.C. name" << endl;
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				
			}
			
			
			
			// diffusion term
			for(int i=str; i<end; ++i){
				SEMO_Face& face = mesh.faces[i];
				
				double area = face.area;
				vector<double> nvec;
				nvec.push_back(face.unitNormals[0]);
				nvec.push_back(face.unitNormals[1]);
				nvec.push_back(face.unitNormals[2]);
				
				vector<double> distanceCells;
				distanceCells.push_back(face.distCells[0]);
				distanceCells.push_back(face.distCells[1]);
				distanceCells.push_back(face.distCells[2]);
				
				double UL = face.varL[controls.fU];
				double VL = face.varL[controls.fV];
				double WL = face.varL[controls.fW];
				
				double UR = face.varR[controls.fU];
				double VR = face.varR[controls.fV];
				double WR = face.varR[controls.fW];
	
				// diffusion term
				double muF = face.varL[controls.fmu];
				double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
				double dPN = sqrt(pow(distanceCells[0],2.0) + 
								  pow(distanceCells[1],2.0) + 
								  pow(distanceCells[2],2.0));
				// double Ef[3];
				// Ef[0] = distanceCells[0]/dPN;
				// Ef[1] = distanceCells[1]/dPN;
				// Ef[2] = distanceCells[2]/dPN;
				 
				// // over-relaxed approach
				// double overRelaxCoeff = 1.0;
				// // double overRelaxCoeff = Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2];
				// // double overRelaxCoeff = 1.0 / (Ef[0]*nvec[0]+Ef[1]*nvec[1]+Ef[2]*nvec[2]);
				
				// // double dPNEff = dPN * overRelaxCoeff;
				
				// // double Tf[3];
				// // Tf[0] = nvec[0] - Ef[0];
				// // Tf[1] = nvec[1] - Ef[1];
				// // Tf[2] = nvec[2] - Ef[2];
					
				double orgUL = mesh.cells[face.owner].var[controls.U];
				double orgVL = mesh.cells[face.owner].var[controls.V];
				double orgWL = mesh.cells[face.owner].var[controls.W];
				
				linA[face.owner] += muF/dPN_e*area;
				
				

				double corrUL = orgUL;
				double corrVL = orgVL;
				double corrWL = orgWL;
				// double tmpVec[3];
				// tmpVec[0] = distanceCells[0]*(1.0-
					// (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
				// tmpVec[1] = distanceCells[1]*(1.0-
					// (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
				// tmpVec[2] = distanceCells[2]*(1.0-
					// (distanceCells[0]*nvec[0]+distanceCells[1]*nvec[1]+distanceCells[2]*nvec[2])/dPN);
					
				// corrUL += gradP[face.owner][0]*tmpVec[0];
				// corrUL += gradP[face.owner][1]*tmpVec[1];
				// corrUL += gradP[face.owner][2]*tmpVec[2];
				
				// corrVL += gradP[face.owner][0]*tmpVec[0];
				// corrVL += gradP[face.owner][1]*tmpVec[1];
				// corrVL += gradP[face.owner][2]*tmpVec[2];
				
				// corrWL += gradP[face.owner][0]*tmpVec[0];
				// corrWL += gradP[face.owner][1]*tmpVec[1];
				// corrWL += gradP[face.owner][2]*tmpVec[2];
				
				
				linB0[face.owner] += muF*(UR-corrUL)/dPN_e*area;
				linB1[face.owner] += muF*(VR-corrVL)/dPN_e*area;
				linB2[face.owner] += muF*(WR-corrWL)/dPN_e*area;

				// // viscous term, properties
				// double delPhiECF;
				
				// double gradUf[3];
				// gradUf[0] = (gradU[face.owner][0]);
				// gradUf[1] = (gradU[face.owner][1]);
				// gradUf[2] = (gradU[face.owner][2]);
				// // delPhiECF = gradUf[0]*Ef[0]+gradUf[1]*Ef[1]+gradUf[2]*Ef[2];
				
				// // gradUf[0] += nvec[0]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				// // gradUf[1] += nvec[1]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				// // gradUf[2] += nvec[2]*overRelaxCoeff*((UR-orgUL)/dPN - delPhiECF);
				
				// double gradVf[3];
				// gradVf[0] = (gradV[face.owner][0]);
				// gradVf[1] = (gradV[face.owner][1]);
				// gradVf[2] = (gradV[face.owner][2]);
				// // delPhiECF = gradVf[0]*Ef[0]+gradVf[1]*Ef[1]+gradVf[2]*Ef[2];
				
				// // gradVf[0] += nvec[0]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				// // gradVf[1] += nvec[1]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				// // gradVf[2] += nvec[2]*overRelaxCoeff*((VR-orgVL)/dPN - delPhiECF);
				
				// double gradWf[3];
				// gradWf[0] = (gradW[face.owner][0]);
				// gradWf[1] = (gradW[face.owner][1]);
				// gradWf[2] = (gradW[face.owner][2]);
				// // delPhiECF = gradWf[0]*Ef[0]+gradWf[1]*Ef[1]+gradWf[2]*Ef[2];
				
				// // gradWf[0] += nvec[0]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
				// // gradWf[1] += nvec[1]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
				// // gradWf[2] += nvec[2]*overRelaxCoeff*((WR-orgWL)/dPN - delPhiECF);
						
						
				// // viscous term, non-orthogonal
				// linB0[face.owner] += muF*
					// (gradUf[0]*Tf[0]+gradUf[1]*Tf[1]+gradUf[2]*Tf[2])*area;
				// linB1[face.owner] += muF*
					// (gradVf[0]*Tf[0]+gradVf[1]*Tf[1]+gradVf[2]*Tf[2])*area;
				// linB2[face.owner] += muF*
					// (gradWf[0]*Tf[0]+gradWf[1]*Tf[1]+gradWf[2]*Tf[2])*area;
				
				
				// double tmpCoeff2 = overRelaxCoeff*(nvec[0]*Tf[0] + nvec[1]*Tf[1] + nvec[2]*Tf[2]);
			
				// linA[face.owner] += muF/dPN*area*tmpCoeff2;
				
				
				
			}
			
			
		}
		
	}
	
	
	
	// linear solver : PETSc library
	vector<double> resiVar0(mesh.cells.size(),0.0);
	vector<double> resiVar1(mesh.cells.size(),0.0);
	vector<double> resiVar2(mesh.cells.size(),0.0);
	
	// cout << "AAAAAAA" << endl;
	
	if(controls.iterPBs != controls.iterPBsMax-1){
	
		solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			controls.solverU, controls.toleranceU, 
			controls.relTolU, controls.preconditionerU,
			controls.maxIterU);
		solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			controls.solverU, controls.toleranceU, 
			controls.relTolU, controls.preconditionerU,
			controls.maxIterU);
		solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			controls.solverU, controls.toleranceU, 
			controls.relTolU, controls.preconditionerU,
			controls.maxIterU);
	
		// solveHYPRE(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solveHYPRE(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solveHYPRE(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		
	}
	else{
	
		solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			controls.solverFinalU, controls.toleranceFinalU, 
			controls.relTolFinalU, controls.preconditionerFinalU,
			controls.maxIterFinalU);
		solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			controls.solverFinalU, controls.toleranceFinalU, 
			controls.relTolFinalU, controls.preconditionerFinalU,
			controls.maxIterFinalU);
		solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			controls.solverFinalU, controls.toleranceFinalU, 
			controls.relTolFinalU, controls.preconditionerFinalU,
			controls.maxIterFinalU);
	
		// solveHYPRE(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solveHYPRE(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solveHYPRE(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
	}
	
	// cout << "BBBBBBB" << endl;
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		// resiVar0[i] = linB0[i]/linA[i];
		// resiVar1[i] = linB1[i]/linA[i];
		// resiVar2[i] = linB2[i]/linA[i];
		
		cell.var[controls.U] += controls.momVelURF * resiVar0[i];
		cell.var[controls.V] += controls.momVelURF * resiVar1[i];
		cell.var[controls.W] += controls.momVelURF * resiVar2[i];
		
		
		
		// if( cell.var[controls.P] <= controls.minP ) 
			// cell.var[controls.P] = controls.minP;
		// if( cell.var[controls.P] >= controls.maxP ) 
			// cell.var[controls.P] = controls.maxP;
		
		if( cell.var[controls.U] <= controls.minU ) 
			cell.var[controls.U] = controls.minU;
		if( cell.var[controls.U] >= controls.maxU ) 
			cell.var[controls.U] = controls.maxU;
		
		if( cell.var[controls.V] <= controls.minV ) 
			cell.var[controls.V] = controls.minV;
		if( cell.var[controls.V] >= controls.maxV ) 
			cell.var[controls.V] = controls.maxV;
		
		if( cell.var[controls.W] <= controls.minW ) 
			cell.var[controls.W] = controls.minW;
		if( cell.var[controls.W] >= controls.maxW ) 
			cell.var[controls.W] = controls.maxW;
		
		
		// cell.var.at(controls.UDV[0]) += gradP[i][0];
		// cell.var.at(controls.UDV[1]) += gradP[i][1];
		// cell.var.at(controls.UDV[2]) += gradP[i][2];
		
		
		residuals[1] += pow(resiVar0[i],2.0)*cell.volume;
		residuals[2] += pow(resiVar1[i],2.0)*cell.volume;
		residuals[3] += pow(resiVar2[i],2.0)*cell.volume;
		
		
		linAD[i] = linA[i] / cell.volume;
		linAOD[i] = 0.0;
		for(auto& j : cell.faces){
			linAOD[i] += linAR[j] / cell.volume;
		}
		
		
		
		gradP[i][0] /= cell.volume;
		gradP[i][1] /= cell.volume;
		gradP[i][2] /= cell.volume;
		
		
	}
	
	
	//========================
	// calc Un 
	// gradient P
	// vector<vector<double>> gradP;
	// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
	// vector<double> limGradP;
	// math.calcLimiterGradient(mesh, controls.P, controls.fP, gradP, limGradP);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] *= limGradP[i]; 
		// }
	// }
	// math.gradientGGBoundaryTreatment(mesh, controls, gradP);
	

	vector<double> gradPx_recv;
	vector<double> gradPy_recv;
	vector<double> gradPz_recv;
	vector<double> resiVar0_recv;
	vector<double> resiVar1_recv;
	vector<double> resiVar2_recv;
	vector<double> linAD_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> gradPx_send;
		vector<double> gradPy_send;
		vector<double> gradPz_send;
		vector<double> resiVar0_send;
		vector<double> resiVar1_send;
		vector<double> resiVar2_send;
		vector<double> linAD_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				gradPx_send.push_back(gradP[face.owner][0]);
				gradPy_send.push_back(gradP[face.owner][1]);
				gradPz_send.push_back(gradP[face.owner][2]);
				resiVar0_send.push_back(resiVar0[face.owner]);
				resiVar1_send.push_back(resiVar1[face.owner]);
				resiVar2_send.push_back(resiVar2[face.owner]);
				linAD_recv.push_back(linAD[face.owner]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatasDouble(
					gradPx_send, gradPx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradPy_send, gradPy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					gradPz_send, gradPz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		
		mpi.setProcsFaceDatasDouble(
					resiVar0_send, resiVar0_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					resiVar1_send, resiVar1_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatasDouble(
					resiVar2_send, resiVar2_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		mpi.setProcsFaceDatasDouble(
					linAD_send, linAD_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		// gradPx_send.clear();
		// gradPy_send.clear();
		// gradPz_send.clear();
	}
	
	
	
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];

		linAFL[i] = linAL[i];
		linAFR[i] = linAR[i];
		
		
		
		double area = face.area;
		vector<double> nvec(3,0.0);
		nvec[0] = face.unitNormals[0];
		nvec[1] = face.unitNormals[1];
		nvec[2] = face.unitNormals[2];
		
		vector<double> distanceCells(3,0.0);
		distanceCells[0] = face.distCells[0];
		distanceCells[1] = face.distCells[1];
		distanceCells[2] = face.distCells[2];
		
		double wCL = face.wC;
		double wCR = 1.0-face.wC;
		
		double RhoL = face.varL[controls.fRho];
		
		double RhoR = face.varR[controls.fRho];
		
		double dPN_e = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		
		double dPN = sqrt(pow(distanceCells[0],2.0) + 
						  pow(distanceCells[1],2.0) + 
						  pow(distanceCells[2],2.0));
				
		double nonOrtholimiter = 1.0;
		double cosAlpha = dPN_e/dPN;
		if( cosAlpha < 0.766 && cosAlpha >= 0.5 ){
			nonOrtholimiter = 0.5;
		}
		else if( cosAlpha < 0.5 && cosAlpha >= 0.342 ){
			nonOrtholimiter = 0.333;
		}
		else if( cosAlpha < 0.342 ){
			nonOrtholimiter = 0.0;
		}
				
		double Ef[3];
		Ef[0] = distanceCells[0]/dPN;
		Ef[1] = distanceCells[1]/dPN;
		Ef[2] = distanceCells[2]/dPN;
				
		double Tf[3];
		Tf[0] = nvec[0] - distanceCells[0]/dPN_e;
		Tf[1] = nvec[1] - distanceCells[1]/dPN_e;
		Tf[2] = nvec[2] - distanceCells[2]/dPN_e;
		
		double tmp1 = (wCL/RhoL + wCR/RhoR)*controls.timeStep;
		

		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				
			face.var[controls.Un] += controls.momVelURF * (
				(wCL*resiVar0[face.owner]+wCR*resiVar0[face.neighbour])*nvec[0] +
				(wCL*resiVar1[face.owner]+wCR*resiVar1[face.neighbour])*nvec[1] +
				(wCL*resiVar2[face.owner]+wCR*resiVar2[face.neighbour])*nvec[2] );
				
				
			// pressure correction
			
			
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = mesh.cells[face.neighbour].var[controls.P];
				
			tmp1 = wCL/linAD[face.owner] + wCR/linAD[face.neighbour];
			face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				wCL*1.0/linAD[face.owner]*( 
					gradP[face.owner][0]*nvec[0] +
					gradP[face.owner][1]*nvec[1] +
					gradP[face.owner][2]*nvec[2] ) +
				wCR*1.0/linAD[face.neighbour]*( 
					gradP[face.neighbour][0]*nvec[0] +
					gradP[face.neighbour][1]*nvec[1] +
					gradP[face.neighbour][2]*nvec[2] ) );
			// face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// wCL*controls.timeStep/RhoL*( 
					// gradP[face.owner][0]*nvec[0] +
					// gradP[face.owner][1]*nvec[1] +
					// gradP[face.owner][2]*nvec[2] ) +
				// wCR*controls.timeStep/RhoR*( 
					// gradP[face.neighbour][0]*nvec[0] +
					// gradP[face.neighbour][1]*nvec[1] +
					// gradP[face.neighbour][2]*nvec[2] ) );
				
			face.var[controls.Un] -= coeffDiss * controls.momVelURF * tmp1*(orgPR-orgPL)/dPN_e;
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradP[face.neighbour][0]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradP[face.neighbour][1]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradP[face.neighbour][2]);
			
			face.var[controls.Un] -= coeffDiss * controls.momVelURF 
				* nonOrtholimiter * tmp1*
				(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				
			face.var[controls.Un] += controls.momVelURF * (
				(wCL*resiVar0[face.owner]+wCR*resiVar0_recv[proc_num])*nvec[0] +
				(wCL*resiVar1[face.owner]+wCR*resiVar1_recv[proc_num])*nvec[1] +
				(wCL*resiVar2[face.owner]+wCR*resiVar2_recv[proc_num])*nvec[2] );
				
				
			// pressure correction
			double orgPL = mesh.cells[face.owner].var[controls.P];
			double orgPR = face.varR[controls.fP];
				
			tmp1 = wCL/linAD[face.owner] + wCR/linAD_recv[proc_num];
			face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				wCL*1.0/linAD[face.owner]*( 
					gradP[face.owner][0]*nvec[0] +
					gradP[face.owner][1]*nvec[1] +
					gradP[face.owner][2]*nvec[2] ) +
				wCR*1.0/linAD_recv[proc_num]*( 
					gradPx_recv[proc_num]*nvec[0] +
					gradPy_recv[proc_num]*nvec[1] +
					gradPz_recv[proc_num]*nvec[2] ) );
			// face.var[controls.Un] += coeffDiss * controls.momVelURF * (
				// wCL*controls.timeStep/RhoL*( 
					// gradP[face.owner][0]*nvec[0] +
					// gradP[face.owner][1]*nvec[1] +
					// gradP[face.owner][2]*nvec[2] ) +
				// wCR*controls.timeStep/RhoR*( 
					// gradPx_recv[proc_num]*nvec[0] +
					// gradPy_recv[proc_num]*nvec[1] +
					// gradPz_recv[proc_num]*nvec[2] ) );
				
				
			face.var[controls.Un] -= coeffDiss * controls.momVelURF * tmp1*(orgPR-orgPL)/dPN_e;
			
			// non-orthogonal, over-relaxed approach
			double gradPf[3];
			gradPf[0] = (wCL*gradP[face.owner][0]+wCR*gradPx_recv[proc_num]);
			gradPf[1] = (wCL*gradP[face.owner][1]+wCR*gradPy_recv[proc_num]);
			gradPf[2] = (wCL*gradP[face.owner][2]+wCR*gradPz_recv[proc_num]);
			
			face.var[controls.Un] -= coeffDiss * controls.momVelURF 
				* nonOrtholimiter * tmp1*
				(gradPf[0]*Tf[0] + gradPf[1]*Tf[1] + gradPf[2]*Tf[2]);
			
				
			++proc_num;
			
		}
		

	}
	
	
	linA.clear();
	linAL.clear();
	linAR.clear();
	linB0.clear();
	linB1.clear();
	linB2.clear();
	resiVar0.clear();
	resiVar1.clear();
	resiVar2.clear();
	
	
}








// // add Rho
// void SEMO_Solvers_Builder::calcMomentumEqs(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& residuals){
		
		

	// SEMO_MPI_Builder mpi;

	// SEMO_Utility_Math math;
	
	// // gradient P
	// vector<vector<double>> gradP(mesh.cells.size(),vector<double>(3,0.0));
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// double PF = 0.5*(face.varL[controls.fP]+face.varR[controls.fP]);
		
		// for(int j=0; j<3; ++j){
			// gradP[face.owner][j] += PF*face.unitNormals[j]*face.area;
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// gradP[face.neighbour][j] -= PF*face.unitNormals[j]*face.area;
			// }
		// }
	// }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// for(int j=0; j<3; ++j){
			// gradP[i][j] /= mesh.cells[i].volume;
		// }
	// }
	// // processor faces
	// // gradP , 
	// vector<double> gradPx_send;
	// vector<double> gradPy_send;
	// vector<double> gradPz_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// gradPx_send.push_back(gradP[face.owner][0]);
			// gradPy_send.push_back(gradP[face.owner][1]);
			// gradPz_send.push_back(gradP[face.owner][2]);
		// }
	// }
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// mpi.setProcsFaceDatasDouble(
				// gradPx_send, gradPx_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPy_send, gradPy_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPz_send, gradPz_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// gradPx_send.clear();
	// gradPy_send.clear();
	// gradPz_send.clear();
	
	
		
	// vector<double> linA;
	// vector<double> linB0;
	// vector<double> linB1;
	// vector<double> linB2;
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		
		// linA.push_back(cell.var[controls.Rho]*cell.volume/controls.timeStep);
		
		// linB0.push_back(-
			// (cell.var[controls.Rho]*cell.var[controls.U]
			// -cell.var[controls.oldRho]*cell.var[controls.oldU])
			// *cell.volume/controls.timeStep);
			
		// linB1.push_back(-
			// (cell.var[controls.Rho]*cell.var[controls.V]
			// -cell.var[controls.oldRho]*cell.var[controls.oldV])
			// *cell.volume/controls.timeStep);
			
		// linB2.push_back(-
			// (cell.var[controls.Rho]*cell.var[controls.W]
			// -cell.var[controls.oldRho]*cell.var[controls.oldW])
			// *cell.volume/controls.timeStep);
		
		// linB0.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[0];
		// linB1.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[1];
		// linB2.back() += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[2];
		
	// }
	
	// // cout << "AAAAAAAAAAAA" << endl;
	
	
	// vector<double> linAL;
	// vector<double> linAR;

	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double area = face.area;
		// vector<double> nvec;
		// nvec.push_back(face.unitNormals[0]);
		// nvec.push_back(face.unitNormals[1]);
		// nvec.push_back(face.unitNormals[2]);
		
		// vector<double> distanceCells;
		// distanceCells.push_back(face.distCells[0]);
		// distanceCells.push_back(face.distCells[1]);
		// distanceCells.push_back(face.distCells[2]);
	
		// double UL = face.varL[controls.fU];
		// double VL = face.varL[controls.fV];
		// double WL = face.varL[controls.fW];
		// double PL = face.varL[controls.fP];
		// double RhoL = face.varL[controls.fRho];
		
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		// double PR = face.varR[controls.fP];
		// double RhoR = face.varR[controls.fRho];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		// double UnF = 0.5*(UnL+UnR);

		// double dPN = nvec[0]*distanceCells[0] + nvec[1]*distanceCells[1] + nvec[2]*distanceCells[2];
		// double coeffP = area/dPN;

		// double coeff1 = 0.5*(1.0/RhoL + 1.0/RhoR)*controls.timeStep;
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // UnF = UnF 
				// // - coeff1
				// // *((PR-PL)/dPN - 0.5*(
				 // // gradP[face.owner][0]*nvec[0]+gradP[face.neighbour][0]*nvec[0]
				// // +gradP[face.owner][1]*nvec[1]+gradP[face.neighbour][1]*nvec[1]
				// // +gradP[face.owner][2]*nvec[2]+gradP[face.neighbour][2]*nvec[2]
				// // ));
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // UnF = UnF 
				// // - coeff1
				// // *((PR-PL)/dPN - 0.5*(
				 // // gradP[face.owner][0]*nvec[0]+gradPx_recv[proc_num]*nvec[0]
				// // +gradP[face.owner][1]*nvec[1]+gradPy_recv[proc_num]*nvec[1]
				// // +gradP[face.owner][2]*nvec[2]+gradPz_recv[proc_num]*nvec[2]
				// // ));
			// ++proc_num;
		// }
		
		
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// double UF = UL*weightL + UR*weightR;
		// double VF = VL*weightL + VR*weightR;
		// double WF = WL*weightL + WR*weightR;
		
		// double RhoF = RhoL*weightL + RhoR*weightR;
		
		// linAL.push_back(-weightL*RhoF*UnF*area);
		// linAR.push_back(+weightR*RhoF*UnF*area);
		
		// // convective term
		// linB0.at(face.owner) -= UF*RhoF*UnF*area;
		// linB1.at(face.owner) -= VF*RhoF*UnF*area;
		// linB2.at(face.owner) -= WF*RhoF*UnF*area;
		
		// // pressure term
		// linB0.at(face.owner) -= 0.5*(PL+PR)*nvec[0]*area;
		// linB1.at(face.owner) -= 0.5*(PL+PR)*nvec[1]*area;
		// linB2.at(face.owner) -= 0.5*(PL+PR)*nvec[2]*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA.at(face.owner)     += (-linAL.back());
			// linA.at(face.neighbour) += (-linAR.back());
			
			// // convective term
			// linB0.at(face.neighbour) += UF*RhoF*UnF*area;
			// linB1.at(face.neighbour) += VF*RhoF*UnF*area;
			// linB2.at(face.neighbour) += WF*RhoF*UnF*area;
			
			// // pressure term
			// linB0.at(face.neighbour) += 0.5*(PL+PR)*nvec[0]*area;
			// linB1.at(face.neighbour) += 0.5*(PL+PR)*nvec[1]*area;
			// linB2.at(face.neighbour) += 0.5*(PL+PR)*nvec[2]*area;
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA.at(face.owner) += (-linAL.back());
			
		// }
		
	// }
	
	
	
	// // // pressure term
	// // proc_num=0;
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // SEMO_Cell& cell = mesh.cells[i];
		
		// // linB0.at(i) -= gradP[i][0]*cell.volume;
		// // linB1.at(i) -= gradP[i][1]*cell.volume;
		// // linB2.at(i) -= gradP[i][2]*cell.volume;
		
		
	// // }
	
	
	
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(
			// boundary.type[1] == "zeroGradient" || 
			// boundary.type[1] == "slip"){
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					// double UnL = face.varL[controls.fU]*face.unitNormals[0]+
					             // face.varL[controls.fV]*face.unitNormals[1]+
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0]+
					             // face.varR[controls.fV]*face.unitNormals[1]+
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// double RhoF = RhoL*weightL + RhoR*weightR;
					
					// linA.at(face.owner) += RhoF*UnL*face.area;
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					// double UnL = face.varL[controls.fU]*face.unitNormals[0]+
					             // face.varL[controls.fV]*face.unitNormals[1]+
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0]+
					             // face.varR[controls.fV]*face.unitNormals[1]+
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// double RhoF = RhoL*weightL + RhoR*weightR;
					
					// linA.at(face.owner) += 0.5*RhoL*UnL*face.area;
				// }
				
			// }
			
		// }
		
	// }
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar0(mesh.cells.size(),0.0);
	// vector<double> resiVar1(mesh.cells.size(),0.0);
	// vector<double> resiVar2(mesh.cells.size(),0.0);
	
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
	// }
	
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		
		// cell.var.at(controls.U) += controls.momVelURF * resiVar0.at(i);
		// cell.var.at(controls.V) += controls.momVelURF * resiVar1.at(i);
		// cell.var.at(controls.W) += controls.momVelURF * resiVar2.at(i);
		
		
		// residuals[1] += pow(resiVar0[i],2.0)*cell.volume;
		// residuals[2] += pow(resiVar1[i],2.0)*cell.volume;
		// residuals[3] += pow(resiVar2[i],2.0)*cell.volume;
		
		
		// linAD[i] = linA[i] / cell.volume;
		// linAOD[i] = 0.0;
		// for(auto& j : cell.faces){
			// linAOD[i] += linAR[j] / cell.volume;
		// }
		
		
	// }
	



	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB0.clear();
	// linB1.clear();
	// linB2.clear();
	// resiVar0.clear();
	// resiVar1.clear();
	// resiVar2.clear();
	
	
// }



















// void SEMO_Solvers_Builder::calcMomentumEqs(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<double>& linAD,
	// vector<double>& linAOD,
	// vector<double>& residuals){
		
		

	// SEMO_MPI_Builder mpi;

	// SEMO_Utility_Math math;
	// vector<vector<double>> gradP;
	// math.calcLeastSquare(mesh, controls.P, controls.fP, gradP);
	
	

	// // processor faces
	// // gradP , 
	// vector<double> gradPx_send;
	// vector<double> gradPy_send;
	// vector<double> gradPz_send;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// gradPx_send.push_back(gradP[face.owner][0]);
			// gradPy_send.push_back(gradP[face.owner][1]);
			// gradPz_send.push_back(gradP[face.owner][2]);
		// }
	// }
	
	// vector<double> gradPx_recv;
	// vector<double> gradPy_recv;
	// vector<double> gradPz_recv;
	// mpi.setProcsFaceDatasDouble(
				// gradPx_send, gradPx_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPy_send, gradPy_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// mpi.setProcsFaceDatasDouble(
				// gradPz_send, gradPz_recv,
				// mesh.countsProcFaces, mesh.countsProcFaces, 
				// mesh.displsProcFaces, mesh.displsProcFaces);
	// gradPx_send.clear();
	// gradPy_send.clear();
	// gradPz_send.clear();
		
		
		
		
	
	// // vector<double> linA(mesh.cells.size(),0.0);
	// // vector<double> linAL(mesh.faces.size(),0.0);
	// // vector<double> linAR(mesh.faces.size(),0.0);
	// // vector<double> linB0(mesh.cells.size(),0.0);
	// // vector<double> linB1(mesh.cells.size(),0.0);
	// // vector<double> linB2(mesh.cells.size(),0.0);
	
	// vector<double> linA;
	// vector<double> linB0;
	// vector<double> linB1;
	// vector<double> linB2;
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// // linA[i] = cell.volume/controls.timeStep;
		
		// // linB0[i] = -(cell.var[controls.U]-cell.var[controls.oldU])*cell.volume/controls.timeStep;
		// // linB1[i] = -(cell.var[controls.V]-cell.var[controls.oldV])*cell.volume/controls.timeStep;
		// // linB2[i] = -(cell.var[controls.W]-cell.var[controls.oldW])*cell.volume/controls.timeStep;
		
		// // linB0[i] += cell.volume*controls.gravityAcceleration[0];
		// // linB1[i] += cell.volume*controls.gravityAcceleration[1];
		// // linB2[i] += cell.volume*controls.gravityAcceleration[2];
		
		// linA.push_back(cell.volume/controls.timeStep);
		
		// linB0.push_back(
			// -(cell.var[controls.U]-cell.var[controls.oldU])*cell.volume/controls.timeStep);
		// linB1.push_back(
			// -(cell.var[controls.V]-cell.var[controls.oldV])*cell.volume/controls.timeStep);
		// linB2.push_back(
			// -(cell.var[controls.W]-cell.var[controls.oldW])*cell.volume/controls.timeStep);
			
			
			
		// linB0.back() += cell.volume*controls.gravityAcceleration[0];
		// linB1.back() += cell.volume*controls.gravityAcceleration[1];
		// linB2.back() += cell.volume*controls.gravityAcceleration[2];
		
		// // cout << linB2.back() << endl;
		
		// // linA[i] = cell.var[controls.Rho]*cell.volume/controls.timeStep;
		
		// // linB0[i] = -cell.var[controls.Rho]*(cell.var[controls.U]-cell.var[controls.oldU])*cell.volume/controls.timeStep;
		// // linB1[i] = -cell.var[controls.Rho]*(cell.var[controls.V]-cell.var[controls.oldV])*cell.volume/controls.timeStep;
		// // linB2[i] = -cell.var[controls.Rho]*(cell.var[controls.W]-cell.var[controls.oldW])*cell.volume/controls.timeStep;
		
		// // linB0[i] += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[0];
		// // linB1[i] += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[1];
		// // linB2[i] += cell.var[controls.Rho]*cell.volume*controls.gravityAcceleration[2];
		
	// }
	
	
	
	
	// vector<double> linAL;
	// vector<double> linAR;

	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// double area = face.area;
		// vector<double> nvec(3,0.0);
		// nvec[0] = face.unitNormals[0];
		// nvec[1] = face.unitNormals[1];
		// nvec[2] = face.unitNormals[2];
		
		// vector<double> distanceCells(3,0.0);
		// distanceCells[0] = face.distCells[0];
		// distanceCells[1] = face.distCells[1];
		// distanceCells[2] = face.distCells[2];
	
		// // if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			// // distanceCells[0] *= 2.0;
			// // distanceCells[1] *= 2.0;
			// // distanceCells[2] *= 2.0;
		// // }
		// double UL = face.varL[controls.fU];
		// double VL = face.varL[controls.fV];
		// double WL = face.varL[controls.fW];
		// double UR = face.varR[controls.fU];
		// double VR = face.varR[controls.fV];
		// double WR = face.varR[controls.fW];
		
		// double PL = face.varL[controls.fP];
		// double PR = face.varR[controls.fP];
		
		// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
		// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
		
		
		// double RhoL = face.varL[controls.fRho];
		// double RhoR = face.varR[controls.fRho];
		
		// // if( RhoL != RhoR){
		// // cout << RhoL << " " << RhoR << endl;
		// // }
	
		// // double dPN = 
			// // sqrt(std::inner_product(std::begin(distanceCells), std::end(distanceCells), 
			                        // // std::begin(distanceCells), 0.0));
		// // double dPNEff = area/dPN/dPN*
			// // (std::inner_product(std::begin(nvec), std::end(nvec), 
			                    // // std::begin(distanceCells), 0.0));
		// double dPN = sqrt(pow(distanceCells[0],2.0)+pow(distanceCells[1],2.0)+pow(distanceCells[2],2.0));
		
		// double UnF = 0.5*(UnL+UnR);
		
		// // double tmp1 = 0.5*(1.0/RhoL+1.0/RhoR)*controls.timeStep;
		// // if( face.getType() == SEMO_Types::INTERNAL_FACE ){
			// // // // tmp1 = 0.5*
				// // // // (mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // // // +mesh.cells[face.neighbour].volume/linAD[face.neighbour]/RhoR);
			
			// // double tmp0 = gradP[face.owner][0]*nvec[0]+gradP[face.neighbour][0]*nvec[0];
			// // tmp0 += gradP[face.owner][1]*nvec[1]+gradP[face.neighbour][1]*nvec[1];
			// // tmp0 += gradP[face.owner][2]*nvec[2]+gradP[face.neighbour][2]*nvec[2];
			// // tmp0 = tmp1*((PR-PL)/dPN - 0.5*tmp0);
			// // UnF -= tmp0;
		// // }
		// // else if( face.getType() == SEMO_Types::PROCESSOR_FACE ){
			// // // // tmp1 = 0.5*
				// // // // (mesh.cells[face.owner].volume/linAD[face.owner]/RhoL
				// // // // +volume_recv[proc_num]/linAD_recv[proc_num]/RhoR);
				
			// // double tmp0 = gradP[face.owner][0]*nvec[0]+gradPx_recv[proc_num]*nvec[0];
			// // tmp0 += gradP[face.owner][1]*nvec[1]+gradPy_recv[proc_num]*nvec[1];
			// // tmp0 += gradP[face.owner][2]*nvec[2]+gradPz_recv[proc_num]*nvec[2];
			// // tmp0 = tmp1*((PR-PL)/dPN - 0.5*tmp0);
			// // UnF -= tmp0;
			
			// // ++proc_num;
		// // }
		// // // else{
			// // // // // tmp1 = mesh.cells[face.owner].volume/linAD[face.owner]/RhoL;
			
		// // // }
		
		// // // UnF -= tmp1*(PR-PL)/dPN;
		
		
		
		
		

			// // tmp1 = 0.5*(PL+PR);

	// // // SLAU
	// // double chat = 10000.0;
	// // double ML = UnL/chat; 
	// // double MR = UnR/chat;
	// // double KLR = sqrt(0.5*(UL*UL+VL*VL+WL*WL+UR*UR+VR*VR+WR*WR));
	
    // // double Mbar = ( RhoL*abs(ML)+RhoR*abs(MR) ) / ( RhoL + RhoR );
	// // double Mcy = min(1.0,KLR/chat);
	// // double phi_c = pow((1.0-Mcy),2.0);
	// // double g_c = 1.0 + max( min(ML,0.0), -1.0 )*min( max(MR,0.0), 1.0 );

    // // double D_L = ML+(1.0-g_c)*abs(ML);
    // // double D_R = MR-(1.0-g_c)*abs(MR);
    // // double D_rho = Mbar*g_c;
	
	// // double MLP_SLAU = 0.5*(D_L+D_rho);
	// // double MRM_SLAU = 0.5*(D_R-D_rho);
	// // double mdot = RhoL*chat*MLP_SLAU + RhoR*chat*MRM_SLAU
	   // // - 0.5
		 // // *(1.0-0.5*(1.0-cos(3.141592*min(1.0,max(abs(ML),abs(MR))))))
		 // // *(1.0-0.5*(1.0+cos(3.141592*min(abs(PL/PR),abs(PR/PL)))))
		 // // *(PR-PL)
		 // // /chat;
		 
	
	// // // double MLP  = M_func(ML,1.0,0.0); 
	// // // double MRM = M_func(MR,-1.0,0.0);
	// // double preP = pre_func(ML,1.0,0.0); 
	// // double preM = pre_func(MR,-1.0,0.0);
	// // double gam = 0.6;
	// // double PLR = 0.5*(PL+PR) 
				// // - 1.0*(KLR/chat)*0.5*preP*preM*0.5*(PL+PR)/chat*(UnR-UnL)
				// // + max(0.2,gam)*(KLR/chat)*0.5*(PL+PR)*(preP+preM-1.0)
				// // - 0.5*(preP-preM)*(PR-PL);
				
				
	// // UnF = mdot / (0.5*(RhoL+RhoR));
	// // tmp1 = PLR;

		// // double signUnF = (UnF > 0.0) ? 1.0 : -1.0;
		// // double weightL = 0.5*(1.0+signUnF);
		// // double weightR = 0.5*(1.0-signUnF);
		// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
		// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
		
		// // weightL = 0.5;
		// // weightR = 0.5;
		
		// double UF = UL*weightL + UR*weightR;
		// double VF = VL*weightL + VR*weightR;
		// double WF = WL*weightL + WR*weightR;
		
		// double RhoF = RhoL*weightL + RhoR*weightR;
		

		// linAL.push_back(-weightL*UnF*area);
		// linAR.push_back(+weightR*UnF*area);
		// // linAL.push_back(0.0);
		// // linAR.push_back(0.0);
		
		// // convective term
		// linB0[face.owner] -= UF*UnF*area;
		// linB1[face.owner] -= VF*UnF*area;
		// linB2[face.owner] -= WF*UnF*area;
		
		// // pressure term
		// linB0[face.owner] -= 0.5*(PL+PR)/RhoL*nvec[0]*area;
		// linB1[face.owner] -= 0.5*(PL+PR)/RhoL*nvec[1]*area;
		// linB2[face.owner] -= 0.5*(PL+PR)/RhoL*nvec[2]*area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			// linA[face.neighbour] += (-linAR[i]);
			
			// // convective term
			// linB0[face.neighbour] += UF*UnF*area;
			// linB1[face.neighbour] += VF*UnF*area;
			// linB2[face.neighbour] += WF*UnF*area;
			
			// // pressure term
			// linB0[face.neighbour] += 0.5*(PL+PR)/RhoR*nvec[0]*area;
			// linB1[face.neighbour] += 0.5*(PL+PR)/RhoR*nvec[1]*area;
			// linB2[face.neighbour] += 0.5*(PL+PR)/RhoR*nvec[2]*area;
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// linA[face.owner] += (-linAL[i]);
			
			
		// }
		
		// // if((RhoL > 1.17766 && RhoL < 1.17764) || (RhoR > 1.17766 && RhoR < 1.17764)){
		// // cout << RhoL << " " << RhoR << endl;
		// // }
		
	// }
	
	// // boundary
	// for(auto& boundary : mesh.boundary){
		
		// if(boundary.neighbProcNo == -1){
			
			// int str = boundary.startFace;
			// int end = str + boundary.nFaces;
			
			// if(
			// boundary.type[1] == "zeroGradient" || 
			// boundary.type[1] == "slip"){
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					// double UnL = face.varL[controls.fU]*face.unitNormals[0]+
					             // face.varL[controls.fV]*face.unitNormals[1]+
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0]+
					             // face.varR[controls.fV]*face.unitNormals[1]+
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// double RhoF = RhoL*weightL + RhoR*weightR;
					
					// // linA[face.owner] += weightL*RhoF*UnF*face.area;
					// // linA[face.owner] += RhoL*UnL*face.area;
					// linA[face.owner] += UnL*face.area;
				// }
			// }
			// else{
				// for(int i=str; i<end; ++i){
					// SEMO_Face& face = mesh.faces[i];
					// double UnL = face.varL[controls.fU]*face.unitNormals[0]+
					             // face.varL[controls.fV]*face.unitNormals[1]+
								 // face.varL[controls.fW]*face.unitNormals[2];
					// double UnR = face.varR[controls.fU]*face.unitNormals[0]+
					             // face.varR[controls.fV]*face.unitNormals[1]+
								 // face.varR[controls.fW]*face.unitNormals[2];
					
					// double UnF = 0.5*(UnL+UnR);
					// double weightL = (UnF > 0.0) ? 1.0 : 0.0;
					// double weightR = (UnF < 0.0) ? 1.0 : 0.0;
					
					// double RhoL = face.varL[controls.fRho];
					// double RhoR = face.varR[controls.fRho];
					
					// double RhoF = RhoL*weightL + RhoR*weightR;
					
					// // linA[face.owner] += weightL*RhoF*UnF*face.area;
					// linA[face.owner] += 0.5*UnL*face.area;
				// }
				
			// }
			
		// }
		
	// }
	
	
	
	
	


	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // cout << linB1[i] << endl;
		
	// // }
	
	
	
	
	
	
	// // linear solver : PETSc library
	// vector<double> resiVar0(mesh.cells.size(),0.0);
	// vector<double> resiVar1(mesh.cells.size(),0.0);
	// vector<double> resiVar2(mesh.cells.size(),0.0);
	
	
	// if(controls.iterPBs != controls.iterPBsMax-1){
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverU, controls.toleranceU, 
			// controls.relTolU, controls.preconditionerU,
			// controls.maxIterU);
		
	// }
	// else{
	
		// solvePETSc(mesh, resiVar0, linA, linAL, linAR, linB0,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar1, linA, linAL, linAR, linB1,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
		// solvePETSc(mesh, resiVar2, linA, linAL, linAR, linB2,
			// controls.solverFinalU, controls.toleranceFinalU, 
			// controls.relTolFinalU, controls.preconditionerFinalU,
			// controls.maxIterFinalU);
	// }
	
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// SEMO_Cell& cell = mesh.cells[i];
		
		
		// // resiVar0[i] = linB0[i]/linA[i];
		// // resiVar1[i] = linB1[i]/linA[i];
		// // resiVar2[i] = linB2[i]/linA[i];
		
		// cell.var[controls.U] += controls.momVelURF * resiVar0[i];
		
		// cell.var[controls.V] += controls.momVelURF * resiVar1[i];
		
		// cell.var[controls.W] += controls.momVelURF * resiVar2[i];
		
		// residuals[1] += pow(resiVar0[i],2.0)*cell.volume;
		// residuals[2] += pow(resiVar1[i],2.0)*cell.volume;
		// residuals[3] += pow(resiVar2[i],2.0)*cell.volume;
		
		
		
		
		// linAD[i] = linA[i];
		// linAOD[i] = 0.0;
		// for(auto& j : cell.faces){
			// linAOD[i] += linAR[j];
		// }
		
		
		// // if(resiVar1[i]>0.0){
		// // cout << resiVar0[i] << " " << resiVar1[i] << " " << resiVar2[i] << endl;
		// // }
		
		
		// // cout << residuals[1] << " " << residuals[2] << " " <<residuals[3] << endl;
		
		
		// // cout << linB0[i] << " " << linB1[i] << " " << linB2[i] << endl;
		
		// // if( cell.var[controls.U] <= controls.minU ) 
			// // cell.var[controls.U] = controls.minU;
		// // if( cell.var[controls.U] >= controls.maxU ) 
			// // cell.var[controls.U] = controls.maxU;
		
		// // if( cell.var[controls.V] <= controls.minV ) 
			// // cell.var[controls.V] = controls.minV;
		// // if( cell.var[controls.V] >= controls.maxV ) 
			// // cell.var[controls.V] = controls.maxV;
		
		// // if( cell.var[controls.W] <= controls.minW ) 
			// // cell.var[controls.W] = controls.minW;
		// // if( cell.var[controls.W] >= controls.maxW ) 
			// // cell.var[controls.W] = controls.maxW;
		
		
	// }
	



	// linA.clear();
	// linAL.clear();
	// linAR.clear();
	// linB0.clear();
	// linB1.clear();
	// linB2.clear();
	// resiVar0.clear();
	// resiVar1.clear();
	// resiVar2.clear();
	
	
// }