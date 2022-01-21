#include "build.h"
#include <cmath>
#include <array>

void SEMO_Solvers_Builder::calcRHS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	vector<vector<double>>& residuals){
	
	
	this->calcFluxes(mesh, controls, species, residuals);
	
	
	this->sourceTerms(mesh, controls, species, residuals);
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		SEMO_Cell& cell = mesh.cells[i];
		
		for(int nq=0; nq<controls.nEq; ++nq){
			residuals[i][nq] *= -cell.var[controls.dtPseudo]/cell.volume;
		}
	}
	
	
}


void SEMO_Solvers_Builder::calcSingleRHS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	vector<vector<double>>& residuals){
	
	// cout << "bbbbb0" << endl;
	
	this->calcFluxes(mesh, controls, species, residuals);
	// cout << "bbbbb1" << endl;
	
	
	// this->sourceTerms(mesh, controls, residuals);
	bool boolGravityOn = true;
	bool boolSurfTensOn = true;
	double magGrav = controls.gravityAcceleration[0]*controls.gravityAcceleration[0];
	magGrav += controls.gravityAcceleration[1]*controls.gravityAcceleration[1];
	magGrav += controls.gravityAcceleration[2]*controls.gravityAcceleration[2];
	if(magGrav < 1.e-200) boolGravityOn = false;
	if(abs(species[0].sigma) < 1.e-200) boolSurfTensOn = false;
	if(boolGravityOn) this->calcSourceGravity(mesh, controls, species);
	if(boolSurfTensOn) this->calcSourceSurfaceTension(mesh, controls, species);

	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		double volume = cell.volume;
		// gravity
		// residuals[i][1] -= cell.var[controls.Rho] * controls.gravityAcceleration[0] * cell.volume; 
		// residuals[i][2] -= cell.var[controls.Rho] * controls.gravityAcceleration[1] * cell.volume;
		// residuals[i][3] -= cell.var[controls.Rho] * controls.gravityAcceleration[2] * cell.volume;
		if(boolGravityOn){
			residuals[i][1] -= volume*cell.var[controls.sourceGravity[0]];
			residuals[i][2] -= volume*cell.var[controls.sourceGravity[1]];
			residuals[i][3] -= volume*cell.var[controls.sourceGravity[2]];
		}
		
		// surface tension force terms
		if(boolSurfTensOn){
			residuals[i][1] -= volume*cell.var[controls.sourceSurfaceTension[0]];
			residuals[i][2] -= volume*cell.var[controls.sourceSurfaceTension[1]];
			residuals[i][3] -= volume*cell.var[controls.sourceSurfaceTension[2]];
		}
	}
	// cout << "bbbbb2" << endl;
	

	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		SEMO_Cell& cell = mesh.cells[i];
		
		for(int nq=0; nq<controls.nEq; ++nq){
			residuals[i][nq] *= (-controls.timeStep/cell.volume);
		// cout << residuals[i][nq] << endl; 
		}
		
	}
	// cout << "bbbbb3" << endl;
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}





void SEMO_Solvers_Builder::calcFluxes(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species,
	vector<vector<double>>& residuals){
		


	// calc gradient
	SEMO_Utility_Math math;
	// vector<vector<double>> gradU;
	// math.calcGaussGreen(mesh, controls.U, controls.fU, gradU);
	// // math.calcLeastSquare2nd(mesh, controls.U, controls.fU, gradU);
	// math.calcGradientFace(mesh, gradU, controls.U, controls.fU, 
					// controls.dUdx, controls.dUdy, controls.dUdz);
	// vector<vector<double>> gradV;
	// math.calcGaussGreen(mesh, controls.V, controls.fV, gradV);
	// // math.calcLeastSquare2nd(mesh, controls.V, controls.fV, gradV);
	// math.calcGradientFace(mesh, gradV, controls.V, controls.fV, 
					// controls.dVdx, controls.dVdy, controls.dVdz);
	// vector<vector<double>> gradW;
	// math.calcGaussGreen(mesh, controls.W, controls.fW, gradW);
	// // math.calcLeastSquare2nd(mesh, controls.W, controls.fW, gradW);
	// math.calcGradientFace(mesh, gradW, controls.W, controls.fW, 
					// controls.dWdx, controls.dWdy, controls.dWdz);
	
	mesh.cellsGradientVar[controls.P].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.U].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.V].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.W].resize(mesh.cells.size(),vector<double>(3,0.0));
	mesh.cellsGradientVar[controls.VF[0]].resize(mesh.cells.size(),vector<double>(3,0.0));
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1){
			int str = boundary.startFace;
			int end = str + boundary.nFaces;
			for(int i=str; i<end; ++i){
				auto& face = mesh.faces[i];
				for(int ii=0; ii<3; ++ii){
					mesh.cellsGradientVar[controls.P][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.U][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.V][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.W][face.owner][ii] = 0.0;
					mesh.cellsGradientVar[controls.VF[0]][face.owner][ii] = 0.0;
				}
			}
		}
	}
	{
		vector<double> dummyVec;
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.P, controls.fP, dummyVec, mesh.cellsGradientVar[controls.P]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.U, controls.fU, dummyVec, mesh.cellsGradientVar[controls.U]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.V, controls.fV, dummyVec, mesh.cellsGradientVar[controls.V]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.W, controls.fW, dummyVec, mesh.cellsGradientVar[controls.W]);
		math.calcLeastSquare(mesh, "cellVertex", "1st", "cell", 
			controls.VF[0], controls.fVF[0], dummyVec, mesh.cellsGradientVar[controls.VF[0]]);
	}
	
	
	
	
	
	

	
	for(auto& face : mesh.faces){
		
		vector<double> cFlux(controls.nEq,0.0);
		vector<double> vFlux(controls.nEq,0.0);
		
		vector<double> MFL;
		vector<double> MFR;
		
		double MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			MFL.push_back(face.varL[controls.fMF[ns]]);
			MFnSp += face.varL[controls.fMF[ns]];
		}
		MFL.push_back(1.0 - MFnSp);
		
		MFnSp = 0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			MFR.push_back(face.varR[controls.fMF[ns]]);
			MFnSp += face.varR[controls.fMF[ns]];
		}
		MFR.push_back(1.0 - MFnSp);
		
		// if(face.unitNormals[2]>0.0 && face.unitNormals[2]<1.0){
		// cout << face.unitNormals[2] << endl;
		// }
		
		
		// if(0
		// // face.getType() == SEMO_Types::INTERNAL_FACE 
		// // ||
		// // face.getType() == SEMO_Types::PROCESSOR_FACE
		// )
		{
			auto& cellLeft = mesh.cells[face.owner];

			double maxPhiLR = -1.e12;
			double minPhiLR = 1.e12;
			// double minPBarLR = 1.e12;
			for(auto j : cellLeft.stencil){
				auto& cellSten = mesh.cells[j];
				
				if(face.neighbour == j) continue;

				double shockDis = cellSten.var[controls.P] 
					+ 0.1 * cellSten.var[controls.Rho]*cellSten.var[controls.C]*cellSten.var[controls.C];

				maxPhiLR = max(maxPhiLR, shockDis);
				minPhiLR = min(minPhiLR, shockDis);
				
			}
			
			if(face.getType() == SEMO_Types::INTERNAL_FACE ){
				auto& cellRight = mesh.cells[face.neighbour];
				
				for(auto j : cellRight.stencil){
					auto& cellSten = mesh.cells[j];
					
					if(face.owner == j) continue;

					double shockDis = cellSten.var[controls.P] 
						+ 0.1 * cellSten.var[controls.Rho]*cellSten.var[controls.C]*cellSten.var[controls.C];

					maxPhiLR = max(maxPhiLR, shockDis);
					minPhiLR = min(minPhiLR, shockDis);
					
				}
			}
			
			double w2 = minPhiLR/maxPhiLR;
			// w2 = 1.0 - w2*w2;
			
			double minRhoC2 = min(
				face.varL[controls.fRho]*face.varL[controls.fC]*face.varL[controls.fC],
				face.varR[controls.fRho]*face.varR[controls.fC]*face.varR[controls.fC]);
			double pBarL = face.varL[controls.fP] + 0.1 * minRhoC2;
			double pBarR = face.varR[controls.fP] + 0.1 * minRhoC2;
			double w1 = min(pBarL/pBarR,pBarR/pBarL);
			// w1 = 1.0 - w1*w1*w1;
			
			double w3 = min(w1, minPhiLR);
			
			
			

			vector<double> MFhat;
			double dummy1;
			double dummy2;
			vector<double> dummy3;
			double RT = sqrt(face.varR[controls.fRho]/face.varL[controls.fRho]);
			double MFnSp = 0.0;
			for(int ns=0; ns<controls.nSp-1; ++ns){
				double MFhatTemp = (face.varL[controls.fMF[ns]]+RT*face.varR[controls.fMF[ns]])/(1.0+RT);
				MFhat.push_back(MFhatTemp);
				MFnSp += face.varL[controls.fMF[ns]];
			}
			MFhat.push_back(1.0 - MFnSp);
			
			double Phat = 0.5*(face.varL[controls.fP]+face.varR[controls.fP]);
			double Uhat = (face.varL[controls.fU]+RT*face.varR[controls.fU])/(1.0+RT);
			double Vhat = (face.varL[controls.fV]+RT*face.varR[controls.fV])/(1.0+RT);
			double What = (face.varL[controls.fW]+RT*face.varR[controls.fW])/(1.0+RT);
			double That = (face.varL[controls.fT]+RT*face.varR[controls.fT])/(1.0+RT);
			double fC = 0.0;
			
			this->getValuesFromEOSMF( species, 
				Phat, Uhat, Vhat, What, That, MFhat, dummy1, fC, dummy2, dummy3 );
				
			// this->calcFluxRoe(
			this->calcFluxSLAU2(
			// this->calcFluxRoeAMTP(
			// this->calcFluxAPRoe(
			// this->calcFluxHAUS(
				face.varL[controls.fP],
				face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW],
				face.varL[controls.fT], MFL,
				face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
				face.varR[controls.fP],
				face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW],
				face.varR[controls.fT], MFR,
				face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
				face.unitNormals,
				cFlux );
				
			// this->calcFluxAUSMPWP_N(
			// this->calcFluxRoeM_N(
				// face.varL[controls.fP],
				// face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW],
				// face.varL[controls.fT], MFL,
				// face.varL[controls.fRho], face.varL[controls.fC], face.varL[controls.fHt],
				// face.varR[controls.fP],
				// face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW],
				// face.varR[controls.fT], MFR,
				// face.varR[controls.fRho], face.varR[controls.fC], face.varR[controls.fHt],
				// face.unitNormals,
				// w1 , w2 , w3 , fC, controls.Uco, controls.timeStep, controls.Lch,
				// cFlux );
			
		}
		// else {

			// vector<double> nvec;
			// nvec.push_back(face.unitNormals[0]);
			// nvec.push_back(face.unitNormals[1]);
			// nvec.push_back(face.unitNormals[2]);
		
			// double UL = face.varL[controls.fU];
			// double VL = face.varL[controls.fV];
			// double WL = face.varL[controls.fW];
			// double PL = face.varL[controls.fP];
			// double rhoL = face.varL[controls.fRho];
			// double HtL = face.varL[controls.fHt];
			
			// double UR = face.varR[controls.fU];
			// double VR = face.varR[controls.fV];
			// double WR = face.varR[controls.fW];
			// double PR = face.varR[controls.fP];
			// double rhoR = face.varR[controls.fRho];
			// double HtR = face.varR[controls.fHt];

			// // properties of Left
			// double UnL = UL*nvec[0] + VL*nvec[1] + WL*nvec[2];
			// // properties of Right
			// double UnR = UR*nvec[0] + VR*nvec[1] + WR*nvec[2];
			
			// double UnF = face.wC*UnL + (1.0-face.wC)*UnR;
			// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// // UnF = face.var[controls.Un];
			// // }
			
			// double f1L;
			// double f1R;
			// if(UnF>=0.0){
				// f1L=rhoL*UnF;
				// f1R=0.0;
			// }
			// else{
				// f1L=0.0;
				// f1R=rhoR*UnF;
			// }

			// double PLR = 0.5*(PL+PR);
			
			// // flux
			// cFlux.clear();
			// cFlux.push_back( f1L    + f1R );
			// cFlux.push_back( f1L*UL + f1R*UR   + PLR*nvec[0] );
			// cFlux.push_back( f1L*VL + f1R*VR   + PLR*nvec[1] );
			// cFlux.push_back( f1L*WL + f1R*WR   + PLR*nvec[2] );
			// cFlux.push_back( f1L*HtL+ f1R*HtR );
			// for(int i=0; i<MFL.size()-1; ++i){
				// cFlux.push_back( f1L*MFL[i] + f1R*MFR[i] );
			// }
			
			
		// }
		


		
		this->calcViscousFlux(
			face.varL[controls.fmu],
			face.varL[controls.fU], face.varL[controls.fV], face.varL[controls.fW],
			face.varL[controls.fRho],
			face.varR[controls.fmu],
			face.varR[controls.fU], face.varR[controls.fV], face.varR[controls.fW],
			face.varR[controls.fRho],
			face.var[controls.dUdx], face.var[controls.dUdy], face.var[controls.dUdz],
			face.var[controls.dVdx], face.var[controls.dVdy], face.var[controls.dVdz],
			face.var[controls.dWdx], face.var[controls.dWdy], face.var[controls.dWdz],
			face.unitNormals,
			vFlux );
			
			
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(int nq=0; nq<controls.nEq; ++nq){
				residuals[face.owner][nq] += (cFlux[nq]-vFlux[nq])*face.area;
				residuals[face.neighbour][nq] -= (cFlux[nq]-vFlux[nq])*face.area;
			}
		}
		else{
			for(int nq=0; nq<controls.nEq; ++nq){
				residuals[face.owner][nq] += (cFlux[nq]-vFlux[nq])*face.area;
			}
		}
		
	}
	
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		
		// auto& cell = mesh.cells[i];
		
		// cell.var[controls.UDV[0]] = residuals[i][0];
		// cell.var[controls.UDV[1]] = residuals[i][1];
		// cell.var[controls.UDV[2]] = residuals[i][2];
		
	// }
		
		
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
}