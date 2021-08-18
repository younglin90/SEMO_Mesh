#include "build.h"
#include <cmath>
#include <array>




// void SEMO_Solvers_Builder::calcLinearSolverPETSc(
	// SEMO_Mesh_Builder& mesh,
	// SEMO_Controls_Builder& controls,
	// vector<vector<double>>& residuals){
	
	// SEMO_Utility_Math math;
	
	// int nEq = controls.nEq;
	
	
	// vector<vector<vector<double>>> linA;
	// // vector<vector<double>> linB;
	// vector<vector<vector<double>>> linAL; // dR/dQL
	// vector<vector<vector<double>>> linAR; // dR/dQR
	
	// vector<vector<double>> resiVar;
	
	// // construct A diagonal matrix
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		
		// double timeFrac = cell.var[controls.dtPseudo] / controls.timeStep;
		// vector<vector<double>> matrixA(nEq,vector<double>(nEq,0.0));
		// // vector<vector<double>> matrixB(nEq,vector<double>(nEq,0.0));
		
		// vector<vector<double>> matrixP(nEq,vector<double>(nEq,0.0));
		// vector<vector<double>> matrixT(nEq,vector<double>(nEq,0.0));
		
		// this->constructMatrix(cell, controls, matrixP, matrixT);
		
		// // A matrix (diagonal terms)
		// for(int in=0; in<nEq; ++in){
			// for(int jn=0; jn<nEq; ++jn){
				// // 2nd-order
				// matrixA[in][jn] = matrixP[in][jn] + 1.5 * timeFrac * matrixT[in][jn];
				
				// // matrixA[in][jn] += DsigmaApS[in][jn];
				
			// }
			
			// // matrixB[jn] = residuals[i][jn];
		// }
		
		// linA.push_back(matrixA);
		// // linB.push_back(matrixB);
		
		// vector<double> tempResi(nEq,0.0);
		// resiVar.push_back(tempResi);
		
	// }
	
	
	
	// // construct A off-diagonal matrix
	// int proc_num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// vector<vector<double>> matrixAL(nEq,vector<double>(nEq,0.0));
		// vector<vector<double>> matrixAR(nEq,vector<double>(nEq,0.0));
		
		// vector<vector<double>> jacALmat;
		// vector<vector<double>> jacARmat;
		
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// auto& cellL = mesh.cells[face.owner];
			// auto& cellR = mesh.cells[face.neighbour];
			
			// double rDtauPerVolL = cellL.var[controls.dtPseudo] / cellL.volume;
			// double rDtauPerVolR = cellR.var[controls.dtPseudo] / cellR.volume;
			
			// this->getNumericalFluxJacobianUpwind(
				// cellL,cellR,controls,face.unitNormals,"L",jacALmat);
				
			// this->getNumericalFluxJacobianUpwind(
				// cellL,cellR,controls,face.unitNormals,"R",jacARmat);
							
			// for(int in=0; in<nEq; ++in){
				// for(int jn=0; jn<nEq; ++jn){
					// matrixAL[in][jn] = (-rDtauPerVolR*jacALmat[in][jn]*face.area);
					// matrixAR[in][jn] = (+rDtauPerVolL*jacARmat[in][jn]*face.area);
				// }
			// }
			// linAL.push_back(matrixAL);
			// linAR.push_back(matrixAR);
		
		
			// for(int in=0; in<nEq; ++in){
				// for(int jn=0; jn<nEq; ++jn){
					// linA[face.owner][in][jn]     += (-(-rDtauPerVolL*jacALmat[in][jn]*face.area));
					// linA[face.neighbour][in][jn] += (-(+rDtauPerVolR*jacARmat[in][jn]*face.area));
				// }
			// }
			
			
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// for(int in=0; in<nEq; ++in){
				// for(int jn=0; jn<nEq; ++jn){
					// linA[face.owner][in][jn]     += (-(-rDtauPerVolL*jacALmat[in][jn]*face.area));
				// }
			// }

		// }
		// else{
			
			// for(int in=0; in<nEq; ++in){
				// for(int jn=0; jn<nEq; ++jn){
					// linA[face.owner][in][jn]     += (-(-rDtauPerVolL*jacALmat[in][jn]*face.area));
				// }
			// }
			
		// }
		
	
	// }
	
	// // linear solver : PETSc library
	
	// solvePETSc(mesh, resiVar, linA, linAL, linAR, residuals,
		// controls.solverU, controls.toleranceU, 
		// controls.relTolU, controls.preconditionerU,
		// controls.maxIterU);
	
	
	
// }





void SEMO_Solvers_Builder::calcLinearSolverLUSGS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	SEMO_Utility_Math math;
	
	int nEq = controls.nEq;
	
	//> forward sweep , 1 -> ncells
	for(int i=0; i<mesh.cells.size(); ++i){
	// for(int i=mesh.cells.size()-1; i>-1; --i){
		
		int own = mesh.cells[i].RCM;
		// int own = i;
		
		SEMO_Cell& cell = mesh.cells[own];
		
		int number_own = i;
		vector<vector<double>> DsigmaApS(nEq,vector<double>(nEq,0.0));
		vector<double> LorUmat(nEq,0.0);
		double rDtauPerVol = cell.var[controls.dtPseudo] / cell.volume;
		
		for(auto& j : cell.faces){
			auto& face = mesh.faces[j];
			
			vector<double> nvec;
			if( face.owner == own ){
				nvec.push_back(face.unitNormals[0]);
				nvec.push_back(face.unitNormals[1]);
				nvec.push_back(face.unitNormals[2]);
			}
			else{
				nvec.push_back(-face.unitNormals[0]);
				nvec.push_back(-face.unitNormals[1]);
				nvec.push_back(-face.unitNormals[2]);
			}
			
			
			vector<vector<double>> Amat;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				int ngb = face.neighbour;
				if( own == ngb ) ngb = face.owner;
				
				int number_ngb = mesh.cells[ngb].invRCM;
				// int number_ngb = ngb;
				
				if(number_ngb < number_own){
					
					this->getNumericalFluxJacobianUpwind(
						cell,mesh.cells[ngb],controls,nvec,"R",Amat);
					
					vector<double> sigmaAmS(nEq,0.0);
					for(int k=0; k<nEq; ++k){
						for(int l=0; l<nEq; ++l){
							sigmaAmS[k] += residuals[ngb][l] * Amat[k][l];
							
							// cout << Amat[k][l] << endl;
						}
						LorUmat[k] += rDtauPerVol * sigmaAmS[k] * face.area;
						
						// cout << sigmaAmS[k] << endl;
					}
				}
				
				this->getNumericalFluxJacobianUpwind(
					cell,mesh.cells[ngb],controls,nvec,"L",Amat);
				
			}
			else{
				this->getNumericalFluxJacobianUpwind(
					cell,cell,controls,nvec,"L",Amat);
				
			}
			
			for(int k=0; k<nEq; ++k){
				for(int l=0; l<nEq; ++l){
					DsigmaApS[k][l] += rDtauPerVol * Amat[k][l] * face.area;
				}
			}
			
			
		}
		
		
		
		
		vector<vector<double>> matrixP(nEq,vector<double>(nEq,0.0));
		vector<vector<double>> matrixT(nEq,vector<double>(nEq,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(nEq,vector<double>(nEq,0.0));
		
		// A matrix (diagonal terms)
		for(int in=0; in<nEq; ++in){
			for(int jn=0; jn<nEq; ++jn){
				double timeFrac = cell.var[controls.dtPseudo] / controls.timeStep;
				// 2nd-order
				matrixA[in][jn] = matrixP[in][jn] + 1.5 * timeFrac * matrixT[in][jn];
				
				matrixA[in][jn] += DsigmaApS[in][jn];
			}
		}
		
		
		// B matrix
		vector<double> vectorB;
		for(int j=0; j<nEq; ++j){
			double LURes = residuals[own][j]-LorUmat[j];
			vectorB.push_back(LURes);
		}
		
		math.GaussSeidelSOR(matrixA);
		
		for(int j=0; j<nEq; ++j){
			double updateResi = 0.0;
			for(int k=0; k<nEq; ++k){
				updateResi += matrixA[j][k]*vectorB[k];
				
				// cout << matrixA[j][k] << endl;
			}
			residuals[own][j] = updateResi;
		}
		
	}
	
	
	
	
	//> backward sweep , ncells -> 1
	for(int i=mesh.cells.size()-1; i>-1; --i){
		int own = mesh.cells[i].RCM;
		// int own = i;
		
		SEMO_Cell& cell = mesh.cells[own];
		
		int number_own = i;
		vector<vector<double>> DsigmaApS(nEq,vector<double>(nEq,0.0));
		vector<double> LorUmat(nEq,0.0);
		double rDtauPerVol = cell.var[controls.dtPseudo] / cell.volume;
		
		for(auto& j : cell.faces){
			auto& face = mesh.faces[j];
			
			vector<double> nvec;
			if( face.owner == own ){
				nvec.push_back(face.unitNormals[0]);
				nvec.push_back(face.unitNormals[1]);
				nvec.push_back(face.unitNormals[2]);
			}
			else{
				nvec.push_back(-face.unitNormals[0]);
				nvec.push_back(-face.unitNormals[1]);
				nvec.push_back(-face.unitNormals[2]);
			}
			
			
			vector<vector<double>> Amat;
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				int ngb = face.neighbour;
				if( own == ngb ) ngb = face.owner;
				
				int number_ngb = mesh.cells[ngb].invRCM;
				// int number_ngb = ngb;
				
				if(number_ngb > number_own){
					
					this->getNumericalFluxJacobianUpwind(
						cell,mesh.cells[ngb],controls,nvec,"R",Amat);
					
					vector<double> sigmaAmS(nEq,0.0);
					for(int k=0; k<nEq; ++k){
						for(int l=0; l<nEq; ++l){
							sigmaAmS[k] += residuals[ngb][l] * Amat[k][l];
						}
						LorUmat[k] += rDtauPerVol * sigmaAmS[k] * face.area;
					}
				}
				
				this->getNumericalFluxJacobianUpwind(
					cell,mesh.cells[ngb],controls,nvec,"L",Amat);
				
			}
			else{
				this->getNumericalFluxJacobianUpwind(
					cell,cell,controls,nvec,"L",Amat);
				
			}
			
			for(int k=0; k<nEq; ++k){
				for(int l=0; l<nEq; ++l){
					DsigmaApS[k][l] += rDtauPerVol * Amat[k][l] * face.area;
				}
			}
			
			
		}
		
		
		vector<vector<double>> matrixP(nEq,vector<double>(nEq,0.0));
		vector<vector<double>> matrixT(nEq,vector<double>(nEq,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(nEq,vector<double>(nEq,0.0));
		
		// A matrix (diagonal terms)
		for(int in=0; in<nEq; ++in){
			for(int jn=0; jn<nEq; ++jn){
				double timeFrac = cell.var[controls.dtPseudo] / controls.timeStep;
				// 2nd-order
				matrixA[in][jn] = matrixP[in][jn] + 1.5 * timeFrac * matrixT[in][jn];
				
				matrixA[in][jn] += DsigmaApS[in][jn];
			}
		}
		
		
		// B matrix
		vector<double> vectorB;
		for(int j=0; j<nEq; ++j){
			double LURes = LorUmat[j];
			vectorB.push_back(LURes);
		}
		
		math.GaussSeidelSOR(matrixA);
		
		for(int j=0; j<nEq; ++j){
			double updateResi = 0.0;
			for(int k=0; k<nEq; ++k){
				updateResi += matrixA[j][k]*vectorB[k];
			}
			residuals[own][j] = residuals[own][j] - updateResi;
		}
		
	}
	
	

	
}




void SEMO_Solvers_Builder::getNumericalFluxJacobianUpwind(
	SEMO_Cell& cellL, SEMO_Cell& cellR, SEMO_Controls_Builder& controls,
	vector<double>& nvec,
	string LorR,
	vector<vector<double>>& fluxJac
	){
	
	// left variables
	double PL=cellL.var[controls.P]; 
	double UL=cellL.var[controls.U]; 
	double VL=cellL.var[controls.V]; 
	double WL=cellL.var[controls.W]; 
	double TL=cellL.var[controls.T]; 
	vector<double> YL;
	for(int ns=0; ns<controls.nSp; ++ns){
		YL.push_back(cellL.var[controls.MF[ns]]);
	}
	double DRDPL = cellL.var[controls.dRhoDP];
	double DHDPL = cellL.var[controls.dHtDP];
	double DRDTL = cellL.var[controls.dRhoDT];
	double DHDTL = cellL.var[controls.dHtDT];
	vector<double> DRDYL, DHDYL;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		DRDYL.push_back(cellL.var[controls.dRhoDMF[ns]]);
		DHDYL.push_back(cellL.var[controls.dHtDMF[ns]]);
	}
	double RhoL = cellL.var[controls.Rho];
	double HtL = cellL.var[controls.Ht];
	double CL = cellL.var[controls.C];
	
	double UnL=UL*nvec[0]+VL*nvec[1]+WL*nvec[2];
	
	// right variables
	double PR=cellR.var[controls.P]; 
	double UR=cellR.var[controls.U]; 
	double VR=cellR.var[controls.V]; 
	double WR=cellR.var[controls.W]; 
	double TR=cellR.var[controls.T]; 
	vector<double> YR;
	for(int ns=0; ns<controls.nSp; ++ns){
		YR.push_back(cellR.var[controls.MF[ns]]);
	}
	double DRDPR = cellR.var[controls.dRhoDP];
	double DHDPR = cellR.var[controls.dHtDP];
	double DRDTR = cellR.var[controls.dRhoDT];
	double DHDTR = cellR.var[controls.dHtDT];
	vector<double> DRDYR, DHDYR;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		DRDYR.push_back(cellR.var[controls.dRhoDMF[ns]]);
		DHDYR.push_back(cellR.var[controls.dHtDMF[ns]]);
	}
	double RhoR = cellR.var[controls.Rho];
	double HtR = cellR.var[controls.Ht];
	double CR = cellR.var[controls.C];
	
	double UnR=UR*nvec[0]+VR*nvec[1]+WR*nvec[2];
	
	// roe-averaged
	double RT = sqrt(RhoL/RhoR);
	double Rhohat = 0.5*(RhoL+RhoR);
    double Chat= 0.5*(CL+CR);
	double DHDThat = (RT*DHDTL+DHDTR)/(RT+1.0);
	double DRDThat = (RT*DRDTL+DRDTR)/(RT+1.0);
	double DHDPhat = (RT*DHDPL+DHDPR)/(RT+1.0);
	double DRDPhat = (RT*DRDPL+DRDPR)/(RT+1.0);
	double Unhat = 0.5*(UnL+UnR);
	double Uhat = (RT*UL+UR)/(RT+1.0);
	double Vhat = (RT*VL+VR)/(RT+1.0);
	double What = (RT*WL+WR)/(RT+1.0);
	double Hthat = (RT*HtL+HtR)/(RT+1.0);
	vector<double> Yhat;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		Yhat.push_back((RT*YL[ns]+YR[ns])/(RT+1.0));
	}
	
	double fgL, fgR;
	if( LorR == "L" ) {
		fgL = 1.0; fgR = 0.0;
	}
	else if( LorR == "R" ) {
		fgL = 0.0; fgR = 1.0;
	}
	else{
		cout << "ERROR!!!!!!" << endl;
	}
	
	double ULR = UL*fgL + UR*fgR; 
	double VLR = VL*fgL + VR*fgR; 
	double WLR = WL*fgL + WR*fgR;
	double HtLR = HtL*fgL + HtR*fgR; 
	double RhoLR = RhoL*fgL + RhoR*fgR; 
	double DRDPLR = DRDPL*fgL + DRDPR*fgR; 
	double DHDPLR = DHDPL*fgL + DHDPR*fgR;
	double DRDTLR = DRDTL*fgL + DRDTR*fgR; 
	double DHDTLR = DHDTL*fgL + DHDTR*fgR;
	double UnLR = UnL*fgL + UnR*fgR;
	vector<double> YLR, DRDYLR, DHDYLR;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		YLR.push_back(YL[ns]*fgL + YR[ns]*fgR);
		DRDYLR.push_back(DRDYL[ns]*fgL + DRDYR[ns]*fgR);
		DHDYLR.push_back(DHDYL[ns]*fgL + DHDYR[ns]*fgR);
	}
	
	
	
	vector<vector<double>> fluxLRdot(controls.nEq,vector<double>(controls.nEq,0.0));
	
    fluxLRdot[0][0] = DRDPLR*UnLR;
	fluxLRdot[1][0] = DRDPLR*UnLR*ULR + nvec[0];
	fluxLRdot[2][0] = DRDPLR*UnLR*VLR + nvec[1];
	fluxLRdot[3][0] = DRDPLR*UnLR*WLR + nvec[2];
	fluxLRdot[4][0] = DRDPLR*UnLR*HtLR + RhoLR*UnLR*DHDPLR;
	
    fluxLRdot[0][1] = RhoLR*nvec[0];
	fluxLRdot[1][1] = RhoLR*nvec[0]*ULR + RhoLR*UnLR;
	fluxLRdot[2][1] = RhoLR*nvec[0]*VLR;
	fluxLRdot[3][1] = RhoLR*nvec[0]*WLR;
	fluxLRdot[4][1] = RhoLR*nvec[0]*HtLR + RhoLR*UnLR*ULR;

    fluxLRdot[0][2] = RhoLR*nvec[1];
	fluxLRdot[1][2] = RhoLR*nvec[1]*ULR;
	fluxLRdot[2][2] = RhoLR*nvec[1]*VLR + RhoLR*UnLR;
	fluxLRdot[3][2] = RhoLR*nvec[1]*WLR;
	fluxLRdot[4][2] = RhoLR*nvec[1]*HtLR + RhoLR*UnLR*VLR;
	
    fluxLRdot[0][3] = RhoLR*nvec[2];
	fluxLRdot[1][3] = RhoLR*nvec[2]*ULR;
	fluxLRdot[2][3] = RhoLR*nvec[2]*VLR;
	fluxLRdot[3][3] = RhoLR*nvec[2]*WLR + RhoLR*UnLR;
	fluxLRdot[4][3] = RhoLR*nvec[2]*HtLR + RhoLR*UnLR*WLR;
	
    fluxLRdot[0][4] = DRDTLR*UnLR;
	fluxLRdot[1][4] = DRDTLR*UnLR*ULR;
	fluxLRdot[2][4] = DRDTLR*UnLR*VLR;
	fluxLRdot[3][4] = DRDTLR*UnLR*WLR;
    fluxLRdot[4][4] = DRDTLR*UnLR*HtLR + RhoLR*UnLR*DHDTLR;


	for(int ns=0; ns<controls.nSp-1; ++ns){
        fluxLRdot[5+ns][0] = DRDPLR*YLR[ns]*UnLR;
        fluxLRdot[5+ns][1] = RhoLR*YLR[ns]*nvec[0];
        fluxLRdot[5+ns][2] = RhoLR*YLR[ns]*nvec[1];
        fluxLRdot[5+ns][3] = RhoLR*YLR[ns]*nvec[2];
        fluxLRdot[5+ns][4] = DRDTLR*YLR[ns]*UnLR;
		
		for(int ns2=0; ns2<controls.nSp-1; ++ns2){
            fluxLRdot[5+ns][5+ns2] = DRDYLR[ns2]*YLR[ns]*UnLR;
		}
		fluxLRdot[5+ns][5+ns] += RhoLR*UnLR;

        fluxLRdot[0][5+ns] = DRDYLR[ns]*UnLR;
        fluxLRdot[1][5+ns] = DRDYLR[ns]*UnLR*ULR;
        fluxLRdot[2][5+ns] = DRDYLR[ns]*UnLR*VLR;
        fluxLRdot[3][5+ns] = DRDYLR[ns]*UnLR*WLR;
        fluxLRdot[4][5+ns] = DRDYLR[ns]*UnLR*HtLR + RhoLR*UnLR*DHDYLR[ns];
	}
	
	
	
	// Mmul = dsqrt(uhat*uhat+vhat*vhat+what*what)/chat
	// cmul = dmin1(dmax1(Mmul,1.d-5),1.d0)*chat
	// theta = dmin1(dmax1(Mmul**2.d0,1.d-5),1.d0)
	// Umul = 0.5d0*(1.d0+theta)*unhat
	// tdash = dmin1(Mmul**2.d0,1.d0)
	// Udash = 0.5d0*(1.d0+tdash)*unhat
	// cdash = 0.5d0*(4.d0*chat*chat*tdash+(1.d0-tdash)**2.d0*unhat**2.d0)
	
	
	//> Roe scheme
	double dissUp = (Chat-abs(Unhat))/Rhohat/Chat/Chat;
	double dissUu = Unhat/Chat;
	double dissPp = Unhat/Chat;
	double dissPu = (Chat-abs(Unhat))*Rhohat;
	double AD = abs(Unhat);
	
	vector<vector<double>> dissimat(controls.nEq,vector<double>(controls.nEq,0.0));
	
	dissimat[0][0] = AD*DRDPLR + dissUp*Rhohat;
	dissimat[1][0] = AD*DRDPLR*ULR + dissUp*Rhohat*Uhat + dissPp*nvec[0];
	dissimat[2][0] = AD*DRDPLR*VLR + dissUp*Rhohat*Vhat + dissPp*nvec[1];
	dissimat[3][0] = AD*DRDPLR*WLR + dissUp*Rhohat*What + dissPp*nvec[2];
	dissimat[4][0] = AD*(DRDPLR*HtLR+RhoLR*DHDPLR-1.0) + dissUp*Rhohat*Hthat + dissPp*Unhat;
	
	dissimat[0][1] = dissUu*nvec[0]*Rhohat;
	dissimat[1][1] = AD*RhoLR + dissUu*nvec[0]*Rhohat*Uhat + dissPu*nvec[0]*nvec[0];
	dissimat[2][1] = dissUu*nvec[0]*Rhohat*Vhat + dissPu*nvec[0]*nvec[1];
	dissimat[3][1] = dissUu*nvec[0]*Rhohat*What + dissPu*nvec[0]*nvec[2];
	dissimat[4][1] = AD*RhoLR*ULR + dissUu*nvec[0]*Rhohat*Hthat + dissPu*nvec[0]*Unhat;
	
	dissimat[0][2] = dissUu*nvec[1]*Rhohat;
	dissimat[1][2] = dissUu*nvec[1]*Rhohat*Uhat + dissPu*nvec[1]*nvec[0];
	dissimat[2][2] = AD*RhoLR + dissUu*nvec[1]*Rhohat*Vhat + dissPu*nvec[1]*nvec[1];
	dissimat[3][2] = dissUu*nvec[1]*Rhohat*What + dissPu*nvec[1]*nvec[2];
	dissimat[4][2] = AD*RhoLR*VLR + dissUu*nvec[1]*Rhohat*Hthat + dissPu*nvec[1]*Unhat;
	
	dissimat[0][3] = dissUu*nvec[2]*Rhohat;
	dissimat[1][3] = dissUu*nvec[2]*Rhohat*Uhat + dissPu*nvec[2]*nvec[0];
	dissimat[2][3] = dissUu*nvec[2]*Rhohat*Vhat + dissPu*nvec[2]*nvec[1];
	dissimat[3][3] = AD*RhoLR + dissUu*nvec[2]*Rhohat*What + dissPu*nvec[2]*nvec[2];
	dissimat[4][3] = AD*RhoLR*WLR + dissUu*nvec[2]*Rhohat*Hthat + dissPu*nvec[2]*Unhat;
	
	dissimat[0][4] = AD*DRDTLR;
	dissimat[1][4] = AD*DRDTLR*ULR;
	dissimat[2][4] = AD*DRDTLR*VLR;
	dissimat[3][4] = AD*DRDTLR*WLR;
	dissimat[4][4] = AD*(DRDTLR*HtLR+RhoLR*DHDTLR);
	

	// for(int ns=0; ns<controls.nSp-1; ++ns){
		
		// dissimat[5+ns][0] = AD*DRDPLR*YLR[ns] + dissUp*Rhohat*Yhat[ns];
		// dissimat[5+ns][1] = dissUu*nvec[0]*Rhohat*Yhat[ns];
		// dissimat[5+ns][2] = dissUu*nvec[1]*Rhohat*Yhat[ns];
		// dissimat[5+ns][3] = dissUu*nvec[2]*Rhohat*Yhat[ns];
		// dissimat[5+ns][4] = AD*DHDTLR*YLR[ns];

		// for(int ns2=0; ns2<controls.nSp-1; ++ns2){
			// dissimat[5+ns][5+ns2] = AD*DRDYLR[ns2]*YLR[ns];
		// }
		// dissimat[5+ns][5+ns] += AD*RhoLR;
		
		// dissimat[0][5+ns] = AD*DRDYLR[ns];
		// dissimat[1][5+ns] = AD*DRDYLR[ns]*ULR;
		// dissimat[2][5+ns] = AD*DRDYLR[ns]*VLR;
		// dissimat[3][5+ns] = AD*DRDYLR[ns]*WLR;
		// dissimat[4][5+ns] = AD*(DRDYLR[ns]*HtLR + RhoLR*DHDYLR[ns]);
	// }
	
	
	
	fluxJac.clear();
	fluxJac.resize(controls.nEq,vector<double>(controls.nEq,0.0));
	if( LorR == "L" ) {
		for(int i=0; i<controls.nEq; ++i){
			for(int j=0; j<controls.nEq; ++j){
				fluxJac[i][j] = 0.5*fluxLRdot[i][j] + 0.5*dissimat[i][j];
			}
				
			// // viscous
			// double viscEigen = 
				// (cellL.var[controls.mu] + cellR.var[controls.mu])/Rhohat/
					// pow(0.5*(cellL.volume+cellR.volume),0.3333);
			// fluxJac[i][i] += viscEigen;
		}
	}
	else if( LorR == "R" ) {
		for(int i=0; i<controls.nEq; ++i){
			for(int j=0; j<controls.nEq; ++j){
				fluxJac[i][j] = 0.5*fluxLRdot[i][j] - 0.5*dissimat[i][j];
				
			}
			
			// // viscous
			// double viscEigen = 
				// (cellL.var[controls.mu] + cellR.var[controls.mu])/Rhohat/
					// pow(0.5*(cellL.volume+cellR.volume),0.3333);
			// fluxJac[i][i] += viscEigen;
		}
	}
	else{
		cout << "ERROR!!!!!!" << endl;
	}
	


}






void SEMO_Solvers_Builder::calcLinearSolver(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	SEMO_Utility_Math math;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		int matSize = controls.nEq;
		
		
		vector<vector<double>> matrixP(matSize,vector<double>(matSize,0.0));
		vector<vector<double>> matrixT(matSize,vector<double>(matSize,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(matSize,vector<double>(matSize,0.0));
		
		for(int in=0; in<matSize; ++in){
			for(int jn=0; jn<matSize; ++jn){
					
				// 1st-order euler
				// matrixA[in][jn] = matrixP[in][jn] + 
					// cell.var[controls.dtPseudo] / controls.timeStep * matrixT[in][jn];
					
				// // 2nd-order
				matrixA[in][jn] = matrixP[in][jn] + 
					1.5 * cell.var[controls.dtPseudo] / controls.timeStep * matrixT[in][jn];
					
			}
		}
		
		vector<double> vectorB;
		for(auto& j : residuals[i]){
			vectorB.push_back(j);
		}
		
		
		math.GaussSeidelSOR(matrixA);
		
		
		for(int j=0; j<residuals[i].size(); ++j){
			double updateResi = 0.0;
			for(int k=0; k<residuals[i].size(); ++k){
				updateResi += matrixA[j][k]*vectorB[k];
			}
			residuals[i][j] = updateResi;
		}
		
		
	}
	

	
}




void SEMO_Solvers_Builder::calcSingleLinearSolver(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	SEMO_Utility_Math math;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		
		int matSize = controls.nEq;
		
		
		vector<vector<double>> matrixP(matSize,vector<double>(matSize,0.0));
		vector<vector<double>> matrixT(matSize,vector<double>(matSize,0.0));
		
		this->constructMatrix(cell, controls, matrixP, matrixT);
		
		
		vector<vector<double>> matrixA(matSize,vector<double>(matSize,0.0));
		
		for(int in=0; in<matSize; ++in){
			for(int jn=0; jn<matSize; ++jn){
				matrixA[in][jn] = matrixT[in][jn];
				// cout << matrixA[in][jn] << endl;
					
			}
				// cout << matrixA[0][0] << " " << matrixA[0][1] << " " << matrixA[0][2] << " " << matrixA[0][3] << " " << matrixA[0][4] << endl;
		}
		
		vector<double> vectorB;
		for(auto& j : residuals[i]){
			vectorB.push_back(j);
		}
		
		
		math.GaussSeidelSOR(matrixA);
		
		
		for(int j=0; j<residuals[i].size(); ++j){
			double updateResi = 0.0;
			for(int k=0; k<residuals[i].size(); ++k){
				updateResi += matrixA[j][k]*vectorB[k];
				
			}
			residuals[i][j] = updateResi;
		}
		
		
	}
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
}






void SEMO_Solvers_Builder::constructMatrix(
	SEMO_Cell& cell, SEMO_Controls_Builder& controls,
	vector<vector<double>>& matrixP,
	vector<vector<double>>& matrixT){


	int N = matrixP.size();

	double DRDP = cell.var[controls.dRhoDP];
	double DHDP = cell.var[controls.dHtDP];
	double DRDT = cell.var[controls.dRhoDT];
	double DHDT = cell.var[controls.dHtDT];
	vector<double> DRDY, DHDY;
	for(int ns=0; ns<controls.nSp-1; ++ns){
		DRDY.push_back(cell.var[controls.dRhoDMF[ns]]);
		DHDY.push_back(cell.var[controls.dHtDMF[ns]]);
	}
	// cell.var[controls.dRhoDMF[ns]]
	double Rho = cell.var[controls.Rho];
	double U = cell.var[controls.U];
	double V = cell.var[controls.V];
	double W = cell.var[controls.W];
	double T = cell.var[controls.T];
	vector<double> Y;
	for(int ns=0; ns<controls.nSp; ++ns){
		Y.push_back(cell.var[controls.MF[ns]]);
	}
	double Ht = cell.var[controls.Ht];
	
	
	matrixT[0][0] = DRDP;
	matrixT[1][0] = DRDP*U;
	matrixT[2][0] = DRDP*V;
	matrixT[3][0] = DRDP*W;
	matrixT[4][0] = DRDP*Ht + Rho*DHDP - 1.0;
	
	// cout << DRDP << " " << Ht << " " << Rho << " " << DHDP << endl;
	
	matrixT[1][1] = Rho;
	matrixT[4][1] = Rho*U;
	
	matrixT[2][2] = Rho;
	matrixT[4][2] = Rho*V;
	
	matrixT[3][3] = Rho;
	matrixT[4][3] = Rho*W;
	
	matrixT[0][4] = DRDT;
	matrixT[1][4] = DRDT*U;
	matrixT[2][4] = DRDT*V;
	matrixT[3][4] = DRDT*W;
	matrixT[4][4] = DRDT*Ht + Rho*DHDT;
	
	// cout << DRDT*Ht + Rho*DHDT << endl;
	
	for(int ns=0; ns<controls.nSp-1; ++ns){
		int nq = 5 + ns;
		matrixT[nq][0] = DRDP*Y[ns];
		matrixT[nq][4] = DRDT*Y[ns];
		
		// cout << DRDT << " " << Y[ns] << endl;
		
		for(int ns2=0; ns2<controls.nSp-1; ++ns2){
			int nq2 = 5 + ns2;
			matrixT[nq][nq2] = Y[ns]*DRDY[ns2];
		}
		matrixT[nq][nq] += Rho;
		
		matrixT[0][nq] = DRDY[ns];
		matrixT[1][nq] = DRDY[ns]*U;
		matrixT[2][nq] = DRDY[ns]*V;
		matrixT[3][nq] = DRDY[ns]*W;
		matrixT[4][nq] = DRDY[ns]*Ht + Rho*DHDY[ns];
	}
	
	
		// cout << matrixT[4][5] << endl;
	

	for(int i=0; i<N; ++i){
		for(int j=0; j<N; ++j){
			matrixP[i][j] = matrixT[i][j];
		}
	}
	
	// double beta = 0.5*DRDP;
	// double beta = 1.0/pow(cell.var[controls.Ur],2.0);
	double beta = 1.0/pow(cell.var[controls.Ur],2.0) 
					+ DRDP - 1.0/pow(cell.var[controls.C],2.0);
	
	matrixP[0][0] = beta;
	matrixP[1][0] = beta*U;
	matrixP[2][0] = beta*V;
	matrixP[3][0] = beta*W;
	matrixP[4][0] = beta*Ht + Rho*DHDP - 1.0;
	
	for(int ns=0; ns<controls.nSp-1; ++ns){
		int nq = 5 + ns;
		matrixP[nq][0] = beta*Y[ns];
	}
	




}