#include "build.h"
#include <cmath>
#include <array>


void SEMO_Solvers_Builder::hybridBased(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	
	
	this->calcIncomCellEOSVF(mesh, controls, species);
	
	
	
	// save old value 
	for(auto& cell : mesh.cells){
		cell.var[controls.oldP] = cell.var[controls.P];
		cell.var[controls.oldU] = cell.var[controls.U];
		cell.var[controls.oldV] = cell.var[controls.V];
		cell.var[controls.oldW] = cell.var[controls.W];
		cell.var[controls.oldT] = cell.var[controls.T];
		for(int i=0; i<controls.nSp-1; ++i){
			cell.var[controls.oldVF[i]] = cell.var[controls.VF[i]];
			cell.var[controls.oldMF[i]] = cell.var[controls.MF[i]];
		}
		cell.var[controls.oldRho] = cell.var[controls.Rho];
		
		
		for(int i=0; i<controls.nEq; ++i){
			cell.var[controls.Qm[i]] = cell.var[controls.Qn[i]];
		}
		
		cell.var[controls.Qn[0]] = cell.var[controls.Rho];
		cell.var[controls.Qn[1]] = cell.var[controls.Rho] * cell.var[controls.U];
		cell.var[controls.Qn[2]] = cell.var[controls.Rho] * cell.var[controls.V];
		cell.var[controls.Qn[3]] = cell.var[controls.Rho] * cell.var[controls.W];
		cell.var[controls.Qn[4]] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
		for(int ns=0; ns<controls.nSp-1; ++ns){
			cell.var[controls.Qn[5+ns]] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
		}
		
	}
	
	
	
	// // incompressible pressure based
	bool boolIncomPressureBased = true;
	
	// if(boolIncomPressureBased){
	
		// vector<double> norm0(controls.nEq,0.0);
		// vector<double> normNM(controls.nEq,0.0);
		
		// controls.iterPBs = 0;
		// while(controls.iterPBs<controls.iterPBsMax){
			
			// // this->calcIncomRealTimeStep(mesh, controls);
			// this->calcUnderRelaxationFactors(controls);
			
			// vector<double> residualsM(controls.nEq,0.0);
			// vector<double> residualsP(controls.nEq,0.0);
			// vector<double> residualsVF(controls.nEq,0.0);
			
			// vector<double> linAD(mesh.cells.size(),0.0);
			// vector<double> linAOD(mesh.cells.size(),0.0);
			
		
			// controls.iterMom = 0;
			// while(controls.iterMom<controls.iterMomMax){
				
				// this->setIncomValuesLeftRightFace(mesh, controls, species);
				// // this->setCompValuesLeftRightFace(mesh, controls, species);
				
				// this->calcMomentumEqs(mesh, controls, linAD, linAOD, residualsM);
				
				// ++controls.iterMom;
				
			// }

			// controls.iterPre = 0;
			// while(controls.iterPre<controls.iterPreMax){
				
				// this->setIncomValuesLeftRightFace(mesh, controls, species);
				// // this->setCompValuesLeftRightFace(mesh, controls, species);
				
				// this->calcPressureEq(mesh, controls, linAD, linAOD, residualsP);
				
				
				
				// ++controls.iterPre;
			// }
			// // MPI_Barrier(MPI_COMM_WORLD);
			// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			
			// // if(controls.iterPBs == controls.iterPBsMax-1){
			// controls.iterVof = 0;
			// while(controls.iterVof<controls.iterVofMax){
				
				// // this->setIncomValuesLeftRightFace(mesh, controls, species);
				// this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
				// // this->setCompValuesLeftRightFaceWithNVD(mesh, controls, species);
				
				// this->calcVolfracEq(mesh, controls, residualsVF);
				
				
				// this->calcIncomCellEOSVF(mesh, controls, species);
				

				// // for(auto& cell : mesh.cells){
					// // cell.var[controls.Rho] = 
						// // 1000.0*cell.var[controls.VF[0]] + 
						// // 1.0*(1.0-cell.var[controls.VF[0]]);
				// // }
				
				
				// ++controls.iterVof;
			// }
			// // }
			
			
			// vector<double> residuals(controls.nEq,0.0);
			// for(int i=0; i<controls.nEq; ++i){
				// residuals[i] += residualsM[i];
				// residuals[i] += residualsP[i];
				// residuals[i] += residualsVF[i];
			// }

			// vector<double> norm;
			// this->calcNormResiduals(mesh, controls, residuals, norm);
			
			// if(rank==0) {
				// cout << "|  └ PIMPLE step = " << controls.iterPBs << " | ";
				// cout.precision(3);
				// for(int i=0; i<controls.nEq; ++i){
					// // if(controls.iterPBs==0) norm0[i] = norm[i];
					// // cout << scientific << abs(norm[i] - normNM[i])/(norm0[i]+1.e-200) << " | ";
					// // normNM[i] = norm[i];
					// if(controls.iterPBs==0) norm0[i] = norm[i];
					// cout << scientific << norm[i]/(norm0[i]+1.e-200) << " | ";
				// }
				// cout.unsetf(ios::scientific);
				// cout << endl;
				
			// }
			
			
			// ++controls.iterPBs;
			
		// }
		
	// }
	
	
	
	
	
	
	
	// compressible density based
	bool boolComDensityBased = true;
	
	if(boolComDensityBased){
		controls.iterPseudo = 0;
		while(controls.iterPseudo<controls.iterPseudoMax){
			
			if(rank==0) cout << "|  └ pseudo step = " << controls.iterPseudo << " | ";
			
			this->calcUnderRelaxationFactorsDualTime(controls);
		
			this->calcPseudoTimeStep(mesh, controls);
			
			this->setCompValuesLeftRightFace(mesh, controls, species);
			
			vector<vector<double>> residuals(
				mesh.cells.size(),vector<double>(controls.nEq,0.0));
		
			this->calcRHS(mesh, controls, residuals);

			this->calcLinearSolver(mesh, controls, residuals);

			vector<double> norm;
			this->calcNormResiduals(mesh, controls, residuals, norm);
			
			if(rank==0) {
				cout.precision(3);
				for(int i=0; i<controls.nEq; ++i){
					cout << scientific << norm[i] << " | ";
				}
				cout.unsetf(ios::scientific);
				cout << endl;
			}
			
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
			this->updateValues(mesh, controls, residuals);
			
			this->calcCellEOSMF(mesh, controls, species);
			
			++controls.iterPseudo;
			
		}
		
	}
	
}