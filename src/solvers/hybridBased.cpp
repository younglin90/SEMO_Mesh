#include "build.h"
#include <cmath>
#include <array>


void SEMO_Solvers_Builder::hybridBased(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	
	vector<vector<double>> tmpQm;
	
	if( controls.iterReal == 0 ){

		// this->calcIncomCellEOSVF(mesh, controls, species);
		this->calcCellEOSVF(mesh, controls, species);
		this->calcCellTransport(mesh, controls, species);
	
		
		for(auto& cell : mesh.cells){
			
			cell.var[controls.Qn[0]] = cell.var[controls.Rho];
			cell.var[controls.Qn[1]] = cell.var[controls.Rho] * cell.var[controls.U];
			cell.var[controls.Qn[2]] = cell.var[controls.Rho] * cell.var[controls.V];
			cell.var[controls.Qn[3]] = cell.var[controls.Rho] * cell.var[controls.W];
			cell.var[controls.Qn[4]] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			for(int ns=0; ns<controls.nSp-1; ++ns){
				cell.var[controls.Qn[5+ns]] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
			}
			
			vector<double> tmpVec;
			for(int i=0; i<controls.nEq; ++i){
				tmpVec.push_back(cell.var[controls.Qn[i]]);
			}
			tmpQm.push_back(tmpVec);
			
			cell.var[controls.Qm[0]] = cell.var[controls.U];
			cell.var[controls.Qm[1]] = cell.var[controls.V];
			cell.var[controls.Qm[2]] = cell.var[controls.W];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qm[3+i]] = cell.var[controls.VF[i]];
			}
			
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
			
		}
		
		
	}
	else{
		for(auto& cell : mesh.cells){
			vector<double> tmpVec;
			for(int i=0; i<controls.nEq; ++i){
				tmpVec.push_back(cell.var[controls.Qn[i]]);
			}
			tmpQm.push_back(tmpVec);
			
			cell.var[controls.Qn[0]] = cell.var[controls.Rho];
			cell.var[controls.Qn[1]] = cell.var[controls.Rho] * cell.var[controls.U];
			cell.var[controls.Qn[2]] = cell.var[controls.Rho] * cell.var[controls.V];
			cell.var[controls.Qn[3]] = cell.var[controls.Rho] * cell.var[controls.W];
			cell.var[controls.Qn[4]] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			for(int ns=0; ns<controls.nSp-1; ++ns){
				cell.var[controls.Qn[5+ns]] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
			}
			
			cell.var[controls.Qm[0]] = cell.var[controls.oldU];
			cell.var[controls.Qm[1]] = cell.var[controls.oldV];
			cell.var[controls.Qm[2]] = cell.var[controls.oldW];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qm[3+i]] = cell.var[controls.oldVF[i]];
			}
			
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
			
		}
	}
	
	
	
	
	
	
	
	//===================================================
	// // incompressible pressure based
	bool boolIncomPressureBased = true;
	
	if(boolIncomPressureBased){


		this->calcIncomRealTimeStep(mesh, controls);
		if(rank==0) {
			cout << " | timeStep = " << controls.timeStep;
		}
		double corantNum = -10.0;
		this->calcCorantNumberForPrint(mesh, controls, corantNum);
		if(rank==0) {
			cout << " | corantNumber = " << corantNum << endl;
		}
		
		
		vector<double> norm0(controls.nEq,0.0);
		vector<double> normNM(controls.nEq,0.0);
		
		controls.iterPBs = 0;
		while(controls.iterPBs<controls.iterPBsMax){
			
			// this->calcIncomRealTimeStep(mesh, controls);
			this->calcUnderRelaxationFactors(controls);
			
			vector<double> residualsM(controls.nEq,0.0);
			vector<double> residualsP(controls.nEq,0.0);
			vector<double> residualsVF(controls.nEq,0.0);
			
			vector<double> linAD(mesh.cells.size(),0.0);
			vector<double> linAOD(mesh.cells.size(),0.0);
			vector<double> linAFL(mesh.faces.size(),0.0);
			vector<double> linAFR(mesh.faces.size(),0.0);

			controls.iterMom = 0;
			while(controls.iterMom<controls.iterMomMax){
				double tmppreVelURF = controls.momVelURF;
				if(
				controls.iterPBsMax == 1 &&
				controls.iterMom == controls.iterMomMax-1
				){
					controls.momVelURF = 1.0;
				}
				
				this->setCompValuesLeftRightFaceWithReconPVForSegregated(mesh, controls, species);
				// this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
				
				this->calcMomentumEqs(mesh, controls, species, linAD, linAOD, linAFL, linAFR, residualsM);
				
				controls.momVelURF = tmppreVelURF;
				
				++controls.iterMom;
				
			}
			


			controls.iterPre = 0;
			while(controls.iterPre<controls.iterPreMax){
				double tmpprePreURF = controls.prePreURF;
				double tmppreVelURF = controls.preVelURF;
				if(
				controls.iterPBsMax == 1 &&
				controls.iterPre == controls.iterPreMax-1
				){
					controls.prePreURF = 1.0;
					controls.preVelURF = 1.0;
				}
				
				this->setCompValuesLeftRightFaceForSegregated(mesh, controls, species);
				// this->setIncomValuesLeftRightFace(mesh, controls, species);
				
				residualsP.clear();
				residualsP.resize(controls.nEq,0.0);
				
				
				this->calcPressureEq(mesh, controls, linAD, linAOD, linAFL, linAFR, residualsP);
				
				
				controls.prePreURF = tmpprePreURF;
				controls.preVelURF = tmppreVelURF;
				
				++controls.iterPre;
			}
			
			
			// if(controls.iterPBs == controls.iterPBsMax-1){
			controls.iterVof = 0;
			while(controls.iterVof<controls.iterVofMax){
				
				this->setCompValuesLeftRightFaceForSegregatedWithVfMSTACS(mesh, controls, species);
				// this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
				
				this->calcVolfracEq(mesh, controls, linAD, linAOD, linAFL, linAFR, residualsVF);
				
				
				// this->calcIncomCellEOSVF(mesh, controls, species);
				
				++controls.iterVof;
			}
			// }
			
			
			

			// this->calcIncomCellEOSVF(mesh, controls, species);
			this->calcCellEOSVF(mesh, controls, species);
			
			this->calcCellTransport(mesh, controls, species);
			
			
			
			vector<double> residuals(controls.nEq,0.0);
			for(int i=0; i<controls.nEq; ++i){
				residuals[i] += residualsM[i];
				residuals[i] += residualsP[i];
				residuals[i] += residualsVF[i];
			}

			vector<double> norm;
			this->calcNormResiduals(mesh, controls, residuals, norm);
			
			if(rank==0) {
				double dClock = clock() - controls.startClock;
				dClock /= CLOCKS_PER_SEC;
				cout << "|  └ PIMPLE step = " << controls.iterPBs << " | ";
				cout.precision(3);
				for(int i=0; i<controls.nEq; ++i){
					if(controls.iterPBs==0) norm0[i] = norm[i];
					
					if(controls.iterPBsMax==1){
						cout << scientific << norm[i] << " | ";
					}
					else{
						cout << scientific << norm[i] << " | ";
					}
				}
				cout.unsetf(ios::scientific);
				cout << dClock << " s | ";
				cout << endl;
				
			}
			
		
			++controls.iterPBs;
			
			
			++controls.iterTotal;
		
		}
		
	}
	

	
	
	//===================================================
	// compressible density based
	bool boolComDensityBased = true;
	
	if(boolComDensityBased){
		
		
		int num=0;
		for(auto& cell : mesh.cells){
			for(int i=0; i<controls.nEq; ++i){
				cell.var[controls.Qm[i]] = tmpQm[num][i];
			}
			++num;
		}
		
		
		controls.iterPseudo = 0;
		while(controls.iterPseudo<controls.iterPseudoMax){
			
			if(rank==0) cout << "|  └ pseudo step = " << controls.iterPseudo << " | ";
			
			
			this->calcUnderRelaxationFactorsDualTime(controls);
		
			this->calcPseudoTimeStep(mesh, controls);
			
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			this->setCompValuesLeftRightFaceWithRecon(mesh, controls, species);
			


			//=========================
			// Un initialization
			if( controls.iterPseudo == 0){
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
			
			
			vector<vector<double>> residuals(
				mesh.cells.size(),vector<double>(controls.nEq,0.0));
		
			this->calcRHS(mesh, controls, species, residuals);

			this->calcLinearSolver(mesh, controls, residuals);

			vector<double> norm;
			this->calcNormResiduals(mesh, controls, residuals, norm);
			
			if(rank==0) {
				double dClock = clock() - controls.startClock;
				dClock /= CLOCKS_PER_SEC;
				cout.precision(3);
				for(int i=0; i<controls.nEq; ++i){
					cout << scientific << norm[i] << " | ";
				}
				cout.unsetf(ios::scientific);
				cout << dClock << " s | ";
				cout << endl;
			}
			
			this->updateValues(mesh, controls, residuals);
			// this->updateValuesSingle(mesh, controls, residuals);
			
			this->calcCellEOSMF(mesh, controls, species);
			
			this->calcCellTransport(mesh, controls, species);
			
			++controls.iterPseudo;
			

			if(controls.saveControl == "pseudoTimeStep"){
				if(controls.iterPseudo % controls.saveInterval == 0){
					SEMO_Mesh_Save save;
					
					string foldername;
					std::ostringstream streamObj;
					streamObj << controls.iterPseudo;
					foldername = "./save/" + streamObj.str() + "/";
					save.vtu(foldername, mesh, controls, species);
				}
			}
		
		
			++controls.iterTotal;
			
		}
		
	}
	
	
	
	
	
	
	

	// if(rank==0) {

		SEMO_Utility_Math math;
		vector<vector<double>> gradP;
		math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
		ofstream outputFile;
 
		string filenamePlot = "pressure_point" + to_string(rank);
		outputFile.open(filenamePlot, ios::app);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	
		outputFile << controls.time << " ";
		outputFile.precision(4);
		for(int i=0; i<mesh.cells.size(); ++i){
			auto& cell = mesh.cells[i];
			if(cell.x<0.008 && cell.y<0.008){
				
				double PF = 
					cell.var[controls.P] +
					 ( (gradP[i][0])*(0.0-cell.x)
					  +(gradP[i][1])*(0.003-cell.y));
				// double PF = cell.var[controls.P];
				
				
				outputFile << scientific << PF << " ";
				cout << endl;
				cout << PF << endl;
			}
		}
		outputFile.unsetf(ios::scientific);
		outputFile << endl;
		
		outputFile.close();
	// }
	
	
	
	
	
	
	
}