#include "build.h"
#include <cmath>
#include <array>


void SEMO_Solvers_Builder::incompressiblePressureBased(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
		
	int rank = MPI::COMM_WORLD.Get_rank();
    int size = MPI::COMM_WORLD.Get_size();
	
	
	
	if( controls.iterReal == 0 ){
		this->calcIncomCellEOSVF(mesh, controls, species);
		this->calcCellTransport(mesh, controls, species);
		// Un
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				face.var[controls.Un] =  0.5*(mesh.cells[face.owner].var[controls.U]*face.unitNormals[0]);
				face.var[controls.Un] += 0.5*(mesh.cells[face.neighbour].var[controls.U]*face.unitNormals[0]);
				face.var[controls.Un] += 0.5*(mesh.cells[face.owner].var[controls.V]*face.unitNormals[1]);
				face.var[controls.Un] += 0.5*(mesh.cells[face.neighbour].var[controls.V]*face.unitNormals[1]);
				face.var[controls.Un] += 0.5*(mesh.cells[face.owner].var[controls.W]*face.unitNormals[2]);
				face.var[controls.Un] += 0.5*(mesh.cells[face.neighbour].var[controls.W]*face.unitNormals[2]);
			}
		}
	}
	
	if( controls.iterReal == 0 ){
		for(auto& cell : mesh.cells){
			cell.var[controls.Qm[0]] = cell.var[controls.U];
			cell.var[controls.Qm[1]] = cell.var[controls.V];
			cell.var[controls.Qm[2]] = cell.var[controls.W];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qm[3+i]] = cell.var[controls.VF[i]];
			}
			
		}
	}
	else{
		for(auto& cell : mesh.cells){
			cell.var[controls.Qm[0]] = cell.var[controls.oldU];
			cell.var[controls.Qm[1]] = cell.var[controls.oldV];
			cell.var[controls.Qm[2]] = cell.var[controls.oldW];
			for(int i=0; i<controls.nSp-1; ++i){
				cell.var[controls.Qm[3+i]] = cell.var[controls.oldVF[i]];
			}
			
		}
	}
	
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
		
		
	}
	
	
	// save old timesteps
	controls.old2TimeStep = controls.oldTimeStep;
	controls.oldTimeStep = controls.timeStep;
	
	
	
	
	double corantNum = -10.0;
	if(controls.adjustTimeStep == "yes"){
		corantNum = this->calcTimeStepFromCorantNumber(mesh, controls, species);
	}
	else{
		this->calcCorantNumberForPrint(mesh, controls, corantNum);
	}
	if(rank==0) cout << " | timeStep = " << controls.timeStep;
	if(rank==0) cout << " | corantNumber = " << corantNum << endl;
	

	// save old timesteps
	if( controls.iterReal == 0 ){
		controls.oldTimeStep = controls.timeStep;
		controls.old2TimeStep = controls.oldTimeStep;
	}
	
	
	
	vector<double> norm0(controls.nEq,0.0);
	vector<double> normNM(controls.nEq,0.0);
	vector<double> norm(2,0.0);
	

	// for(int iterVFs=0; iterVFs<controls.iterPBsMax; ++iterVFs){
		controls.iterVof = 0;
		while(controls.iterVof<controls.iterVofMax){
			
			// this->setIncomValuesLeftRightFace(mesh, controls, species);
			// this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
			this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
			// this->setCompValuesLeftRightFaceWithNVD(mesh, controls, species);
			
			norm[1] = this->calcVolfracEq(mesh, controls, species);
			
			
			this->calcIncomCellEOSVF(mesh, controls, species);
			this->calcCellTransport(mesh, controls, species);

			vector<double> norm_global(2,0.0);
			MPI_Allreduce(norm.data(), norm_global.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			if(rank==0) {
				double dClock = clock() - controls.startClock;
				dClock /= CLOCKS_PER_SEC;
				cout << "|  └ VF step = " << controls.iterVof << " | ";
				cout.precision(3);
				cout << scientific << norm_global[1] << " | ";
				cout.unsetf(ios::scientific);
				cout << dClock << " s | ";
				cout << endl;
				
			}
			
			++controls.iterVof;
			
		}
		
	// }
	
	
	this->setIncomValuesLeftRightFace(mesh, controls, species);
	this->setCellVarMinMax(mesh, controls.VF[0], controls.maximumVF[0], controls.minimumVF[0]);
	this->calcSourceGravity(mesh, controls, species);
	this->calcSourceSurfaceTension(mesh, controls, species);

	
	controls.iterPBs = 0;
	while(controls.iterPBs<controls.iterPBsMax){
		
		
		this->calcUnderRelaxationFactors(controls);
		
		
		//####################################################
		// // volume fraction
		// for(int iterVFs=0; iterVFs<1; ++iterVFs){
			// controls.iterVof = 0;
			// while(controls.iterVof<controls.iterVofMax){
				
				// // if(rank==0) cout << "|========== volume fraction" << endl;
				// // this->setIncomValuesLeftRightFace(mesh, controls, species);
				// // this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
				// this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
				// // this->setCompValuesLeftRightFaceWithNVD(mesh, controls, species);
				
				// norm[1] = this->calcVolfracEq(mesh, controls);
				
				
				// this->calcIncomCellEOSVF(mesh, controls, species);
				// this->calcCellTransport(mesh, controls, species);
				
				// ++controls.iterVof;
			
				// vector<double> norm_global(2,0.0);
				// MPI_Allreduce(norm.data(), norm_global.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				// if(rank==0) {
					// double dClock = clock() - controls.startClock;
					// dClock /= CLOCKS_PER_SEC;
					// cout << "|  └ VF step = " << iterVFs << " | ";
					// cout.precision(3);
					// cout << scientific << norm_global[1] << " | ";
					// cout.unsetf(ios::scientific);
					// cout << dClock << " s | ";
					// cout << endl;
					
				// }
			// }
			
		// }
		
 
 
 
		// this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
		// norm[1] = this->calcVolfracEq(mesh, controls);
		// this->calcIncomCellEOSVF(mesh, controls, species);
		// this->calcCellTransport(mesh, controls, species);
		
		
		
		// momentum
		// if(rank==0) cout << "|========== momentum" << endl;
		this->setIncomValuesLeftRightFace(mesh, controls, species);
		norm[0] = this->calcMomentumEqs(mesh, controls, species);
		
		
		
		// pressure
		controls.iterPre = 0;
		while(controls.iterPre<controls.iterPreMax){
			// if(rank==0) cout << "|========== pressure" << endl;
			
			this->setIncomValuesLeftRightFace(mesh, controls, species);
			norm[1] = this->calcPressureEq(mesh, controls, species, 1);
			
			++controls.iterPre;
		}
		
		
		// // coupled
		// // if(rank==0) cout << "|========== coupled" << endl;
		// this->setIncomValuesLeftRightFace(mesh, controls, species);
		// norm[0] = this->calcCoupledEq(mesh, controls, species);
		



		vector<double> norm_global(2,0.0);
		MPI_Allreduce(norm.data(), norm_global.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		if(rank==0) {
			double dClock = clock() - controls.startClock;
			dClock /= CLOCKS_PER_SEC;
			cout << "|  └ Coupled step = " << controls.iterPBs << " | ";
			cout.precision(3);
			cout << scientific << norm_global[0] << " | ";
			cout << scientific << norm_global[1] << " | ";
			// for(int i=0; i<norm.size(); ++i){
				// // // if(controls.iterPBs==0) norm0[i] = norm[i];
				// // // cout << scientific << abs(norm[i] - normNM[i])/(norm0[i]+1.e-200) << " | ";
				// // // normNM[i] = norm[i];
				// // if(controls.iterPBs==0) norm0[i] = norm[i];
				
				// // if(controls.iterPBsMax==1){
					// // cout << scientific << norm[i] << " | ";
				// // }
				// // else{
					// // // cout << scientific << norm[i]/(norm0[i]+1.e-200) << " | ";
					// // cout << scientific << norm[i] << " | ";
				// // }
				
				// cout << scientific << norm_global[i] << " | ";
			// }
			cout.unsetf(ios::scientific);
			cout << dClock << " s | ";
			cout << endl;
			
		}
		
		//####################################################
		
		

		// //####################################################
		// vector<double> residualsM(controls.nEq,0.0);
		// vector<double> residualsP(controls.nEq,0.0);
		// vector<double> residualsVF(controls.nEq,0.0);
		
		// vector<double> linAD(mesh.cells.size(),0.0);
		// vector<double> linAOD(mesh.cells.size(),0.0);
		// vector<double> linAFL(mesh.faces.size(),0.0);
		// vector<double> linAFR(mesh.faces.size(),0.0);
		
	
		// controls.iterMom = 0;
		// while(controls.iterMom<controls.iterMomMax){
			// double tmppreVelURF = controls.momVelURF;
			// if(
			// controls.iterPBsMax == 1 &&
			// controls.iterMom == controls.iterMomMax-1
			// ){
				// controls.momVelURF = 1.0;
			// }
			
			// // this->setIncomValuesLeftRightFace(mesh, controls, species);
			// this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
			// // this->setCompValuesLeftRightFace(mesh, controls, species);
			
			// norm[1] = this->calcMomentumEqs(mesh, controls, species, linAD, linAOD, linAFL, linAFR, residualsM);
			
			// controls.momVelURF = tmppreVelURF;
			
			// ++controls.iterMom;
			
		// }
		


		// controls.iterPre = 0;
		// while(controls.iterPre<controls.iterPreMax){
			// double tmpprePreURF = controls.prePreURF;
			// double tmppreVelURF = controls.preVelURF;
			// if(
			// controls.iterPBsMax == 1 &&
			// controls.iterPre == controls.iterPreMax-1
			// ){
				// controls.prePreURF = 1.0;
				// controls.preVelURF = 1.0;
			// }
			
			// this->setIncomValuesLeftRightFace(mesh, controls, species);
			// // this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
			// // this->setCompValuesLeftRightFace(mesh, controls, species);
			
			// residualsP.clear();
			// residualsP.resize(controls.nEq,0.0);
			
			
			// norm[0] = this->calcPressureEq(mesh, controls, linAD, linAOD, linAFL, linAFR, residualsP);
			
			
			// // vector<double> norm;
			// // this->calcNormResiduals(mesh, controls, residualsP, norm);
			// // for(int i=0; i<controls.nEq; ++i){
				// // cout << scientific << norm[i] << " | ";
			// // }
			// // cout.unsetf(ios::scientific);
			// // cout << endl;
			
			
			// controls.prePreURF = tmpprePreURF;
			// controls.preVelURF = tmppreVelURF;
			
			// ++controls.iterPre;
		// }
		
		
		// // if(controls.iterPBs == controls.iterPBsMax-1){
		// controls.iterVof = 0;
		// while(controls.iterVof<controls.iterVofMax){
			
			// // this->setIncomValuesLeftRightFace(mesh, controls, species);
			// // this->setIncomValuesLeftRightFaceWithReconPV(mesh, controls, species);
			// this->setIncomValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
			// // this->setCompValuesLeftRightFaceWithNVD(mesh, controls, species);
			
			// this->calcVolfracEq(mesh, controls, linAD, linAOD, linAFL, linAFR, residualsVF);
			
			
			// this->calcIncomCellEOSVF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
			
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
			// double dClock = clock() - controls.startClock;
			// dClock /= CLOCKS_PER_SEC;
			// cout << "|  └ PIMPLE step = " << controls.iterPBs << " | ";
			// cout.precision(3);
			// for(int i=0; i<norm.size(); ++i){
				// // if(controls.iterPBs==0) norm0[i] = norm[i];
				// // cout << scientific << abs(norm[i] - normNM[i])/(norm0[i]+1.e-200) << " | ";
				// // normNM[i] = norm[i];
				// if(controls.iterPBs==0) norm0[i] = norm[i];
				
				// if(controls.iterPBsMax==1){
					// cout << scientific << norm[i] << " | ";
				// }
				// else{
					// // cout << scientific << norm[i]/(norm0[i]+1.e-200) << " | ";
					// cout << scientific << norm[i] << " | ";
				// }
				
			// }
			// cout.unsetf(ios::scientific);
			// cout << dClock << " s | ";
			// cout << endl;
			
		// }
		// //####################################################
		
		
		
		
		
		++controls.iterPBs;
		

		// if(controls.saveControl == "pseudoTimeStep"){
			// if(controls.iterPseudo % controls.saveInterval == 0){
				// SEMO_Mesh_Save save;
				
				// string foldername;
				// std::ostringstream streamObj;
				// streamObj << controls.iterPseudo;
				// foldername = "./save/" + streamObj.str() + "/";
				// save.vtu(foldername, mesh, controls, species);
			// }
		// }
		
		
		++controls.iterTotal;
		
	}
	
	
	
	// this->saveSurfaceNormalVelocity(mesh, controls, species);
	
	
	
	
	
	
	
	// double maxPressure = 0.0;
	// // int iiiiisave = 0;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(maxPressure<mesh.cells[i].var[controls.P]){
			// maxPressure = mesh.cells[i].var[controls.P];
			// // iiiiisave = i;
		// }
	// }
	// double maxPressureReduced;
    // MPI_Allreduce(&maxPressure, &maxPressureReduced, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	// if(rank==0) cout << maxPressureReduced << endl;
	
	

	//==================================
	// extract for point data
	// // if(rank==0) {

		// SEMO_Utility_Math math;
		// vector<vector<double>> gradP;
		// math.calcLeastSquare2nd(mesh, controls.P, controls.fP, gradP);
	
		// ofstream outputFile;
 
		// string filenamePlot = "pressure_point" + to_string(rank);
		// outputFile.open(filenamePlot, ios::app);
		// if(outputFile.fail()){
			// cerr << "Unable to write file for writing." << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
	
		// outputFile << controls.time << " ";
		// outputFile.precision(4);
		// for(int i=0; i<mesh.cells.size(); ++i){
			// auto& cell = mesh.cells[i];
			// if(cell.x<0.008 && cell.y<0.008){
				
				// double PF = 
					// cell.var[controls.P] +
					 // ( (gradP[i][0])*(0.0-cell.x)
					  // +(gradP[i][1])*(0.003-cell.y));
				// // double PF = cell.var[controls.P];
				
				
				// outputFile << scientific << PF << " ";
				// cout << endl;
				// cout << PF << endl;
			// }
		// }
		// outputFile.unsetf(ios::scientific);
		// outputFile << endl;
		
		// outputFile.close();
	// // }
	//==================================
		
	
	
}




void SEMO_Solvers_Builder::calcUnderRelaxationFactors(
	SEMO_Controls_Builder& controls){
	
	
	// relaxationFactors
	if(controls.momVelAdjustRF == "yes"){
		for(int i=1; i<controls.momVelAdjustSteps.size(); ++i){
			if( 
			controls.iterPBs >= controls.momVelAdjustSteps[i-1] &&
			controls.iterPBs < controls.momVelAdjustSteps[i] ){
				double stri = (double)controls.momVelAdjustSteps[i-1];
				double endi = (double)controls.momVelAdjustSteps[i];
				double str = (double)controls.momVelAdjustValues[i-1];
				double end = (double)controls.momVelAdjustValues[i];
				double prnt = (double)controls.iterPBs;
				
				controls.momVelURF = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	
	if(controls.prePreAdjustRF == "yes"){
		for(int i=1; i<controls.prePreAdjustSteps.size(); ++i){
			if( 
			controls.iterPBs >= controls.prePreAdjustSteps[i-1] &&
			controls.iterPBs < controls.prePreAdjustSteps[i] ){
				double stri = (double)controls.prePreAdjustSteps[i-1];
				double endi = (double)controls.prePreAdjustSteps[i];
				double str = (double)controls.prePreAdjustValues[i-1];
				double end = (double)controls.prePreAdjustValues[i];
				double prnt = (double)controls.iterPBs;
				
				controls.prePreURF = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	if(controls.preVelAdjustRF == "yes"){
		for(int i=1; i<controls.preVelAdjustSteps.size(); ++i){
			if( 
			controls.iterPBs >= controls.preVelAdjustSteps[i-1] &&
			controls.iterPBs < controls.preVelAdjustSteps[i] ){
				double stri = (double)controls.preVelAdjustSteps[i-1];
				double endi = (double)controls.preVelAdjustSteps[i];
				double str = (double)controls.preVelAdjustValues[i-1];
				double end = (double)controls.preVelAdjustValues[i];
				double prnt = (double)controls.iterPBs;
				
				controls.preVelURF = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	
	if(controls.vofVofAdjustRF == "yes"){
		for(int i=1; i<controls.vofVofAdjustSteps.size(); ++i){
			if( 
			controls.iterPBs >= controls.vofVofAdjustSteps[i-1] &&
			controls.iterPBs < controls.vofVofAdjustSteps[i] ){
				double stri = (double)controls.vofVofAdjustSteps[i-1];
				double endi = (double)controls.vofVofAdjustSteps[i];
				double str = (double)controls.vofVofAdjustValues[i-1];
				double end = (double)controls.vofVofAdjustValues[i];
				double prnt = (double)controls.iterPBs;
				
				controls.vofVofURF = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	
	
	
	
}