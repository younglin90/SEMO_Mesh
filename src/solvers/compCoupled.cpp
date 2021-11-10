#include "build.h"
#include <cmath>
#include <array>


void SEMO_Solvers_Builder::compressibleCoupled(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	
	
	// for(auto& cell : mesh.cells){
		// for(int i=0; i<controls.nSp-1; ++i){
			// cell.var[controls.oldMF[i]] = max(1.e-8,min(1.0-1.e-8,cell.var[controls.oldMF[i]]));
			// cell.var[controls.MF[i]] = max(1.e-8,min(1.0-1.e-8,cell.var[controls.MF[i]]));
		// }
	// } 
	
	
	// save Q^(n-1)
	for(auto& cell : mesh.cells){
		cell.var[controls.Qm[0]] = cell.var[controls.oldP];
		cell.var[controls.Qm[1]] = cell.var[controls.oldU];
		cell.var[controls.Qm[2]] = cell.var[controls.oldV];
		cell.var[controls.Qm[3]] = cell.var[controls.oldW];
		cell.var[controls.Qm[4]] = cell.var[controls.oldT];
		cell.var[controls.Qm[5]] = cell.var[controls.oldRho];
		cell.var[controls.Qm[6]] = cell.var[controls.oldHt];
		for(int i=0; i<controls.nSp-1; ++i){
			cell.var[controls.Qm[7+i]] = cell.var[controls.oldMF[i]];
		}
	}
	
	// for(auto& face : mesh.faces){
		// face.var[controls.old2Un] = face.var[controls.oldUn];
	// }
	
	// save Q^(n)
	for(auto& cell : mesh.cells){
		cell.var[controls.oldP] = cell.var[controls.P];
		cell.var[controls.oldU] = cell.var[controls.U];
		cell.var[controls.oldV] = cell.var[controls.V];
		cell.var[controls.oldW] = cell.var[controls.W];
		cell.var[controls.oldT] = cell.var[controls.T];
		cell.var[controls.oldRho] = cell.var[controls.Rho];
		cell.var[controls.oldHt] = cell.var[controls.Ht];
		for(int i=0; i<controls.nSp-1; ++i){
			cell.var[controls.oldMF[i]] = cell.var[controls.MF[i]];
		}
	}
	
	for(auto& face : mesh.faces){
		face.var[controls.oldUn] = face.var[controls.Un];
	}
	
	
	
	//===================================================
	double corantNum = -10.0;
	if(controls.adjustTimeStep == "yes"){
		corantNum = this->calcTimeStepFromCorantNumber(mesh, controls, species);
	}
	else{
		this->calcCorantNumberForPrint(mesh, controls, corantNum);
	}
	if(rank==0) cout << " | timeStep = " << controls.timeStep;
	if(rank==0) cout << " | corantNumber = " << corantNum << endl;
	
	vector<double> norm(2,0.0);
	
	controls.iterPBs = 0;
	while(controls.iterPBs<controls.iterPBsMax){
	// while(controls.iterPBs<1){
		
		this->calcUnderRelaxationFactors(controls);
		
		
		
		// // mass fraction
		// for(int iter=0; iter<1; ++iter)
		// {
			// // reconstruction
			// // this->setCompValuesLeftRightFace(mesh, controls, species);
			// this->setCompValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
			
			// // solve
			// norm[0] = this->calcCompMassfracEq(mesh, controls, species);
			
			// // EOS, Transport
			// this->calcCellEOSMF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
		// }
		

		
		// // flows equations
		// // for(int iter=0; iter<20; ++iter)
		// {
			// // reconstruction
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			// // this->setCompValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
			
			// // solve
			// norm[0] = this->calcCompFlowsEq(mesh, controls, species);
			
			// // EOS, Transport
			// this->calcCellEOSMF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
		// }
			
		
		
		
		// // momentum equations
		// {
			// // reconstruction
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			
			// // solve
			// norm[0] = this->calcCompMomentumEq(mesh, controls, species);
			
			// // EOS, Transport
			// this->calcCellEOSMF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
		// }
		
		
		
		// // pressure equation
		// for(int iter=0; iter<5; ++iter)
		// {
			// // reconstruction
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			
			// // solve
			// norm[0] = this->calcCompPressureEq(mesh, controls, species);
			
			// // EOS, Transport
			// this->calcCellEOSMF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
		// }
		
		
		
		// // energy equations
		// {
			// // reconstruction
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			
			// // solve
			// norm[0] = this->calcCompEnergyEq(mesh, controls, species);
			
			// // EOS, Transport
			// this->calcCellEOSMF(mesh, controls, species);
			// this->calcCellTransport(mesh, controls, species);
		// }
		
		
		
		
		// coupled equations
		// for(int iter=0; iter<5; ++iter)
		{
			// reconstruction
			// this->setCompValuesLeftRightFace(mesh, controls, species);
			this->setCompValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
			
			// solve
			norm[0] = this->calcCompCoupledEq(mesh, controls, species);
			
			// EOS, Transport
			this->calcCellEOSMF(mesh, controls, species);
			this->calcCellTransport(mesh, controls, species);
		}
			
			
			
			
			
		// // print residual norm
		// vector<double> norm_global(2,0.0);
		// MPI_Allreduce(norm.data(), norm_global.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		// if(rank==0) {
			// double dClock = clock() - controls.startClock;
			// dClock /= CLOCKS_PER_SEC;
			// cout << "|  └ Coupled step = " << controls.iterPBs << " | ";
			// cout.precision(3);
			// cout << scientific << norm_global[0] << " | ";
			// cout << scientific << norm_global[1] << " | ";
			// cout.unsetf(ios::scientific);
			// cout << dClock << " s | ";
			// cout << endl;
			
		// }
		
		

		
		++controls.iterPBs;
		++controls.iterTotal;
		
	}
	
	

		
		// // for(auto& cell : mesh.cells){
			// // cell.var[controls.oldP] = cell.var[controls.P];
			// // cell.var[controls.oldU] = cell.var[controls.U];
			// // cell.var[controls.oldV] = cell.var[controls.V];
			// // cell.var[controls.oldW] = cell.var[controls.W];
			// // cell.var[controls.oldT] = cell.var[controls.T];
			// // for(int i=0; i<controls.nSp-1; ++i){
				// // cell.var[controls.oldVF[i]] = cell.var[controls.VF[i]];
				// // cell.var[controls.oldMF[i]] = cell.var[controls.MF[i]];
			// // }
		// // }
		

		// if(rank==0) cout << " | timeStep = " << controls.timeStep << endl;
		// if(rank==0) cout << "|  └ residuals" << " | ";
		
		// this->setCompValuesLeftRightFace(mesh, controls, species);
		// // this->setCompValuesLeftRightFaceWithRecon(mesh, controls, species);
		
		// vector<vector<double>> residuals(
			// mesh.cells.size(),vector<double>(controls.nEq,0.0));

		// this->calcSingleRHS(mesh, controls, species, residuals);

		// this->calcSingleLinearSolver(mesh, controls, residuals);

		// vector<double> norm;
		// this->calcNormResiduals(mesh, controls, residuals, norm);
		
		// if(rank==0) {
			// double dClock = clock() - controls.startClock;
			// dClock /= CLOCKS_PER_SEC;
			
			// cout.precision(3);
			// for(int i=0; i<controls.nEq; ++i){
				// cout << scientific << norm[i] << " | ";
			// }
			// cout.unsetf(ios::scientific);
			// cout << dClock << " s | ";
			// cout << endl;
		// }
		

		// this->updateValuesSingle(mesh, controls, residuals);
		
		// this->calcCellEOSMF(mesh, controls, species);
			
		// this->calcCellTransport(mesh, controls, species);
		
		// ++controls.iterPseudo;
		
		
		// ++controls.iterTotal;
		
		
	
	
	//####################################################
	
	
}