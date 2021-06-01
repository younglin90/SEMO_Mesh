#include "build.h"
#include <cmath>
#include <array>





void SEMO_Solvers_Builder::compressibleDensityBasedSingleTime(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
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
	}
	
	
	
	
	if(rank==0) cout << "|  └ residuals" << " | ";

	this->calcRealTimeStep(mesh, controls);
	
	this->setCompValuesLeftRightFace(mesh, controls, species);
	
	vector<vector<double>> residuals(
		mesh.cells.size(),vector<double>(controls.nEq,0.0));

	this->calcSingleRHS(mesh, controls, residuals);

	this->calcSingleLinearSolver(mesh, controls, residuals);

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
	

	this->updateValuesSingle(mesh, controls, residuals);
	
	this->calcCellEOSMF(mesh, controls, species);
		
	this->calcCellTransport(mesh, controls, species);
	
	++controls.iterPseudo;
	
	
	
	
	
	
}




void SEMO_Solvers_Builder::compressibleDensityBasedDualTime(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
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
		
		
		
	}
	
}






void SEMO_Solvers_Builder::updateValuesSingle(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		cell.var[controls.P] += residuals[i][0];
		
		cell.var[controls.U] += residuals[i][1];
		cell.var[controls.V] += residuals[i][2];
		cell.var[controls.W] += residuals[i][3];
		
		cell.var[controls.T] += residuals[i][4];
		
		for(int ns=0; ns<controls.nSp-1; ++ns){
			cell.var[controls.MF[ns]] += residuals[i][5+ns];
		}
		
		double tmp=0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			cell.var[controls.MF[ns]] = max(0.0,min(1.0,cell.var[controls.MF[ns]]));
			tmp += cell.var[controls.MF[ns]];
		}
		cell.var[controls.MF[controls.nSp-1]] = 1.0 - tmp;
		
	}
	
	
	
}







void SEMO_Solvers_Builder::updateValues(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<vector<double>>& residuals){
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		// for(int nq=0; nq<controls.nEq; ++nq){
			// cell.var[nq] += residuals[i][nq];
		// }
		cell.var[controls.P] += controls.dualTimeURF_P * residuals[i][0];
		
		cell.var[controls.U] += controls.dualTimeURF_U * residuals[i][1];
		cell.var[controls.V] += controls.dualTimeURF_U * residuals[i][2];
		cell.var[controls.W] += controls.dualTimeURF_U * residuals[i][3];

		cell.var[controls.T] += controls.dualTimeURF_T * residuals[i][4];
		
		for(int ns=0; ns<controls.nSp-1; ++ns){
			cell.var[controls.MF[ns]] += controls.dualTimeURF_MF * residuals[i][5+ns];
		}
		
		
		if( cell.var[controls.P] <= controls.minP ) 
			cell.var[controls.P] = controls.minP;
		if( cell.var[controls.P] >= controls.maxP ) 
			cell.var[controls.P] = controls.maxP;
		
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
		
		if( cell.var[controls.T] <= controls.minT ) 
			cell.var[controls.T] = controls.minT;
		if( cell.var[controls.T] >= controls.maxT ) 
			cell.var[controls.T] = controls.maxT;
		
		double tmp=0.0;
		for(int ns=0; ns<controls.nSp-1; ++ns){
			cell.var[controls.MF[ns]] = max(0.0,min(1.0,cell.var[controls.MF[ns]]));
			tmp += cell.var[controls.MF[ns]];
		}
		cell.var[controls.MF[controls.nSp-1]] = 1.0 - tmp;
		
	}
	
	
	
}








void SEMO_Solvers_Builder::calcUnderRelaxationFactorsDualTime(
	SEMO_Controls_Builder& controls){
	
	
	// relaxationFactors
	if(controls.dualTimeAdjustRF_P == "yes"){
		for(int i=1; i<controls.dualTimeAdjustSteps_P.size(); ++i){
			if( 
			controls.iterPseudo >= controls.dualTimeAdjustSteps_P[i-1] &&
			controls.iterPseudo < controls.dualTimeAdjustSteps_P[i] ){
				double stri = (double)controls.dualTimeAdjustSteps_P[i-1];
				double endi = (double)controls.dualTimeAdjustSteps_P[i];
				double str = (double)controls.dualTimeAdjustValues_P[i-1];
				double end = (double)controls.dualTimeAdjustValues_P[i];
				double prnt = (double)controls.iterPseudo;
				
				controls.dualTimeURF_P = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	if(controls.dualTimeAdjustRF_U == "yes"){
		for(int i=1; i<controls.dualTimeAdjustSteps_U.size(); ++i){
			if( 
			controls.iterPseudo >= controls.dualTimeAdjustSteps_U[i-1] &&
			controls.iterPseudo < controls.dualTimeAdjustSteps_U[i] ){
				double stri = (double)controls.dualTimeAdjustSteps_U[i-1];
				double endi = (double)controls.dualTimeAdjustSteps_U[i];
				double str = (double)controls.dualTimeAdjustValues_U[i-1];
				double end = (double)controls.dualTimeAdjustValues_U[i];
				double prnt = (double)controls.iterPseudo;
				
				controls.dualTimeURF_U = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	if(controls.dualTimeAdjustRF_T == "yes"){
		for(int i=1; i<controls.dualTimeAdjustSteps_T.size(); ++i){
			if( 
			controls.iterPseudo >= controls.dualTimeAdjustSteps_T[i-1] &&
			controls.iterPseudo < controls.dualTimeAdjustSteps_T[i] ){
				double stri = (double)controls.dualTimeAdjustSteps_T[i-1];
				double endi = (double)controls.dualTimeAdjustSteps_T[i];
				double str = (double)controls.dualTimeAdjustValues_T[i-1];
				double end = (double)controls.dualTimeAdjustValues_T[i];
				double prnt = (double)controls.iterPseudo;
				
				controls.dualTimeURF_T = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
	
	
	
	if(controls.dualTimeAdjustRF_MF == "yes"){
		for(int i=1; i<controls.dualTimeAdjustSteps_MF.size(); ++i){
			if( 
			controls.iterPseudo >= controls.dualTimeAdjustSteps_MF[i-1] &&
			controls.iterPseudo < controls.dualTimeAdjustSteps_MF[i] ){
				double stri = (double)controls.dualTimeAdjustSteps_MF[i-1];
				double endi = (double)controls.dualTimeAdjustSteps_MF[i];
				double str = (double)controls.dualTimeAdjustValues_MF[i-1];
				double end = (double)controls.dualTimeAdjustValues_MF[i];
				double prnt = (double)controls.iterPseudo;
				
				controls.dualTimeURF_MF = 
					str + (end-str)/(endi-stri)*(prnt-stri);
				
				break;
			}	
		}
	}
	
}