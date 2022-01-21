#include "build.h"
#include <cmath>
#include <array>
#include <time.h>


void SEMO_Solvers_Builder::compressibleDensityBasedSingleTime(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	double strSolverClock = clock();
	

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
	
	
	

	// save old timesteps
	controls.old2TimeStep = controls.oldTimeStep;
	controls.oldTimeStep = controls.timeStep;
	
	

	this->calcRealTimeStep(mesh, controls);
	// if(rank==0) cout << " | timeStep = " << controls.timeStep << endl;
	// if(rank==0) cout << "|  └ residuals" << " | ";
	
	//===================================================
	double corantNum = -10.0;
	if(controls.adjustTimeStep == "yes"){
		// corantNum = this->calcTimeStepFromCorantNumber(mesh, controls, species);
		corantNum = this->calcExplicitTimeStepFromCorantNumber(mesh, controls, species);
	}
	else{
		// this->calcCorantNumberForPrint(mesh, controls, corantNum);
		this->calcExplicitCorantNumberForPrint(mesh, controls, corantNum);
	}
	if(rank==0) cout << " | timeStep = " << controls.timeStep;
	if(rank==0) cout << " | corantNumber = " << corantNum << endl;
	if(rank==0) cout << "|  └ residuals" << " | ";
	
	
	
	this->setCompValuesLeftRightFace(mesh, controls, species);
	// this->setCompValuesLeftRightFaceWithVfMSTACS(mesh, controls, species);
	// this->setCompValuesLeftRightFaceWithRecon(mesh, controls, species);
	vector<vector<double>> residuals(
		mesh.cells.size(),vector<double>(controls.nEq,0.0));
	// cout << "aaaaa1" << endl;

	this->calcSingleRHS(mesh, controls, species, residuals);
	// cout << "aaaaa2" << endl;

	this->calcSingleLinearSolver(mesh, controls, residuals);
	// cout << "aaaaa3" << endl;

	vector<double> norm;
	this->calcNormResiduals(mesh, controls, residuals, norm);
	// cout << "aaaaa4" << endl;
	
	if(rank==0) {
		// double dClock = clock() - controls.startClock;
		double dClock = clock() - strSolverClock;
		dClock /= CLOCKS_PER_SEC;
		
		cout.precision(3);
		for(int i=0; i<controls.nEq; ++i){
			cout << scientific << norm[i] << " | ";
		}
		cout.unsetf(ios::scientific);
		cout << dClock << " s | ";
		cout << endl;
	}
	

	this->updateValuesSingle(mesh, controls, residuals);
	
	this->calcCellEOSMF(mesh, controls, species);
		
	this->calcCellTransport(mesh, controls, species);
	
	++controls.iterPseudo;
	
	
	++controls.iterTotal;
	
	
	
	
}




void SEMO_Solvers_Builder::compressibleDensityBasedDualTime(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	

	if( controls.iterReal == 0 ){
		// this->calcIncomCellEOSVF(mesh, controls, species);
		this->calcCellEOSVF(mesh, controls, species);
		this->calcCellTransport(mesh, controls, species);
	}
	
	
	
	// if( controls.iterReal == 0 ){
		// for(auto& cell : mesh.cells){
			// cell.var[controls.oldP] = cell.var[controls.P];
			// cell.var[controls.oldU] = cell.var[controls.U];
			// cell.var[controls.oldV] = cell.var[controls.V];
			// cell.var[controls.oldW] = cell.var[controls.W];
			// cell.var[controls.oldT] = cell.var[controls.T];
			// for(int i=0; i<controls.nSp-1; ++i){
				// cell.var[controls.oldVF[i]] = cell.var[controls.VF[i]];
				// cell.var[controls.oldMF[i]] = cell.var[controls.MF[i]];
			// }
			
			// cell.var[controls.Qn[0]] = cell.var[controls.Rho];
			// cell.var[controls.Qn[1]] = cell.var[controls.Rho] * cell.var[controls.U];
			// cell.var[controls.Qn[2]] = cell.var[controls.Rho] * cell.var[controls.V];
			// cell.var[controls.Qn[3]] = cell.var[controls.Rho] * cell.var[controls.W];
			// cell.var[controls.Qn[4]] = cell.var[controls.Rho] * cell.var[controls.Ht] - cell.var[controls.P];
			// for(int ns=0; ns<controls.nSp-1; ++ns){
				// cell.var[controls.Qn[5+ns]] = cell.var[controls.Rho] * cell.var[controls.MF[ns]];
			// }
			
			// for(int i=0; i<controls.nEq; ++i){
				// cell.var[controls.Qm[i]] = cell.var[controls.Qn[i]];
			// }
			
			
		// }
	// }
	// else{
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
	// }
	
	
	//===========================================
	// SER
	controls.allowableCFL = controls.specifiedCFL;
	vector<double> norm_ago;
	//===========================================
	
	
	controls.iterPseudo = 0;
	while(controls.iterPseudo<controls.iterPseudoMax){
		
		// if(rank==0) cout << " | timeStep = " << controls.timeStep << endl;
		if(rank==0) cout << "|  └ pseudo step = " << controls.iterPseudo << " | ";
		
		this->calcUnderRelaxationFactorsDualTime(controls);
	
		this->calcPseudoTimeStep(mesh, controls);
		
		// this->setCompValuesLeftRightFace(mesh, controls, species);
		this->setCompValuesLeftRightFaceWithRecon(mesh, controls, species);
		
		vector<vector<double>> residuals(
			mesh.cells.size(),vector<double>(controls.nEq,0.0));
	
		this->calcRHS(mesh, controls, species, residuals);




		this->calcLinearSolver(mesh, controls, residuals);
		// this->calcLinearSolverLUSGS(mesh, controls, residuals);
		
		
		//===========================================
		// SER
		int myReturnSER = 0;
		for(int i=0; i<mesh.cells.size(); ++i){
			if(
			abs(residuals[i][0])>0.01*mesh.cells[i].var[controls.P] ||
			abs(residuals[i][4])>0.01*mesh.cells[i].var[controls.T] 
			){
				myReturnSER = 1;
			}
		}
		
		int returnSER;
		MPI_Allreduce(&myReturnSER, &returnSER, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		
		if(returnSER>0){
				if(rank==0) cout << "| Retunrs calculation for SER method" << endl;
				controls.allowableCFL *= 0.5; 
				continue;
		}
		//===========================================
		

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
		
		
		//===========================================
		// SER
		if(returnSER==0 && controls.iterPseudo!=0){
			double tmp_cfl = controls.allowableCFL*norm_ago[0]/norm[0];
			tmp_cfl = min(tmp_cfl, controls.allowableCFL*norm_ago[4]/norm[4]);
			controls.allowableCFL = min(tmp_cfl, controls.specifiedCFL); 
		}
		
		norm_ago = norm;
		//===========================================
		
		
		
		this->updateValues(mesh, controls, residuals);
		// this->updateValuesSingle(mesh, controls, residuals);
		
		this->calcCellEOSMF(mesh, controls, species);
		
		this->calcCellTransport(mesh, controls, species);
		
		++controls.iterPseudo;
		

		if(controls.saveControl == "pseudoTimeStep"){
			if(controls.iterPseudo % (int)controls.saveInterval == 0){
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
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
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
	
	
	
	
	
	

	vector<double> residualsx_recv;
	vector<double> residualsy_recv;
	vector<double> residualsz_recv;
	if(size>1){
		// processor faces
		// gradP , 
		vector<double> residualsx_send;
		vector<double> residualsy_send;
		vector<double> residualsz_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				residualsx_send.push_back(residuals[face.owner][1]);
				residualsy_send.push_back(residuals[face.owner][2]);
				residualsz_send.push_back(residuals[face.owner][3]);
			}
		}
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					residualsx_send, residualsx_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					residualsy_send, residualsy_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		mpi.setProcsFaceDatas(
					residualsz_send, residualsz_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
	}
	
	
	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];

		double wCL = face.wC;
		double wCR = 1.0-face.wC;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
				
				
			face.var[controls.Un] += controls.dualTimeURF_U * (
				(wCL*residuals[face.owner][1]+wCR*residuals[face.neighbour][1])*face.unitNormals[0] +
				(wCL*residuals[face.owner][2]+wCR*residuals[face.neighbour][2])*face.unitNormals[1] +
				(wCL*residuals[face.owner][3]+wCR*residuals[face.neighbour][3])*face.unitNormals[2] );
				
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				
				
			face.var[controls.Un] += controls.dualTimeURF_U * (
				(wCL*residuals[face.owner][1]+wCR*residualsx_recv[proc_num])*face.unitNormals[0] +
				(wCL*residuals[face.owner][2]+wCR*residualsy_recv[proc_num])*face.unitNormals[1] +
				(wCL*residuals[face.owner][3]+wCR*residualsz_recv[proc_num])*face.unitNormals[2] );
				
			++proc_num;
			
		}
		

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