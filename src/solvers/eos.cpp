#include "build.h"
#include <cmath>
#include <array>
#include <tuple>





void SEMO_Solvers_Builder::calcIncomCellEOSVF(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int nSp = controls.nSp;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		vector<double> rhoi(nSp,0.0);
		vector<double> Ci(nSp,0.0);
		vector<double> Hti(nSp,0.0);
		vector<double> drhodPi(nSp,0.0);
		vector<double> drhodTi(nSp,0.0);
		vector<double> dhdPi(nSp,0.0);
		vector<double> dhdTi(nSp,0.0);
		
		double P = 101325.0;
		double T = 300.0;

		for(int ns=0; ns<nSp; ++ns){
			if( species[ns].rhoType == "const" ){
				rhoi[ns] = species[ns].rho;
				Ci[ns] = 340.0;
				Hti[ns] = 10000.0;
				drhodPi[ns] = 1.0e+8;
				drhodTi[ns] = 1.0e+8;
				dhdPi[ns] = 1.0e+8;
				dhdTi[ns] = 1.0e+8;
			}
			else if( species[ns].rhoType == "Ideal" ){
				this->eosIdeal(
						species[ns],
						P, 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						T,
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
			else if( species[ns].rhoType == "NASG" ){
				this->eosNASG(
						species[ns],
						P, 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						T,
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
				
		}


		vector<double> VF;
		double tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			double MF_tmp = max(0.0,min(1.0,cell.var[controls.VF[ns]]));
			VF.push_back( max(0.0,min(1.0,MF_tmp)) );
			tmp1 += MF_tmp;
		}
		VF.push_back( max(0.0,min(1.0,1.0 - tmp1)) );
		
		
		double rho = 0.0;
		for(int ns=0; ns<nSp; ++ns){
			rho += VF[ns]*rhoi[ns];
		}
		
		// if(rho<10.0) cout << rho << endl;
		
		
		cell.var[controls.Rho] = rho;
		
		
		// mass fraction
		vector<double> MF(nSp,0.0);
		tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			MF[ns] = VF[ns]*rhoi[ns]/rho;
			
			MF[ns] = max(0.0,min(1.0,MF[ns]));
			
			tmp1 += MF[ns];
		}
		MF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
		
		// double testrho = 0.0;
		// for(int ns=0; ns<nSp; ++ns){
			// testrho += MF[ns]/rhoi[ns];
		// }
		// testrho = 1.0 / testrho;
		
		// if(rho != testrho){
			// cout << rho << " " << testrho << endl;
		// }
		
		
		// double Ht=0.0;
		// double drdp=0.0;
		// double drdt=0.0;
		// double dhdp=0.0;
		// double dhdt=0.0;
		// for(int ns=0; ns<nSp; ++ns){
			// Ht += MF[ns]*Hti[ns];
			// drdp += VF[ns]*drhodPi[ns];
			// drdt += VF[ns]*drhodTi[ns];
			// dhdp += MF[ns]*dhdPi[ns];
			// dhdt += MF[ns]*dhdTi[ns];
		// }



		for(int ns=0; ns<nSp; ++ns){
			cell.var[controls.MF[ns]] = MF[ns];
			cell.var[controls.VF[ns]] = VF[ns];
		}
		
		// cell.var[controls.Ht] = Ht;
		// cell.var[controls.dRhoDP] = drdp;
		// cell.var[controls.dRhoDT] = drdt;
		// cell.var[controls.dHtDP] = dhdp;
		// cell.var[controls.dHtDT] = dhdt;
		
		// for(int ns=0; ns<nSp-1; ++ns){
			// cell.var[controls.dRhoDMF[ns]] = -rho*rho*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]);
			// // cout << cell.var[controls.dRhoDMF[ns]] << endl;
			// cell.var[controls.dHtDMF[ns]] = dhdTi[ns] - dhdTi[nSp-1];
		// }
		
		// // cout << drdp << " " << rho << " " << drdt << " " << dhdt << " " << dhdp << endl;
		
		// cell.var[controls.C] = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
		// cell.var[controls.C] = sqrt( 1.0 / cell.var[controls.C] );
	
		// // cout << cell.var[controls.C] << endl;
		
	}
	
	
	
	
}











void SEMO_Solvers_Builder::calcCellEOSVF(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int nSp = controls.nSp;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		vector<double> rhoi(nSp,0.0);
		vector<double> Ci(nSp,0.0);
		vector<double> Hti(nSp,0.0);
		vector<double> drhodPi(nSp,0.0);
		vector<double> drhodTi(nSp,0.0);
		vector<double> dhdPi(nSp,0.0);
		vector<double> dhdTi(nSp,0.0);

		for(int ns=0; ns<nSp; ++ns){
			if( species[ns].rhoType == "const" ){
				rhoi[ns] = species[ns].rho;
				Ci[ns] = 340.0;
				Hti[ns] = 10000.0;
				drhodPi[ns] = 1.0e+8;
				drhodTi[ns] = 1.0e+8;
				dhdPi[ns] = 1.0e+8;
				dhdTi[ns] = 1.0e+8;
			}
			else if( species[ns].rhoType == "Ideal" ){
				this->eosIdeal(
						species[ns],
						cell.var[controls.P], 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						cell.var[controls.T],
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
			else if( species[ns].rhoType == "NASG" ){
				this->eosNASG(
						species[ns],
						cell.var[controls.P], 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						cell.var[controls.T],
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
				
		}


		vector<double> VF;
		double tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			double MF_tmp = max(0.0,min(1.0,cell.var[controls.VF[ns]]));
			VF.push_back( max(0.0,min(1.0,MF_tmp)) );
			tmp1 += MF_tmp;
		}
		VF.push_back( max(0.0,min(1.0,1.0 - tmp1)) );
		
		// cout << nSp << " " << VF[0] << " " << VF[1] << endl;
		
		double rho = 0.0;
		for(int ns=0; ns<nSp; ++ns){
			rho += VF[ns]*rhoi[ns];
		}
		
		// if(rho<10.0) cout << rho << endl;
		
		
		cell.var[controls.Rho] = rho;
		
		
		// mass fraction
		vector<double> MF(nSp,0.0);
		tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			MF[ns] = VF[ns]*rhoi[ns]/rho;
			
			MF[ns] = max(0.0,min(1.0,MF[ns]));
			
			tmp1 += MF[ns];
		}
		MF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
		
		// double testrho = 0.0;
		// for(int ns=0; ns<nSp; ++ns){
			// testrho += MF[ns]/rhoi[ns];
		// }
		// testrho = 1.0 / testrho;
		
		// if(rho != testrho){
			// cout << rho << " " << testrho << endl;
		// }
		
		
		double Ht=0.0;
		double drdp=0.0;
		double drdt=0.0;
		double dhdp=0.0;
		double dhdt=0.0;
		for(int ns=0; ns<nSp; ++ns){
			Ht += MF[ns]*Hti[ns];
			drdp += VF[ns]*drhodPi[ns];
			drdt += VF[ns]*drhodTi[ns];
			dhdp += MF[ns]*dhdPi[ns];
			dhdt += MF[ns]*dhdTi[ns];
		}



		for(int ns=0; ns<nSp; ++ns){
			cell.var[controls.MF[ns]] = MF[ns];
			cell.var[controls.VF[ns]] = VF[ns];
		}
		
		cell.var[controls.Ht] = Ht;
		cell.var[controls.dRhoDP] = drdp;
		cell.var[controls.dRhoDT] = drdt;
		cell.var[controls.dHtDP] = dhdp;
		cell.var[controls.dHtDT] = dhdt;
		
		for(int ns=0; ns<nSp-1; ++ns){
			cell.var[controls.dRhoDMF[ns]] = -rho*rho*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]);
			// cout << cell.var[controls.dRhoDMF[ns]] << endl;
			cell.var[controls.dHtDMF[ns]] = dhdTi[ns] - dhdTi[nSp-1];
		}
		
		// cout << drdp << " " << rho << " " << drdt << " " << dhdt << " " << dhdp << endl;
		
		cell.var[controls.C] = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
		cell.var[controls.C] = sqrt( 1.0 / cell.var[controls.C] );
	
		// cout << cell.var[controls.C] << endl;
		
	}
	
	
	
	
}







void SEMO_Solvers_Builder::calcCellEOSMF(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species){
	
	int nSp = controls.nSp;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		SEMO_Cell& cell = mesh.cells[i];
		
		vector<double> rhoi(nSp,0.0);
		vector<double> Ci(nSp,0.0);
		vector<double> Hti(nSp,0.0);
		vector<double> drhodPi(nSp,0.0);
		vector<double> drhodTi(nSp,0.0);
		vector<double> dhdPi(nSp,0.0);
		vector<double> dhdTi(nSp,0.0);

		for(int ns=0; ns<nSp; ++ns){
			if( species[ns].rhoType == "const" ){
				rhoi[ns] = species[ns].rho;
				Ci[ns] = 340.0;
				Hti[ns] = 10000.0;
				drhodPi[ns] = 1.0e+8;
				drhodTi[ns] = 1.0e+8;
				dhdPi[ns] = 1.0e+8;
				dhdTi[ns] = 1.0e+8;
			}
			else if( species[ns].rhoType == "Ideal" ){
				this->eosIdeal(
						species[ns],
						cell.var[controls.P], 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						cell.var[controls.T],
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
			else if( species[ns].rhoType == "NASG" ){
				this->eosNASG(
						species[ns],
						cell.var[controls.P], 
						cell.var[controls.U], cell.var[controls.V], cell.var[controls.W],
						cell.var[controls.T],
						rhoi[ns], Ci[ns], Hti[ns],
						drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
			}
				
		}
		
		vector<double> MF;
		double tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			double MF_tmp = max(0.0,min(1.0,cell.var[controls.MF[ns]]));
			MF.push_back( max(0.0,min(1.0,MF_tmp)) );
			tmp1 += MF_tmp;
		}
		MF.push_back( max(0.0,min(1.0,1.0 - tmp1)) );
		
		
		double rho = 0.0;
		for(int ns=0; ns<nSp; ++ns){
			rho += MF[ns]/rhoi[ns];
		}
		rho = 1.0 / rho;
		
		cell.var[controls.Rho] = rho;
		
		
		// volume fraction
		vector<double> VF(nSp,0.0);
		tmp1=0.0;
		for(int ns=0; ns<nSp-1; ++ns){
			VF[ns] = rho*MF[ns]/rhoi[ns];
			
			VF[ns] = max(0.0,min(1.0,VF[ns]));
			
			tmp1 += VF[ns];
		}
		VF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
		
		double Ht=0.0;
		double drdp=0.0;
		double drdt=0.0;
		double dhdp=0.0;
		double dhdt=0.0;
		for(int ns=0; ns<nSp; ++ns){
			Ht += MF[ns]*Hti[ns];
			drdp += rho*rho*(MF[ns]/rhoi[ns]/rhoi[ns]*drhodPi[ns]);
			drdt += rho*rho*(MF[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
			dhdp += MF[ns]*dhdPi[ns];
			dhdt += MF[ns]*dhdTi[ns];
		}


		for(int ns=0; ns<nSp; ++ns){
			cell.var[controls.MF[ns]] = MF[ns];
			cell.var[controls.VF[ns]] = VF[ns];
		}
		

		cell.var[controls.Ht] = Ht;
		cell.var[controls.dRhoDP] = drdp;
		cell.var[controls.dRhoDT] = drdt;
		cell.var[controls.dHtDP] = dhdp;
		cell.var[controls.dHtDT] = dhdt;
		
		for(int ns=0; ns<nSp-1; ++ns){
			cell.var[controls.dRhoDMF[ns]] = -rho*rho*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]);
			cell.var[controls.dHtDMF[ns]] = Hti[ns] - Hti[nSp-1];
		}
		
		cell.var[controls.C] = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
		cell.var[controls.C] = sqrt( 1.0 / cell.var[controls.C] );
	
		
	}
	
	
	
	
}








void SEMO_Solvers_Builder::getValuesFromEOSVF(
	vector<SEMO_Species>& species,
	double& P, double& U, double& V, double& W, double& T, vector<double>& VF,
	double& rho, double& C, double& Ht,
	vector<double>& MF){

	int nSp = VF.size();

	for(int ns=0; ns<nSp; ++ns){
		VF[ns] = max(0.0,min(1.0,VF[ns]));
	}
	
	vector<double> rhoi(nSp,0.0);
	vector<double> Ci(nSp,0.0);
	vector<double> Hti(nSp,0.0);
	vector<double> drhodPi(nSp,0.0);
	vector<double> drhodTi(nSp,0.0);
	vector<double> dhdPi(nSp,0.0);
	vector<double> dhdTi(nSp,0.0);

	for(int ns=0; ns<nSp; ++ns){
		if( species[ns].rhoType == "const" ){
			rhoi[ns] = species[ns].rho;
			Ci[ns] = 340.0;
			Hti[ns] = 10000.0;
			drhodPi[ns] = 1.0e+8;
			drhodTi[ns] = 1.0e+8;
			dhdPi[ns] = 1.0e+8;
			dhdTi[ns] = 1.0e+8;
		}
		else if( species[ns].rhoType == "Ideal" ){
			this->eosIdeal(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
		else if( species[ns].rhoType == "NASG" ){
			this->eosNASG(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
			
	}
	
	rho = 0.0;
	for(int ns=0; ns<nSp; ++ns){
		rho += VF[ns]*rhoi[ns];
	}
	
	
	// mass fraction
	MF.resize(nSp,0.0);
	double tmp1=0.0;
	for(int ns=0; ns<nSp-1; ++ns){
		MF[ns] = VF[ns]*rhoi[ns]/rho;
		
		MF[ns] = max(0.0,min(1.0,MF[ns]));
		
		tmp1 += MF[ns];
	}
	MF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
	
	Ht=0.0;
	double drdp=0.0;
	double drdt=0.0;
	double dhdp=0.0;
	double dhdt=0.0;
	for(int ns=0; ns<nSp; ++ns){
		Ht += MF[ns]*Hti[ns];
		drdp += VF[ns]*drhodPi[ns];
		drdt += VF[ns]*drhodTi[ns];
		dhdp += MF[ns]*dhdPi[ns];
		dhdt += MF[ns]*dhdTi[ns];
	}
	
	C = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
	C = sqrt( 1.0 / C );
	
	// forall(mc=1:nspec-1) dhdy(i,mc) = hti(mc) - hti(nspec)
	// forall(mc=1:nspec-1) drdy(i,mc) = -rho(i)**r2*(r1/rho_i(mc)-r1/rho_i(nspec))


}








void SEMO_Solvers_Builder::getValuesFromEOSMF(
	vector<SEMO_Species>& species,
	double& P, double& U, double& V, double& W, double& T, vector<double>& MF,
	double& rho, double& C, double& Ht,
	vector<double>& VF){

	int nSp = MF.size();
	
	for(int ns=0; ns<nSp; ++ns){
		MF[ns] = max(0.0,min(1.0,MF[ns]));
	}
	
	vector<double> rhoi(nSp,0.0);
	vector<double> Ci(nSp,0.0);
	vector<double> Hti(nSp,0.0);
	vector<double> drhodPi(nSp,0.0);
	vector<double> drhodTi(nSp,0.0);
	vector<double> dhdPi(nSp,0.0);
	vector<double> dhdTi(nSp,0.0);

	for(int ns=0; ns<nSp; ++ns){
		if( species[ns].rhoType == "const" ){
			rhoi[ns] = species[ns].rho;
			Ci[ns] = 340.0;
			Hti[ns] = 10000.0;
			drhodPi[ns] = 1.0e+8;
			drhodTi[ns] = 1.0e+8;
			dhdPi[ns] = 1.0e+8;
			dhdTi[ns] = 1.0e+8;
		}
		else if( species[ns].rhoType == "Ideal" ){
			this->eosIdeal(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
		else if( species[ns].rhoType == "NASG" ){
			this->eosNASG(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
			
	}
	
	rho = 0.0;
	for(int ns=0; ns<nSp; ++ns){
		rho += MF[ns]/rhoi[ns];
	}
	rho = 1.0/rho;
	
	
	// mass fraction
	VF.resize(nSp,0.0);
	double tmp1=0.0;
	for(int ns=0; ns<nSp-1; ++ns){
		VF[ns] = MF[ns]*rho/rhoi[ns];
		
		VF[ns] = max(0.0,min(1.0,VF[ns]));
		
		tmp1 += VF[ns];
	}
	VF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
	
	Ht=0.0;
	double drdp=0.0;
	double drdt=0.0;
	double dhdp=0.0;
	double dhdt=0.0;
	for(int ns=0; ns<nSp; ++ns){
		Ht += MF[ns]*Hti[ns];
		drdp += VF[ns]*drhodPi[ns];
		drdt += VF[ns]*drhodTi[ns];
		dhdp += MF[ns]*dhdPi[ns];
		dhdt += MF[ns]*dhdTi[ns];
	}
	
	C = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
	C = sqrt( 1.0 / C );
	
	// forall(mc=1:nspec-1) dhdy(i,mc) = hti(mc) - hti(nspec)
	// forall(mc=1:nspec-1) drdy(i,mc) = -rho(i)**r2*(r1/rho_i(mc)-r1/rho_i(nspec))


}




void SEMO_Solvers_Builder::getValuesFromEOSMF(
	vector<SEMO_Species>& species,
	double& P, double& U, double& V, double& W, double& T, vector<double>& MF, 
	vector<double>& VF, double& rho, double& C, double& Ht,
	double& drdp, double& dhdp,
	double& drdt, double& dhdt,
	vector<double>& drdMF, vector<double>& dhdMF){

	int nSp = MF.size();
	
	for(int ns=0; ns<nSp; ++ns){
		MF[ns] = max(0.0,min(1.0,MF[ns]));
	}
	
	vector<double> rhoi(nSp,0.0);
	vector<double> Ci(nSp,0.0);
	vector<double> Hti(nSp,0.0);
	vector<double> drhodPi(nSp,0.0);
	vector<double> drhodTi(nSp,0.0);
	vector<double> dhdPi(nSp,0.0);
	vector<double> dhdTi(nSp,0.0);

	for(int ns=0; ns<nSp; ++ns){
		if( species[ns].rhoType == "const" ){
			rhoi[ns] = species[ns].rho;
			Ci[ns] = 340.0;
			Hti[ns] = 10000.0;
			drhodPi[ns] = 1.0e+8;
			drhodTi[ns] = 1.0e+8;
			dhdPi[ns] = 1.0e+8;
			dhdTi[ns] = 1.0e+8;
		}
		else if( species[ns].rhoType == "Ideal" ){
			this->eosIdeal(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
		else if( species[ns].rhoType == "NASG" ){
			this->eosNASG(
				species[ns],
				P, U, V, W, T,
				rhoi[ns], Ci[ns], Hti[ns],
				drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		}
			
	}
	
	rho = 0.0;
	for(int ns=0; ns<nSp; ++ns){
		rho += MF[ns]/rhoi[ns];
	}
	rho = 1.0/rho;
	
	
	// mass fraction
	VF.resize(nSp,0.0);
	double tmp1=0.0;
	for(int ns=0; ns<nSp-1; ++ns){
		VF[ns] = MF[ns]*rho/rhoi[ns];
		
		VF[ns] = max(0.0,min(1.0,VF[ns]));
		
		tmp1 += VF[ns];
	}
	VF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
	
	Ht=0.0;
	drdp=0.0;
	drdt=0.0;
	dhdp=0.0;
	dhdt=0.0;
	drdMF.clear();
	dhdMF.clear();
	for(int ns=0; ns<nSp; ++ns){
		Ht += MF[ns]*Hti[ns];
		drdp += rho*rho*(MF[ns]/rhoi[ns]/rhoi[ns]*drhodPi[ns]);
		drdt += rho*rho*(MF[ns]/rhoi[ns]/rhoi[ns]*drhodTi[ns]);
		dhdp += MF[ns]*dhdPi[ns];
		dhdt += MF[ns]*dhdTi[ns];
		
		drdMF.push_back( -rho*rho*(1.0/rhoi[ns]-1.0/rhoi[nSp-1]) );
		dhdMF.push_back( Hti[ns]-Hti[nSp-1] );
	}
	
	C = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
	C = sqrt( 1.0 / C );
	
	// return tuple(VF, rho, C, Ht, drdp, dhdp, drdt, dhdt, drdMF, dhdMF);
	

}





// void SEMO_Solvers_Builder::getValuesFromEOSMF(
	// vector<SEMO_Species>& species,
	// double& P, double& U, double& V, double& W, double& T, vector<double>& MF,
	// double& rho, double& C, double& Ht){

	// int nSp = MF.size();
	
	// vector<double> rhoi(nSp,0.0);
	// vector<double> Ci(nSp,0.0);
	// vector<double> Hti(nSp,0.0);
	// vector<double> drhodPi(nSp,0.0);
	// vector<double> drhodTi(nSp,0.0);
	// vector<double> dhdPi(nSp,0.0);
	// vector<double> dhdTi(nSp,0.0);

	// for(int ns=0; ns<nSp; ++ns){
		// if( species[ns].rhoType == "const" ){
			// rhoi[ns] = species[ns].rho;
			// Ci[ns] = 340.0;
			// Hti[ns] = 10000.0;
			// drhodPi[ns] = 1.0e+8;
			// drhodTi[ns] = 1.0e+8;
			// dhdPi[ns] = 1.0e+8;
			// dhdTi[ns] = 1.0e+8;
		// }
		// else if( species[ns].rhoType == "Ideal" ){
			// this->eosIdeal(
				// species[ns],
				// P, U, V, W, T,
				// rhoi[ns], Ci[ns], Hti[ns],
				// drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		// }
		// else if( species[ns].rhoType == "NASG" ){
			// this->eosNASG(
				// species[ns],
				// P, U, V, W, T,
				// rhoi[ns], Ci[ns], Hti[ns],
				// drhodPi[ns], drhodTi[ns], dhdPi[ns], dhdTi[ns]);
		// }
			
	// }
	
	// rho = 0.0;
	// for(int ns=0; ns<nSp; ++ns){
		// rho += MF[ns]/rhoi[ns];
	// }
	// rho = 1.0 / rho;
	
	
	// // volume fraction
	// vector<double> VF(nSp,0.0);
	// double tmp1=0.0;
	// for(int ns=0; ns<nSp-1; ++ns){
		// VF[ns] = rho*MF[ns]/rhoi[ns];
		// tmp1 += VF[ns];
	// }
	// VF[nSp-1] = max(0.0,min(1.0,1.0 - tmp1));
	
	// Ht=0.0;
	// double drdp=0.0;
	// double drdt=0.0;
	// double dhdp=0.0;
	// double dhdt=0.0;
	// for(int ns=0; ns<nSp; ++ns){
		// Ht += MF[ns]*Hti[ns];
		// drdp += VF[ns]*drhodPi[ns];
		// drdt += VF[ns]*drhodTi[ns];
		// dhdp += MF[ns]*dhdPi[ns];
		// dhdt += MF[ns]*dhdTi[ns];
	// }
	
	// C = drdp + 1.0/rho*drdt/dhdt*(1.0-rho*dhdp);
	// C = sqrt( 1.0 / C );
	
	// // forall(mc=1:nspec-1) dhdy(i,mc) = hti(mc) - hti(nspec)
	// // forall(mc=1:nspec-1) drdy(i,mc) = -rho(i)**r2*(r1/rho_i(mc)-r1/rho_i(nspec))


// }




void SEMO_Solvers_Builder::eosIdeal(
	SEMO_Species& species,
	double& P, double& U, double& V, double& W, double& T,
	double& rho, double& C, double& Ht,
	double& drhodP, double& drhodT, double& dhdP, double& dhdT){
		
		
	double usqrt = U*U+V*V+W*W;
	
		
	double cv = species.cv;
    double gam = species.gamma;
	
	double cp = gam*cv;
		
	// density of each phase
	rho = 1.0/( (gam-1.0)*cv*T/(P) );
	C = sqrt( gam/(rho*rho)*(P)/(1.0/rho) );

	// d(rho)/d(p)
	drhodP = rho*rho*(1.0/rho)/(P); 
	
	// d(rho)/d(T)
	drhodT = -rho*rho*(1.0/rho)/T;

	// d(h)/d(p)
	dhdP = 0.0;
	// d(h)/d(T)
	dhdT = gam*cv;

	// internal energy of each phase
	// internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	// enthalpy of each phase
	double enthalpy = gam*cv*T;

		   
	// eti = internal_energy + 0.5*usqrt
	Ht = enthalpy + 0.5*usqrt;

	// cvi = cv
	// cpi = cp	


		
		
}


void SEMO_Solvers_Builder::eosNASG(
	SEMO_Species& species,
	double& P, double& U, double& V, double& W, double& T,
	double& rho, double& C, double& Ht,
	double& drhodP, double& drhodT, double& dhdP, double& dhdT){
		
		
	double usqrt = U*U+V*V+W*W;
	
		
    double pinf = species.Pinf; 
	double cv = species.cv;
    double gam = species.gamma;
    double bNASG = species.b; 
	double q = species.q;
	
	double cp = gam*cv;
		
	// density of each phase
	rho = 1.0/( (gam-1.0)*cv*T/(P+pinf)+bNASG );
	C = sqrt( gam/(rho*rho)*(P+pinf)/(1.0/rho-bNASG) );

	// d(rho)/d(p)
	drhodP = rho*rho*(1.0/rho-bNASG)/(P+pinf); 
	
	// d(rho)/d(T)
	drhodT = -rho*rho*(1.0/rho-bNASG)/T;

	// d(h)/d(p)
	dhdP = bNASG;
	// d(h)/d(T)
	dhdT = gam*cv;

	// internal energy of each phase
	// internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	// enthalpy of each phase
	double enthalpy = gam*cv*T + bNASG*P + q;

		   
	// eti = internal_energy + 0.5*usqrt
	Ht = enthalpy + 0.5*usqrt;

	// cvi = cv
	// cpi = cp	


		
		
}