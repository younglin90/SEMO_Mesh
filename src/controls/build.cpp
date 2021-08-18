#include "build.h"
#include <cmath>
#include <array>
#include <time.h>


void SEMO_Controls_Builder::readSpecies(vector<SEMO_Species>& species){
	
	
	SEMO_Utility_Read read;


	// thermophysicalProperties
	map<string,string> thermophysicalProperties;
	read.file("./constant/thermophysicalProperties",thermophysicalProperties);
	
	vector<string> tmpspecies;
	read.vecters(thermophysicalProperties["phases"], tmpspecies);
	
	
	for(int i=0; i<tmpspecies.size(); ++i){
		SEMO_Species tmp;
		species.push_back(tmp);
	}
	
	for(int i=0; i<tmpspecies.size(); ++i){
		species[i].name = read.trim(tmpspecies[i]);
	}
	
	for(int i=0; i<tmpspecies.size(); ++i){
		map<string,string> tmp;
		read.file("./constant/thermophysicalProperties",
			tmpspecies[i], "thermodynamics", "rho", tmp);
		
		// cout << tmpspecies[i] << " " << tmp["value"] << endl;

		species[i].rhoType = read.trim(tmp["type"]);
		if( tmp["type"] == "const" ){
			species[i].rho = stod(tmp["value"]);
		}
		else if( tmp["type"] == "NASG" ){
			species[i].Pinf = stod(tmp["Pinf"]);
			species[i].cv = stod(tmp["cv"]);
			species[i].gamma = stod(tmp["gamma"]);
			species[i].b = stod(tmp["b"]);
			species[i].q = stod(tmp["q"]);
		}
		else if( tmp["type"] == "Ideal" ){
			species[i].cv = stod(tmp["cv"]);
			species[i].gamma = stod(tmp["gamma"]);
		}
	}
	
	for(int i=0; i<tmpspecies.size(); ++i){
		
		// cout << " AAA " << endl;
		
		map<string,string> tmp;
		read.file("./constant/thermophysicalProperties",
			tmpspecies[i], "transport", "mu", tmp);
		
		// cout << species.name[i] << " " << endl;
		
		if( tmp["type"] == "const" ){
			species[i].mu = stod(tmp["value"]);
		}
		
		
		map<string,string> tmp2;
		read.file("./constant/thermophysicalProperties",
			tmpspecies[i], "transport", "sigma", tmp2);
		
		if( tmp2["type"] == "const" ){
			species[i].sigma = stod(tmp2["value"]);
		}
	}
	
	
	
	
	
}



	
void SEMO_Controls_Builder::readConfigures(){
	
	
	SEMO_Utility_Read read;

	vector<string> tmp;
	
	// turbulenceProperties
	map<string,string> turbulenceProperties;
	read.file("./constant/turbulenceProperties",turbulenceProperties);
	
	this->turbType = turbulenceProperties["simulationType"];
	
	if(this->turbType == "LES"){
		map<string,string> turbulenceProperties_LES;
		read.file("./constant/turbulenceProperties", "LES", turbulenceProperties_LES);
		this->turbLESModel = turbulenceProperties_LES["LESModel"];
	}
	else if(this->turbType == "RANS"){
		map<string,string> turbulenceProperties_RANS;
		read.file("./constant/turbulenceProperties", "RANS", turbulenceProperties_RANS);
		this->turbRANSModel = turbulenceProperties_RANS["RANSModel"];
	}
	else if(this->turbType == "DES"){
		map<string,string> turbulenceProperties_LES;
		read.file("./constant/turbulenceProperties", "LES", turbulenceProperties_LES);
		map<string,string> turbulenceProperties_RANS;
		read.file("./constant/turbulenceProperties", "RANS", turbulenceProperties_RANS);
		this->turbLESModel = turbulenceProperties_LES["LESModel"];
		this->turbRANSModel = turbulenceProperties_RANS["RANSModel"];
	}
	
	
	
	// controlDict
	map<string,string> controlDict;
	read.file("./system/controlDict",controlDict);

	this->application = controlDict["application"];
	this->startFrom = controlDict["startFrom"];
	this->stopAt = controlDict["stopAt"];
	this->timeStep = stod(controlDict["timeStep"]);
	this->orgTimeStep = this->timeStep;
	this->saveInterval = stoi(controlDict["saveInterval"]);
	this->saveControl = controlDict["saveControl"];
	this->saveFormat = controlDict["saveFormat"];
	this->saveCompression = controlDict["saveCompression"];
	this->timeFormat = controlDict["timeFormat"];
	this->adjustTimeStep = controlDict["adjustTimeStep"];
	this->maxCo = stod(controlDict["maxCo"]);
	this->maxVfCo = stod(controlDict["maxVfCo"]);
	this->maxTimeStep = stod(controlDict["maxTimeStep"]);
	this->writePrecision = stoi(controlDict["writePrecision"]);


	// gravity acceleration
	map<string,string> gravity;
	read.file("./constant/g",gravity);
	read.vecters(gravity["value"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->gravityAcceleration.push_back(stod(tmp[i]));
		// cout << gravityAcceleration.back() << endl;
	}
	
	


	
	// fvSchemes
	map<string,string> fvSchemes_fluxScheme;
	read.file("./system/fvSchemes", "fluxScheme", fvSchemes_fluxScheme);
	
	map<string,string> fvSchemes_ddtSchemes;
	read.file("./system/fvSchemes", "ddtSchemes", fvSchemes_ddtSchemes);
	
	map<string,string> fvSchemes_gradSchemes;
	read.file("./system/fvSchemes", "gradSchemes", fvSchemes_gradSchemes);
	
	map<string,string> fvSchemes_divSchemes;
	read.file("./system/fvSchemes", "divSchemes", fvSchemes_divSchemes);
	
	map<string,string> fvSchemes_laplacianSchemes;
	read.file("./system/fvSchemes", "laplacianSchemes", fvSchemes_laplacianSchemes);
	
	map<string,string> fvSchemes_interpolationSchemes;
	read.file("./system/fvSchemes", "interpolationSchemes", fvSchemes_interpolationSchemes);
	
	map<string,string> fvSchemes_snGradSchemes;
	read.file("./system/fvSchemes", "snGradSchemes", fvSchemes_snGradSchemes);
	

	this->fluxScheme = fvSchemes_fluxScheme["default"];
	
	this->gradScheme = fvSchemes_gradSchemes["default"];
	
	this->divScheme = fvSchemes_divSchemes["default"];
	
	// cout << this->fluxScheme << this->divScheme << " " << this->gradScheme << " " << endl;
	
	
	
	
	
	
	
	// fvSolution
	map<string,string> fvSolution_solvers_P;
	read.file("./system/fvSolution", 
		"solvers", "P", 
		fvSolution_solvers_P);

	this->solverP = read.trim(fvSolution_solvers_P["solver"]);
	this->toleranceP = stod(fvSolution_solvers_P["tolerance"]);
	this->relTolP = stod(fvSolution_solvers_P["relTol"]);
	this->preconditionerP = read.trim(fvSolution_solvers_P["preconditioner"]);
	this->maxIterP = stoi(fvSolution_solvers_P["maxIter"]);
	
	this->solverFinalP = read.trim(fvSolution_solvers_P["solverFinal"]);
	this->toleranceFinalP = stod(fvSolution_solvers_P["toleranceFinal"]);
	this->relTolFinalP = stod(fvSolution_solvers_P["relTolFinal"]);
	this->preconditionerFinalP = read.trim(fvSolution_solvers_P["preconditionerFinal"]);
	this->maxIterFinalP = stoi(fvSolution_solvers_P["maxIterFinal"]);
		
		
	map<string,string> fvSolution_solvers_U;
	read.file("./system/fvSolution", 
		"solvers", "U", 
		fvSolution_solvers_U);

	this->solverU = read.trim(fvSolution_solvers_U["solver"]);
	this->toleranceU = stod(fvSolution_solvers_U["tolerance"]);
	this->relTolU = stod(fvSolution_solvers_U["relTol"]);
	this->preconditionerU = read.trim(fvSolution_solvers_U["preconditioner"]);
	this->maxIterU = stoi(fvSolution_solvers_U["maxIter"]);
	
	this->solverFinalU = read.trim(fvSolution_solvers_U["solverFinal"]);
	this->toleranceFinalU = stod(fvSolution_solvers_U["toleranceFinal"]);
	this->relTolFinalU = stod(fvSolution_solvers_U["relTolFinal"]);
	this->preconditionerFinalU = read.trim(fvSolution_solvers_U["preconditionerFinal"]);
	this->maxIterFinalU = stoi(fvSolution_solvers_U["maxIterFinal"]);
		
		
		
	map<string,string> fvSolution_solvers_VF;
	read.file("./system/fvSolution", 
		"solvers", "VF", 
		fvSolution_solvers_VF);

	this->solverVF.push_back(read.trim(fvSolution_solvers_VF["solver"]));
	this->toleranceVF.push_back(stod(fvSolution_solvers_VF["tolerance"]));
	this->relTolVF.push_back(stod(fvSolution_solvers_VF["relTol"]));
	this->preconditionerVF.push_back(read.trim(fvSolution_solvers_VF["preconditioner"]));
	this->maxIterVF.push_back(stoi(fvSolution_solvers_VF["maxIter"]));
	
	this->solverFinalVF.push_back(read.trim(fvSolution_solvers_VF["solverFinal"]));
	this->toleranceFinalVF.push_back(stod(fvSolution_solvers_VF["toleranceFinal"]));
	this->relTolFinalVF.push_back(stod(fvSolution_solvers_VF["relTolFinal"]));
	this->preconditionerFinalVF.push_back(read.trim(fvSolution_solvers_VF["preconditionerFinal"]));
	this->maxIterFinalVF.push_back(stoi(fvSolution_solvers_VF["maxIterFinal"]));
	
	
		
	map<string,string> fvSolution_PIMPLE;
	read.file("./system/fvSolution", 
		"PIMPLE",  
		fvSolution_PIMPLE);
		
	map<string,string> fvSolution_relaxationFactors_momentumEq_U;
	read.file("./system/fvSolution", 
		"relaxationFactors", "momentumEq", "U",
		fvSolution_relaxationFactors_momentumEq_U);
		
	map<string,string> fvSolution_relaxationFactors_pressureEq_P;
	read.file("./system/fvSolution", 
		"relaxationFactors", "pressureEq", "P",
		fvSolution_relaxationFactors_pressureEq_P);
		
	map<string,string> fvSolution_relaxationFactors_pressureEq_U;
	read.file("./system/fvSolution", 
		"relaxationFactors", "pressureEq", "U",
		fvSolution_relaxationFactors_pressureEq_U);
		
	map<string,string> fvSolution_relaxationFactors_speciesEq_VF;
	read.file("./system/fvSolution", 
		"relaxationFactors", "speciesEq", "VF",
		fvSolution_relaxationFactors_speciesEq_VF);






	map<string,string> fvSolution_dualTime;
	read.file("./system/fvSolution", 
		"dualTime",  
		fvSolution_dualTime);
		
	map<string,string> fvSolution_relaxationFactorsDualTime_P;
	read.file("./system/fvSolution", 
		"relaxationFactorsDualTime", "P", 
		fvSolution_relaxationFactorsDualTime_P);
		
	map<string,string> fvSolution_relaxationFactorsDualTime_U;
	read.file("./system/fvSolution", 
		"relaxationFactorsDualTime", "P", 
		fvSolution_relaxationFactorsDualTime_U);
		
	map<string,string> fvSolution_relaxationFactorsDualTime_T;
	read.file("./system/fvSolution", 
		"relaxationFactorsDualTime", "P", 
		fvSolution_relaxationFactorsDualTime_T);
		
	map<string,string> fvSolution_relaxationFactorsDualTime_MF;
	read.file("./system/fvSolution", 
		"relaxationFactorsDualTime", "P", 
		fvSolution_relaxationFactorsDualTime_MF);






	map<string,string> fvSolution_limiter;
	read.file("./system/fvSolution", 
		"limiter",  
		fvSolution_limiter);
		

	
	// URF
	this->iterPBsMax = stoi(fvSolution_PIMPLE["nCorrectors"]);
	this->iterMomMax = stoi(fvSolution_PIMPLE["nMomentumEq"]);
	this->iterPreMax = stoi(fvSolution_PIMPLE["nPressureEq"]);
	this->iterVofMax = stoi(fvSolution_PIMPLE["nSpeciesEq"]);
	
	
	this->iterPseudoMax = stoi(fvSolution_dualTime["nPseudo"]);
	// cout << fvSolution_dualTime["nPseudo"] << " " << fvSolution_dualTime["pseudoCo"] << endl;
	this->pseudoCo = stod(fvSolution_dualTime["pseudoCo"]);
	
	this->Uco = stod(fvSolution_dualTime["Uco"]);
	this->Lch = stod(fvSolution_dualTime["Lch"]);

	// relaxationFactors
	this->momVelURF = stod(fvSolution_relaxationFactors_momentumEq_U["default"]);
	this->momVelAdjustRF = fvSolution_relaxationFactors_momentumEq_U["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_momentumEq_U["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->momVelAdjustSteps.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_momentumEq_U["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->momVelAdjustValues.push_back(stod(tmp[i]));
	}
	
	
	this->prePreURF = stod(fvSolution_relaxationFactors_pressureEq_P["default"]);
	this->prePreAdjustRF = fvSolution_relaxationFactors_pressureEq_P["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_pressureEq_P["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->prePreAdjustSteps.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_pressureEq_P["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->prePreAdjustValues.push_back(stod(tmp[i]));
	}
	
	
	this->preVelURF = stod(fvSolution_relaxationFactors_pressureEq_U["default"]);
	this->preVelAdjustRF = fvSolution_relaxationFactors_pressureEq_U["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_pressureEq_U["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->preVelAdjustSteps.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_pressureEq_U["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->preVelAdjustValues.push_back(stod(tmp[i]));
	}
	
	
	this->vofVofURF = stod(fvSolution_relaxationFactors_speciesEq_VF["default"]);
	this->vofVofAdjustRF = fvSolution_relaxationFactors_speciesEq_VF["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_speciesEq_VF["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->vofVofAdjustSteps.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactors_speciesEq_VF["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->vofVofAdjustValues.push_back(stod(tmp[i]));
	}
	




	//===================
	this->dualTimeURF_P = stod(fvSolution_relaxationFactorsDualTime_P["default"]);
	this->dualTimeAdjustRF_P = fvSolution_relaxationFactorsDualTime_P["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_P["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustSteps_P.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_P["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustValues_P.push_back(stod(tmp[i]));
	}
	
	
	this->dualTimeURF_U = stod(fvSolution_relaxationFactorsDualTime_U["default"]);
	this->dualTimeAdjustRF_U = fvSolution_relaxationFactorsDualTime_U["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_U["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustSteps_U.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_U["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustValues_U.push_back(stod(tmp[i]));
	}
	
	
	this->dualTimeURF_T = stod(fvSolution_relaxationFactorsDualTime_T["default"]);
	this->dualTimeAdjustRF_T = fvSolution_relaxationFactorsDualTime_T["adjustRF"];
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_T["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustSteps_T.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_T["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustValues_T.push_back(stod(tmp[i]));
	}
	
	this->dualTimeURF_MF = stod(fvSolution_relaxationFactorsDualTime_MF["default"]);
	this->dualTimeAdjustRF_MF = fvSolution_relaxationFactorsDualTime_MF["adjustRF"]; 
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_MF["adjustSteps"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustSteps_MF.push_back(stoi(tmp[i]));
	}
	
	tmp.clear();
	read.vecters(fvSolution_relaxationFactorsDualTime_MF["adjustValues"], tmp);
	for(int i=0; i<tmp.size(); ++i){
		this->dualTimeAdjustValues_MF.push_back(stod(tmp[i]));
	}
	
	//===========================






	// limiter
	this->maxP = stod(fvSolution_limiter["maxP"]);
	this->minP = stod(fvSolution_limiter["minP"]);
	this->maxU = stod(fvSolution_limiter["maxU"]);
	this->minU = stod(fvSolution_limiter["minU"]);
	this->maxV = stod(fvSolution_limiter["maxV"]);
	this->minV = stod(fvSolution_limiter["minV"]);
	this->maxW = stod(fvSolution_limiter["maxW"]);
	this->minW = stod(fvSolution_limiter["minW"]);
	this->maxT = stod(fvSolution_limiter["maxT"]);
	this->minT = stod(fvSolution_limiter["minT"]);
	
	// this->timestep = 1.e-5;
		
	// momVelURF = 0.1;
	// prePreURF = 0.1;
	// preVelURF = 0.1;
	// vofVofURF = 0.1;
	
	// iterator
	this->iterRealMax = std::numeric_limits<int>::max();
	this->iterReal = 0;
	
	this->iterTotal = 0;
	this->startClock = clock();
	
	
	// this->iterPseudoMax = 15;

}









void SEMO_Controls_Builder::setValues(vector<SEMO_Species>& species){
	



	this->nSp = species.size();
	
	this->nEq = 5 + this->nSp-1;
	
	
	
	
	// cell var
	this->nTotalCellVar = 0;
	
	this->P = this->nTotalCellVar++; 
	this->name.push_back("pressure");
	// ++this->nEq;
	this->U = this->nTotalCellVar++; 
	this->name.push_back("x-velocity");
	// ++this->nEq;
	this->V = this->nTotalCellVar++;
	this->name.push_back("y-velocity");
	// ++this->nEq;
	this->W = this->nTotalCellVar++;
	this->name.push_back("z-velocity");
	// ++this->nEq;
	this->T = this->nTotalCellVar++;
	this->name.push_back("temperature");
	// ++this->nEq;
	
	// int nvar = this->nTotalCellVar;
	
	for(int i=0; i<this->nSp; ++i){
		this->VF.push_back(this->nTotalCellVar++);
		string name = "volume-fraction-" + species[i].name;
		this->name.push_back(name);
	}
	
	for(int i=0; i<this->nSp; ++i){
		this->MF.push_back(this->nTotalCellVar++);
		string name = "mass-fraction-" + species[i].name;
		this->name.push_back(name);
	}


	this->Rho = this->nTotalCellVar++;
	this->name.push_back("density");
	
	this->C = this->nTotalCellVar++;
	this->name.push_back("speed-of-sound");
	
	this->Ht = this->nTotalCellVar++;
	this->name.push_back("total-enthalpy");
	
	this->dRhoDP = this->nTotalCellVar++;
	this->name.push_back("dRhoDP");
	
	this->dHtDP = this->nTotalCellVar++;
	this->name.push_back("dHtDP");
	
	this->dRhoDT = this->nTotalCellVar++;
	this->name.push_back("dRhoDT");
	
	this->dHtDT = this->nTotalCellVar++;
	this->name.push_back("dHtDT");
	
	for(int i=0; i<this->nSp; ++i){
		this->dRhoDVF.push_back(this->nTotalCellVar++);
		this->name.push_back("dRhoDVF");
	
		this->dHtDVF.push_back(this->nTotalCellVar++);
		this->name.push_back("dHtDVF");
	
	}
	for(int i=0; i<this->nSp; ++i){
		this->dRhoDMF.push_back(this->nTotalCellVar++);
		this->name.push_back("dRhoDMF");
	
		this->dHtDMF.push_back(this->nTotalCellVar++);
		this->name.push_back("dHtDMF");
	
	}
	
	// for(int i=0; i<this->nEq; ++i){
		// this->Q.push_back(this->nTotalCellVar++);
	// }
	for(int i=0; i<this->nEq; ++i){
		this->Qn.push_back(this->nTotalCellVar++);
	}
	for(int i=0; i<this->nEq; ++i){
		this->Qm.push_back(this->nTotalCellVar++);
	}

	this->oldP = this->nTotalCellVar++;
	this->name.push_back("pressure-old");
	
	this->oldU = this->nTotalCellVar++;
	this->name.push_back("x-velocity-old");
	
	this->oldV = this->nTotalCellVar++;
	this->name.push_back("y-velocity-old");
	
	this->oldW = this->nTotalCellVar++;
	this->name.push_back("z-velocity-old");
	
	this->oldT = this->nTotalCellVar++;
	this->name.push_back("temperature-old");
	
	for(int i=0; i<this->nSp-1; ++i){
		this->oldVF.push_back(this->nTotalCellVar++);
		string name = "volume-fraction-" + species[i].name;
		name.append("-old");
		this->name.push_back(name);
	}
	
	for(int i=0; i<this->nSp-1; ++i){
		this->oldMF.push_back(this->nTotalCellVar++);
		string name = "mass-fraction-" + species[i].name;
		name.append("-old");
		this->name.push_back(name);
	}
	
	// cout << "AAAAAAAAAA" << endl;
	
	this->oldRho = this->nTotalCellVar++;
	this->name.push_back("density-old");
	
	
	// cout << "BBBBBBBB" << endl;
	
	
	this->dtPseudo = this->nTotalCellVar++;
	this->name.push_back("pseudo-time-step");
	
	
	this->Ur = this->nTotalCellVar++;
	this->name.push_back("preconditioner");
	
	
	this->mu = this->nTotalCellVar++;
	this->name.push_back("viscosity");
	
	this->muT = this->nTotalCellVar++;
	this->name.push_back("turbulence-viscosity");
	
	this->muEffective = this->nTotalCellVar++;
	this->name.push_back("effective-viscosity");
	
	this->kappa = this->nTotalCellVar++;
	this->name.push_back("thermal-conductivity");
	
	this->kappaEffective = this->nTotalCellVar++;
	this->name.push_back("conductivity");
	
	this->D = this->nTotalCellVar++;
	this->name.push_back("mass-diffusivity");
	
	this->DEffective = this->nTotalCellVar++;
	this->name.push_back("effective-mass-diffusivity");
	
	this->kSGS = this->nTotalCellVar++;
	this->name.push_back("subgrid-tubulent-kinetic-energy");
	
	this->cv = this->nTotalCellVar++;
	this->name.push_back("constant-vloume-heat-capacity");
	
	this->cp = this->nTotalCellVar++;
	this->name.push_back("constant-pressure-heat-capacity");
	
	// this->RCM = this->nTotalCellVar++;
	// this->name.push_back("RCM");
	
	// this->invRCM = this->nTotalCellVar++;
	// this->name.push_back("invRCM");
	
	
	this->indicatorAMR = this->nTotalCellVar++;
	this->name.push_back("indicatorAMR");
	
	this->UDV.push_back(this->nTotalCellVar++);
	this->name.push_back("UDV0");
	this->UDV.push_back(this->nTotalCellVar++);
	this->name.push_back("UDV1");
	this->UDV.push_back(this->nTotalCellVar++);
	this->name.push_back("UDV2");
	
	


	// face LR var
	this->nTotalFaceLRVar = 0;
	this->fP = this->nTotalFaceLRVar++;
	this->fU = this->nTotalFaceLRVar++;
	this->fV = this->nTotalFaceLRVar++;
	this->fW = this->nTotalFaceLRVar++;
	this->fT = this->nTotalFaceLRVar++;
	for(int i=0; i<this->nSp-1; ++i){
		this->fVF.push_back(this->nTotalFaceLRVar++);
	}
	for(int i=0; i<this->nSp-1; ++i){
		this->fMF.push_back(this->nTotalFaceLRVar++);
	}
	this->fRho = this->nTotalFaceLRVar++;
	this->fC = this->nTotalFaceLRVar++;
	this->fHt = this->nTotalFaceLRVar++;

	this->fmu = this->nTotalFaceLRVar++;
	
	
	
	// face var
	this->nTotalFaceVar = 0;
	
	this->Un = this->nTotalFaceVar++;
	
	this->dUdx = this->nTotalFaceVar++;
	this->dUdy = this->nTotalFaceVar++;
	this->dUdz = this->nTotalFaceVar++;
	this->dVdx = this->nTotalFaceVar++;
	this->dVdy = this->nTotalFaceVar++;
	this->dVdz = this->nTotalFaceVar++;
	this->dWdx = this->nTotalFaceVar++;
	this->dWdy = this->nTotalFaceVar++;
	this->dWdz = this->nTotalFaceVar++;
	
	// this->wC = this->nTotalFaceVar++;
	
	// this->fDistCells.push_back(this->nTotalFaceVar++);
	// this->fDistCells.push_back(this->nTotalFaceVar++);
	// this->fDistCells.push_back(this->nTotalFaceVar++);
	
	// // bc
	// this->subinletU = 0.0;
	// this->subinletV = 0.0;
	// this->subinletW = -10.0;
	// this->subinletVF.push_back(1.0);
	
	// this->suboutletP = 101325.0;
	
	// this->supinletP = 101325.0;
	// this->supinletU = 0.0;
	// this->supinletV = 0.0;
	// this->supinletW = 0.0;
	// this->supinletVF.push_back(1.0);


	this->time = 0.0;
	
	this->PrT = 0.85;
	
	this->ScT = 0.35;
	
}