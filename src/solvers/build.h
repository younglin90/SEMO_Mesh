#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

#include "../species/species.h"
#include "../mesh/build.h"
#include "../controls/build.h"

class SEMO_Solvers_Builder{
	public:
		void calcUnderRelaxationFactors(
			SEMO_Controls_Builder& controls);
			
		void calcUnderRelaxationFactorsDualTime(
			SEMO_Controls_Builder& controls);
	
		void setInitValues(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls);
			
		// solvers
		void incompressiblePressureBased(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
		
		void compressibleDensityBasedSingleTime(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void compressibleDensityBasedDualTime(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void hybridBased(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		// reconstructions
		void setIncomValuesLeftRightFace(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFaceForSegregated(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFaceWithReconPVForSegregated(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFace(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setIncomValuesLeftRightFaceWithReconPV(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFaceWithRecon(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFaceForSegregatedWithVfMSTACS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setCompValuesLeftRightFaceWithVfMSTACS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void setIncomValuesLeftRightFaceWithVelQUICK(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
			
		void setIncomValuesLeftRightFaceWithVfMSTACS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
			
		void reconZeroOrder(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void reconIncomZeroOrder(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void calcOtherDataFromEOS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void calcOtherDataFromEOSMF(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void calcOtherDataIncomFromEOS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		// pressure based
		void calcMomentumEqs(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species,
			vector<double>& linAD,
			vector<double>& linAOD,
			vector<double>& linAFL,
			vector<double>& linAFR,
			vector<double>& residuals);
			
		void calcPressureEq(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<double>& linAD,
			vector<double>& linAOD,
			vector<double>& linAFL,
			vector<double>& linAFR,
			vector<double>& residuals);
			
		void calcVolfracEq(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<double>& linAD,
			vector<double>& linAOD,
			vector<double>& linAFL,
			vector<double>& linAFR,
			vector<double>& residuals);
			
		void solvePETSc(
			SEMO_Mesh_Builder& mesh,
			vector<double>& resiVar, 
			vector<double>& linA, vector<double>& linAL, vector<double>& linAR, 
			vector<double>& linB, 
			string solver, double tolerance, double relTol, string preconditioner,
			int maxIter);
			
		void solveHYPRE(
			SEMO_Mesh_Builder& mesh,
			vector<double>& resiVar, 
			vector<double>& linA, vector<double>& linAL, vector<double>& linAR, 
			vector<double>& linB, 
			string solver, double tolerance, double relTol, string preconditioner,
			int maxIter);

		void calcPseudoTimeStep(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls);

		void calcRealTimeStep(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls);

		void calcIncomRealTimeStep(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls);

		void calcCorantNumberForPrint(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			double& corantNum);
			
		void calcRHS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species,
			vector<vector<double>>& residuals);
			
		void calcSingleRHS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species,
			vector<vector<double>>& residuals);
			
		void calcFluxes(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species,
			vector<vector<double>>& residuals);
			
		void calcFluxHAUS(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& RhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& RhoR, double& CR, double& HtR,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		void calcFluxRoe(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& RhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& RhoR, double& CR, double& HtR,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		void calcFluxAPRoe(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& RhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& RhoR, double& CR, double& HtR,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		void calcFluxRoeAMTP(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& RhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& RhoR, double& CR, double& HtR,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		void calcFluxUpwind(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& RhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& RhoR, double& CR, double& HtR,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		void calcFluxAUSMPWP_N(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& rhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& rhoR, double& CR, double& HtR,
			vector<double>& nvec,
			double& w1 , double& w2 , double& w3 , double& fC, double Uco, double dtstep, double Lch,
			vector<double>& flux
			);
			
		void calcFluxRoeM_N(
			double& PL, double& UL, double& VL, double& WL, double& TL, 
			vector<double>& YL, double& rhoL, double& CL, double& HtL, 
			double& PR, double& UR, double& VR, double& WR, double& TR, 
			vector<double>& YR, double& rhoR, double& CR, double& HtR,
			vector<double>& nvec,
			double& w1 , double& w2 , double& w3 , double& fC, double Uco, double dtstep, double Lch,
			vector<double>& flux
			);
			
		void calcViscousFlux(
			double& muL, double& UL, double& VL, double& WL, double& RhoL, 
			double& muR, double& UR, double& VR, double& WR, double& RhoR, 
			double& dUdxF, double& dUdyF, double& dUdzF,
			double& dVdxF, double& dVdyF, double& dVdzF,
			double& dWdxF, double& dWdyF, double& dWdzF,
			vector<double>& nvec,
			vector<double>& flux
			);
			
		double M_func(double M, double op, double alp);
		double pre_func(double M, double op, double alp);
		double f_func(double M);
		
		void sourceTerms(
			SEMO_Mesh_Builder& mesh, 
			SEMO_Controls_Builder& controls, 
			vector<SEMO_Species>& species,
			vector<vector<double>>& residuals);
			
		void calcLinearSolver(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals);
			
		void calcSingleLinearSolver(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals);
			
		void constructMatrix(
			SEMO_Cell& cell, SEMO_Controls_Builder& controls,
			vector<vector<double>>& matrixP,
			vector<vector<double>>& matrixT);
			
		void updateValues(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals);
			
		void updateValuesSingle(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals);
			
			
		void calcIncomCellEOSVF(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
			
		void calcCellTransport(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void calcCellEOSMF(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
		void calcCellEOSVF(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<SEMO_Species>& species);
			
			
		// void getValuesFromEOSMF(
			// vector<SEMO_Species>& species,
			// double& P, double& U, double& V, double& W, double& T, vector<double>& MF,
			// double& rho, double& C, double& Ht);
			
		void getValuesFromEOSVF(
			vector<SEMO_Species>& species,
			double& P, double& U, double& V, double& W, double& T, vector<double>& VF,
			double& rho, double& C, double& Ht,
			vector<double>& MF);
			
		void getValuesFromEOSMF(
			vector<SEMO_Species>& species,
			double& P, double& U, double& V, double& W, double& T, vector<double>& MF,
			double& rho, double& C, double& Ht,
			vector<double>& VF);
			
		void eosIdeal(
			SEMO_Species& species,
			double& P, double& U, double& V, double& W, double& T,
			double& rho, double& C, double& Ht,
			double& drhodP, double& drhodT, double& dhdP, double& dhdT);
			
		void eosNASG(
			SEMO_Species& species,
			double& P, double& U, double& V, double& W, double& T,
			double& rho, double& C, double& Ht,
			double& drhodP, double& drhodT, double& dhdP, double& dhdT);
			
			
			
		void calcNormResiduals(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals,
			vector<double>& norm);
			
		void calcNormResiduals(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<double>& residuals,
			vector<double>& norm);
			
		// NVD
		void calcVanLeer(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcQUICK(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcMINMOD(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcBoundedCD(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcOSHER(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcSMART(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcModifiedSMART(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcSTOIC(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcModifiedSTOIC(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcMUSCL(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcSUPERBEE(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcModifiedSUPERBEE(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
		void calcMSTACS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			int cn, int fn,
			vector<vector<double>>& inpDX);
			
			
			
		// curvature
		void calcCurvature(
			SEMO_Mesh_Builder& mesh,
			int cn,
			vector<double>& kappa);
			
			
			
			
			
			
		void calcLinearSolverLUSGS(
			SEMO_Mesh_Builder& mesh,
			SEMO_Controls_Builder& controls,
			vector<vector<double>>& residuals);
	
		void getNumericalFluxJacobianUpwind(
			SEMO_Cell& cellL, SEMO_Cell& cellR, SEMO_Controls_Builder& controls,
			vector<double>& nvec,
			string LorR,
			vector<vector<double>>& fluxJac
			);
			
			
			
};

