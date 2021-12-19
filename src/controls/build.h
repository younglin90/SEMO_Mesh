#pragma once
#include <iostream>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
using namespace std;

#include "../mesh/build.h"
// #include "../utility/read.h"
#include "../load/load.h"
#include "../species/species.h"


class SEMO_Controls_Builder{
	public:
		void readSpecies(vector<SEMO_Species>& species);
		void readConfigures();
		void setValues(vector<SEMO_Species>& species);
		
		string application;
		string startFrom;
		string stopAt;
		double saveInterval;
		string saveControl;
		string saveFormat;
		string timeFormat;
		string adjustTimeStep;
		int saveCompression;
		double maxCo;
		double maxVfCo;
		double maxTimeStep;
		double orgTimeStep;
		int writePrecision;
		
		// user save datas
		map<string,bool> saveMeshData;
		map<string,bool> saveGradientData;
		map<string,bool> saveThermodynamicData;
		map<string,bool> saveBodyForceData;
		
		// extract 관련
		double residual;
		vector<string> extractFieldDatas;
		vector<string> extractAverageDatas;
		vector<string> extractNames;
		vector<string> extractCellValueTargets;
		vector<vector<double>> extractCenterPoints;
		vector<double> extractRadii;
		vector<vector<double>> extractDatas;
		
		string fluxScheme;
		string gradScheme;
		string divScheme;
	
	
		
		vector<double> gravityAcceleration;
		
		int nTotalCellVar;
		int nEq;
		int nSp;
		int P;
		int U;
		int V;
		int W;
		int T;
		vector<int> VF;
		vector<int> MF;
		int Rho;
		int C;
		int Ht;
		vector<int> indicatorAMR;
		

		// vector<int> Q;
		vector<int> Qn;
		vector<int> Qm;
		
		int oldP;
		int oldU;
		int oldV;
		int oldW;
		int oldT;
		vector<int> oldVF;
		vector<int> oldMF;
		int oldRho;
		int oldHt;
		
		int dRhoDP;
		int dHtDP;
		int dRhoDT;
		int dHtDT;
		vector<int> dRhoDVF;
		vector<int> dHtDVF;
		vector<int> dRhoDMF;
		vector<int> dHtDMF;
		
		int dtPseudo;
		int Ur;
		
		vector<int> UDV;
		
		int Un;
		int oldUn;
		// int wC;
		
		
		// int RCM;
		// int invRCM;
		
		double Uco;
		double Lch;
		
		// // bc
		// double subinletU;
		// double subinletV;
		// double subinletW;
		// vector<double> subinletVF;
		
		// double suboutletP;
		
		// double supinletP;
		// double supinletU;
		// double supinletV;
		// double supinletW;
		// vector<double> supinletVF;
		
		
		vector<string> name;
		
		// face LR var
		int nTotalFaceLRVar;
		int fP;
		int fU;
		int fV;
		int fW;
		int fT;
		vector<int> fVF;
		vector<int> fMF;
		int fRho;
		int fC;
		int fHt;
		int fmu;
		
		int fdRhoDP;
		int fdRhoDT;
		int fdHtDP;
		int fdHtDT;
		vector<int> fdRhoDMF;
		vector<int> fdHtDMF;
		
		vector<int> fMF_HO;
		int fRho_HO;
		int fHt_HO;
		
		int fVF_NVD;
	
		
		
		// face var
		int nTotalFaceVar;
		// vector<int> fDistCells;
		
		int dUdx;
		int dUdy;
		int dUdz;
		int dVdx;
		int dVdy;
		int dVdz;
		int dWdx;
		int dWdy;
		int dWdz;
		
		// 최대 최소 값
		int maximumP;
		int minimumP;
		int maximumU;
		int minimumU;
		int maximumV;
		int minimumV;
		int maximumW;
		int minimumW;
		int maximumT;
		int minimumT;
		vector<int> maximumVF;
		vector<int> minimumVF;
		vector<int> maximumMF;
		vector<int> minimumMF;
		
		// turbulence
		string turbType;
		string turbLESModel;
		string turbRANSModel;
		
		
		// solvers
		string solverP;
		double toleranceP;
		double relTolP;
		string preconditionerP;
		int maxIterP;
		
		string solverFinalP;
		double toleranceFinalP;
		double relTolFinalP;
		string preconditionerFinalP;
		int maxIterFinalP;
		
		string solverU;
		double toleranceU;
		double relTolU;
		string preconditionerU;
		int maxIterU;
		
		string solverFinalU;
		double toleranceFinalU;
		double relTolFinalU;
		string preconditionerFinalU;
		int maxIterFinalU;
		
		vector<string> solverVF;
		vector<double> toleranceVF;
		vector<double> relTolVF;
		vector<string> preconditionerVF;
		vector<int> maxIterVF;
		
		vector<string> solverFinalVF;
		vector<double> toleranceFinalVF;
		vector<double> relTolFinalVF;
		vector<string> preconditionerFinalVF;
		vector<int> maxIterFinalVF;
		
		
		// URF
		double momVelURF;
		string momVelAdjustRF;
		vector<int> momVelAdjustSteps;
		vector<double> momVelAdjustValues;
		
		double prePreURF;
		string prePreAdjustRF;
		vector<int> prePreAdjustSteps;
		vector<double> prePreAdjustValues;
		
		double preVelURF;
		string preVelAdjustRF;
		vector<int> preVelAdjustSteps;
		vector<double> preVelAdjustValues;
		
		double vofVofURF;
		string vofVofAdjustRF;
		vector<int> vofVofAdjustSteps;
		vector<double> vofVofAdjustValues;
		
		// URF
		double dualTimeURF_P;
		string dualTimeAdjustRF_P;
		vector<int> dualTimeAdjustSteps_P;
		vector<double> dualTimeAdjustValues_P;
		
		double dualTimeURF_U;
		string dualTimeAdjustRF_U;
		vector<int> dualTimeAdjustSteps_U;
		vector<double> dualTimeAdjustValues_U;
		
		double dualTimeURF_T;
		string dualTimeAdjustRF_T;
		vector<int> dualTimeAdjustSteps_T;
		vector<double> dualTimeAdjustValues_T;
		
		double dualTimeURF_MF;
		string dualTimeAdjustRF_MF;
		vector<int> dualTimeAdjustSteps_MF;
		vector<double> dualTimeAdjustValues_MF;
		
		// limiter
		double maxP;
		double minP;
		double maxU;
		double minU;
		double maxV;
		double minV;
		double maxW;
		double minW;
		double maxT;
		double minT;
		
		
		// scheme
		double time;
		double timeStep;
		double oldTimeStep;
		double old2TimeStep;
		double pseudoCo;
		double specifiedCFL;
		double allowableCFL;
		
		// iterator
		int iterRealMax;
		int iterReal;
		
		int iterPBs;
		int iterPBsMax;
		int iterMom;
		int iterMomMax;
		int iterPre;
		int iterPreMax;
		int iterVof;
		int iterVofMax;
	
		
		int iterPseudo;
		int iterPseudoMax;
		
		int iterTotal;
		double startClock;
		
		// transport & turbulence
		int mu;
		int muT;
		int muEffective;
		int k;
		int kEffective;
		int D;
		int DEffective;
		int cv;
		int cp;
		int kSGS;
		
		// source terms
		vector<int> sourceGravity;
		vector<int> sourceSurfaceTension;
		

		double PrT;
		double ScT;
		
		int kappa;
		
		// dynamic Mesh
		vector<double> indicatorCriterion;
		int intervalRefine;
		// double indicatorRefine;
		int maxLevelRefine;
		int maxCellsRefine;
		double minVolumeRefine;
		
		int bufferLayer;
		
		int intervalUnrefine;
		// double indicatorUnrefine;
		
};

