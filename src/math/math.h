#pragma once
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

#include "../mesh/build.h"

class SEMO_Utility_Math{
	public:
		double *crossProduct(double v_A[], double v_B[]);
		
		
		int CompareDoubleAbsoulteAndUlps(double x,
										 double y,
										 double absTolerance,
										 int ulpsTolerance);
		bool approximatelyEqualAbsRel(double a, 
			double b, double absEpsilon, double relEpsilon);
		
		template<typename TReal>
		bool isApproximatelyEqual(TReal a, TReal b, TReal tolerance);
		
		void GaussSeidelSOR(vector<vector<double>>& A);
		
		double determinant(vector<vector<double>>& matrix, int n);
		
		void getCofactor(vector<vector<double>>& mat, vector<vector<double>>& temp, int p,
                                    int q, int n);
		bool isSingular(vector<vector<double>>& mat, int n);
		
		void initLeastSquare(
			SEMO_Mesh_Builder& mesh);
			
		void calcLeastSquare(
			SEMO_Mesh_Builder& mesh,
			string sStencil, string sOrder, string sCellorInp,
			int cn, int fn,
			vector<double>& phi,
			vector<vector<double>>& gradient
			);
			

		void printMatrix(vector<vector<double>>& a);
		void printVector(vector<double>& v);
		double pythag(double a, double b);
		void svdcmp(vector<vector<double>>& a, vector<double>& w, vector<vector<double>>& v);
		void matmul(vector<vector<double>>& a, vector<vector<double>>& b, vector<vector<double>>& out);
		void transpose(vector<vector<double>>& a, vector<vector<double>>& out);
			
		// void calcLeastSquare(
			// SEMO_Mesh_Builder& mesh,
			// int cn, int fn,
			// vector<vector<double>>& gradient);
			
		// void calcLeastSquare(
			// SEMO_Mesh_Builder& mesh,
			// int cn,
			// vector<double> phi,
			// vector<vector<double>>& gradient);
			
		// void initLeastSquare2nd(
			// SEMO_Mesh_Builder& mesh);
			
		// void calcLeastSquare2nd(
			// SEMO_Mesh_Builder& mesh,
			// int cn, int fn,
			// vector<vector<double>>& gradient);
			
		// void calcLeastSquare2nd(
			// SEMO_Mesh_Builder& mesh,
			// int cn,
			// vector<double> phi,
			// vector<vector<double>>& gradient
			// );
			
		void calcGaussGreen(
			SEMO_Mesh_Builder& mesh,
			int cn, int fn,
			vector<vector<double>>& gradient);
			
		void calcGaussGreen(
			SEMO_Mesh_Builder& mesh,
			int cn,
			vector<double> phi,
			vector<vector<double>>& gradient
			);
		void calcGaussGreen(
			SEMO_Mesh_Builder& mesh,
			int cn,
			vector<double> phi0, vector<double> phi1, vector<double> phi2,
			vector<vector<double>>& gradient
			);
			
		// void calcGGLSQ(
			// SEMO_Mesh_Builder& mesh,
			// int cn, int fn,
			// vector<vector<double>>& gradient);
			
		// void calcMGG(
			// SEMO_Mesh_Builder& mesh,
			// int cn, int fn,
			// int iterMax, double toler, 
			// vector<vector<double>>& gradient);
			
			
		// void gradientGGBoundaryTreatment(
			// SEMO_Mesh_Builder& mesh,
			// SEMO_Controls_Builder& controls,
			// vector<vector<double>>& gradient);
			
		// void gradientLSQBoundaryTreatment(
			// SEMO_Mesh_Builder& mesh,
			// SEMO_Controls_Builder& controls,
			// vector<vector<double>>& gradient);
			
		void calcLimiterGradient(
			SEMO_Mesh_Builder& mesh,
			int cn, int fn,
			vector<vector<double>>& gradient,
			vector<double>& limGrad);
			
		void calcGradientFace(
			SEMO_Mesh_Builder& mesh,
			vector<vector<double>>& gradient,
			int cn, int fn,
			int inX, int inY, int inZ);
		
		void lib_RCM(SEMO_Mesh_Builder& mesh);
	
};

