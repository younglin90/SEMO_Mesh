#include <algorithm>
#include <cmath>
#include <limits>
#include "math.h"

using namespace std;

double* SEMO_Utility_Math::crossProduct(double v_A[], double v_B[]) {
	static double c_P[3];
	c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
	c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
	return c_P;
}


int SEMO_Utility_Math::CompareDoubleAbsoulteAndUlps(double x,
								 double y,
								 double absTolerance = (1.0e-8),
								 int ulpsTolerance = 4)
{
	double diff = x - y;
	if (fabs(diff) <= absTolerance)
		return 0;

	__int64 nx = *((__int64*)&x);
	__int64 ny = *((__int64*)&y);

	if ((nx & 0x8000000000000000) != (ny & 0x8000000000000000))
		return (diff > 0) ? 1 : -1;

	__int64 ulpsDiff = nx - ny;
	if ((ulpsDiff >= 0 ? ulpsDiff : -ulpsDiff) <= ulpsTolerance)
		return 0;

	return (diff > 0) ? 1 : -1;
}


bool SEMO_Utility_Math::approximatelyEqualAbsRel(double a, double b, double absEpsilon, double relEpsilon) { 
	// Check if the numbers are really close -- needed when comparing numbers near zero. 
	
	// absEpsilon = 1.e-5;
	// relEpsilon = 1.e-1;
	
	double diff = fabs(a - b); 
	if (diff <= absEpsilon) 
		return true; 
	// Otherwise fall back to Knuth's algorithm 
	return diff <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * relEpsilon); 
}



// bool SEMO_Utility_Math::AreSame(double a, double b) {
    // return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
// }


// bool SEMO_Utility_Math::AreSame(double a, double b)
// {
	// double VERYSMALL = 1.E-8;
    // double absDiff = fabs(a - b);
    // if (absDiff < VERYSMALL)
    // {
        // return true;
    // }

    // double maxAbs  = max(fabs(a) - fabs(b));
    // return (absDiff/maxAbs) < EPSILON;
// }


//use this most of the time, tolerance needs to be meaningful in your context
template<typename TReal>
bool SEMO_Utility_Math::isApproximatelyEqual(
TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon()){
    TReal diff = std::fabs(a - b);
    if (diff <= tolerance)
        return true;

    if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance)
        return true;

    return false;
}

// //supply tolerance that is meaningful in your context
// //for example, default tolerance may not work if you are comparing double with float
// template<typename TReal>
// static bool SEMO_Utility_Math::isApproximatelyZero(TReal a, TReal tolerance = std::numeric_limits<TReal>::epsilon())
// {
    // if (std::fabs(a) <= tolerance)
        // return true;
    // return false;
// }

// //use this when you want to be on safe side
// //for example, don&#39;t start rover unless signal is above 1
// template<typename TReal>
// static bool SEMO_Utility_Math::isDefinitelyLessThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
// {
    // TReal diff = a - b;
    // if (diff < tolerance)
        // return true;

    // if (diff < std::fmax(std::fabs(a), std::fabs(b)) * tolerance)
        // return true;

    // return false;
// }
// template<typename TReal>
// static bool SEMO_Utility_Math::isDefinitelyGreaterThan(TReal a, TReal b, TReal tolerance = std::numeric_limits<TReal>::epsilon())
// {
    // TReal diff = a - b;
    // if (diff > tolerance)
        // return true;

    // if (diff > std::fmax(std::fabs(a), std::fabs(b)) * tolerance)
        // return true;

    // return false;
// }

// //implements ULP method
// //use this when you are only concerned about floating point precision issue
// //for example, if you want to see if a is 1.0 by checking if its within
// //10 closest representable floating point numbers around 1.0.
// template<typename TReal>
// static bool SEMO_Utility_Math::isWithinPrecisionInterval(TReal a, TReal b, unsigned int interval_size = 1)
// {
    // TReal min_a = a - (a - std::nextafter(a, std::numeric_limits<TReal>::lowest())) * interval_size;
    // TReal max_a = a + (std::nextafter(a, std::numeric_limits<TReal>::max()) - a) * interval_size;

    // return min_a <= b && max_a >= b;
// }




void SEMO_Utility_Math::GaussSeidelSOR(vector<vector<double>>& A) { 


	double N = A.size();
	
	vector<vector<double>> L(N,vector<double>(N,0.0));
	vector<vector<double>> U(N,vector<double>(N,0.0));
	vector<double> B(N,0.0);
	vector<double> D(N,0.0);
	vector<double> X(N,0.0);
	
	for(int k=1; k<=N-1; ++k){
		for(int i=k+1; i<=N; ++i){
			double coeff = A[i-1][k-1]/A[k-1][k-1];
			L[i-1][k-1] = coeff;
			for(int j=k+1; j<=N; ++j){
				A[i-1][j-1] = A[i-1][j-1] - coeff*A[k-1][j-1];
			}
		}
	}
	
	for(int i=1; i<=N; ++i){
		L[i-1][i-1] = 1.0;
	}
	
	for(int j=1; j<=N; ++j){
		for(int i=1; i<=j; ++i){
			U[i-1][j-1] = A[i-1][j-1];
		}
	}
	
	for(int k=1; k<=N; ++k){
		B[k-1] = 1.0;
		D[0] = B[0];
		for(int i=2; i<=N; ++i){
			D[i-1] = B[i-1];
			for(int j=1; j<=i-1; ++j){
				D[i-1] = D[i-1] - L[i-1][j-1]*D[j-1];
			}
		}
		
		X[N-1] = D[N-1] / U[N-1][N-1];
		
		// cout << U[N-1][N-1] << endl;
		
		for(int i=N-1; i>=1; --i){
			X[i-1] = D[i-1];
			for(int j=N; j>=i+1; --j){
				X[i-1] = X[i-1] - U[i-1][j-1]*X[j-1];
			}
			X[i-1] = X[i-1]/U[i-1][i-1];
		}
		
		for(int i=1; i<=N; ++i){
			A[i-1][k-1] = X[i-1];
		}
		
		B[k-1] = 0.0;
	}


}
