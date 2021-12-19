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


double SEMO_Utility_Math::determinant(vector<vector<double>>& matrix, int n) { 

	// int n = matrix.size();

	double det = 0.0;
	vector<vector<double>> submatrix(n,vector<double>(n,0.0));
	if (n == 2) return ((matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]));
	else {
		for (int x = 0; x < n; x++) {
			int subi = 0;
			for (int i = 1; i < n; i++) {
				int subj = 0;
				for (int j = 0; j < n; j++) {
					if (j == x)
					continue;
					submatrix[subi][subj] = matrix[i][j];
					subj++;
				}
				subi++;
			}
			det = det + (pow(-1, x) * matrix[0][x] * determinant( submatrix, n - 1 ));
		}
	}
	return det;
	
  // int det=0, p, h, k, i, j;
  // vector<vector<double>> temp(9,vector<double>(9,0.0));
  // if(n==1) {
    // return matrix[0][0];
  // } else if(n==2) {
    // det=(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]);
    // return det;
  // } else {
    // for(p=0;p<n;p++) {
      // h = 0;
      // k = 0;
      // for(i=1;i<n;i++) {
        // for( j=0;j<n;j++) {
          // if(j==p) {
            // continue;
          // }
          // temp[h][k] = matrix[i][j];
          // k++;
          // if(k==n-1) {
            // h++;
            // k = 0;
          // }
        // }
      // }
      // det=det+matrix[0][p]*pow(-1,p)*determinant(temp,n-1);
    // }
    // return det;
  // }
	
	

}




// Function to get cofactor of mat[p][q] in temp[][].
// n is current dimension of mat[][]
void SEMO_Utility_Math::getCofactor(vector<vector<double>>& mat, vector<vector<double>>& temp, int p,
                                    int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
 
            // Copying into temporary matrix only
            // those element which are not in given
            // row and column
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
 
                // Row is filled, so increase row
                // index and reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
 
/* Recursive function to check if mat[][] is
   singular or not. */
bool SEMO_Utility_Math::isSingular(vector<vector<double>>& mat, int n)
{
    int D = 0; // Initialize result
 
    // Base case : if matrix contains single element
    if (n == 1)
        return mat[0][0];
 
    vector<vector<double>> temp(9,vector<double>(9,0.0)); // To store cofactors
 
    double sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++) {
 
        // Getting Cofactor of mat[0][f]
        getCofactor(mat, temp, 0, f, n);
        D += sign * mat[0][f] * isSingular(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}







// prints an arbitrary size matrix to the standard output
void SEMO_Utility_Math::printMatrix(vector<vector<double>>& a) {
	int i, j;

	for (i = 0; i < a.size(); i++) {
		for (j = 0; j < a[0].size(); j++) {
			printf("%.4lf ", a[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

// prints an arbitrary size vector to the standard output
void SEMO_Utility_Math::printVector(vector<double>& v) {
	int i;

	for (i = 0; i < v.size(); i++) {
		printf("%.4lf ", v[i]);
	}
	printf("\n\n");
}

// calculates sqrt( a^2 + b^2 ) with decent precision
double SEMO_Utility_Math::pythag(double a, double b) {
	double absa, absb;

	absa = abs(a);
	absb = abs(b);

	if (absa > absb)
		return (absa * sqrt(1.0 + (absb / absa)*(absb / absa)));
	else
		return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb)*(absa / absb)));
}

/*
 Modified from Numerical Recipes in C
 Given a matrix a[nRows][nCols], svdcmp() computes its singular value
 decomposition, A = U * W * Vt.  A is replaced by U when svdcmp
 returns.  The diagonal matrix W is output as a vector w[nCols].
 V (not V transpose) is output as the matrix V[nCols][nCols].
 */
// ref : numerical recip
void SEMO_Utility_Math::svdcmp(vector<vector<double>>& a, vector<double>& w, vector<vector<double>>& v) {
	int flag, i, its, j, jj, k, l, nm;
	double anorm, c, f, g, h, s, scale, x, y, z;

	int nRows = a.size();
	int nCols = a[0].size();

	vector<double> rv1(nCols, 0.0);

	g = scale = anorm = 0.0;
	for (i = 0; i < nCols; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i < nRows) {
			for (k = i; k < nRows; k++)
				scale += abs(a[k][i]);
			if (scale) {
				for (k = i; k < nRows; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -(f > 0.0 ? abs(sqrt(s)) : -abs(sqrt(s)));
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j < nCols; j++) {
					for (s = 0.0, k = i; k < nRows; k++)
						s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k < nRows; k++)
						a[k][j] += f * a[k][i];
				}
				for (k = i; k < nRows; k++)
					a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i < nRows && i != nCols - 1) {
			for (k = l; k < nCols; k++)
				scale += abs(a[i][k]);
			if (scale) {
				for (k = l; k < nCols; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -(f > 0.0 ? abs(sqrt(s)) : -abs(sqrt(s)));
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k < nCols; k++)
					rv1[k] = a[i][k] / h;
				for (j = l; j < nRows; j++) {
					for (s = 0.0, k = l; k < nCols; k++)
						s += a[j][k] * a[i][k];
					for (k = l; k < nCols; k++)
						a[j][k] += s * rv1[k];
				}
				for (k = l; k < nCols; k++)
					a[i][k] *= scale;
			}
		}
		anorm = max(anorm, (abs(w[i]) + abs(rv1[i])));

		// printf(".");
		// fflush(stdout);
	}

	for (i = nCols - 1; i >= 0; i--) {
		if (i < nCols - 1) {
			if (g) {
				for (j = l; j < nCols; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j < nCols; j++) {
					for (s = 0.0, k = l; k < nCols; k++)
						s += a[i][k] * v[k][j];
					for (k = l; k < nCols; k++)
						v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j < nCols; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
		// printf(".");
		// fflush(stdout);
	}

	for (i = min(nRows, nCols) - 1; i >= 0; i--) {
		l = i + 1;
		g = w[i];
		for (j = l; j < nCols; j++)
			a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j < nCols; j++) {
				for (s = 0.0, k = l; k < nRows; k++)
					s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k = i; k < nRows; k++)
					a[k][j] += f * a[k][i];
			}
			for (j = i; j < nRows; j++)
				a[j][i] *= g;
		}
		else
			for (j = i; j < nRows; j++)
				a[j][i] = 0.0;
		++a[i][i];
		// printf(".");
		// fflush(stdout);
	}

	for (k = nCols - 1; k >= 0; k--) {
		for (its = 0; its < 30; its++) {
			flag = 1;
			for (l = k; l >= 0; l--) {
				nm = l - 1;
				if ((abs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((abs(w[nm]) + anorm) == anorm)
					break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((abs(f) + anorm) == anorm)
						break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 0; j < nRows; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 0; j < nCols; j++)
						v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29)
				printf("no convergence in 30 svdcmp iterations\n");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + (f > 0.0 ? abs(g) : -abs(g)))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 0; jj < nCols; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 0; jj < nRows; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
		// printf(".");
		// fflush(stdout);
	}
	// printf("\n");
}


void SEMO_Utility_Math::matmul(vector<vector<double>>& a, vector<vector<double>>& b, vector<vector<double>>& out) {
	out.clear();
	out.resize(a.size(), vector<double>(b[0].size(), 0.0));

	for (int i = 0; i < a.size(); ++i) {
		for (int j = 0; j < b[0].size(); ++j) {
			double tmp = 0.0;
			for (int m = 0; m < b.size(); ++m) {
				tmp += a[i][m] * b[m][j];
			}
			out[i][j] = tmp;
		}
	}

}

void SEMO_Utility_Math::transpose(vector<vector<double>>& a, vector<vector<double>>& out) {
	out.clear();
	out.resize(a[0].size(), vector<double>(a.size(), 0.0));

	for (int i = 0; i < a[0].size(); ++i) {
		for (int j = 0; j < a.size(); ++j) {
			out[i][j] = a[j][i];
		}
	}
}
