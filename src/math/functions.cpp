
class SEMO_Utility_Math{
	public:
		double *crossProduct(double v_A[], double v_B[]) {
			static double c_P[3];
			c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
			c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
			c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
			return c_P;
		}
};


