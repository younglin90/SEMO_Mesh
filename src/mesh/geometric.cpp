#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
// #include <numeric>
using namespace std;

#include "geometric.h"


/*
new, 2021-10-15
http://paulbourke.net/geometry/polygonmesh/
*/
void SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
	int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	vector<double>& unitNormals, double& area,
	double& x, double& y, double& z,
	double& VSn, vector<double>& cellCentroid ){
		
	double x_tmp = 0.0;  double y_tmp = 0.0; double z_tmp = 0.0;
	double Vx0_old = 0.0;  double Vy0_old = 0.0; double Vz0_old = 0.0;
	vector<double> unitNormals_tmp(3,0.0);
	vector<double> cellCentroid_tmp(3,0.0);
	double area_tmp = 0.0;
	double VSn_tmp = 0.0;
	
	double x_avg, y_avg, z_avg;
	x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)Vx.size();
	y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)Vy.size();
	z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)Vz.size();
	
	for(int two=0; two<5; ++two){
		
		unitNormals_tmp[0] = 0.0;
		unitNormals_tmp[1] = 0.0;
		unitNormals_tmp[2] = 0.0;
		
		cellCentroid_tmp[0] = 0.0;
		cellCentroid_tmp[1] = 0.0;
		cellCentroid_tmp[2] = 0.0;
		
		area_tmp = 0.0;
		VSn_tmp = 0.0;
		double VSn_x = 0.0; double VSn_y = 0.0; double VSn_z = 0.0;
		double Vx0 = 0.0; double Vy0 = 0.0; double Vz0 = 0.0;
		
		if(two==0){
			// Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			
			Vx0 = x_avg; Vy0 = y_avg; Vz0 = z_avg;
		}
		else{
			Vx0 = Vx0_old + 1.0*(x_tmp-Vx0_old); 
			Vy0 = Vy0_old + 1.0*(y_tmp-Vy0_old); 
			Vz0 = Vz0_old + 1.0*(z_tmp-Vz0_old); 
		}
		Vx0_old = Vx0; Vy0_old = Vy0; Vz0_old = Vz0;
			// Vx0 = 0; Vy0 = 0; Vz0 = 0;
			// Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			// Vx0 = Vx.back(); Vy0 = Vy.back(); Vz0 = Vz.back();
			
		x_tmp = 0.0;  y_tmp = 0.0; z_tmp = 0.0;
		
		double totArea = 0.0;
		
		for(int i = 0; i < n; ++i ) {
			int b=i;
			int c=i+1;
			if(i==n-1) c=0;
			
			double vect_A[3];
			double vect_B[3];
			
			vect_A[0] = Vx[b]-Vx0; vect_A[1] = Vy[b]-Vy0; vect_A[2] = Vz[b]-Vz0;
			vect_B[0] = Vx[c]-Vx0; vect_B[1] = Vy[c]-Vy0; vect_B[2] = Vz[c]-Vz0;
			// vect_A[0] = Vx[b]; vect_A[1] = Vy[b]; vect_A[2] = Vz[b];
			// vect_B[0] = Vx[c]; vect_B[1] = Vy[c]; vect_B[2] = Vz[c];
		  
			double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
			double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
			double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

			unitNormals_tmp[0] += Nx_tmp;
			unitNormals_tmp[1] += Ny_tmp;
			unitNormals_tmp[2] += Nz_tmp;
			
			double tmp_area = sqrtl(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp);
			// double tmp_area = pow(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp,0.5);
			double rx = (Vx0+Vx[b]+Vx[c]) / 3.0;
			double ry = (Vy0+Vy[b]+Vy[c]) / 3.0;
			double rz = (Vz0+Vz[b]+Vz[c]) / 3.0;
			x_tmp += rx*tmp_area; y_tmp += ry*tmp_area; z_tmp += rz*tmp_area;
			
			totArea += tmp_area;
			
			VSn_tmp += Vx0*Nx_tmp; VSn_tmp += Vy0*Ny_tmp; VSn_tmp += Vz0*Nz_tmp;
			
			cellCentroid_tmp[0] += 2.0*Nx_tmp*(
				(Vx0+Vx[b])*(Vx0+Vx[b]) + (Vx[b]+Vx[c])*(Vx[b]+Vx[c]) + (Vx[c]+Vx0)*(Vx[c]+Vx0));
			cellCentroid_tmp[1] += 2.0*Ny_tmp*(
				(Vy0+Vy[b])*(Vy0+Vy[b]) + (Vy[b]+Vy[c])*(Vy[b]+Vy[c]) + (Vy[c]+Vy0)*(Vy[c]+Vy0));
			cellCentroid_tmp[2] += 2.0*Nz_tmp*(
				(Vz0+Vz[b])*(Vz0+Vz[b]) + (Vz[b]+Vz[c])*(Vz[b]+Vz[c]) + (Vz[c]+Vz0)*(Vz[c]+Vz0));
			
			
			// VSn_tmp += Vx0*2.0*Nx_tmp / 3.0; 
			// VSn_tmp += Vy0*2.0*Ny_tmp / 3.0; 
			// VSn_tmp += Vz0*2.0*Nz_tmp / 3.0;
			
			// cellCentroid_tmp[0] += Vx0*2.0*Nx_tmp * 0.25*(Vx0+Vx[b]+Vx[c]) / 3.0;
			// cellCentroid_tmp[1] += Vy0*2.0*Ny_tmp * 0.25*(Vy0+Vy[b]+Vy[c]) / 3.0;
			// cellCentroid_tmp[2] += Vz0*2.0*Nz_tmp * 0.25*(Vz0+Vz[b]+Vz[c]) / 3.0;
			
		}
		
		
		

		// Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			
		// x_tmp = 0.0;  y_tmp = 0.0; z_tmp = 0.0;
		
		// for(int i = 2; i < n; ++i ) {
			// int b=i-1;
			// int c=i;
			
			// double vect_A[3];
			// double vect_B[3];
			
			// vect_A[0] = Vx[b]-Vx0; vect_A[1] = Vy[b]-Vy0; vect_A[2] = Vz[b]-Vz0;
			// vect_B[0] = Vx[c]-Vx0; vect_B[1] = Vy[c]-Vy0; vect_B[2] = Vz[c]-Vz0;
		  
			// double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
			// double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
			// double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

			// unitNormals_tmp[0] += Nx_tmp;
			// unitNormals_tmp[1] += Ny_tmp;
			// unitNormals_tmp[2] += Nz_tmp;
			
			// double tmp_area = sqrtl(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp);
			// // double tmp_area = pow(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp,0.5);
			// double rx = (Vx0+Vx[b]+Vx[c]) / 3.0;
			// double ry = (Vy0+Vy[b]+Vy[c]) / 3.0;
			// double rz = (Vz0+Vz[b]+Vz[c]) / 3.0;
			// x_tmp += rx*tmp_area; y_tmp += ry*tmp_area; z_tmp += rz*tmp_area;
			
			// VSn_tmp += Vx0*Nx_tmp; VSn_tmp += Vy0*Ny_tmp; VSn_tmp += Vz0*Nz_tmp;
			
			// cellCentroid_tmp[0] += 2.0*Nx_tmp*(
				// (Vx0+Vx[b])*(Vx0+Vx[b]) + (Vx[b]+Vx[c])*(Vx[b]+Vx[c]) + (Vx[c]+Vx0)*(Vx[c]+Vx0));
			// cellCentroid_tmp[1] += 2.0*Ny_tmp*(
				// (Vy0+Vy[b])*(Vy0+Vy[b]) + (Vy[b]+Vy[c])*(Vy[b]+Vy[c]) + (Vy[c]+Vy0)*(Vy[c]+Vy0));
			// cellCentroid_tmp[2] += 2.0*Nz_tmp*(
				// (Vz0+Vz[b])*(Vz0+Vz[b]) + (Vz[b]+Vz[c])*(Vz[b]+Vz[c]) + (Vz[c]+Vz0)*(Vz[c]+Vz0));
			
		// }
		
		
		
		// VSn = VSn_x + VSn_y + VSn_z;
		
		double mag_unitNormals = sqrtl(
			unitNormals_tmp[0]*unitNormals_tmp[0]+
			unitNormals_tmp[1]*unitNormals_tmp[1]+
			unitNormals_tmp[2]*unitNormals_tmp[2]);
		unitNormals_tmp[0] /= mag_unitNormals;
		unitNormals_tmp[1] /= mag_unitNormals;
		unitNormals_tmp[2] /= mag_unitNormals;
		
		area_tmp = mag_unitNormals;
		
		x_tmp /= totArea;
		y_tmp /= totArea;
		z_tmp /= totArea;
	
	}
	
	unitNormals.clear();
	unitNormals.resize(3,0.0);
	cellCentroid.clear();
	cellCentroid.resize(3,0.0);
	
	for(int i=0; i<3; ++i){
		unitNormals[i] = unitNormals_tmp[i];
		cellCentroid[i] = cellCentroid_tmp[i];
	}
	
	area = area_tmp;
	
	VSn = VSn_tmp;
		
	x = x_tmp; y = y_tmp; z = z_tmp;
	
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
}





// void SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
	// int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	// vector<double>& unitNormals, double& area,
	// double& x, double& y, double& z,
	// double& VSn, vector<double>& cellCentroid ){
		
	// for(int two=0; two<1; ++two){
		
		// unitNormals.clear();
		// unitNormals.resize(3,0.0);
		// cellCentroid.clear();
		// cellCentroid.resize(3,0.0);
		// area = 0.0;
		// double VSn_x = 0.0; double VSn_y = 0.0; double VSn_z = 0.0;
		// double Vx0 = 0.0;
		// double Vy0 = 0.0;
		// double Vz0 = 0.0;
		// if(two==0){
			// // Vx0 = Vx[0]; Vy0 = Vy[0]; Vz0 = Vz[0];
			// double x_avg, y_avg, z_avg;
			// x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)Vx.size();
			// y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)Vx.size();
			// z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)Vx.size();
			
			// Vx0 = x_avg; Vy0 = y_avg; Vz0 = z_avg;
		// }
		// else{
			// Vx0 = x; Vy0 = y; Vz0 = z;
		// }
		// VSn = 0.0;
		// x = 0.0;  y = 0.0; z = 0.0;
		
		// for(int i = 0; i < n; ++i ) {
			// int b=i;
			// int c=i+1;
			// if(i==n-1) c=0;
			
			// double vect_A[3];
			// double vect_B[3];
			
			// vect_A[0] = Vx[b]-Vx0; vect_A[1] = Vy[b]-Vy0; vect_A[2] = Vz[b]-Vz0;
			// vect_B[0] = Vx[c]-Vx0; vect_B[1] = Vy[c]-Vy0; vect_B[2] = Vz[c]-Vz0;
			// // vect_A[0] = Vx[b]; vect_A[1] = Vy[b]; vect_A[2] = Vz[b];
			// // vect_B[0] = Vx[c]; vect_B[1] = Vy[c]; vect_B[2] = Vz[c];
		  
			// double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
			// double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
			// double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

			// unitNormals[0] += Nx_tmp;
			// unitNormals[1] += Ny_tmp;
			// unitNormals[2] += Nz_tmp;
			
			// double tmp_area = sqrt(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp);
			// double rx = (Vx0+Vx[b]+Vx[c])/ 3.0;
			// double ry = (Vy0+Vy[b]+Vy[c])/ 3.0;
			// double rz = (Vz0+Vz[b]+Vz[c])/ 3.0;
			// x += rx*tmp_area; y += ry*tmp_area; z += rz*tmp_area;
			
			// VSn += Vx0*Nx_tmp;
			// VSn += Vy0*Ny_tmp;
			// VSn += Vz0*Nz_tmp;
			// // VSn_x += Vx[a]*Nx_tmp; VSn_y += Vy[a]*Ny_tmp; VSn_z += Vz[a]*Nz_tmp;
			// // area += tmp_area;
			
			// cellCentroid[0] += 2.0*Nx_tmp*(
				// pow(Vx0+Vx[b],2.0) + pow(Vx[b]+Vx[c],2.0) + pow(Vx[c]+Vx0,2.0));
			// cellCentroid[1] += 2.0*Ny_tmp*(
				// pow(Vy0+Vy[b],2.0) + pow(Vy[b]+Vy[c],2.0) + pow(Vy[c]+Vy0,2.0));
			// cellCentroid[2] += 2.0*Nz_tmp*(
				// pow(Vz0+Vz[b],2.0) + pow(Vz[b]+Vz[c],2.0) + pow(Vz[c]+Vz0,2.0));
			
		// }
		
		// // VSn = VSn_x + VSn_y + VSn_z;
		
		// double mag_unitNormals = sqrt(
			// unitNormals[0]*unitNormals[0]+
			// unitNormals[1]*unitNormals[1]+
			// unitNormals[2]*unitNormals[2]);
		// unitNormals[0] /= mag_unitNormals;
		// unitNormals[1] /= mag_unitNormals;
		// unitNormals[2] /= mag_unitNormals;
		
		// area = mag_unitNormals;
		
		// x /= area;
		// y /= area;
		// z /= area;
	
	// }
	
	// if(area < std::numeric_limits<double>::min()) {
		// cerr << endl;
		// cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		// cerr << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
// }

  

// /*
// https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
// */
// void SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
	// int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	// vector<double>& unitNormals, double& area,
	// double& x, double& y, double& z ){
		
	// double sumArea = 0.0;
	// unitNormals.clear();
	// unitNormals.resize(3,0.0);
	// x = 0.0;  y = 0.0; z = 0.0;
	// for(int i = 0; i < n; ++i ) {
		// double vect_A[3];
		// double vect_B[3];
		
		// vect_A[0] = Vx[i]-Vx[0];
		// vect_A[1] = Vy[i]-Vy[0];
		// vect_A[2] = Vz[i]-Vz[0];
		
		// int ne = i+1;
		// if(i==n-1) ne = 0;
		// vect_B[0] = Vx[ne]-Vx[0];
		// vect_B[1] = Vy[ne]-Vy[0];
		// vect_B[2] = Vz[ne]-Vz[0];
	  
		// double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
		// double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
		// double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

		// unitNormals[0] += Nx_tmp;
		// unitNormals[1] += Ny_tmp;
		// unitNormals[2] += Nz_tmp;

		// double area_tmp = sqrt(Nx_tmp*Nx_tmp + Ny_tmp*Ny_tmp + Nz_tmp*Nz_tmp);
		// sumArea += area_tmp;
		// x += area_tmp*(Vx[0]+Vx[i]+Vx[ne])/3.0;
		// y += area_tmp*(Vy[0]+Vy[i]+Vy[ne])/3.0;
		// z += area_tmp*(Vz[0]+Vz[i]+Vz[ne])/3.0;
	// }
	
	// x /= sumArea;
	// y /= sumArea;
	// z /= sumArea;
	
	// unitNormals[0] /= sumArea;
	// unitNormals[1] /= sumArea;
	// unitNormals[2] /= sumArea;
	
	// area = sumArea;
	
	// if(area < std::numeric_limits<double>::min()) {
		// cerr << endl;
		// cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		// cerr << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	
// }

/* 
Gradient Calculation Methods on ArbitraryPolyhedral Unstructured Meshes for Cell-CenteredCFD Solvers
*/
void SEMO_Mesh_Geometric::calcUnitNormals_ArbitraryPolyhedral(
	int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	vector<double>& unitNormals, double& area,
	double& x, double& y, double& z){
	
	// int nPoint = Vx.size();

	// n = nPoint;

	// midpoint of face
	vector<double> rf(3,0.0);
	double dbN = (double)n;
	for(int i=0; i<n; ++i){
		rf[0] += Vx[i] / dbN;
		rf[1] += Vy[i] / dbN;
		rf[2] += Vz[i] / dbN;
	}
	
	int num=0;
	
	vector<double> A_f(3,0.0);
	while(1){
	
		vector<vector<double>> A_tri(n,vector<double>(3,0.0));
		for(int i=0; i<n; ++i){
			double vect_A[3];
			double vect_B[3];
			
			vect_A[0] = Vx[i]-rf[0];
			vect_A[1] = Vy[i]-rf[1];
			vect_A[2] = Vz[i]-rf[2];
			
			if(i==n-1){
				vect_B[0] = Vx[0]-rf[0];
				vect_B[1] = Vy[0]-rf[1];
				vect_B[2] = Vz[0]-rf[2];
			}
			else{
				vect_B[0] = Vx[i+1]-rf[0];
				vect_B[1] = Vy[i+1]-rf[1];
				vect_B[2] = Vz[i+1]-rf[2];
			}
			
			// cross product
			A_tri[i][0] = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
			A_tri[i][1] = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
			A_tri[i][2] = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );
		}
		
		
		A_f[0] = 0.0;
		A_f[1] = 0.0;
		A_f[2] = 0.0;
		for(int i=0; i<n; ++i){
			A_f[0] += A_tri[i][0];
			A_f[1] += A_tri[i][1];
			A_f[2] += A_tri[i][2];
		}
		double mag_A_f = sqrt(A_f[0]*A_f[0] + A_f[1]*A_f[1] + A_f[2]*A_f[2]);
		
		// new midpoint of face
		vector<double> rf_new(3,0.0);
		for(int i=0; i<n; ++i){
			double mag_A_tri = 
				sqrt(A_tri[i][0]*A_tri[i][0] + A_tri[i][1]*A_tri[i][1] + A_tri[i][2]*A_tri[i][2]);
				
			if(i==n-1){
				rf_new[0] += 1.0/mag_A_f * mag_A_tri*(Vx[i]+Vx[0]+rf[0])/3.0;
				rf_new[1] += 1.0/mag_A_f * mag_A_tri*(Vy[i]+Vy[0]+rf[1])/3.0;
				rf_new[2] += 1.0/mag_A_f * mag_A_tri*(Vz[i]+Vz[0]+rf[2])/3.0;
			}
			else{
				rf_new[0] += 1.0/mag_A_f * mag_A_tri*(Vx[i]+Vx[i+1]+rf[0])/3.0;
				rf_new[1] += 1.0/mag_A_f * mag_A_tri*(Vy[i]+Vy[i+1]+rf[1])/3.0;
				rf_new[2] += 1.0/mag_A_f * mag_A_tri*(Vz[i]+Vz[i+1]+rf[2])/3.0;
			}
		}
		
		double error = 
			abs(rf_new[0]-rf[0]) +
			abs(rf_new[1]-rf[1]) +
			abs(rf_new[2]-rf[2]);
			
			
			

			break;
			// rf[0] = rf_new[0];
			// rf[1] = rf_new[1];
			// rf[2] = rf_new[2];
			// break;
			
		cout << endl;
		cout << error << " " << mag_A_f << endl;
		cout << num << endl;
		cout << rf[0] << " " << rf[1] << " " << rf[2] << endl;
		cout << rf_new[0] << " " << rf_new[1] << " " << rf_new[2] << endl;
			
		// if( isnan(error) ){
			
			// for(int i=0; i<n; ++i){
				// cout << Vx[i] << " " << Vy[i] << " " << Vz[i] << endl;
			// }
			
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }

		rf[0] = rf[0] + 1.0*(rf_new[0]-rf[0]);
		rf[1] = rf[1] + 1.0*(rf_new[1]-rf[1]);
		rf[2] = rf[2] + 1.0*(rf_new[2]-rf[2]);
			
		if( error < 1.e-8 ) break;
		
		++num;
	
	}
	
	
	x = rf[0];
	y = rf[1];
	z = rf[2];
	
	
	area = sqrt(A_f[0]*A_f[0] + A_f[1]*A_f[1] + A_f[2]*A_f[2]);
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	unitNormals.resize(3,0.0);
	unitNormals[0] = A_f[0]/area; 
	unitNormals[1] = A_f[1]/area; 
	unitNormals[2] = A_f[2]/area;
	
	
	// vector<double> tmpA(3,0.0);
	// double vect_A[3];
	// double vect_B[3];
	
	// vect_A[0] = Vx[1]-rf[0];
	// vect_A[1] = Vy[1]-rf[1];
	// vect_A[2] = Vz[1]-rf[2];
	// vect_B[0] = Vx[0]-rf[0];
	// vect_B[1] = Vy[0]-rf[1];
	// vect_B[2] = Vz[0]-rf[2];
	// tmpA[0] = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
	// tmpA[1] = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
	// tmpA[2] = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );
	
	// double length = sqrt(pow(tmpA[0],2.0)+pow(tmpA[1],2.0)+pow(tmpA[2],2.0));
	// // if(length < std::numeric_limits<double>::min()) {
		// // cerr << endl;
		// // cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		// // cerr << endl;
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // }
	
	// unitNormals.resize(3,0.0);
	// unitNormals[0] = tmpA[0]/length; 
	// unitNormals[1] = tmpA[1]/length; 
	// unitNormals[2] = tmpA[2]/length;
	
	
	// area = sqrt(A_f[0]*A_f[0]*unitNormals[0]*unitNormals[0] + 
				// A_f[1]*A_f[1]*unitNormals[1]*unitNormals[1] + 
				// A_f[2]*A_f[2]*unitNormals[2]*unitNormals[2]);
	
	
	
	
}



void SEMO_Mesh_Geometric::calcUnitNormals(
	vector<double>& unitNormals,
	double x1,double y1,double z1,
	double x2,double y2,double z2,
	double x3,double y3,double z3){
	
	double v1[3] = {x2-x1,y2-y1,z2-z1};
	double v2[3] = {x3-x1,y3-y1,z3-z1};
	
	unitNormals.resize(3,0.0);
	
	unitNormals[0] = v1[1] * v2[2] - v1[2] * v2[1];
	unitNormals[1] = -(v1[0] * v2[2] - v1[2] * v2[0]);
	unitNormals[2] = v1[0] * v2[1] - v1[1] * v2[0];
	double length = sqrt(pow(unitNormals[0], 2) + pow(unitNormals[1], 2) + pow(unitNormals[2], 2));
	if(length < std::numeric_limits<double>::min()) {
		cerr << "from calcUnitNormals, length ~= 0.0" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
		

	unitNormals[0] = unitNormals[0]/length;
	unitNormals[1] = unitNormals[1]/length;
	unitNormals[2] = unitNormals[2]/length;
}


double SEMO_Mesh_Geometric::calcArea3dPolygon( 
	int n, 
	vector<double> Vx, vector<double> Vy, vector<double> Vz,
	double Nx, double Ny, double Nz ) {
		
	double area = 0;
	double an, ax, ay, az; // abs value of normal and its coords
	int  coord;           // coord to ignore: 1=x, 2=y, 3=z
	int  i, j, k;         // loop indices
	
	printf("%d\n",n);

	for (auto item : Vx) {
		printf("%lf\n",item);
	}
	for (auto item : Vy) {
		printf("%lf\n",item);
	}
	for (auto item : Vz) {
		printf("%lf\n",item);
	}
	printf("%lf %lf %lf\n",Vx[0],Vy[0],Vz[0]);
	printf("%lf %lf %lf\n",Nx,Ny,Nz);

	if (n < 3) return 0;  // a degenerate polygon

	// select largest abs coordinate to ignore for projection
	ax = (Nx>0.0 ? Nx : -Nx);    // abs x-coord
	ay = (Ny>0.0 ? Ny : -Ny);    // abs y-coord
	az = (Nz>0.0 ? Nz : -Nz);    // abs z-coord

	coord = 3;                    // ignore z-coord
	if (ax > ay) {
		if (ax > az) coord = 1;   // ignore x-coord
	}
	else if (ay > az) coord = 2;  // ignore y-coord

	// compute area of the 2D projection
	switch (coord) {
	  case 1:
		for (i=1, j=2, k=0; i<n; i++, j++, k++)
			area += (Vy[i] * (Vz[j] - Vz[k]));
		break;
	  case 2:
		for (i=1, j=2, k=0; i<n; i++, j++, k++)
			area += (Vz[i] * (Vx[j] - Vx[k]));
		break;
	  case 3:
		for (i=1, j=2, k=0; i<n; i++, j++, k++)
			area += (Vx[i] * (Vy[j] - Vy[k]));
		break;
	}
	switch (coord) {    // wrap-around term
	  case 1:
		area += (Vy[n] * (Vz[1] - Vz[n-1]));
		break;
	  case 2:
		area += (Vz[n] * (Vx[1] - Vx[n-1]));
		break;
	  case 3:
		area += (Vx[n] * (Vy[1] - Vy[n-1]));
		break;
	}

	// scale to get area before projection
	an = sqrt( ax*ax + ay*ay + az*az); // length of normal vector
	switch (coord) {
	  case 1:
		area *= (an / (2.0 * Nx));
		break;
	  case 2:
		area *= (an / (2.0 * Ny));
		break;
	  case 3:
		area *= (an / (2.0 * Nz));
	}
	
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	return area;
}






bool SEMO_Mesh_Geometric::isInsidePolygon(
vector<double> Vx, vector<double> Vy, vector<double> Vz,
double Vx_tar, double Vy_tar, double Vz_tar,
double maxAng, double eps
) {

	int n = Vx.size();
	double degree = 0.0;
	bool clockwise_x = true;
	bool clockwise_y = true;
	bool clockwise_z = true;
	for (int i = 0; i < n; ++i) {
		int a = i;
		int b = i + 1;
		if (i == n - 1) b = 0;

		vector<double> ta(3, 0.0);
		vector<double> tb(3, 0.0);
		ta[0] = Vx[a] - Vx_tar;
		ta[1] = Vy[a] - Vy_tar;
		ta[2] = Vz[a] - Vz_tar;
		double mag_ta = sqrt(ta[0] * ta[0] + ta[1] * ta[1] + ta[2] * ta[2]);

		tb[0] = Vx[b] - Vx_tar;
		tb[1] = Vy[b] - Vy_tar;
		tb[2] = Vz[b] - Vz_tar;
		double mag_tb = sqrt(tb[0] * tb[0] + tb[1] * tb[1] + tb[2] * tb[2]);

		vector<double> cross(3, 0.0);
		cross[0] = (ta[1] * tb[2] - ta[2] * tb[1]) / mag_ta / mag_tb;
		cross[1] = (ta[2] * tb[0] - ta[0] * tb[2]) / mag_ta / mag_tb;
		cross[2] = (ta[0] * tb[1] - ta[1] * tb[0]) / mag_ta / mag_tb;
		if (cross[0] < eps && cross[0] > 0.0) cross[0] = eps;
		if (cross[1] < eps && cross[1] > 0.0) cross[1] = eps;
		if (cross[2] < eps && cross[2] > 0.0) cross[2] = eps;

		if (cross[0] < 0.0 && cross[0] > -eps) cross[0] = eps;
		if (cross[1] < 0.0 && cross[1] > -eps) cross[1] = eps;
		if (cross[2] < 0.0 && cross[2] > -eps) cross[2] = eps;

		double dist_a_b = sqrt(pow(Vx[a]- Vx[b],2.0)+ pow(Vy[a] - Vy[b], 2.0)+ pow(Vz[a] - Vz[b], 2.0));
		double dist_a_t = sqrt(pow(ta[0], 2.0) + pow(ta[1], 2.0) + pow(ta[2], 2.0));
		double dist_b_t = sqrt(pow(tb[0], 2.0) + pow(tb[1], 2.0) + pow(tb[2], 2.0));

		bool clockwise = false;
		if (i == 0) {
			clockwise_x = cross[0] < 0;
			clockwise_y = cross[1] < 0;
			clockwise_z = cross[2] < 0;
			clockwise = true;
		}
		else {
			if (
				clockwise_x == (cross[0] < 0) &&
				clockwise_y == (cross[1] < 0) &&
				clockwise_z == (cross[2] < 0)) {
				clockwise = true;
			}
			else {
				clockwise = false;
			}
		}
		
		if (clockwise) {
			degree += acos(
				(dist_a_t * dist_a_t + dist_b_t * dist_b_t - dist_a_b * dist_a_b) 
				/ (2.0 * dist_a_t * dist_b_t));
		}
		else {
			degree -= acos(
				(dist_a_t * dist_a_t + dist_b_t * dist_b_t - dist_a_b * dist_a_b) 
				/ (2.0 * dist_a_t * dist_b_t));
		}

	}

	degree = degree * 180.0/3.141592;
	if (abs(degree - 360.0) <= maxAng) {
		return true;
	}
	else {
		return false;
	}
}




bool SEMO_Mesh_Geometric::isInsidePolyhedron(
	vector<double> Fx, vector<double> Fy, vector<double> Fz,
	vector<double> Nx, vector<double> Ny, vector<double> Nz,
	double Vx_tar, double Vy_tar, double Vz_tar
) {

	int n = Fx.size();
	for (int i = 0; i < n; ++i) {

		vector<double> Vec(3, 0.0);

		Vec[0] = Fx[i] - Vx_tar;
		Vec[1] = Fy[i] - Vy_tar;
		Vec[2] = Fz[i] - Vz_tar;
		double mag = sqrt(Vec[0] * Vec[0] + Vec[1] * Vec[1] + Vec[2] * Vec[2]);

		double direction = Vec[0] * Nx[i] + Vec[1] * Ny[i] + Vec[2] * Nz[i];
		direction /= mag;
		if (direction < -0.9) {
			// cout << direction << endl;
			return false;
		}
	}

	return true;
}







void SEMO_Mesh_Geometric::init(SEMO_Mesh_Builder& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0){
		// cout << "┌────────────────────────────────────────────────────" << endl;
		// cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
	}
	
	
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(auto& cell : mesh.cells){
		cell.volume = 0.0;
		cell.x = 0.0;
		cell.y = 0.0;
		cell.z = 0.0;
	}
	
	// polygon face normal vectors & polygon face area
	// polygon face center x,y,z
	// 3D Polygon Areas
	for(auto& face : mesh.faces){
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : face.points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		
		double VSn=0.0;
		vector<double> cellCentroid;
		
		this->calcUnitNormals_Area3dPolygon(
			face.points.size(), Vx,Vy,Vz,
			face.unitNormals, face.area,
			face.x, face.y, face.z,
			VSn, cellCentroid);
			
		mesh.cells[face.owner].volume += VSn / 3.0;
		// mesh.cells[face.owner].volume += VSn;
		mesh.cells[face.owner].x += cellCentroid[0];
		mesh.cells[face.owner].y += cellCentroid[1];
		mesh.cells[face.owner].z += cellCentroid[2];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[face.neighbour].volume -= VSn / 3.0;
			// mesh.cells[face.neighbour].volume -= VSn;
			mesh.cells[face.neighbour].x -= cellCentroid[0];
			mesh.cells[face.neighbour].y -= cellCentroid[1];
			mesh.cells[face.neighbour].z -= cellCentroid[2];
		}
			
		// double x_avg, y_avg, z_avg;
		// x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)face.points.size();
		// y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)face.points.size();
		// z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)face.points.size();
		
		// face.x = x_avg; face.y = y_avg; face.z = z_avg;
		
		// if(abs(x_avg-face.x)+abs(y_avg-face.y)+abs(z_avg-face.z)>1.e-3){
			// cout << endl;
			// cout << Vx.size() << endl;
			// cout << x_avg << " " << y_avg << " " << z_avg << " " << endl;
			// cout << face.x << " " << face.y << " " << face.z << " " << endl;
			// // abs(x_avg-face.x)+abs(y_avg-face.y)+abs(z_avg-face.z) << endl;
		// }
		
		// this->calcUnitNormals_ArbitraryPolyhedral(
			// face.points.size(), Vx,Vy,Vz,
			// face.unitNormals, face.area, 
			// face.x, face.y, face.z );
		// face.x = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)face.points.size();
		// face.y = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)face.points.size();
		// face.z = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)face.points.size();
		
		// if(face.unitNormals[2]>0.0001 && face.unitNormals[2]<1.0){
		// cout << face.unitNormals[2] << endl;
		// }
		
	}
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	for(auto& cell : mesh.cells){
		cell.x *= 1.0/(24.0*2.0*cell.volume);
		cell.y *= 1.0/(24.0*2.0*cell.volume);
		cell.z *= 1.0/(24.0*2.0*cell.volume);
		
		// cell.x /= cell.volume;
		// cell.y /= cell.volume;
		// cell.z /= cell.volume;
		
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : cell.points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		
		// double volume = 0.0;
		// vector<double> ABv(3,0.0);
		// vector<double> ACv(3,0.0);
		// vector<double> ADv(3,0.0);
		// ABv[0] = Vx[1] - Vx[0];
		// ABv[1] = Vy[1] - Vy[0];
		// ABv[2] = Vz[1] - Vz[0];
		
		// ACv[0] = Vx[2] - Vx[0];
		// ACv[1] = Vy[2] - Vy[0];
		// ACv[2] = Vz[2] - Vz[0];
		
		// ADv[0] = Vx[3] - Vx[0];
		// ADv[1] = Vy[3] - Vy[0];
		// ADv[2] = Vz[3] - Vz[0];
		
		// double Nx_tmp = ( ABv[1] * ACv[2] - ABv[2] * ACv[1] );
		// double Ny_tmp = ( ABv[2] * ACv[0] - ABv[0] * ACv[2] );
		// double Nz_tmp = ( ABv[0] * ACv[1] - ABv[1] * ACv[0] );
		
		// double aaaa = Nx_tmp*ADv[0] + Ny_tmp*ADv[1] + Nz_tmp*ADv[2];
		// aaaa = 1.0/6.0*abs(aaaa);
		
		
		// if(Vx.size() == 4 && abs(aaaa-cell.volume) > 1.e-12*aaaa){
			// cout << endl;
			// cout << "volume error : " << aaaa << " " << cell.volume << endl;
		// }
		
		// double x_avg, y_avg, z_avg;
		// x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)cell.points.size();
		// y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)cell.points.size();
		// z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)cell.points.size();
		
		// cell.x = x_avg; cell.y = y_avg; cell.z = z_avg;
		
		// if(abs(x_avg-cell.x)+abs(y_avg-cell.y)+abs(z_avg-cell.z)>1.e-5){
			// cout << endl;
			// cout << x_avg << " " << y_avg << " " << z_avg << " " << endl;
			// cout << cell.x << " " << cell.y << " " << cell.z << " " << endl;
		// }
	}
	
	
	
	
	
	// for(auto& cell : mesh.cells){
		
		// vector<double> Vx, Vy, Vz;
		// for(auto iPoint : cell.points){
			// Vx.push_back(mesh.points[iPoint].x);
			// Vy.push_back(mesh.points[iPoint].y);
			// Vz.push_back(mesh.points[iPoint].z);
		// }
		
		// double x_avg, y_avg, z_avg;
		// x_avg = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)cell.points.size();
		// y_avg = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)cell.points.size();
		// z_avg = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)cell.points.size();
		
		// double cellVol = 0.0;
		// double cellX = 0.0;
		// double cellY = 0.0;
		// double cellZ = 0.0;
		
		// for(auto& iFace : cell.faces){
			// auto& face = mesh.faces[iFace];
			
			// vector<double> fVx, fVy, fVz;
			// for(auto& iPoint : face.points){
				// fVx.push_back(mesh.points[iPoint].x);
				// fVy.push_back(mesh.points[iPoint].y);
				// fVz.push_back(mesh.points[iPoint].z);
			// }
			
			// double fx_avg, fy_avg, fz_avg;
			// fx_avg = accumulate(fVx.begin(), fVx.end(), 0.0) / (double)fVx.size();
			// fy_avg = accumulate(fVy.begin(), fVy.end(), 0.0) / (double)fVy.size();
			// fz_avg = accumulate(fVz.begin(), fVz.end(), 0.0) / (double)fVz.size();
			
			// vector<double> faceNormals(3,0.0);
			// double faceArea = 0.0;
			// double faceX = 0.0;
			// double faceY = 0.0;
			// double faceZ = 0.0;
			
			// int n = fVx.size();
			// for(int i = 0; i < n; ++i ) {
				// int b=i;
				// int c=i+1;
				// if(i==n-1) c=0;
				
				// double vect_A[3];
				// double vect_B[3];
				
				// vect_A[0] = fVx[b]-fx_avg; vect_A[1] = fVy[b]-fy_avg; vect_A[2] = fVz[b]-fz_avg;
				// vect_B[0] = fVx[c]-fx_avg; vect_B[1] = fVy[c]-fy_avg; vect_B[2] = fVz[c]-fz_avg;
			  
				// double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
				// double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
				// double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );
				
				// double dV = 0.0;
				// dV += Nx_tmp*(fx_avg-x_avg)/3.0;
				// dV += Ny_tmp*(fy_avg-y_avg)/3.0;
				// dV += Nz_tmp*(fz_avg-z_avg)/3.0;
				// dV = abs(dV);
				
				// double dX = dV*0.25*(fVx[b]+fVx[c]+fx_avg+x_avg);
				// double dY = dV*0.25*(fVy[b]+fVy[c]+fy_avg+y_avg);
				// double dZ = dV*0.25*(fVz[b]+fVz[c]+fz_avg+z_avg);
				
				// cellVol += dV;
				// cellX += dX;
				// cellY += dY;
				// cellZ += dZ;
				
				// double dArea = sqrt(Nx_tmp*Nx_tmp+Ny_tmp*Ny_tmp+Nz_tmp*Nz_tmp);
				
				// faceNormals[0] += Nx_tmp;
				// faceNormals[1] += Ny_tmp;
				// faceNormals[2] += Nz_tmp;
				// faceArea += dArea;
				// faceX += dArea*(fVx[b]+fVx[c]+fx_avg)/3.0;
				// faceY += dArea*(fVy[b]+fVy[c]+fy_avg)/3.0;
				// faceZ += dArea*(fVz[b]+fVz[c]+fz_avg)/3.0;
				
			// }
		
			// faceX /= faceArea;
			// faceY /= faceArea;
			// faceZ /= faceArea;
			
			// double mag_unitNormal = sqrt(
				// faceNormals[0]*faceNormals[0] +
				// faceNormals[1]*faceNormals[1] +
				// faceNormals[2]*faceNormals[2]);
				
			// face.unitNormals[0] = faceNormals[0]/mag_unitNormal;
			// face.unitNormals[1] = faceNormals[1]/mag_unitNormal;
			// face.unitNormals[2] = faceNormals[2]/mag_unitNormal;
			
			// face.area = mag_unitNormal;
			// face.x = faceX;
			// face.y = faceY;
			// face.z = faceZ;
		
		// }
		
		// cellX /= cellVol;
		// cellY /= cellVol;
		// cellZ /= cellVol;
		
		// cell.volume = cellVol;
		// cell.x = cellX;
		// cell.y = cellY;
		// cell.z = cellZ;
		
		
		// // cell.x = x_avg; cell.y = y_avg; cell.z = z_avg;
	// }
	
	
	
	
	
	
	
	int num_wrong_face_centroid = 0;
	int num_wrong_cell_centroid = 0;
	int id_ = 0;
	for(auto& face : mesh.faces){
		
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : face.points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}

		if( isInsidePolygon(Vx, Vy, Vz, face.x, face.y,face.z, 180.0, 1.e-8) == false ){
			++num_wrong_face_centroid;
			// cout << endl;
			// cout << " Face Centroid is NOT inside polygon face, id = " << id_ << endl;
			// for(int j=0; j<Vx.size(); ++j){
				// cout << " face points = " << Vx[j] << " " << Vy[j] << " " << Vz[j] << " " << endl;
			// }
			// cout << " face centroid = " << face.x << " " << face.y << " " << face.z << " " << endl;
		}
		++id_;
	}
	
	id_ = 0;
	for(auto& cell : mesh.cells){

		vector<double> Fx, Fy, Fz;
		vector<double> Nx, Ny, Nz;
		for(auto iFace : cell.faces){
			auto& face = mesh.faces[iFace];
			Fx.push_back(face.x);
			Fy.push_back(face.y);
			Fz.push_back(face.z);
			
			if(face.owner == id_){
				Nx.push_back(face.unitNormals[0]);
				Ny.push_back(face.unitNormals[1]);
				Nz.push_back(face.unitNormals[2]);
			}
			else{
				Nx.push_back(-face.unitNormals[0]);
				Ny.push_back(-face.unitNormals[1]);
				Nz.push_back(-face.unitNormals[2]);
			}
		}

		if( isInsidePolyhedron(Fx, Fy, Fz, Nx, Ny, Nz, cell.x, cell.y,cell.z) == false ){
			++num_wrong_cell_centroid;
			// cout << endl;
			// cout << " Cell Centroid is NOT inside polyhedron cell, id = " << id_ << endl;
			// for(int j=0; j<Fx.size(); ++j){
				// cout << " face centroid = " << Fx[j] << " " << Fy[j] << " " << Fz[j] << " " << endl;
			// }
			// for(int j=0; j<Nx.size(); ++j){
				// cout << " face normal vector = " << Nx[j] << " " << Ny[j] << " " << Nz[j] << " " << endl;
			// }
			// cout << " cell centroid = " << cell.x << " " << cell.y << " " << cell.z << " " << endl;
		}
		++id_;
	}

	int num_wrong_face_centroid_global;
	int num_wrong_cell_centroid_global;
	MPI_Allreduce(&num_wrong_face_centroid, &num_wrong_face_centroid_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&num_wrong_cell_centroid, &num_wrong_cell_centroid_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	// if(rank==0) cout << endl;
	// if(rank==0) cout << "| # of Wrong face centroid = " << num_wrong_face_centroid_global << endl;
	// if(rank==0) cout << "| # of Wrong cell centroid = " << num_wrong_cell_centroid_global << endl;
	
	
	
	
	int tmp_icell = 0;
	double myMaxVolume = 0.0;
	double myMinVolume = 1.e10;
	for(auto& cell : mesh.cells){
		
		// cout << "aaa" << endl;
		// cout << cell.volume << endl;
		// cout << cell.x << " " << cell.y << " " << cell.z << endl;
		// for(auto& i : cell.faces){
			// cout << " face : " << mesh.faces[i].x << " " << mesh.faces[i].y << " " << mesh.faces[i].z << endl;
			
			// for(auto& j : mesh.faces[i].points){
				// cout << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << endl;
			// }
		// }
		
		myMaxVolume = max(myMaxVolume,cell.volume);
		myMinVolume = min(myMinVolume,cell.volume);
		
		
		if(cell.volume < std::numeric_limits<double>::min()) {
			cerr << endl;
			cerr << endl;
			cerr << "  #error, from calc cell volume, cell volume = " << cell.volume << " < cpu_min_val " << endl;
			cerr << tmp_icell << endl;
			cerr << "Cell level = " << cell.level << endl;
			cerr << endl;
			int ii = 0;
			int jj = 0;
			for(auto& i : cell.faces){
				if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
					cout << i << " " << "INTERNAL_FACE" << endl;
				}
				else if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
					cout << i << " " << "PROCESSOR_FACE" << endl;
				}
				else if(mesh.faces[i].getType() == SEMO_Types::BOUNDARY_FACE){
					cout << i << " " << "BOUNDARY_FACE" << endl;
				}
				
				// cout << " level = " << mesh.faces[i].level << endl;
				// cout << " owner = " << mesh.faces[i].owner << " , neighbour = " << mesh.faces[i].neighbour << endl;
				// cout << " owner lv = " << mesh.cells[mesh.faces[i].owner].level << " , neighbour lv = " << mesh.cells[mesh.faces[i].neighbour].level << endl;
				// cout << " owner->face vec = " << mesh.faces[i].x-mesh.cells[mesh.faces[i].owner].x << " " << mesh.faces[i].y-mesh.cells[mesh.faces[i].owner].y << " " << mesh.faces[i].z-mesh.cells[mesh.faces[i].owner].z << endl;
				// cout << " neighbour->face vec = " << mesh.faces[i].x-mesh.cells[mesh.faces[i].neighbour].x << " " << mesh.faces[i].y-mesh.cells[mesh.faces[i].neighbour].y << " " << mesh.faces[i].z-mesh.cells[mesh.faces[i].neighbour].z << endl;
				// cout << " face normal vec = " << mesh.faces[i].unitNormals[0] << " " << mesh.faces[i].unitNormals[1] << " " << mesh.faces[i].unitNormals[2] << endl;
				
				for(auto& j : mesh.faces[i].points){
					cout << ii << " " << jj << " " << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << endl;
					++jj;
				}
				++ii;
			}
			cerr << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		++tmp_icell;
	}
	
	
	
	double maxVolume;
	double minVolume;
	MPI_Allreduce(&myMaxVolume, &maxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&myMinVolume, &minVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	// if(rank==0) cout << endl;
	// if(rank==0) cout << "| Max Cell Volume = " << maxVolume << endl;
	// if(rank==0) cout << "| Min Cell Volume = " << minVolume << endl;
	

	
	
	for(auto& face : mesh.faces){
		
		face.distCells.resize(3,0.0);

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			face.distCells[0] = mesh.cells[face.neighbour].x - mesh.cells[face.owner].x;
			face.distCells[1] = mesh.cells[face.neighbour].y - mesh.cells[face.owner].y;
			face.distCells[2] = mesh.cells[face.neighbour].z - mesh.cells[face.owner].z;
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			face.distCells[0] = face.x - mesh.cells[face.owner].x;
			face.distCells[1] = face.y - mesh.cells[face.owner].y;
			face.distCells[2] = face.z - mesh.cells[face.owner].z;
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.distCells[0] = face.x - mesh.cells[face.owner].x;
			face.distCells[1] = face.y - mesh.cells[face.owner].y;
			face.distCells[2] = face.z - mesh.cells[face.owner].z;
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// cout<<mesh.displsProcFaces[size-1]<< " " <<mesh.countsProcFaces[size-1]<<endl;
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// cout << size << endl;
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	if(size>1){
		
		vector<double> sendValues0;
		vector<double> recvValues0(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> sendValues1;
		vector<double> recvValues1(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> sendValues2;
		vector<double> recvValues2(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
			
		
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues0.push_back(face.distCells[0]);
				sendValues1.push_back(face.distCells[1]);
				sendValues2.push_back(face.distCells[2]);
			}
		}
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		MPI_Alltoallv( sendValues0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   recvValues0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( sendValues1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   recvValues1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
		MPI_Alltoallv( sendValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   recvValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
					   MPI_COMM_WORLD);
					   
		vector<double>::iterator iter0 = recvValues0.begin();
		vector<double>::iterator iter1 = recvValues1.begin();
		vector<double>::iterator iter2 = recvValues2.begin();
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				face.distCells[0] = face.distCells[0] - *iter0;
				++iter0;
				face.distCells[1] = face.distCells[1] - *iter1;
				++iter1;
				face.distCells[2] = face.distCells[2] - *iter2;
				++iter2;
			}
		}
		
	}
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << endl;
	// if(rank==0) cout << "| 5 finished" << endl;
	// if(rank==0) cout << endl;
	

	
	
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			double ownVol = mesh.cells[face.owner].volume;
			double ngbVol = mesh.cells[face.neighbour].volume;
			
			face.wVC = ngbVol/(ownVol+ngbVol);
			
			double dxFP = face.x-mesh.cells[face.owner].x;
			double dyFP = face.y-mesh.cells[face.owner].y;
			double dzFP = face.z-mesh.cells[face.owner].z;
			double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
			
			double dxFN = face.x-mesh.cells[face.neighbour].x;
			double dyFN = face.y-mesh.cells[face.neighbour].y;
			double dzFN = face.z-mesh.cells[face.neighbour].z;
			double dFN = sqrt(pow(dxFN,2.0)+pow(dyFN,2.0)+pow(dzFN,2.0));
			
			double dxPN = mesh.cells[face.neighbour].x-mesh.cells[face.owner].x;
			double dyPN = mesh.cells[face.neighbour].y-mesh.cells[face.owner].y;
			double dzPN = mesh.cells[face.neighbour].z-mesh.cells[face.owner].z;
			double dPN = sqrt(pow(dxPN,2.0)+pow(dyPN,2.0)+pow(dzPN,2.0));
			
			face.vecPF.resize(3,0.0);
			face.vecPF[0] = dxFP;
			face.vecPF[1] = dyFP;
			face.vecPF[2] = dzFP;
			
			face.vecNF.resize(3,0.0);
			face.vecNF[0] = dxFN;
			face.vecNF[1] = dyFN;
			face.vecNF[2] = dzFN;
			
			// at openfoam 
			double dFP_of = abs(
				dxFP*face.unitNormals[0] +
				dyFP*face.unitNormals[1] +
				dzFP*face.unitNormals[2]);
			double dFN_of = abs(
				dxFN*face.unitNormals[0] +
				dyFN*face.unitNormals[1] +
				dzFN*face.unitNormals[2]);
				
			// if( dFN_of > 0.3*(dFP_of+dFN_of) && dFN_of < 0.7*(dFP_of+dFN_of) ){
				// dFP = dFP_of;
				// dFN = dFN_of;
			// }
			
			
			
			// face.wC = dFN/(dFP+dFN);
			
			// at openfoam
			face.wC = dFN_of/(dFP_of+dFN_of);
			// 절단오차 걸러내기 위한 과정
			face.wC = (float)(round(face.wC * 10000000) / 10000000);
			
			
			
			
			// face.wC = abs( dxFN*dxPN/dPN + dyFN*dyPN/dPN + dzFN*dzPN/dPN ) / ( dPN );
			// if( face.wC < 0.1 || face.wC > 0.9 ){
				// // cout << face.wC << endl;
				// face.wC = dFN_of/(dFP_of+dFN_of);
			// }
			
			
			
			
			// // CFD book
			// double mag_a = min(dFP_of, dFN_of);
			// face.magPN = 2.0 * mag_a;
			
			// original
			double mag_a = dPN;
			face.magPN = dPN;
			face.unitNomalsPN.resize(3,0.0);
			face.unitNomalsPN[0] = dxPN/dPN;
			face.unitNomalsPN[1] = dyPN/dPN;
			face.unitNomalsPN[2] = dzPN/dPN;
			double alphaF = 0.0;
			alphaF += face.unitNormals[0]*face.unitNomalsPN[0];
			alphaF += face.unitNormals[1]*face.unitNomalsPN[1];
			alphaF += face.unitNormals[2]*face.unitNomalsPN[2];
			face.alphaF = 1.0/abs(alphaF);
			
			// skewness
			double D_plane = -(
				face.unitNormals[0]*face.x+
				face.unitNormals[1]*face.y+
				face.unitNormals[2]*face.z);
			double u_line = 
				(face.unitNormals[0]*mesh.cells[face.owner].x+
				 face.unitNormals[1]*mesh.cells[face.owner].y+
				 face.unitNormals[2]*mesh.cells[face.owner].z+
				 D_plane) /
				(face.unitNormals[0]*(mesh.cells[face.owner].x-mesh.cells[face.neighbour].x)+
				 face.unitNormals[1]*(mesh.cells[face.owner].y-mesh.cells[face.neighbour].y)+
				 face.unitNormals[2]*(mesh.cells[face.owner].z-mesh.cells[face.neighbour].z));
			face.vecSkewness.resize(3,0.0);
			face.vecSkewness[0] = mesh.cells[face.owner].x-u_line*(mesh.cells[face.owner].x-mesh.cells[face.neighbour].x);
			face.vecSkewness[1] = mesh.cells[face.owner].y-u_line*(mesh.cells[face.owner].y-mesh.cells[face.neighbour].y);
			face.vecSkewness[2] = mesh.cells[face.owner].z-u_line*(mesh.cells[face.owner].z-mesh.cells[face.neighbour].z);
			face.vecSkewness[0] = face.x-face.vecSkewness[0];
			face.vecSkewness[1] = face.y-face.vecSkewness[1];
			face.vecSkewness[2] = face.z-face.vecSkewness[2];
			
			
			
			vector<double> Pd(3,0.0);
			Pd[0] = face.x - mag_a * face.unitNormals[0];
			Pd[1] = face.y - mag_a * face.unitNormals[1];
			Pd[2] = face.z - mag_a * face.unitNormals[2];
			vector<double> Nd(3,0.0);
			Nd[0] = face.x + mag_a * face.unitNormals[0];
			Nd[1] = face.y + mag_a * face.unitNormals[1];
			Nd[2] = face.z + mag_a * face.unitNormals[2];
			
			face.vecPdP.resize(3,0.0);
			face.vecPdP[0] = Pd[0] - mesh.cells[face.owner].x;
			face.vecPdP[1] = Pd[1] - mesh.cells[face.owner].y;
			face.vecPdP[2] = Pd[2] - mesh.cells[face.owner].z;
			face.vecNdN.resize(3,0.0);
			face.vecNdN[0] = Nd[0] - mesh.cells[face.neighbour].x;
			face.vecNdN[1] = Nd[1] - mesh.cells[face.neighbour].y;
			face.vecNdN[2] = Nd[2] - mesh.cells[face.neighbour].z;
			
			
			face.vecPN.resize(3,0.0);
			face.vecPN[0] = mesh.cells[face.neighbour].x - mesh.cells[face.owner].x;
			face.vecPN[1] = mesh.cells[face.neighbour].y - mesh.cells[face.owner].y;
			face.vecPN[2] = mesh.cells[face.neighbour].z - mesh.cells[face.owner].z;
			
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			face.wVC = 0.5;
			
			face.wC = 0.5;
			// 절단오차 걸러내기 위한 과정
			face.wC = (float)(round(face.wC * 10000000) / 10000000);
			
			
			double dxFP = face.x-mesh.cells[face.owner].x;
			double dyFP = face.y-mesh.cells[face.owner].y;
			double dzFP = face.z-mesh.cells[face.owner].z;
			double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
			
			face.vecPF.resize(3,0.0);
			face.vecPF[0] = dxFP;
			face.vecPF[1] = dyFP;
			face.vecPF[2] = dzFP;
			
			face.vecNF.resize(3,0.0);
			face.vecNF[0] = dxFP;
			face.vecNF[1] = dyFP;
			face.vecNF[2] = dzFP;
			
			double dFP_of = abs(
				dxFP*face.unitNormals[0] +
				dyFP*face.unitNormals[1] +
				dzFP*face.unitNormals[2]);
				
			double mag_a = dFP_of;
			face.magPN = 2.0 * mag_a;
			
			face.unitNomalsPN.resize(3,0.0);
			face.unitNomalsPN[0] = dxFP/dFP;
			face.unitNomalsPN[1] = dyFP/dFP;
			face.unitNomalsPN[2] = dzFP/dFP;
			double alphaF = 0.0;
			alphaF += face.unitNormals[0]*face.unitNomalsPN[0];
			alphaF += face.unitNormals[1]*face.unitNomalsPN[1];
			alphaF += face.unitNormals[2]*face.unitNomalsPN[2];
			face.alphaF = 1.0/abs(alphaF);
			
			// skewness
			face.vecSkewness.resize(3,0.0);
			
			vector<double> dP(3,0.0);
			dP[0] = face.x - mag_a * face.unitNormals[0];
			dP[1] = face.y - mag_a * face.unitNormals[1];
			dP[2] = face.z - mag_a * face.unitNormals[2];
			
			face.vecPdP.resize(3,0.0);
			face.vecPdP[0] = dP[0] - mesh.cells[face.owner].x;
			face.vecPdP[1] = dP[1] - mesh.cells[face.owner].y;
			face.vecPdP[2] = dP[2] - mesh.cells[face.owner].z;
			face.vecNdN.resize(3,0.0);
			face.vecNdN[0] = -face.vecPdP[0];
			face.vecNdN[1] = -face.vecPdP[1];
			face.vecNdN[2] = -face.vecPdP[2];
			
			
			face.vecPN.resize(3,0.0);
			face.vecPN[0] = 2.0 * dxFP;
			face.vecPN[1] = 2.0 * dyFP;
			face.vecPN[2] = 2.0 * dzFP;
			
				
			
		}
		
		// face.wC = 0.5;
	}
	if(size>1){
		
		vector<double> sendValues;
		vector<double> recvValues(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		
		vector<double> sendValues0;
		vector<double> sendValues1;
		vector<double> sendValues2;
		vector<double> recvValues0(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> recvValues1(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> recvValues2(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
			
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues.push_back(mesh.cells[face.owner].volume);
				
				sendValues0.push_back(mesh.cells[face.owner].x);
				sendValues1.push_back(mesh.cells[face.owner].y);
				sendValues2.push_back(mesh.cells[face.owner].z);
			}
		}
		
		
		MPI_Alltoallv( 
          sendValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          recvValues.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          MPI_COMM_WORLD);
		
		MPI_Alltoallv( 
          sendValues0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          recvValues0.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          MPI_COMM_WORLD);
		MPI_Alltoallv( 
          sendValues1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          recvValues1.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          MPI_COMM_WORLD);
		MPI_Alltoallv( 
          sendValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          recvValues2.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_DOUBLE, 
          MPI_COMM_WORLD);
					   
		vector<double>::iterator iter = recvValues.begin();
		vector<double>::iterator iter0 = recvValues0.begin();
		vector<double>::iterator iter1 = recvValues1.begin();
		vector<double>::iterator iter2 = recvValues2.begin();
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				double ownVol = mesh.cells[face.owner].volume;
				double ngbVol = *iter;
				
				face.wVC = ngbVol/(ownVol+ngbVol);
				
				++iter;
				
				double dxFP = face.x-mesh.cells[face.owner].x;
				double dyFP = face.y-mesh.cells[face.owner].y;
				double dzFP = face.z-mesh.cells[face.owner].z;
				double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
				
				double dxFN = face.x-*iter0;
				double dyFN = face.y-*iter1;
				double dzFN = face.z-*iter2;
				double dFN = sqrt(pow(dxFN,2.0)+pow(dyFN,2.0)+pow(dzFN,2.0));
			
				double dxPN = *iter0-mesh.cells[face.owner].x;
				double dyPN = *iter1-mesh.cells[face.owner].y;
				double dzPN = *iter2-mesh.cells[face.owner].z;
				double dPN = sqrt(pow(dxPN,2.0)+pow(dyPN,2.0)+pow(dzPN,2.0));
			
				face.vecPF.resize(3,0.0);
				face.vecPF[0] = dxFP;
				face.vecPF[1] = dyFP;
				face.vecPF[2] = dzFP;
				
				face.vecNF.resize(3,0.0);
				face.vecNF[0] = dxFN;
				face.vecNF[1] = dyFN;
				face.vecNF[2] = dzFN;
				
				// at openfoam 
				double dFP_of = abs(
					dxFP*face.unitNormals[0] +
					dyFP*face.unitNormals[1] +
					dzFP*face.unitNormals[2]);
				double dFN_of = abs(
					dxFN*face.unitNormals[0] +
					dyFN*face.unitNormals[1] +
					dzFN*face.unitNormals[2]);
					
				// if( dFN_of > 0.3*(dFP_of+dFN_of) && dFN_of < 0.7*(dFP_of+dFN_of) ){
					// dFP = dFP_of;
					// dFN = dFN_of;
				// }
				
				
				
				// face.wC = dFN/(dFP+dFN);
				
				// at openfoam
				// 절단오차 걸러내기 위해 float 붙임
				face.wC = dFN_of/(dFP_of+dFN_of);
				// 절단오차 걸러내기 위한 과정
				face.wC = (float)(round(face.wC * 10000000) / 10000000);
					
					
					
				// face.wC = abs( dxFN*dxPN + dyFN*dyPN + dzFN*dzPN ) / ( dxPN*dxPN + dyPN*dyPN + dzPN*dzPN );
				// if( face.wC < 0.1 || face.wC > 0.9 ){
					// face.wC = dFN_of/(dFP_of+dFN_of);
				// }
				
				// if( face.wC < 0.4 ) face.wC = 0.4;
				// if( face.wC > 0.6 ) face.wC = 0.6;
				
				
				
			
				// // CFD book
				// double mag_a = min(dFP_of, dFN_of);
				// face.magPN = 2.0 * mag_a;
				
				// original
				// face.magPN = dPN;
				double mag_a = dPN;
				face.magPN = dPN;
				face.unitNomalsPN.resize(3,0.0);
				face.unitNomalsPN[0] = dxPN/dPN;
				face.unitNomalsPN[1] = dyPN/dPN;
				face.unitNomalsPN[2] = dzPN/dPN;
				double alphaF = 0.0;
				alphaF += face.unitNormals[0]*face.unitNomalsPN[0];
				alphaF += face.unitNormals[1]*face.unitNomalsPN[1];
				alphaF += face.unitNormals[2]*face.unitNomalsPN[2];
				face.alphaF = 1.0/abs(alphaF);


				double D_plane = -(
					face.unitNormals[0]*face.x+
					face.unitNormals[1]*face.y+
					face.unitNormals[2]*face.z);
				double u_line = 
					(face.unitNormals[0]*mesh.cells[face.owner].x+
					 face.unitNormals[1]*mesh.cells[face.owner].y+
					 face.unitNormals[2]*mesh.cells[face.owner].z+
					 D_plane) /
					(face.unitNormals[0]*(mesh.cells[face.owner].x-*iter0)+
					 face.unitNormals[1]*(mesh.cells[face.owner].y-*iter1)+
					 face.unitNormals[2]*(mesh.cells[face.owner].z-*iter2));
				face.vecSkewness.resize(3,0.0);
				face.vecSkewness[0] = mesh.cells[face.owner].x-u_line*(mesh.cells[face.owner].x-*iter0);
				face.vecSkewness[1] = mesh.cells[face.owner].y-u_line*(mesh.cells[face.owner].y-*iter1);
				face.vecSkewness[2] = mesh.cells[face.owner].z-u_line*(mesh.cells[face.owner].z-*iter2);
				face.vecSkewness[0] = face.x-face.vecSkewness[0];
				face.vecSkewness[1] = face.y-face.vecSkewness[1];
				face.vecSkewness[2] = face.z-face.vecSkewness[2];
					
				
				vector<double> dP(3,0.0);
				dP[0] = face.x - mag_a * face.unitNormals[0];
				dP[1] = face.y - mag_a * face.unitNormals[1];
				dP[2] = face.z - mag_a * face.unitNormals[2];
				vector<double> dN(3,0.0);
				dN[0] = face.x + mag_a * face.unitNormals[0];
				dN[1] = face.y + mag_a * face.unitNormals[1];
				dN[2] = face.z + mag_a * face.unitNormals[2];
				
				face.vecPdP.resize(3,0.0);
				face.vecPdP[0] = dP[0] - mesh.cells[face.owner].x;
				face.vecPdP[1] = dP[1] - mesh.cells[face.owner].y;
				face.vecPdP[2] = dP[2] - mesh.cells[face.owner].z;
				face.vecNdN.resize(3,0.0);
				face.vecNdN[0] = dN[0] - *iter0;
				face.vecNdN[1] = dN[1] - *iter1;
				face.vecNdN[2] = dN[2] - *iter2;
				
				
				face.vecPN.resize(3,0.0);
				face.vecPN[0] = *iter0 - mesh.cells[face.owner].x;
				face.vecPN[1] = *iter1 - mesh.cells[face.owner].y;
				face.vecPN[2] = *iter2 - mesh.cells[face.owner].z;
				
				
				++iter0;
				++iter1;
				++iter2;
			}
				
				
				// face.wC = 0.5;
			
		}
		
		
		
	}
	
	
	
	
	
	
	// // wC limiter
	// for(auto& face : mesh.faces){
		
		// double wC = face.wC;
		// wC = 2.0*wC;
		// wC = pow(wC,0.1);
		// face.wC = 0.5*wC;
		
		// // face.wC = max(0.5-0.2,min(0.5+0.2,face.wC));
		// // face.wC = max(0.5-0.4,min(0.5+0.4,face.wC));
	// }
	

	// for(auto& face : mesh.faces){
		// face.wC = 0.5;
	// }
	
	
	
	
	
	
	
	

	// for(auto& cell : mesh.cells){
		// if( cell.x < 0.000465 && cell.x > 0.000462 &&
		// cell.y < -0.0680 && cell.y > -0.0683  &&
		// cell.z < 0.00179 && cell.z > 0.00176
		// ){
			
			// cout << cell.x << " " << cell.y << " " << cell.z << " " << endl;
		// }
		
	// }
	
	// for(auto& face : mesh.faces){
		// if(face.wC > 0.6)
		// cout << face.wC << endl;
	// }
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	// for(auto& face : mesh.faces){
		// // face.wC = 0.5;
		// face.wC = face.wVC;
	// }
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==0) cout << endl;
	// if(rank==0) cout << "| Geometric finished" << endl;
	// if(rank==0) cout << endl;
	
	
	this->setStencil(mesh);
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank==0){
		// cout << "-> completed" << endl;
		// cout << "└────────────────────────────────────────────────────" << endl;
	}
		
}












void SEMO_Mesh_Geometric::setStencil(SEMO_Mesh_Builder& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	for(auto& point : mesh.points){
		point.stencil.clear();
	}
	for(auto& face : mesh.faces){
		face.stencil.clear();
	}
	for(auto& cell : mesh.cells){
		cell.stencil.clear();
	}
	
	
	
	
	// point stencils
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		
		for(auto j : cell.points){
			mesh.points[j].stencil.push_back(i);
		}
		
	}
	

	// face stencils
	for(int i=0; i<mesh.faces.size(); ++i){
		
		auto& face = mesh.faces[i];
		
		for(auto j : face.points){
			for(auto k : mesh.points[j].stencil){
				if ( std::find( face.stencil.begin(), face.stencil.end(), k ) 
					  == face.stencil.end() ) {
					face.stencil.push_back(k);
				}
			}
		}
	}
	
	
	// cell stencils
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		
		for(auto j : cell.points){
			for(auto k : mesh.points[j].stencil){
				if(k==i) continue;
				if ( std::find( cell.stencil.begin(), cell.stencil.end(), k ) 
					  == cell.stencil.end() ) {
					cell.stencil.push_back(k);
				}
			}
		}
	}
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		
		// for(auto j : cell.stencil){
			// cout << j << endl;
		// }
		// cout << cell.stencil.size() << endl;
		
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }


}