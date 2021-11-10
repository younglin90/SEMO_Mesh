#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include "build.h"
using namespace std;

class SEMO_Mesh_Geometric{
	public:
		// void calcUnitNormals_Area3dPolygon(
			// int n, vector<long double> Vx, vector<long double> Vy, vector<long double> Vz,
			// vector<double>& unitNormals, double& area,
			// double& x, double& y, double& z,
			// double& VSn, vector<double>& cellCentroid );
			
		void calcUnitNormals_Area3dPolygon(
			int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
			vector<double>& unitNormals, double& area,
			double& x, double& y, double& z,
			double& VSn, vector<double>& cellCentroid );
		
		void calcUnitNormals_ArbitraryPolyhedral(
			int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
			vector<double>& unitNormals, double& area,
			double& x, double& y, double& z );
		
		void calcUnitNormals(
			vector<double>& unitNormals,
			double x1,double y1,double z1,
			double x2,double y2,double z2,
			double x3,double y3,double z3);
		
		
		double calcArea3dPolygon( 
			int n, 
			vector<double> Vx, vector<double> Vy, vector<double> Vz,
			double Nx, double Ny, double Nz );
			

		bool isInsidePolygon(
			vector<double> Vx, vector<double> Vy, vector<double> Vz,
			double Vx_tar, double Vy_tar, double Vz_tar,
			double maxAng, double eps);

		bool isInsidePolyhedron(
			vector<double> Fx, vector<double> Fy, vector<double> Fz,
			vector<double> Nx, vector<double> Ny, vector<double> Nz,
			double Vx_tar, double Vy_tar, double Vz_tar);
		
		double calcVolumePolyhedron(){
			
		}
		
		void init(SEMO_Mesh_Builder& mesh);
		void setStencil(SEMO_Mesh_Builder& mesh);
		
};
