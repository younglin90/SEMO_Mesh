#pragma once
#include <iostream>
#include <vector>
using namespace std;

class SEMO_Mesh_Geometric{
	public:
		void calcUnitNormals_Area3dPolygon(
			int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
			double unitNormals[], double& area );
		
		
		void calcUnitNormals(
			double unitNormals[],
			double x1,double y1,double z1,
			double x2,double y2,double z2,
			double x3,double y3,double z3);
		
		
		double calcArea3dPolygon( 
			int n, 
			vector<double> Vx, vector<double> Vy, vector<double> Vz,
			double Nx, double Ny, double Nz );
		
		double calcVolumePolyhedron(){
			
		}
		
};
