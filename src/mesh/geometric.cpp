#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
// #include <numeric>
using namespace std;

#include "geometric.h"  


void SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
	int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	double unitNormals[], double& area ){
		
	double Nx, Ny, Nz;
	int a,b,c,s;
	
	b=n-2; c=n-1;
	Nx=0.0; Ny=0.0; Nz=0.0;
	for(int i = 0; i < n; ++i ) {
	  a = b; b = c; c = i;
 
	  Nx += Vy[b] * ( Vz[c] - Vz[a] );
	  Ny += Vz[b] * ( Vx[c] - Vx[a] );
	  Nz += Vx[b] * ( Vy[c] - Vy[a] );
	}
	double length = sqrt(pow(Nx, 2) + pow(Ny, 2) + pow(Nz, 2));
	Nx /= length; Ny /= length; Nz /= length;
	
	area = 0.5*length;
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	unitNormals[0] = Nx; unitNormals[1] = Ny; unitNormals[2] = Nz;
	
}

void SEMO_Mesh_Geometric::calcUnitNormals(
	double unitNormals[],
	double x1,double y1,double z1,
	double x2,double y2,double z2,
	double x3,double y3,double z3){
	
	double v1[3] = {x2-x1,y2-y1,z2-z1};
	double v2[3] = {x3-x1,y3-y1,z3-z1};
	
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

