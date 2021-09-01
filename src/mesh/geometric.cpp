#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mpi.h>
// #include <numeric>
using namespace std;

#include "geometric.h"  

/*
https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
*/
void SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
	int n, vector<double> Vx, vector<double> Vy, vector<double> Vz,
	vector<double>& unitNormals, double& area,
	double& x, double& y, double& z ){
		
	double sumArea = 0.0;
	unitNormals.clear();
	unitNormals.resize(3,0.0);
	x = 0.0;  y = 0.0; z = 0.0;
	for(int i = 0; i < n; ++i ) {
		double vect_A[3];
		double vect_B[3];
		
		vect_A[0] = Vx[i]-Vx[0];
		vect_A[1] = Vy[i]-Vy[0];
		vect_A[2] = Vz[i]-Vz[0];
		
		int ne = i+1;
		if(i==n-1) ne = 0;
		vect_B[0] = Vx[ne]-Vx[0];
		vect_B[1] = Vy[ne]-Vy[0];
		vect_B[2] = Vz[ne]-Vz[0];
	  
		double Nx_tmp = 0.5*( vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1] );
		double Ny_tmp = 0.5*( vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2] );
		double Nz_tmp = 0.5*( vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0] );

		unitNormals[0] += Nx_tmp;
		unitNormals[1] += Ny_tmp;
		unitNormals[2] += Nz_tmp;

		double area_tmp = sqrt(Nx_tmp*Nx_tmp + Ny_tmp*Ny_tmp + Nz_tmp*Nz_tmp);
		sumArea += area_tmp;
		x += area_tmp*(Vx[0]+Vx[i]+Vx[ne])/3.0;
		y += area_tmp*(Vy[0]+Vy[i]+Vy[ne])/3.0;
		z += area_tmp*(Vz[0]+Vz[i]+Vz[ne])/3.0;
	}
	
	x /= sumArea;
	y /= sumArea;
	z /= sumArea;
	
	unitNormals[0] /= sumArea;
	unitNormals[1] /= sumArea;
	unitNormals[2] /= sumArea;
	
	area = sumArea;
	
	if(area < std::numeric_limits<double>::min()) {
		cerr << endl;
		cerr << "#error, from calcArea3dPolygon, area = " << area <<  " < cpu_min_val " << endl;
		cerr << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
}

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





void SEMO_Mesh_Geometric::init(SEMO_Mesh_Builder& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute geometric (face normal vectors, face area, face center, cell volume) ... ";
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
		
		this->calcUnitNormals_Area3dPolygon(
			face.points.size(), Vx,Vy,Vz,
			face.unitNormals, face.area,
			face.x, face.y, face.z );
			
		// face.x = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)face.points.size();
		// face.y = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)face.points.size();
		// face.z = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)face.points.size();
		
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
	
	
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(auto& cell : mesh.cells){
		cell.volume = 0.0;
		cell.x = 0.0;
		cell.y = 0.0;
		cell.z = 0.0;
	}
	for(auto& face : mesh.faces){
		double rn = face.x*face.unitNormals[0] + face.y*face.unitNormals[1] + face.z*face.unitNormals[2];
		double vol_tmp = rn*face.area/3.0;
		
		mesh.cells[face.owner].volume += vol_tmp;
		
		mesh.cells[face.owner].x += face.x*3.0/4.0 * vol_tmp;
		mesh.cells[face.owner].y += face.y*3.0/4.0 * vol_tmp;
		mesh.cells[face.owner].z += face.z*3.0/4.0 * vol_tmp;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[face.neighbour].volume -= vol_tmp;
			
			mesh.cells[face.neighbour].x -= face.x*3.0/4.0 * vol_tmp;
			mesh.cells[face.neighbour].y -= face.y*3.0/4.0 * vol_tmp;
			mesh.cells[face.neighbour].z -= face.z*3.0/4.0 * vol_tmp;
		}
		
	}
	for(auto& cell : mesh.cells){
		cell.x /= cell.volume;
		cell.y /= cell.volume;
		cell.z /= cell.volume;
	}
	
	
	int tmp_icell = 0;
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
			
			face.wC = dFN/(dFP+dFN);
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			// double dxFP = face.x-mesh.cells[face.owner].x;
			// double dyFP = face.y-mesh.cells[face.owner].y;
			// double dzFP = face.z-mesh.cells[face.owner].z;
			// double dFP = sqrt(pow(dxFP,2.0)+pow(dyFP,2.0)+pow(dzFP,2.0));
			
			// face.wC = dFP;
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			face.wVC = 0.5;
			
			face.wC = 0.5;
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
				
				// face.wC = dFP;
				
				double dxFN = face.x-*iter0;
				double dyFN = face.y-*iter1;
				double dzFN = face.z-*iter2;
				double dFN = sqrt(pow(dxFN,2.0)+pow(dyFN,2.0)+pow(dzFN,2.0));
				
				face.wC = dFN/(dFP+dFN);
				
				++iter0;
				++iter1;
				++iter2;
			}
		}
		
	}
	
	
	
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
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
		
}












void SEMO_Mesh_Geometric::setStencil(SEMO_Mesh_Builder& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	for(auto& point : mesh.points){
		point.stencil.clear();
	}
	for(auto& cell : mesh.cells){
		cell.stencil.clear();
	}
	
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		
		auto& cell = mesh.cells[i];
		
		for(auto j : cell.points){
			mesh.points[j].stencil.push_back(i);
		}
		
	}
	
	
	
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