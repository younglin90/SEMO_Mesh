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
	vector<double>& unitNormals, double& area ){
		
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
	unitNormals.resize(3,0.0);
	unitNormals[0] = Nx; unitNormals[1] = Ny; unitNormals[2] = Nz;
	
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
	// 3D Polygon Areas : https://thebuildingcoder.typepad.com/blog/2008/12/3d-polygon-areas.html
	for(auto& face : mesh.faces){
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : face.points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		
		SEMO_Mesh_Geometric::calcUnitNormals_Area3dPolygon(
			face.points.size(), Vx,Vy,Vz,
			face.unitNormals, face.area );
		
		face.x = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)face.points.size();
		face.y = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)face.points.size();
		face.z = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)face.points.size();
		
	}
	
	
	// polyhedron cell volume (Green-Gauss Theorem.)
	for(auto& cell : mesh.cells)
		cell.volume = 0.0;
	
	for(auto& face : mesh.faces){
		
		// cout << iter->owner << endl;
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			mesh.cells[face.owner].volume += 
				(face.x*face.unitNormals[0]+
				 face.y*face.unitNormals[1]+
				 face.z*face.unitNormals[2])*face.area;
			
			mesh.cells[face.neighbour].volume -= 
				(face.x*face.unitNormals[0]+
				 face.y*face.unitNormals[1]+
				 face.z*face.unitNormals[2])*face.area;
		}
		else if(
		face.getType() == SEMO_Types::BOUNDARY_FACE ||
		face.getType() == SEMO_Types::PROCESSOR_FACE
		){
			mesh.cells[face.owner].volume += 
				(face.x*face.unitNormals[0]+
				 face.y*face.unitNormals[1]+
				 face.z*face.unitNormals[2])*face.area;
		}
		
	}
	for(auto& cell : mesh.cells){
		cell.volume /= 3.0;
		
		if(cell.volume < std::numeric_limits<double>::min()) {
			cerr << endl;
			cerr << endl;
			cerr << "  #error, from calc cell volume, cell volume = " << cell.volume << " < cpu_min_val " << endl;
			cerr << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		vector<double> Vx, Vy, Vz;
		for(auto iPoint : cell.points){
			Vx.push_back(mesh.points[iPoint].x);
			Vy.push_back(mesh.points[iPoint].y);
			Vz.push_back(mesh.points[iPoint].z);
		}
		cell.x = accumulate(Vx.begin(), Vx.end(), 0.0) / (double)cell.points.size();
		cell.y = accumulate(Vy.begin(), Vy.end(), 0.0) / (double)cell.points.size();
		cell.z = accumulate(Vz.begin(), Vz.end(), 0.0) / (double)cell.points.size();
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
	
	
	
	
	
	
	

	
	
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			double ownVol = mesh.cells[face.owner].volume;
			double ngbVol = mesh.cells[face.neighbour].volume;
			
			// face.wC = ngbVol/(ownVol+ngbVol);
			
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
			face.wC = 0.5;
		}
		
		// face.wC = 0.5;
	}
	if(size>1){
		
		vector<double> sendValues0;
		vector<double> sendValues1;
		vector<double> sendValues2;
		vector<double> recvValues0(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> recvValues1(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
		vector<double> recvValues2(mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1],0.0);
			
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				sendValues0.push_back(mesh.cells[face.owner].x);
				sendValues1.push_back(mesh.cells[face.owner].y);
				sendValues2.push_back(mesh.cells[face.owner].z);
			}
		}
		
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
					   
		vector<double>::iterator iter0 = recvValues0.begin();
		vector<double>::iterator iter1 = recvValues1.begin();
		vector<double>::iterator iter2 = recvValues2.begin();
		for(auto& face : mesh.faces){
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// double ownVol = mesh.cells[face.owner].volume;
				// double ngbVol = *iter;
				
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
		// face.wC = 0.5;
	// }
	
	
	
	this->setStencil(mesh);
	
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
}












void SEMO_Mesh_Geometric::setStencil(SEMO_Mesh_Builder& mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	
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