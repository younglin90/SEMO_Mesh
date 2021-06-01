#include "build.h"
#include <cmath>
#include <array>


void SEMO_Solvers_Builder::calcViscousFlux(
	double& muL, double& UL, double& VL, double& WL, double& RhoL, 
	double& muR, double& UR, double& VR, double& WR, double& RhoR, 
	double& dUdxF, double& dUdyF, double& dUdzF,
	double& dVdxF, double& dVdyF, double& dVdzF,
	double& dWdxF, double& dWdyF, double& dWdzF,
	vector<double>& nvec,
	vector<double>& flux
	){
	
	double nx = nvec[0];
	double ny = nvec[1];
	double nz = nvec[2];
	
	// double dPN = face.distCells[0]*nx+face.distCells[1]*ny+face.distCells[2]*nz;
	// double dPF = (face.x-cell.x)*nx+(face.y-cell.y)*ny+(face.z-cell.z)*nz;
	
	double gc = 0.5; //(dPN-dPF)/dPN;
	
	double muF = gc*muL + (1.0-gc)*muR;
	double rhoF = gc*RhoL + (1.0-gc)*RhoR;
	double UF = gc*UL + (1.0-gc)*UR;
	double VF = gc*VL + (1.0-gc)*VR;
	double WF = gc*WL + (1.0-gc)*WR;
	
	// Ref : Blazek's book, pp.20-21
    double txx = 4.0/3.0 * dUdxF - 2.0/3.0 * (dVdyF + dWdzF);
    double txy = dVdxF + dUdyF;
    double txz = dWdxF + dUdzF;
    double tyx = txy;
    double tyy = 4.0/3.0 * dVdyF - 2.0/3.0 * (dUdxF + dWdzF);
    double tyz = dVdzF + dWdyF;
    double tzx = txz;
    double tzy = tyz;
    double tzz = 4.0/3.0 * dWdzF - 2.0/3.0 * (dUdxF + dVdyF);
			
	// flux.clear();
	// flux.resize(controls.nEq,0.0);
	
	flux[0] = 0.0;
	flux[1] = muF * (txx*nx + txy*ny + txz*nz);
	// flux[1] -= nx * 2.0/3.0 * rhoF * tkei;
	flux[2] = muF * (tyx*nx + tyy*ny + tyz*nz);
	// flux[2] -= ny * 2.0/3.0 * rhoF * tkei;
	flux[3] = muF * (tzx*nx + tzy*ny + tzz*nz);
	// flux[3] -= nz * 2.0/3.0 * rhoF * tkei;
	flux[4] = muF * (txx*UF + txy*VF + txz*WF) * nx
            + muF * (tyx*UF + tyy*VF + tyz*WF) * ny
            + muF * (tzx*UF + tzy*VF + tzz*WF) * nz;
			
	// heat transfer
	// flux.back() += cnti * (dTdxFace(i,1)*nx + dTdxFace(i,2)*ny + dTdxFace(i,3)*nz);
	
    // turbulent terms
	// flux.back() -= 2.0/3.0*prim(1)*tkei*(prim(2)*nx + prim(3)*ny + prim(4)*nz);
	
    // TODO consider more terms
	// flux.back() += mutprti * (dHtdxFace(i,1)*nx + dHtdxFace(i,2)*ny + dHtdxFace(i,3)*nz);
	

}

