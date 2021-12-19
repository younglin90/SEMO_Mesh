
#include "../build.h"
#include <cmath>
#include <array>
#include <numeric>

double SEMO_Solvers_Builder::calcSourceGravity(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	vector<SEMO_Species>& species
	){
		
	// ===============================
	// 중력 항 저장
	for(int i=0; i<mesh.cells.size(); ++i){
		mesh.cells[i].var[controls.sourceGravity[0]] = 0.0;
		mesh.cells[i].var[controls.sourceGravity[1]] = 0.0;
		mesh.cells[i].var[controls.sourceGravity[2]] = 0.0;
		mesh.cells[i].var[controls.sourceGravity[3]] = 0.0;
	}
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		double wCL = face.wC; double wCR = 1.0-wCL;
		double UL = face.varL[controls.fU];
		double UR = face.varR[controls.fU];
		double VL = face.varL[controls.fV];
		double VR = face.varR[controls.fV];
		double WL = face.varL[controls.fW];
		double WR = face.varR[controls.fW];
		double UnL = UL*face.unitNormals[0] + 
					 VL*face.unitNormals[1] + 
					 WL*face.unitNormals[2];
		double UnR = UR*face.unitNormals[0] + 
					 VR*face.unitNormals[1] + 
					 WR*face.unitNormals[2];
		double UnF = wCL * UnL + wCR * UnR;
		double RhoL = face.varL[controls.fRho];
		double RhoR = face.varR[controls.fRho];
		
		// interfacial gravity force
		// harmonic 평균
		double RhoF = 1.0/(wCL/RhoL + wCR/RhoR);
		double gAlphaL = 0.0; double gAlphaR = 0.0;
		for(int ii=0; ii<3; ++ii){
			gAlphaL += controls.gravityAcceleration[ii]*face.vecPF[ii];
			gAlphaR += controls.gravityAcceleration[ii]*face.vecNF[ii];
		}
		
        // ----------------------------
		for(int ii=0; ii<3; ++ii){
			mesh.cells[face.owner].var[controls.sourceGravity[ii]] += 
				( (RhoF * gAlphaL * face.unitNormals[ii]) * face.area / 
				mesh.cells[face.owner].volume );
		}
		mesh.cells[face.owner].var[controls.sourceGravity[3]] += 
			( (RhoF * gAlphaL * UnF) * face.area / 
			mesh.cells[face.owner].volume );
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			for(int ii=0; ii<3; ++ii){
				mesh.cells[face.neighbour].var[controls.sourceGravity[ii]] -= 
					( (RhoF * gAlphaR * face.unitNormals[ii]) * face.area / 
					mesh.cells[face.neighbour].volume );
			}
			mesh.cells[face.neighbour].var[controls.sourceGravity[3]] -= 
				( (RhoF * gAlphaR * UnF) * face.area / 
				mesh.cells[face.neighbour].volume );
		}
	}
	// ===============================
	
	return 0.0;
	
}