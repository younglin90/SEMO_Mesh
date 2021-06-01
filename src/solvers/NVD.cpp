#include "build.h"
#include <cmath>
#include <array>



void SEMO_Solvers_Builder::calcVanLeer(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			// tildeCd = rf (at CFD book)
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			gamma_f = (tildeCd+abs(tildeCd))/(1.0+abs(tildeCd));
			double valueL = alpD + 0.5*gamma_f*(alpA-alpD);
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			gamma_f = (tildeCd+abs(tildeCd))/(1.0+abs(tildeCd));
			double valueR = alpD + 0.5*gamma_f*(alpA-alpD);
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}





void SEMO_Solvers_Builder::calcQUICK(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	
	
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			// double minAlp = 0.0;
			// double maxAlp = 1.0;
			
			// L2 = max( minAlp, min( maxAlp, L2 ));
			// R2 = max( minAlp, min( maxAlp, R2 ));

			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			// QUICK
			if(tildeCd>=0.0 && tildeCd<1.0){
				double tildeCf = 0.375 + 0.75*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<1.0){
				double tildeCf = 0.375 + 0.75*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
			// if(isnan(valueL) || isnan(valueR)){
				// cerr << "| #Error : at QUICK scheme" << endl;
				// cout << valueL << " " << valueR  << endl;
				// cout << L2 << " " << L1  << " " << R1 << " " << R2 << endl;
				// cout << tildeCd << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
			
			
		}
	}
}



void SEMO_Solvers_Builder::calcMINMOD(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<1.0){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<1.0){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}


void SEMO_Solvers_Builder::calcBoundedCD(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<1.0){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<1.0){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}




void SEMO_Solvers_Builder::calcOSHER(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}




void SEMO_Solvers_Builder::calcSMART(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.83333){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.83333 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<0.83333){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.83333 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}





void SEMO_Solvers_Builder::calcModifiedSMART(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.166666){
				double tildeCf = 3.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.166666 && tildeCd<0.7){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.7 && tildeCd<1.0){
				double tildeCf = 0.3333*tildeCd + 0.6666;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			if(tildeCd>=0.0 && tildeCd<0.166666){
				double tildeCf = 3.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.166666 && tildeCd<0.7){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.7 && tildeCd<1.0){
				double tildeCf = 0.3333*tildeCd + 0.6666;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}






void SEMO_Solvers_Builder::calcSTOIC(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.83333){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.83333 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.83333){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.83333 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}






void SEMO_Solvers_Builder::calcModifiedSTOIC(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.2){
				double tildeCf = 3.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.2 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.7){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.7 && tildeCd<1.0){
				double tildeCf = 0.3333*tildeCd + 0.6666;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.2){
				double tildeCf = 3.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.2 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.7){
				double tildeCf = 0.75*tildeCd + 0.375;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.7 && tildeCd<1.0){
				double tildeCf = 0.3333*tildeCd + 0.6666;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}







void SEMO_Solvers_Builder::calcMUSCL(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.25){
				double tildeCf = 2.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.25 && tildeCd<0.75){
				double tildeCf = tildeCd + 0.25;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.75 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.25){
				double tildeCf = 2.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.25 && tildeCd<0.75){
				double tildeCf = tildeCd + 0.25;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.75 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}









void SEMO_Solvers_Builder::calcSUPERBEE(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}








void SEMO_Solvers_Builder::calcModifiedSUPERBEE(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double gamma_f;

			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.3333){
				double tildeCf = 2.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.3333 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueL = gamma_f * alpA + (1.0-gamma_f) * alpD;
			


			// right value
			alpU = R2; 
			alpD = R1; 
			alpA = L1;
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			if(tildeCd>=0.0 && tildeCd<0.3333){
				double tildeCf = 2.0*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.3333 && tildeCd<0.5){
				double tildeCf = 0.5*tildeCd + 0.5;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.6666){
				double tildeCf = 1.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.6666 && tildeCd<1.0){
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double valueR = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// set velues
			face.varL[fn] = valueL;
			face.varR[fn] = valueR;
			
		}
	}
}







// 
void SEMO_Solvers_Builder::calcMSTACS(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	
	
	vector<double> corantNo;

	for(auto& cell : mesh.cells){
		
		double magVel = 
			sqrt(cell.var[controls.U]*cell.var[controls.U] +
			     cell.var[controls.V]*cell.var[controls.V] +
				 cell.var[controls.W]*cell.var[controls.W]);
		
		// double C = cell.var[controls.C];
		
		// double co = pow(cell.volume,0.333) / ( magVel + abs(C) );
		
		double maxA=0.0;
		for(auto& i : cell.faces){
			maxA = max(maxA,mesh.faces[i].area);
		}
		
		double co = magVel * controls.timeStep * maxA / cell.volume;
		
		corantNo.push_back(co);
	}
	
	
	for(auto& face : mesh.faces){

		// SEMO_Cell& own = mesh.cells[face.owner];
		// SEMO_Cell& ngb = mesh.cells[face.neighbour];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			double L1 = mesh.cells[face.owner].var[cn];
			double R1 = mesh.cells[face.neighbour].var[cn];

			double L2=0.0;
			for(int i=0; i<3; ++i){
				L2 += inpDX[face.owner][i]*face.distCells[i];
			}
			L2 = L1 - L2;
			
			double R2=0.0;
			for(int i=0; i<3; ++i){
				R2 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			R2 = R1 + R2;
			
			double minAlp = 0.0;
			double maxAlp = 1.0;
			
			L2 = max( minAlp, min( maxAlp, L2 ));
			R2 = max( minAlp, min( maxAlp, R2 ));



			// left value
			double alpU = L2; 
			double alpD = L1; 
			double alpA = R1; 
			double coDD = corantNo[face.owner];
			
			double cosTheta = 0.0;
			double tmp1 = 0.0;
			for(int i=0; i<3; ++i){
				cosTheta += pow(inpDX[face.owner][i],2.0);
				tmp1 += pow(face.distCells[i],2.0);
			}
			cosTheta = sqrt(cosTheta)*sqrt(tmp1);
			
			tmp1 = 0.0;
			for(int i=0; i<3; ++i){
				tmp1 += inpDX[face.owner][i]*face.distCells[i];
			}
			if(cosTheta>0.0) cosTheta = abs(tmp1) / cosTheta;
			
			double gamF = min(pow(cosTheta,4.0),1.0);
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			double gamma_f;
			

			// CDS-MSTACS
			if(tildeCd>=0.0 && tildeCd<1.0 && coDD>0.0 && coDD<=0.33){
				double tildeCf = min(1.0,tildeCd/coDD);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.0 && tildeCd<1.0 && coDD>0.33 && coDD<=1.0) {
				double tildeCf = min(1.0,3.0*tildeCd);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double vfCompressive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			
			// HR-STOIC
			if(tildeCd>=0.0 && tildeCd<0.2) {
				double tildeCf = min(1.0,3.0*tildeCd);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.2 && tildeCd<0.5) {
				double tildeCf = 0.5 + 0.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.8333) {
				double tildeCf = 0.375 + 0.75*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.8333 && tildeCd<1.0) {
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			double vfDiffusive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			double vf = gamF*vfCompressive + (1.0-gamF)*vfDiffusive;
			
			double vfL = max( 0.0, min( 1.0, vf ));
			


			// right value
			alpU = R2; alpD = R1; alpA = L1; coDD = corantNo[face.neighbour];

			cosTheta = 0.0;
			tmp1 = 0.0;
			for(int i=0; i<3; ++i){
				cosTheta += pow(inpDX[face.neighbour][i],2.0);
				tmp1 += pow(face.distCells[i],2.0);
			}
			cosTheta = sqrt(cosTheta)*sqrt(tmp1);
			
			tmp1 = 0.0;
			for(int i=0; i<3; ++i){
				tmp1 += inpDX[face.neighbour][i]*face.distCells[i];
			}
			if(cosTheta>0.0) cosTheta = abs(tmp1) / cosTheta;
			
			gamF = min(pow(cosTheta,4.0),1.0);
			
			tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			// CDS-MSTACS
			if(tildeCd>=0.0 && tildeCd<1.0 && coDD>0.0 && coDD<=0.33){
				double tildeCf = min(1.0,tildeCd/coDD);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.0 && tildeCd<1.0 && coDD>0.33 && coDD<=1.0) {
				double tildeCf = min(1.0,3.0*tildeCd);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			vfCompressive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			// HR-STOIC
			if(tildeCd>=0.0 && tildeCd<0.2) {
				double tildeCf = min(1.0,3.0*tildeCd);
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.2 && tildeCd<0.5) {
				double tildeCf = 0.5 + 0.5*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.5 && tildeCd<0.8333) {
				double tildeCf = 0.375 + 0.75*tildeCd;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else if(tildeCd>=0.8333 && tildeCd<1.0) {
				double tildeCf = 1.0;
				gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			}
			else{
				gamma_f = 0.0;
			}
			vfDiffusive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			vf = gamF*vfCompressive + (1.0-gamF)*vfDiffusive;
			
			double vfR = max( 0.0, min( 1.0, vf ));
			
			
			
			// set velues
			face.varL[fn] = vfL;
			face.varR[fn] = vfR;
			
			// if(vfL>0) cout << vfL << " " << vfR << " " << L1 << " " << R1 << endl;
			
			
			
		}
		
	
	
	}
	
	
			
			

			
			
			
			
			
}