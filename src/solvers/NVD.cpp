#include "build.h"
#include <cmath>
#include <array>



void SEMO_Solvers_Builder::calcVanLeer(
	SEMO_Mesh_Builder& mesh,
	SEMO_Controls_Builder& controls,
	int cn, int fn,
	vector<vector<double>>& inpDX){
	
	for(auto& face : mesh.faces){

		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			// double L1 = mesh.cells[face.owner].var[cn];
			// double R1 = mesh.cells[face.neighbour].var[cn];

			// double L2=0.0;
			// for(int i=0; i<3; ++i){
				// L2 += inpDX[face.owner][i]*face.distCells[i];
			// }
			// L2 = L1 - L2;
			
			// double R2=0.0;
			// for(int i=0; i<3; ++i){
				// R2 += inpDX[face.neighbour][i]*face.distCells[i];
			// }
			// R2 = R1 + R2;
			
			// double gamma_f;

			// // left value
			// double alpU = L2; 
			// double alpD = L1; 
			// double alpA = R1; 
			
			// // tildeCd = rf (at CFD book)
			// double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			
			// gamma_f = (tildeCd+abs(tildeCd))/(1.0+abs(tildeCd));
			// double valueL = alpD + 0.5*gamma_f*(alpA-alpD);
			


			// // right value
			// alpU = R2; 
			// alpD = R1; 
			// alpA = L1;
			
			// tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );

			// gamma_f = (tildeCd+abs(tildeCd))/(1.0+abs(tildeCd));
			// double valueR = alpD + 0.5*gamma_f*(alpA-alpD);
			
			// // set velues
			// face.varL[fn] = valueL;
			// face.varR[fn] = valueR;
			
		// }
		

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
				double tildeCf = tildeCd + (tildeCd*(1.0-tildeCd))
								/(2.0-4.0*tildeCd+4.0*tildeCd*tildeCd);
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
				double tildeCf = tildeCd + (tildeCd*(1.0-tildeCd))
								/(2.0-4.0*tildeCd+4.0*tildeCd*tildeCd);
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
	vector<vector<double>>& inpDX,
	int fn_out){
	
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();
	
	SEMO_MPI_Builder mpi;
	
	vector<double> corantNo;
	double tmpCFL = 1.0;
	double eps = 1.e-8; 

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
		
		double co = 1.0/tmpCFL * magVel * controls.timeStep * maxA / cell.volume;
		
		corantNo.push_back(co);
	}
	
	
	
	vector<double> corantNo_recv;
	if(size>1){
		// processor faces
		vector<double> corantNo_send;
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				corantNo_send.push_back(corantNo[face.owner]);
			}
		}
		// SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					corantNo_send, corantNo_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		corantNo_send.clear();
	}
	
	
	int proc_num = 0;
	for(auto& face : mesh.faces){

		if(face.getType() == SEMO_Types::BOUNDARY_FACE) continue;
		
		int cellsize = mesh.cells.size();
		int procN = cellsize+proc_num;
		
		double wCL = face.wC;
		double wCR = 1.0-wCL;
		
		double L1 = mesh.cells[face.owner].var[cn];
		double R1 = 0.0;
		
		vector<double> inpDXL(3,0.0);
		vector<double> inpDXR(3,0.0);
		double coDDL=0.0;
		double coDDR=0.0;
		for(int i=0; i<3; ++i){
			inpDXL[i] = inpDX[face.owner][i];
		}
		coDDL = corantNo[face.owner];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			R1 = mesh.cells[face.neighbour].var[cn];
			for(int i=0; i<3; ++i){
				inpDXR[i] = inpDX[face.neighbour][i];
			}
			coDDR = corantNo[face.neighbour];
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			R1 = face.varR[fn];
			for(int i=0; i<3; ++i){
				inpDXR[i] = inpDX[procN][i];
			}
			coDDR = corantNo_recv[proc_num];
		}
		
		
		double L2=0.0;
		for(int i=0; i<3; ++i){
			L2 += inpDXL[i]*face.distCells[i];
		}
		L2 = L1 - L2;
		
		double R2=0.0;
		for(int i=0; i<3; ++i){
			R2 += inpDXR[i]*face.distCells[i];
		}
		R2 = R1 + R2;
		
		double minAlp = 0.0;
		double maxAlp = 1.0;
		
		L2 = max( minAlp, min( maxAlp, L2 ));
		R2 = max( minAlp, min( maxAlp, R2 ));

		double vfL = L1;
		double vfR = R1;

		for(int LR=0; LR<2; ++LR){
			
			if(LR==0 && ( vfL > 1.0-eps || vfL < eps )) continue;
			if(LR==1 && ( vfR > 1.0-eps || vfR < eps )) continue;
			
			double alpU=0;
			double alpD=0;
			double alpA=0;
			double coDD=0;
			double cosTheta = 0.0;
			double tmp1 = 0.0;
			if(LR==0) {
				// left value
				alpU = L2; alpD = L1; alpA = R1; 
				coDD = coDDL;
				
				for(int i=0; i<3; ++i){
					cosTheta += pow(inpDXL[i],2.0);
					tmp1 += pow(face.distCells[i],2.0);
				}
				cosTheta = sqrt(cosTheta)*sqrt(tmp1);
				
				tmp1 = 0.0;
				for(int i=0; i<3; ++i){
					tmp1 += inpDXL[i]*face.distCells[i];
				}
			}
			else{
				// right value
				alpU = R2; alpD = R1; alpA = L1; 
				coDD = coDDR;
				
				for(int i=0; i<3; ++i){
					cosTheta += pow(inpDXR[i],2.0);
					tmp1 += pow(face.distCells[i],2.0);
				}
				cosTheta = sqrt(cosTheta)*sqrt(tmp1);
				
				tmp1 = 0.0;
				for(int i=0; i<3; ++i){
					tmp1 += inpDXR[i]*face.distCells[i];
				}
				
			}
			
			// double gamma_f_skewness = 0.0;
			// gamma_f_skewness += wCL*
				// (inpDX[face.owner][0]*face.vecSkewness[0]+
				 // inpDX[face.owner][1]*face.vecSkewness[1]+
				 // inpDX[face.owner][2]*face.vecSkewness[2]) / (alpA-alpD+1.e-200);

			
			// if(face.getType() == SEMO_Types::INTERNAL_FACE){
				// gamma_f_skewness += wCR*
					// (inpDX[face.neighbour][0]*face.vecSkewness[0]+
					 // inpDX[face.neighbour][1]*face.vecSkewness[1]+
					 // inpDX[face.neighbour][2]*face.vecSkewness[2]) / (alpA-alpD+1.e-200);
			// }
			// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// gamma_f_skewness += wCR*
					// (inpDX[procN][0]*face.vecSkewness[0]+
					 // inpDX[procN][1]*face.vecSkewness[1]+
					 // inpDX[procN][2]*face.vecSkewness[2]) / (alpA-alpD+1.e-200);
			// }
			
			
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
			
			// // skewness
			// gamma_f += gamma_f_skewness;
			// if(gamma_f>1.0) gamma_f=1.0;
			// if(gamma_f<0.0) gamma_f=0.0;
			
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
			
			// // skewness
			// gamma_f += gamma_f_skewness;
			// if(gamma_f>1.0) gamma_f=1.0;
			// if(gamma_f<0.0) gamma_f=0.0;
			
			double vfDiffusive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			
			double vf = gamF*vfCompressive + (1.0-gamF)*vfDiffusive;
			
			if(LR==0) vfL = max( 0.0, min( 1.0, vf ));
			if(LR==1) vfR = max( 0.0, min( 1.0, vf ));
			
		}

		// set velues
		face.varL[fn_out] = vfL;
		face.varR[fn_out] = vfR;
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
	}
	
	
			
			

			
			
			
			
			
}