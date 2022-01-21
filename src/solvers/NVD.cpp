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
	
	// vector<double> corantNo;
	double tmpCFL = 1.0;
	double eps = 1.e-12; 

	vector<double> corantNo(mesh.cells.size(),-1.e80);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		
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
		
		corantNo[i] = co;
	}
	
	// for(auto& face : mesh.faces){
	
		// double normalMagVel = mesh.cells[face.owner].var[controls.U]*face.unitNormals[0];
		// normalMagVel += mesh.cells[face.owner].var[controls.V]*face.unitNormals[1];
		// normalMagVel += mesh.cells[face.owner].var[controls.W]*face.unitNormals[2];
		// normalMagVel = abs(normalMagVel);
		
		// double normalFluxOwn = normalMagVel*face.area;
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// normalMagVel = mesh.cells[face.neighbour].var[controls.U]*face.unitNormals[0];
			// normalMagVel += mesh.cells[face.neighbour].var[controls.V]*face.unitNormals[1];
			// normalMagVel += mesh.cells[face.neighbour].var[controls.W]*face.unitNormals[2];
			// normalMagVel = abs(normalMagVel);
			// double normalFluxNgb = normalMagVel*face.area;
			
			// double normalFlux = max(normalFluxOwn,normalFluxNgb);
			
			// double corantOwn = controls.timeStep / ( mesh.cells[face.owner].volume / (normalFlux + 1.e-100) );
			// double corantNgb = controls.timeStep / ( mesh.cells[face.neighbour].volume / (normalFlux + 1.e-100) );
			
			// corantNo[face.owner] = max(corantNo[face.owner], corantOwn );
			// corantNo[face.neighbour] = max(corantNo[face.neighbour], corantNgb );
		// }
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// normalMagVel = face.varR[controls.fU]*face.unitNormals[0];
			// normalMagVel += face.varR[controls.fV]*face.unitNormals[1];
			// normalMagVel += face.varR[controls.fW]*face.unitNormals[2];
			// normalMagVel = abs(normalMagVel);
			// double normalFluxNgb = normalMagVel*face.area;
			
			// double normalFlux = max(normalFluxOwn,normalFluxNgb);
			
			// double corantOwn = controls.timeStep / ( mesh.cells[face.owner].volume / (normalFlux + 1.e-100) );
			
			// corantNo[face.owner] = max(corantNo[face.owner], corantOwn );
			
			// ++proc_num;
		// }
		// else{
			// double corantOwn = controls.timeStep / ( mesh.cells[face.owner].volume / (normalFluxOwn + 1.e-100) );
			
			// corantNo[face.owner] = max(corantNo[face.owner], corantOwn );
		// }
	// }
	
	
		// cout << "AAAAAAAAA" << endl;
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
	
	
		// // cout << "BBBBBBBBB" << endl;
	// vector<double> maxPhi(mesh.cells.size(),-1.e15);
	// vector<double> minPhi(mesh.cells.size(),+1.e15);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// maxPhi[i] = cell.var[cn];
		// minPhi[i] = cell.var[cn];
		// for(auto j : cell.stencil){
			// auto& cellSten = mesh.cells[j];
			// maxPhi[i] = max(maxPhi[i],cellSten.var[cn]);
			// minPhi[i] = min(minPhi[i],cellSten.var[cn]);
		// }
	// }
	// if(size>1){
		// // processor faces
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// maxPhi[face.owner] = max(maxPhi[face.owner],face.varR[fn]);
				// minPhi[face.owner] = min(minPhi[face.owner],face.varR[fn]);
			// }
		// }
	// }
	// vector<double> maxPhi_recv;
	// vector<double> minPhi_recv;
	// if(size>1){
		// // processor faces
		// vector<double> maxPhi_send;
		// vector<double> minPhi_send;
		// for(int i=0; i<mesh.faces.size(); ++i){
			// auto& face = mesh.faces[i];
			
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// maxPhi_send.push_back(maxPhi[face.owner]);
				// minPhi_send.push_back(minPhi[face.owner]);
			// }
		// }
		
		// mpi.setProcsFaceDatas(
					// maxPhi_send, maxPhi_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// mpi.setProcsFaceDatas(
					// minPhi_send, minPhi_recv,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
		// maxPhi_send.clear();
		// minPhi_send.clear();
	// }
	
	
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
		// double maxPhiL = maxPhi[face.owner];
		// double maxPhiR = 0.0;
		// double minPhiL = minPhi[face.owner];
		// double minPhiR = 0.0;
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
			// maxPhiR = maxPhi[face.neighbour];
			// minPhiR = minPhi[face.neighbour];
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			R1 = face.varR[fn];
			for(int i=0; i<3; ++i){
				inpDXR[i] = inpDX[procN][i];
			}
			coDDR = corantNo_recv[proc_num];
			// maxPhiR = maxPhi_recv[proc_num];
			// minPhiR = minPhi_recv[proc_num];
		}
		
		
		// // skewness
		// for(int i=0; i<3; ++i){
			// L1 += inpDXL[i]*face.vecSkewness[i];
			// R1 += inpDXR[i]*face.vecSkewness[i];
		// }
		
		
		double L2=0.0;
		for(int i=0; i<3; ++i){
			// L2 += inpDXL[i]*face.distCells[i];
			// L2 += (wCL*inpDXL[i]+wCR*inpDXR[i])*face.distCells[i];
			// L2 += inpDXL[i]*face.distCells[i];
			L2 += inpDXL[i]*face.vecPN[i];
		}
		// L2 = L1 - L2;
		L2 = R1 - 2.0*L2;
		
		double R2=0.0;
		for(int i=0; i<3; ++i){
			// R2 += inpDXR[i]*face.distCells[i];
			// R2 += (wCL*inpDXL[i]+wCR*inpDXR[i])*face.distCells[i];
			// R2 += inpDXR[i]*face.distCells[i];
			R2 += inpDXR[i]*face.vecPN[i];
		}
		// R2 = R1 + R2;
		R2 = L1 + 2.0*R2;
		
		L2 = max( 0.0, min( 1.0, L2 ));
		R2 = max( 0.0, min( 1.0, R2 ));
		
		// L2 = max( min(minPhiL,minPhiR), min( max(maxPhiL,maxPhiR), L2 ));
		// R2 = max( min(minPhiL,minPhiR), min( max(maxPhiL,maxPhiR), R2 ));

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
			double xf_tilde = 0.0;
			double xD_tilde = 0.0;
			double wCLR = 0.0;
			if(LR==0) {
				// left value
				alpU = L2; alpD = L1; alpA = R1; 
				// coDD = coDDL;
				coDD = max(coDDL,coDDR);
				
				// for(int i=0; i<3; ++i){
					// cosTheta += pow(inpDXL[i],2.0);
					// tmp1 += pow(face.distCells[i],2.0);
				// }
				// cosTheta = sqrt(cosTheta)*sqrt(tmp1);
				
				// tmp1 = 0.0;
				// for(int i=0; i<3; ++i){
					// tmp1 += inpDXL[i]*face.distCells[i];
				// }
				
				wCLR = 1.0;//1.0/(2.0*wCR);
			}
			else{
				// right value
				alpU = R2; alpD = R1; alpA = L1; 
				// coDD = coDDR;
				coDD = max(coDDL,coDDR);
				
				// for(int i=0; i<3; ++i){
					// cosTheta += pow(inpDXR[i],2.0);
					// tmp1 += pow(face.distCells[i],2.0);
				// }
				// cosTheta = sqrt(cosTheta)*sqrt(tmp1);
				
				// tmp1 = 0.0;
				// for(int i=0; i<3; ++i){
					// tmp1 += inpDXR[i]*face.distCells[i];
				// }
				
				wCLR = 1.0;//1.0/(2.0*wCL);
			}
			
			

			double mfL[3];
			mfL[0] = inpDXL[0];
			mfL[1] = inpDXL[1];
			mfL[2] = inpDXL[2];
			double magMfL = mfL[0]*mfL[0];
			magMfL += mfL[1]*mfL[1];
			magMfL += mfL[2]*mfL[2];
			magMfL = sqrt(magMfL);
			mfL[0] = mfL[0]/(magMfL+1.e-200);
			mfL[1] = mfL[1]/(magMfL+1.e-200);
			mfL[2] = mfL[2]/(magMfL+1.e-200);
			
			double mfR[3];
			mfR[0] = inpDXR[0];
			mfR[1] = inpDXR[1];
			mfR[2] = inpDXR[2];
			double magMfR = mfR[0]*mfR[0];
			magMfR += mfR[1]*mfR[1];
			magMfR += mfR[2]*mfR[2];
			magMfR = sqrt(magMfR);
			mfR[0] = mfR[0]/(magMfR+1.e-200);
			mfR[1] = mfR[1]/(magMfR+1.e-200);
			mfR[2] = mfR[2]/(magMfR+1.e-200);
			
			double mfLR[3];
			mfLR[0] = wCL*inpDXL[0]+wCR*inpDXR[0];
			mfLR[1] = wCL*inpDXL[1]+wCR*inpDXR[1];
			mfLR[2] = wCL*inpDXL[2]+wCR*inpDXR[2];
			double magMfLR = mfLR[0]*mfLR[0];
			magMfLR += mfLR[1]*mfLR[1];
			magMfLR += mfLR[2]*mfLR[2];
			magMfLR = sqrt(magMfLR);
			mfLR[0] = mfLR[0]/(magMfLR+1.e-200);
			mfLR[1] = mfLR[1]/(magMfLR+1.e-200);
			mfLR[2] = mfLR[2]/(magMfLR+1.e-200);
			
			// cosTheta = 0.5*(mfL[0]+mfR[0])*face.unitNomalsPN[0];
			// cosTheta += 0.5*(mfL[1]+mfR[1])*face.unitNomalsPN[1];
			// cosTheta += 0.5*(mfL[2]+mfR[2])*face.unitNomalsPN[2];
			// cosTheta = (wCL*mfL[0]+wCR*mfR[0])*face.unitNormals[0];
			// cosTheta += (wCL*mfL[1]+wCR*mfR[1])*face.unitNormals[1];
			// cosTheta += (wCL*mfL[2]+wCR*mfR[2])*face.unitNormals[2];
			cosTheta = mfLR[0]*face.unitNormals[0];
			cosTheta += mfLR[1]*face.unitNormals[1];
			cosTheta += mfLR[2]*face.unitNormals[2];
			// if(LR==0){
				// // cosTheta = mfL[0]*face.unitNormals[0];
				// // cosTheta += mfL[1]*face.unitNormals[1];
				// // cosTheta += mfL[2]*face.unitNormals[2];
				// cosTheta = mfL[0]*face.unitNomalsPN[0];
				// cosTheta += mfL[1]*face.unitNomalsPN[1];
				// cosTheta += mfL[2]*face.unitNomalsPN[2];
			// }
			// else{
				// // cosTheta = mfR[0]*face.unitNormals[0];
				// // cosTheta += mfR[1]*face.unitNormals[1];
				// // cosTheta += mfR[2]*face.unitNormals[2];
				// cosTheta = mfR[0]*face.unitNomalsPN[0];
				// cosTheta += mfR[1]*face.unitNomalsPN[1];
				// cosTheta += mfR[2]*face.unitNomalsPN[2];
			// }
			cosTheta = abs(cosTheta);
			
			
			
			
			
			// skewness correction
			double gamma_f_skewness = 0.0;
			for(int i=0; i<3; ++i){
				gamma_f_skewness += wCL * inpDXL[i]*face.vecSkewness[i] /
						(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
				gamma_f_skewness += wCR * inpDXR[i]*face.vecSkewness[i] /
						(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			}
			// if(LR==0) gamma_f_skewness *= 1.0/wCR;
			// if(LR==1) gamma_f_skewness *= 1.0/wCL;
			// gamma_f_skewness *= 1.0/wCLR;
			
			
			// if(cosTheta>0.0) cosTheta = abs(tmp1) / cosTheta;
			
			double gamF = min(pow(cosTheta,4.0),1.0);
			// double gamF = 1.0;
			
			double tildeCd = (alpD-alpU)/(abs(alpA-alpU)+1.e-200)*( alpA>alpU ? 1.0 : -1.0 );
			double tildeCr = (alpD-alpU)/(abs(alpA-alpD)+1.e-200)*( alpA>alpD ? 1.0 : -1.0 );
			
			double gamma_f;
			




			// // orthogonal
			// double cosNonOrtho = 0.0;
			// cosNonOrtho += face.unitNomalsPN[0]*face.unitNormals[0];
			// cosNonOrtho += face.unitNomalsPN[1]*face.unitNormals[1];
			// cosNonOrtho += face.unitNomalsPN[2]*face.unitNormals[2];
			// double gamF_nonOrtho = min(pow(cosNonOrtho,2.0),1.0);
			// gamF = min(gamF,gamF_nonOrtho);





			// // SUPERBEE
			// double gamma_f_SUPERBEE;
			// double tildeCf_SUPERBEE;
			// if(tildeCd>=0.0 && tildeCd<0.3333){
				// tildeCf_SUPERBEE = 2.0*tildeCd;
				// gamma_f_SUPERBEE = (tildeCf_SUPERBEE - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.3333 && tildeCd<0.5){
				// tildeCf_SUPERBEE = 0.5*tildeCd + 0.5;
				// gamma_f_SUPERBEE = (tildeCf_SUPERBEE - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.5 && tildeCd<0.6666){
				// tildeCf_SUPERBEE = 1.5*tildeCd;
				// gamma_f_SUPERBEE = (tildeCf_SUPERBEE - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.6666 && tildeCd<1.0){
				// tildeCf_SUPERBEE = 1.0;
				// gamma_f_SUPERBEE = (tildeCf_SUPERBEE - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// tildeCf_SUPERBEE = tildeCd;
				// gamma_f_SUPERBEE = 0.0;
			// }

			// // CN-CBC
			// if(tildeCd>=0.0 && tildeCd<1.0){
				// if(coDD>0.0 && coDD<=0.3){
					// double tildeCf = min(1.0,tildeCd/coDD);
					// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
				// }
				// else if(coDD>0.3 && coDD<=0.6){
					// double tildeCf = min(1.0,tildeCd/0.3);
					// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
				// }
				// else if(coDD>0.6 && coDD<=0.7){
					// double tildeCf = (0.7-coDD)/0.1*min(1.0,tildeCd/0.3)+(coDD-0.6)/0.1*tildeCf_SUPERBEE;
					// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
				// }
				// else if(coDD>0.7 && coDD<=1.0){
					// gamma_f = gamma_f_SUPERBEE;
				// }
				// else{
					// gamma_f = 0.0;
				// }
			// }
			// else{
				// gamma_f = 0.0;
			// }


			// // Hyper-C
			// double gamma_f_Hyper_C = 0.0;
			// if(tildeCd>=0.0 && tildeCd<1.0){
				// double tildeCf = min(1.0,tildeCd/coDD);
				// gamma_f_Hyper_C = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f_Hyper_C = 0.0;
			// }
			// // ULTIMATE QUICKEST (UQ)
			// if(tildeCd>=0.0 && tildeCd<1.0){
				// double tildeCf = min((8.0*coDD*tildeCd+(1.0-coDD)*(6.0*tildeCd+3.0))/8.0,gamma_f_Hyper_C);
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }


			// // BD-FBICS
			// if(tildeCd>=0.0 && tildeCd<0.3333333333333*wCLR){
				// double tildeCf = min(1.0,3.0*tildeCd);
				// // tildeCf *= wCLR;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.3333333333333*wCLR && tildeCd<1.0*wCLR) {
				// double tildeCf = 1.0;
				// // tildeCf *= wCLR;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }


			// // BD-SAISH
			// if(tildeCd>=0.0 && tildeCd<0.25){
				// double tildeCf = min(1.0,4.0*tildeCd);
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.25 && tildeCd<1.0) {
				// double tildeCf = 1.0;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }
			
			

			// CDS-MSTACS
			if(tildeCd>=0.0 && tildeCd<1.0){
				if(coDD>0.0 && coDD<=0.33){
					double tildeCf = min(1.0,tildeCd/coDD);
					gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
				}
				else if(coDD>0.33 && coDD<=1.0){
					double tildeCf = min(1.0,3.0*tildeCd);
					gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
				}
				else{
					gamma_f = 0.0;
				}
			}
			else{
				gamma_f = 0.0;
			}
			
			
			
			
			
			// uniformity correction
			if(LR==0) gamma_f *= 2.0*wCR;
			if(LR==1) gamma_f *= 2.0*wCL;
			gamma_f = max(0.0,min(1.0,gamma_f));
			
			
			
			// // skewness correction
			// gamma_f += gamma_f_skewness;
			// // // TVD condition (1st order)
			// // double minTVD = max(0.0,min(tildeCd,0.5));
			// // double maxTVD = 0.0;
			// // if(LR==0) maxTVD = max(0.0,2.0*tildeCd/wCR);
			// // if(LR==1) maxTVD = max(0.0,2.0*tildeCd/wCL);
			// // gamma_f = min(maxTVD,max(minTVD,gamma_f));
			// // // gamma_f = max(0.0,min(gamma_f,min(tildeCr*(1.0-coDD)/coDD,1.0)));
			// gamma_f = max(0.0,min(1.0,gamma_f));
			
			
			
			// // uniformity correction
			// gamma_f *= 1.0/wCLR;
			
			
			// double vfCompressive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			double vfCompressive = alpD + gamma_f * (alpA - alpD);
			
			// vfCompressive += (1.0-gamma_f)*inpDXL[0]*face.vecPF[0];
			// vfCompressive += (1.0-gamma_f)*inpDXL[1]*face.vecPF[1];
			// vfCompressive += (1.0-gamma_f)*inpDXL[2]*face.vecPF[2];
			// vfCompressive += gamma_f*inpDXR[0]*face.vecNF[0];
			// vfCompressive += gamma_f*inpDXR[1]*face.vecNF[1];
			// vfCompressive += gamma_f*inpDXR[2]*face.vecNF[2];
			
			// vfCompressive += (1.0-gamma_f)*inpDXL[0]*face.vecSkewness[0];
			// vfCompressive += (1.0-gamma_f)*inpDXL[1]*face.vecSkewness[1];
			// vfCompressive += (1.0-gamma_f)*inpDXL[2]*face.vecSkewness[2];
			// vfCompressive += gamma_f*inpDXR[0]*face.vecSkewness[0];
			// vfCompressive += gamma_f*inpDXR[1]*face.vecSkewness[1];
			// vfCompressive += gamma_f*inpDXR[2]*face.vecSkewness[2];
			
			vfCompressive += wCL*inpDXL[0]*face.vecSkewness[0];
			vfCompressive += wCL*inpDXL[1]*face.vecSkewness[1];
			vfCompressive += wCL*inpDXL[2]*face.vecSkewness[2];
			vfCompressive += wCR*inpDXR[0]*face.vecSkewness[0];
			vfCompressive += wCR*inpDXR[1]*face.vecSkewness[1];
			vfCompressive += wCR*inpDXR[2]*face.vecSkewness[2];
			
			vfCompressive = max(0.0,min(1.0,vfCompressive));
			// vfCompressive = max(min(minPhiL,minPhiR),min(max(maxPhiL,maxPhiR),vfCompressive));
			
			double gamma_comp_f = gamma_f;


			// // MUSCL
			// if(tildeCd>=0.0 && tildeCd<0.25){
				// double tildeCf = 2.0*tildeCd;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.25 && tildeCd<0.75){
				// double tildeCf = tildeCd + 0.25;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.75 && tildeCd<1.0){
				// double tildeCf = 1.0;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }
			
			
			// // Hyper-C
			// if(tildeCd>=0.0 && tildeCd<1.0){
				// double tildeCf = min(1.0,tildeCd/coDD);
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }
			
			// // HR-FBICS
			// if(tildeCd>=0.0 && tildeCd<0.125) {
				// double tildeCf = min(1.0,3.0*tildeCd);
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.125 && tildeCd<0.75) {
				// double tildeCf = tildeCd + 0.25;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.75 && tildeCd<1.0) {
				// double tildeCf = 1.0;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }
			
			
			// // HR-SAISH
			// if(tildeCd>=0.0 && tildeCd<0.5) {
				// double tildeCf = min(1.0,tildeCd*(2.0-tildeCd));
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.5 && tildeCd<0.75) {
				// double tildeCf = (tildeCd + 0.25);
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else if(tildeCd>=0.75 && tildeCd<1.0) {
				// double tildeCf = 1.0;
				// gamma_f = (tildeCf - tildeCd) / (1.0 - tildeCd);
			// }
			// else{
				// gamma_f = 0.0;
			// }
			
			
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
			
			
			
			
			// uniformity correction
			if(LR==0) gamma_f *= 2.0*wCR;
			if(LR==1) gamma_f *= 2.0*wCL;
			gamma_f = max(0.0,min(1.0,gamma_f));
			
			
			
			// // skewness correction
			// gamma_f += gamma_f_skewness;
			// // // TVD condition (1st order)
			// // double minTVD = max(0.0,min(tildeCd,0.5));
			// // double maxTVD = 0.0;
			// // if(LR==0) maxTVD = max(0.0,2.0*tildeCd/wCR);
			// // if(LR==1) maxTVD = max(0.0,2.0*tildeCd/wCL);
			// // gamma_f = min(maxTVD,max(minTVD,gamma_f));
			// // // gamma_f = max(0.0,min(gamma_f,min(tildeCr*(1.0-coDD)/coDD,1.0)));
			// gamma_f = max(0.0,min(1.0,gamma_f));
			
			
			
			
			// double vfDiffusive = gamma_f * alpA + (1.0-gamma_f) * alpD;
			double vfDiffusive = alpD + gamma_f * (alpA - alpD);
			
			// vfDiffusive += (1.0-gamma_f)*inpDXL[0]*face.vecPF[0];
			// vfDiffusive += (1.0-gamma_f)*inpDXL[1]*face.vecPF[1];
			// vfDiffusive += (1.0-gamma_f)*inpDXL[2]*face.vecPF[2];
			// vfDiffusive += gamma_f*inpDXR[0]*face.vecNF[0];
			// vfDiffusive += gamma_f*inpDXR[1]*face.vecNF[1];
			// vfDiffusive += gamma_f*inpDXR[2]*face.vecNF[2];
			
			// vfDiffusive += (1.0-gamma_f)*inpDXL[0]*face.vecSkewness[0];
			// vfDiffusive += (1.0-gamma_f)*inpDXL[1]*face.vecSkewness[1];
			// vfDiffusive += (1.0-gamma_f)*inpDXL[2]*face.vecSkewness[2];
			// vfDiffusive += gamma_f*inpDXR[0]*face.vecSkewness[0];
			// vfDiffusive += gamma_f*inpDXR[1]*face.vecSkewness[1];
			// vfDiffusive += gamma_f*inpDXR[2]*face.vecSkewness[2];
			
			vfDiffusive += wCL*inpDXL[0]*face.vecSkewness[0];
			vfDiffusive += wCL*inpDXL[1]*face.vecSkewness[1];
			vfDiffusive += wCL*inpDXL[2]*face.vecSkewness[2];
			vfDiffusive += wCR*inpDXR[0]*face.vecSkewness[0];
			vfDiffusive += wCR*inpDXR[1]*face.vecSkewness[1];
			vfDiffusive += wCR*inpDXR[2]*face.vecSkewness[2];
			
			vfDiffusive = max(0.0,min(1.0,vfDiffusive));
			// vfDiffusive = max(min(minPhiL,minPhiR),min(max(maxPhiL,maxPhiR),vfDiffusive));
			
			
			double gamma_diff_f = gamma_f;
			
			
			// gamF = 1.0;
			
			double vf = gamF*vfCompressive + (1.0-gamF)*vfDiffusive;
			
			
			
			// // skwness correction
			// if(LR==0) {
				// vf += inpDXL[0]*face.vecSkewness[0];
				// vf += inpDXL[1]*face.vecSkewness[1];
				// vf += inpDXL[2]*face.vecSkewness[2];
			// }
			// else{
				// vf += inpDXR[0]*face.vecSkewness[0];
				// vf += inpDXR[1]*face.vecSkewness[1];
				// vf += inpDXR[2]*face.vecSkewness[2];
			// }
			
			
			if(LR==0) vfL = vf;
			if(LR==1) vfR = vf;
			// if(LR==0) vfL = max(min(minPhiL,minPhiR),min(max(maxPhiL,maxPhiR),vf));
			// if(LR==1) vfR = max(min(minPhiL,minPhiR),min(max(maxPhiL,maxPhiR),vf));
			
			if(LR==0) face.varL[controls.fVF_NVD] = gamF*(1.0-gamma_comp_f) + (1.0-gamF)*(1.0-gamma_diff_f);
			if(LR==1) face.varR[controls.fVF_NVD] = gamF*(1.0-gamma_comp_f) + (1.0-gamF)*(1.0-gamma_diff_f);
			
		}

		// set velues
		face.varL[fn_out] = vfL;
		face.varR[fn_out] = vfR;
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) ++proc_num;
		
	}
	
	
			
			

			
			
			
			
			
}