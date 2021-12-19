#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>
#include <unistd.h>
using namespace std;

#include "save.h" 
#include "build.h" 
#include "../controls/build.h" 

void SEMO_Mesh_Save::gnuplot(int iter, double dtime, vector<double>& norm){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	if(rank==0) {
		ofstream outputFile;
		

		string gnu = "residuals.gnu";
		if(access(gnu.c_str(),F_OK) != 0){

			outputFile.open(gnu);
		
			// outputFile << "set logscale y" << endl;
			outputFile << "set title \"Residuals\" #textcolor rgb \"white\"" << endl;
			outputFile << "set ylabel \"Residual\" #textcolor rgb \"white\"" << endl;
			outputFile << "set xlabel \"Iteration\" #textcolor rgb \"white\"" << endl;
			outputFile << "set grid" << endl;
			outputFile << "plot \"<cat residuals\" using 1:2 title 'continuity' with lines lw 1,\\" << endl;
			outputFile << "\"<cat residuals\" using 1:3 title 'x-mom' with lines lw 1,\\" << endl;
			outputFile << "\"<cat residuals\" using 1:4 title 'y-mom' with lines lw 1,\\" << endl;
			outputFile << "\"<cat residuals\" using 1:5 title 'z-mom' with lines lw 1,\\" << endl;
			outputFile << "\"<cat residuals\" using 1:6 title 'energy' with lines lw 1" << endl;
			outputFile << "pause 1" << endl;
			outputFile << "reread" << endl;
			outputFile.close();
		
		}
 
		
		string filenamePlot = "residuals";
		outputFile.open(filenamePlot, ios::app);
		if(outputFile.fail()){
			cerr << "Unable to write file for writing." << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
	
		outputFile << iter << " ";
		outputFile.precision(4);
		for(int i=0; i<norm.size(); ++i){
			outputFile << scientific << norm[i] << " ";
		}
		outputFile << scientific << dtime;
		outputFile.unsetf(ios::scientific);
		outputFile << endl;
		
		outputFile.close();
	}
		
	
}

