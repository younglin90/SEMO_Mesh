#include <fstream>
using namespace std;

#include "save.h"  


template void SEMO_Mesh_Save::writeAscii<int>(ofstream& out, vector<int>& vecInp);
template void SEMO_Mesh_Save::writeAscii<double>(ofstream& out, vector<double>& vecInp);
template<typename T>
void SEMO_Mesh_Save::writeAscii(ofstream& out, vector<T>& vecInp)
{
	for(auto& inp : vecInp){
		out << inp << " ";
	}
	out << endl;
}
