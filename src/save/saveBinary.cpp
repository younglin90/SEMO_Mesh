#include <fstream>
using namespace std;

#include "save.h"  

template void SEMO_Mesh_Save::writeBinary<int>(ofstream& out, vector<int>& vecInp);
template void SEMO_Mesh_Save::writeBinary<double>(ofstream& out, vector<double>& vecInp);
template<typename T>
void SEMO_Mesh_Save::writeBinary(ofstream& out, vector<T>& vecInp)
{
	int datasize = 8;
	int dataByteSize = sizeof(T);
	
	// cout << dataByteSize << endl;
	
	int data_length = datasize + dataByteSize * vecInp.size();
	int encoded_data_length = Base64encode_len(data_length);
	char* base64_string = (char*)malloc(encoded_data_length);
	char* data = (char*)malloc(data_length);
	{
		long long byte_size = datasize * vecInp.size();
		char* data1 = (char*)&byte_size;
		std::copy(data1, data1 + datasize, data);
	}
	int numm = 0;
	for(auto& inp : vecInp) {
		char* data1 = (char*)&inp;
		std::copy(data1, data1 + dataByteSize, data + datasize + numm*dataByteSize);
		++numm;
	}
	Base64encode(base64_string, data, data_length);
	out << base64_string << endl;
	free(data);
	free(base64_string);
	
}
