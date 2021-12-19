#include <cstring>
using namespace std;

#include "load.h" 
#include "../mesh/build.h"  
#include "../controls/build.h" 
#include "../utility/read.h" 
#include "../solvers/build.h" 

template void SEMO_Mesh_Load::readBinary<int>(
	string& word, vector<int>& outData);
template void SEMO_Mesh_Load::readBinary<double>(
	string& word, vector<double>& outData);
template void SEMO_Mesh_Load::readBinary<string>(
	string& word, vector<string>& outData);
template<typename T>
void SEMO_Mesh_Load::readBinary(
	string& word, vector<T>& outData){
		
	outData.clear();
		
	int datasize = 8;
	int dataByteSize = sizeof(T);
		
	const char *data_in = word.c_str();
		
	int decoded_data_length = Base64decode_len(data_in);
	char* data_out = (char*)malloc(decoded_data_length);
	
	Base64decode(data_out, data_in);
	// printf("The string\n[%s]\ndecodes from base64 as:\n[%s]\n", data_in, data_out);
	
		// cout << sizeof(data_in) << endl;
		// cout << decoded_data_length << " " << sizeof(data_out) << endl;
	
	char buffer1[datasize];
	
	int pointer_end = 0;
	pointer_end += datasize;
	
	std::copy(data_out, data_out + pointer_end, buffer1);
	long long total_byte;
	std::memcpy( &total_byte, buffer1, datasize );
	int data_length = total_byte/8;
	
	// cout << data_length << endl;
	
	for(int i=0; i<data_length; ++i){
		int pointer_str = pointer_end;
		pointer_end += dataByteSize;
		
		char buffer2[dataByteSize];
		std::copy(data_out + pointer_str, data_out + pointer_end, buffer2);
		T data;
		std::memcpy( &data, buffer2, dataByteSize );
		
		// cout << data << " ";
		
		outData.push_back(data);
	}
			
			
	free(data_out);
	
	
}