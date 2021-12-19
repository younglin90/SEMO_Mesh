#include <fstream>
using namespace std;

#include "save.h" 

#include <zlib.h>

template void SEMO_Mesh_Save::writeCompress<int>(ofstream& out, vector<int>& vecInp, int compressSize);
template void SEMO_Mesh_Save::writeCompress<double>(ofstream& out, vector<double>& vecInp, int compressSize);
template<typename T>
void SEMO_Mesh_Save::writeCompress(ofstream& out, vector<T>& vecInp, int compressSize)
{
	int value_data_size = vecInp.size();
	
	int headerByteSize = 8;
	int dataByteSize = sizeof(T);
	
	
	// cout << dataByteSize << endl;
	
	int header_size = headerByteSize*4;
	int data_size = dataByteSize*value_data_size;
	long long header[4] = {0};
	header[0] = 1;
	header[1] = data_size;
	header[2] = data_size;
	header[3] = 2*data_size;
	Bytef *deflate_data = (Bytef*)malloc(header[3]);
	uLong deflate_size = header[3];
	Bytef* raw_data = (Bytef*)malloc(data_size);
	uLong raw_size = data_size;
	
	int numm = 0;
	for(auto& inp : vecInp){
		char* header_char0 = (char*)&inp;
		std::copy(header_char0, header_char0 + dataByteSize, raw_data + numm*dataByteSize);
		++numm;
	}
	if(compress2(deflate_data, &deflate_size, raw_data, raw_size, compressSize) != Z_OK){
		cout << "zlib error !!!!!!!!!!!!" << endl;
	}
	// cout << endl;
	// cout << deflate_size <<  endl;
	header[3] = deflate_size;
	int header_encoded_length = Base64encode_len(header_size);
	char* base64_string_header = (char*)malloc(header_encoded_length);
	char* header_data = (char*)malloc(header_size);
	
	char* header_char1 = (char*)&header[0];
	char* header_char2 = (char*)&header[1];
	char* header_char3 = (char*)&header[2];
	char* header_char4 = (char*)&header[3];
	std::copy(header_char1, header_char1 + headerByteSize, header_data+headerByteSize*0);
	std::copy(header_char2, header_char2 + headerByteSize, header_data+headerByteSize*1);
	std::copy(header_char3, header_char3 + headerByteSize, header_data+headerByteSize*2);
	std::copy(header_char4, header_char4 + headerByteSize, header_data+headerByteSize*3);
	Base64encode(base64_string_header, header_data, header_size);
	
	int data_encoded_length = Base64encode_len(deflate_size);
	char* base64_string_data = (char*)malloc(data_encoded_length);
	Base64encode(base64_string_data, (char*)deflate_data, deflate_size);
	
	// cout << base64_string_header << " " << header_encoded_length << endl;
	out << base64_string_header << base64_string_data << endl; 
	
	free(base64_string_header);
	free(base64_string_data);
	free(deflate_data);
	free(header_data);
	free(raw_data); 
}

