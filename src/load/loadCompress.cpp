#include <fstream>
#include <cstring>
using namespace std;

#include "load.h" 
#include "../mesh/build.h"  
#include "../controls/build.h" 
#include "../utility/read.h" 
#include "../solvers/build.h" 

#include <zlib.h>

template void SEMO_Mesh_Load::readCompress<int>(
	string& word, vector<int>& outData);
template void SEMO_Mesh_Load::readCompress<double>(
	string& word, vector<double>& outData);
template void SEMO_Mesh_Load::readCompress<string>(
	string& word, vector<string>& outData);
template<typename T>
void SEMO_Mesh_Load::readCompress(
	string& word, vector<T>& outData){
		
	// 데이터 클리어
	outData.clear();
	
	// string 나누기, 44 는 header 사이즈
	string header_string = word.substr(0, 44);
	string value_string = word.substr(44, word.size()-44);
	
	// 인풋 데이터 캐릭터형으로 복사
	const char *header_in = header_string.c_str();
	const char *value_in = value_string.c_str();
	
	// 헤더 base64 디코딩 후 저장
	int header_byte_size = 8;
	int header_encode_length = header_string.size();
	char header_in_char[header_encode_length];
	std::copy(header_in, header_in + header_encode_length, header_in_char);
	
	int header_decoded_length = Base64decode_len(header_in_char);
	char* header_decoded = (char*)malloc(header_decoded_length);
	Base64decode(header_decoded, header_in_char);
	
	// 헤더 부분 저장
	int pointer_end = 0;
	vector<long long> header;
	for(int i=0; i<4; ++i){
		int pointer_str = pointer_end;
		pointer_end += header_byte_size;
		
		char header_buffer[header_byte_size];
		std::copy(header_decoded + pointer_str, header_decoded + pointer_end, header_buffer);
		long long header_tmp;
		std::memcpy( &header_tmp, header_buffer, header_byte_size );
		header.push_back(header_tmp);
	}
	free(header_decoded);
	
	// 데이터 부분만 base64 디코딩 후 저장
	int value_encode_length = value_string.size();
	char* value_in_char = (char*)malloc(value_encode_length);
	std::copy(value_in, value_in + value_encode_length, value_in_char);
	int value_decoded_length = Base64decode_len(value_in_char);
	char* value_decoded = (char*)malloc(value_decoded_length);
	Base64decode(value_decoded, value_in_char);
	free(value_in_char);

	// zlib 로 압축해제
	Bytef* value_decoded_byte = (Bytef*)malloc(value_decoded_length);
	uLong value_sizee = value_decoded_length;
	char* value_buffer = (char*)malloc(value_decoded_length);
	std::copy(value_decoded, value_decoded + value_decoded_length, value_buffer);
	std::memcpy( value_decoded_byte, value_buffer, value_decoded_length );
	uLong inflate_size = header[1];
	Bytef *inflate_data = (Bytef*)malloc(inflate_size);
	
	int error = uncompress2(inflate_data, &inflate_size, value_decoded_byte, &value_sizee);
	if(error != Z_OK){
		cout << "zlib error !!!!!!!!!!!!" << endl;
		cout << error << endl;
	}
	free(value_decoded);
	free(value_buffer);
	
	// 데이터 부분 저장
	int size_value_byte = sizeof(T);
	pointer_end = 0;
	for(int i=0; i<inflate_size/size_value_byte; ++i){
		int pointer_str = pointer_end;
		pointer_end += size_value_byte;
		
		char header_buffer[size_value_byte];
		std::copy(inflate_data + pointer_str, inflate_data + pointer_end, header_buffer);
		T header_tmp;
		std::memcpy( &header_tmp, header_buffer, size_value_byte );
		outData.push_back(header_tmp);
	}
	free(inflate_data);
	
	
	
}