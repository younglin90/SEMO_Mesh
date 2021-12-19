#pragma once
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <mpi.h>
#include "../species/species.h" 

// #include "build.h"
class SEMO_Mesh_Builder;
class SEMO_Controls_Builder;

class SEMO_Mesh_Load {
private:

public:
    SEMO_Mesh_Load() {
    }
    ~SEMO_Mesh_Load() {
    }
	
	string &ltrim(std::string &s);
	string &rtrim(std::string &s);
	string &trim(std::string &s);
	
	template<typename T>
	void loadDatasAtVTU(
		ifstream& inputFile, string dataName, vector<T>& outData);
		
	template<typename T>
	void readBinary(
		string& word, vector<T>& outData);
		
	template<typename T>
	void readCompress(
		string& word, vector<T>& outData);
		
	int Base64decode_len(const char *bufcoded);
	int Base64decode(char *bufplain, const char *bufcoded);
	
	void file
		(string fileName, map<string,string> &store);
	void file
		(string fileName, string c1, map<string,string> &store);
	void file
		(string fileName, string c1, string c2, map<string,string> &store);
	void file
		(string fileName, string c1, string c2, string c3, map<string,string> &store);

	void vecters
		(string in, vector<string>& out);
	
	void OpenFoam(SEMO_Mesh_Builder &in, string folder);
	
	// void vtu(SEMO_Mesh_Builder &mesh, string folder);
	
	void vtu(string folder, 
		SEMO_Mesh_Builder &mesh, 
		SEMO_Controls_Builder &controls,
		vector<SEMO_Species>& species);
	
	void boundaryProperties(string file, 
			SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<double>& value);
	void boundaryVelocities(string file, 
			SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<vector<double>>& value);
			
	// bool boolBinary = false;
	bool boolCompress = false;
	


};


