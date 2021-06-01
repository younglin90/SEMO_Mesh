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
using namespace std;

class SEMO_Utility_Read {
public:

	string &ltrim(std::string &s);
	string &rtrim(std::string &s);
	string &trim(std::string &s);
	
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
private:

};


