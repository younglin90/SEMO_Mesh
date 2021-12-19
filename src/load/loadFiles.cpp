
using namespace std;

#include "load.h" 
#include "../mesh/build.h"  
// #include "../controls/build.h" 
// #include "../utility/read.h" 
// #include "../solvers/build.h" 



void SEMO_Mesh_Load::boundaryProperties(string file, 
SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<double>& value){
	
	name.clear();
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;
	bool startReading2 = false;
	

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading || startReading2){
				if(nextToken.find("{") != string::npos) startReading=true;
				if(!startReading && nextToken.find("}") != string::npos) startReading2=false;
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("boundaryField") != string::npos ){
					startReading=true;
					startReading2 = true;
					// cout << nextToken << endl;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			for(int i=0; i<mesh.boundary.size(); ++i){
				if(read[line].find(mesh.boundary[i].name) != string::npos){
					
					name.push_back(mesh.boundary[i].name);
					
					for(int line2=line; line2<read.size(); ++line2){
						if(read[line2].find("type") != string::npos){
							istringstream iss2(read[line2]);
							string tempstring2;
							string tempstring3;
							iss2 >> tempstring2 >> tempstring3;
							tempstring3.erase(tempstring3.find(";"),1); 
							type.push_back(tempstring3);
							if(
							type.back() == "fixedValue" ||
							type.back() == "inletOutlet" ||
							type.back() == "switch"){
								istringstream iss2(read[line2+1]);
								string tempstring2;
								string tempstring3;
								iss2 >> tempstring2 >> tempstring3;
								tempstring3.erase(tempstring3.find(";"),1); 
								value.push_back(stod(tempstring3));
							}
							else{
								value.push_back(0.0);
							}
							break;
						}
					}
					
				}
			}
		}
	}
	inputFile.close();
	
}


void SEMO_Mesh_Load::boundaryVelocities(string file, 
SEMO_Mesh_Builder& mesh, vector<string>& name, vector<string>& type, vector<vector<double>>& value){
	
	type.clear();
	value.clear();
	
	vector<string> read;
	bool startReading = false;
	bool startReading2 = false;
	

	ifstream inputFile;
	string openFileName;
	openFileName = file;
	inputFile.open(openFileName);
	string nextToken;
	if(!inputFile.fail()){
		while(getline(inputFile, nextToken)){
			if(startReading || startReading2){
				if(nextToken.find("{") != string::npos) startReading=true;
				if(!startReading && nextToken.find("}") != string::npos) startReading2=false;
				if(nextToken.find("}") != string::npos) startReading=false;
				read.push_back(nextToken);
			}
			else{
				if( nextToken.find("boundaryField") != string::npos ){
					startReading=true;
					startReading2 = true;
					// cout << nextToken << endl;
				}
			}
		}
		
		for(int line=0; line<read.size(); ++line){
			for(int i=0; i<mesh.boundary.size(); ++i){
				if(read[line].find(mesh.boundary[i].name) != string::npos){
					
					name.push_back(mesh.boundary[i].name);
					
					for(int line2=line; line2<read.size(); ++line2){
						if(read[line2].find("type") != string::npos){
							istringstream iss2(read[line2]);
							string tempstring2;
							string tempstring3;
							iss2 >> tempstring2 >> tempstring3;
							tempstring3.erase(tempstring3.find(";"),1); 
							type.push_back(tempstring3);
							if(type.back() == "fixedValue"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								lineee.erase(lineee.find("("),1); 
								lineee.erase(lineee.find(")"),1); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0] >> tmpD[1] >> tmpD[2];
								value.push_back(tmpD);
								break;
							}
							else if(type.back() == "surfaceNormalFixedValue"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0];
								value.push_back(tmpD);
								break;
							}
							else if(type.back() == "inletOutlet"){
								string lineee = read[line2+1];
								lineee.erase(lineee.find("value"),5); 
								
								istringstream iss2(lineee);
								vector<double> tmpD(3,0.0);
								iss2 >> tmpD[0];
								value.push_back(tmpD);
								break;
							}
							else{
								vector<double> tmpD(3,0.0);
								value.push_back(tmpD);
								break;
							}
						}
					}
					
				}
			}
		}
	}
	else{
		cerr << "| #Warning : Unable to open file for reading : " << openFileName << endl;
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	inputFile.close();
	
}






	

void SEMO_Mesh_Load::file(string fileName, map<string,string> &store){

	std::ifstream inputFile;
	inputFile.open(fileName);

	std::string line;
	bool noRead=false;
	while( std::getline(inputFile, line) )
	{
		if(noRead){
			if(line.find("}") != string::npos){
				noRead=false;
			}
		}
		else{
			std::istringstream is_line(line);
			std::string key;
			if( std::getline(is_line, key, '=') )
			{
				std::string value;
				if( std::getline(is_line, value) ) 
					store.insert(make_pair(trim(key), trim(value)));
			}
		
			if( line.find("{") != string::npos ){
				noRead=true;
			}
		}
	}
	
	inputFile.close();
	
}


void SEMO_Mesh_Load::file(
	string fileName, 
	string c1, 
	map<string,string> &store){

	std::ifstream inputFile;
	inputFile.open(fileName);

	std::string line;
	bool start1=false;
	while( std::getline(inputFile, line) )
	{
		if(start1){
			if(line.find("}") != string::npos){
				start1=false;
			}
			else{
				std::istringstream is_line(line);
				std::string key;
				if( std::getline(is_line, key, '=') )
				{
					std::string value;
					if( std::getline(is_line, value) ) 
						store.insert(make_pair(trim(key), trim(value)));
				}
			}
		}
		else{
			if( line.find(c1) != string::npos ){
				if( line.find("{") != string::npos )
					start1=true;
					
				std::getline(inputFile, line);
				if( line.find("{") != string::npos )
					start1=true;
			}
		}
		
	}
	
	inputFile.close();
	
}



void SEMO_Mesh_Load::file(
	string fileName, 
	string c1, string c2, 
	map<string,string> &store){

	std::ifstream inputFile;
	inputFile.open(fileName);

	std::string line;
	bool start1=false;
	bool start2=false;
	while( std::getline(inputFile, line) )
	{
		if(start1){
			if(start2){
				if(line.find("}") != string::npos){
					start1=false;
					start2=false;
				}
				else{
					std::istringstream is_line(line);
					std::string key;
					if( std::getline(is_line, key, '=') )
					{
						std::string value;
						if( std::getline(is_line, value) ) 
							store.insert(make_pair(trim(key), trim(value)));
					}
				}
			}
			else{
				if( line.find(c2) != string::npos ){
					if( line.find("{") != string::npos )
						start2=true;
						
					std::getline(inputFile, line);
					if( line.find("{") != string::npos )
						start2=true;
				}
			}
		}
		else{
			if( line.find(c1) != string::npos ){
				if( line.find("{") != string::npos )
					start1=true;
					
				std::getline(inputFile, line);
				if( line.find("{") != string::npos )
					start1=true;
			}
		}
		
	}
	
	inputFile.close();
	
}



void SEMO_Mesh_Load::file(
	string fileName, 
	string c1, string c2, string c3,
	map<string,string> &store){

	std::ifstream inputFile;
	inputFile.open(fileName);
	
	// cout << "AAAAAAAAAAAA " << endl;

	std::string line;
	bool start1=false;
	bool start2=false;
	bool start3=false;
	while( std::getline(inputFile, line) )
	{
		// cout << line << endl;
		if(start1){
			if(start2){
				if(start3){
					if(line.find("}") != string::npos){
						start1=false;
						start2=false;
						start3=false;
					}
					else{
						std::istringstream is_line(line);
						std::string key;
						if( std::getline(is_line, key, '=') )
						{
							std::string value;
							if( std::getline(is_line, value) ) 
								store.insert(make_pair(trim(key), trim(value)));
						}
					}
				}
				else{
					if( line.find(c3) != string::npos ){
						// start3=true;
						if( line.find("{") != string::npos )
							start3=true;
							
						std::getline(inputFile, line);
						if( line.find("{") != string::npos )
							start3=true;
					}
				}
			}
			else{
				if( line.find(c2) != string::npos ){
					// start2=true;
					if( line.find("{") != string::npos )
						start2=true;
						
					std::getline(inputFile, line);
					if( line.find("{") != string::npos )
						start2=true;
				}
			}
		}
		else{
			if( line.find(c1) != string::npos ){
				if( line.find("{") != string::npos )
					start1=true;
					
				std::getline(inputFile, line);
				if( line.find("{") != string::npos )
					start1=true;
			}
		}
		
	}
	
	inputFile.close();
	
}



void SEMO_Mesh_Load::vecters(string in, vector<string>& out){

	std::istringstream is_line(in);
	std::string key;
	// while(!feof(in)){
	// }

	while( std::getline(is_line, key, ' ') )
	{
		if( key.find("(") != string::npos ){
			key.erase(key.find("("),1);
		}
		if( key.find(")") != string::npos ){
			key.erase(key.find(")"),1);
		}
		trim(key);
		out.push_back(key);
	}
	
}


