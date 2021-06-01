#include "read.h" 


//앞에 있는 개행 문자 제거 
string &SEMO_Utility_Read::ltrim(std::string &s) { 
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace)))); 
	return s; 
}

//뒤에 있는 개행 문자 제거 
string &SEMO_Utility_Read::rtrim(std::string &s) { 
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end()); 
	return s; 
}

//양쪽 끝의 개행 문자 제거 
string &SEMO_Utility_Read::trim(std::string &s) { 
	return ltrim(rtrim(s)); 
}
	
	

void SEMO_Utility_Read::file(string fileName, map<string,string> &store){

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


void SEMO_Utility_Read::file(
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



void SEMO_Utility_Read::file(
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



void SEMO_Utility_Read::file(
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



void SEMO_Utility_Read::vecters(string in, vector<string>& out){

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






// int main(){
	
	// SEMO_Utility_Read read;
	

	// // thermophysicalProperties
	// map<string,string> thermophysicalProperties;
	// read.file("./constant/thermophysicalProperties",thermophysicalProperties);
	
	// vector<string> phases;
	// read.vecters(thermophysicalProperties["phases"], phases);

	
	
	// vector<map<string,string>> rho;
	// for(int i=0; i<phases.size(); ++i){
		// map<string,string> tmp;
		// read.file("./constant/thermophysicalProperties",
			// phases[i], "thermodynamics", "rho", tmp);
		
		// rho.push_back(tmp);
	// }
	
	
	// vector<map<string,string>> mu;
	// for(int i=0; i<phases.size(); ++i){
		// map<string,string> tmp;
		// read.file("./constant/thermophysicalProperties",
			// phases[i], "transport", "mu", tmp);
		
		// mu.push_back(tmp);
	// }
	
	
	
	// // turbulenceProperties
	// map<string,string> turbulenceProperties;
	// read.file("./constant/turbulenceProperties",turbulenceProperties);
	
	// string simulationType = turbulenceProperties["simulationType"];
	
	
	
	// // controlDict
	// map<string,string> controlDict;
	// read.file("./system/controlDict",controlDict);


	// // for(auto& it : controlDict){
		// // cout << "key: " << it.first << " " << "value: " << it.second << endl;
	// // }
	
	
	// // fvSchemes
	// map<string,string> fvSchemes_fluxScheme;
	// read.file("./system/fvSchemes", "fluxScheme", fvSchemes_fluxScheme);
	
	// map<string,string> fvSchemes_ddtSchemes;
	// read.file("./system/fvSchemes", "ddtSchemes", fvSchemes_ddtSchemes);
	
	// map<string,string> fvSchemes_gradSchemes;
	// read.file("./system/fvSchemes", "gradSchemes", fvSchemes_gradSchemes);
	
	// map<string,string> fvSchemes_divSchemes;
	// read.file("./system/fvSchemes", "divSchemes", fvSchemes_divSchemes);
	
	// map<string,string> fvSchemes_laplacianSchemes;
	// read.file("./system/fvSchemes", "laplacianSchemes", fvSchemes_laplacianSchemes);
	
	// map<string,string> fvSchemes_interpolationSchemes;
	// read.file("./system/fvSchemes", "interpolationSchemes", fvSchemes_interpolationSchemes);
	
	// map<string,string> fvSchemes_snGradSchemes;
	// read.file("./system/fvSchemes", "snGradSchemes", fvSchemes_snGradSchemes);
	
	
	
	// // fvSolution
	// map<string,string> fvSolution_solvers_P;
	// read.file("./system/fvSolution", 
		// "solvers", "P", 
		// fvSolution_solvers_P);
		
	// map<string,string> fvSolution_solvers_U;
	// read.file("./system/fvSolution", 
		// "solvers", "U", 
		// fvSolution_solvers_U);
		
	// map<string,string> fvSolution_solvers_VF;
	// read.file("./system/fvSolution", 
		// "solvers", "VF", 
		// fvSolution_solvers_VF);
		
	// map<string,string> fvSolution_PIMPLE;
	// read.file("./system/fvSolution", 
		// "PIMPLE",  
		// fvSolution_PIMPLE);
		
	// map<string,string> fvSolution_relaxationFactors_momentumEq_U;
	// read.file("./system/fvSolution", 
		// "relaxationFactors", "momentumEq", "U",
		// fvSolution_relaxationFactors_momentumEq_U);
		
	// map<string,string> fvSolution_relaxationFactors_pressureEq_P;
	// read.file("./system/fvSolution", 
		// "relaxationFactors", "pressureEq", "P",
		// fvSolution_relaxationFactors_pressureEq_P);
		
	// map<string,string> fvSolution_relaxationFactors_pressureEq_U;
	// read.file("./system/fvSolution", 
		// "relaxationFactors", "pressureEq", "U",
		// fvSolution_relaxationFactors_pressureEq_U);
		
	// map<string,string> fvSolution_relaxationFactors_speciesEq_VF;
	// read.file("./system/fvSolution", 
		// "relaxationFactors", "speciesEq", "VF",
		// fvSolution_relaxationFactors_speciesEq_VF);

	// map<string,string> fvSolution_limiter;
	// read.file("./system/fvSolution", 
		// "limiter",  
		// fvSolution_limiter);
		
	
	
// }