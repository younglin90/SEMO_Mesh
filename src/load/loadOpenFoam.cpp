
#include "load.h" 
#include "../mesh/build.h"  
#include "../controls/build.h" 
#include "../utility/read.h" 
#include "../solvers/build.h" 

void SEMO_Mesh_Load::OpenFoam(SEMO_Mesh_Builder &mesh, string folder){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute load OpenFoam files ... ";
	
	if(rank == 0){
		
		string gridFolderName = folder;
		string pointsName = "points";
		string facesName = "faces";
		string ownerName = "owner";
		string neighbourName = "neighbour";
		string boundaryName = "boundary";
		
		ifstream inputFile;
		string openFileName;
		
		// points
		openFileName = gridFolderName + "/" + pointsName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		string nextToken;
		bool startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					vector<double> xyz(3,0.0);
					nextToken.erase(nextToken.find("("),1); 
					nextToken.erase(nextToken.find(")"),1); 
					stringstream sstream(nextToken);
					string word;
					char del = ' ';
					int num=0;
						// cout << "A" << endl;
					while (getline(sstream, word, del)){
						// cout << word << endl;
						xyz[num] = stold(word);
						++num;
					}
						// cout << "B" << endl;
					mesh.addPoint();
					mesh.points.back().x = xyz[0];
					mesh.points.back().y = xyz[1];
					mesh.points.back().z = xyz[2];
						// cout << "C" << endl;
					
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		// cout << "------------------------------------" << endl;
		// cout << "point x,y,z size : " << mesh.points.size() << endl;
		inputFile.close();
		
		

		// faces
		openFileName = gridFolderName + "/" + facesName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		startInput=false;
		bool continueInput=false;
		string saveToken;
		// saveToken.append(" ");
		// saveToken.append("asdfadsf ");
		// saveToken.append("asdfadsf ");
		// saveToken.append("asdfadsf ");
		// saveToken.append("asdfadsf ");
		// cout<<saveToken<< endl;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" && !continueInput ){
					break;
				}
				else{
					if(nextToken.size()==1) continue;
					
					if(nextToken.find(")") == string::npos){
						saveToken.append(" ");
						rtrim(nextToken);
						saveToken.append(nextToken);
						// cout<<saveToken<< endl;
						continueInput = true;
						continue;
					}
					saveToken.append(" ");
					rtrim(nextToken);
					saveToken.append(nextToken);
					
					saveToken.replace(saveToken.find("("),1," ");
					saveToken.replace(saveToken.find(")"),1," ");
					// cout << saveToken << endl;
					istringstream iss(saveToken);
					int tempint;
					iss >> tempint;
					
					
					mesh.addFace();
					
					while(iss >> tempint){
						mesh.faces.back().points.push_back(tempint);
					}
					saveToken.clear();
					continueInput = false;
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		
		
		
		
		// cout << "face size : " << mesh.faces.size() << endl;
		inputFile.close();
		
		
		
		// owner
		openFileName = gridFolderName + "/" + ownerName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		int temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					
					
					while(iss >> tempint){
						mesh.faces[temp_num].owner = tempint;
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		// cout << "owner size : " << temp_num << endl;
		inputFile.close();
		
		
		
		
		// neighbour
		openFileName = gridFolderName + "/" + neighbourName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		temp_num = 0;
		startInput=false;
		while(getline(inputFile, nextToken)){
			string asignToken;
			
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				else{
					istringstream iss(nextToken);
					int tempint;
					while(iss >> tempint){
						if(tempint<0) break;
						mesh.faces[temp_num].neighbour = tempint;
						// cout << tempint << endl;
						
						++temp_num;
					}
				}
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}			
		}
		// cout << "neighbour size : " << temp_num << endl;
		inputFile.close();
		
		
		
		// boundary
		openFileName = gridFolderName + "/" + boundaryName;
		inputFile.open(openFileName);
		if(inputFile.fail()){
			cerr << "Unable to open file for reading : " << openFileName << endl;
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		vector<string> boundary_name;
		vector<char> boundary_type;
		vector<int> boundary_nFaces;
		vector<int> boundary_startFace;
		vector<int> boundary_myProcNo;
		vector<int> boundary_neighbProcNo;
		
		string backToken;
		startInput=false;
		vector<string> setToken;
		while(getline(inputFile, nextToken)){
			string asignToken;
			if(startInput){
				if( asignToken.assign(nextToken, 0, 1) == ")" ){
					break;
				}
				setToken.push_back(nextToken.c_str());
			}
			else{
				if( asignToken.assign(nextToken, 0, 1) == "(" ){
					startInput=true;
				}
			}
			backToken = nextToken;
		}
		
		string names;
		vector<string> setToken2;
		startInput=false;
		for (auto item : setToken) {
			string asignToken;
			if(startInput){
				if( item.find("}") != string::npos ){
					trim(names);
					
			// cout << names.length() << endl;
					
					boundary_name.push_back(names);
					int l=0;
					for (auto item2 : setToken2) {
						if( item2.find("nFaces") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_nFaces.push_back(temptempint);
						}
						if( item2.find("startFace") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_startFace.push_back(temptempint);
						}
						if( item2.find("myProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_myProcNo.push_back(temptempint);
							++l;
						}
						if( item2.find("neighbProcNo") != string::npos ){
							istringstream iss(item2);
							string temptemp;
							int temptempint;
							iss >> temptemp >> temptempint;
							boundary_neighbProcNo.push_back(temptempint);
							++l;
						}
					}
					if(l==0) {
						boundary_myProcNo.push_back(-1);
						boundary_neighbProcNo.push_back(-1);
					}
					
					startInput=false;
					setToken2.clear();
				}
				setToken2.push_back(item.c_str());
			}
			else{ 
				if( item.find("{") != string::npos ){
					names.clear();
					names = backToken;
					startInput=true;
				}
			}
			backToken = item;
		}
		// cout << mesh.faces.size() << " " << mesh.cells.size() << " " << mesh.points.size() << " " << mesh.boundary.size() << endl;
		
		inputFile.close();
		mesh.boundary.clear();
		int nbcs=boundary_name.size();
		for (int i=0; i<boundary_name.size(); ++i) {
			// cout << endl;
			// cout << endl;
			// cout << trim(boundary_name[i]) << endl;
			// cout << boundary_nFaces[i] << endl;
			// cout << boundary_startFace[i] << endl;
			// cout << boundary_myProcNo[i] << endl;
			// cout << boundary_neighbProcNo[i] << endl;
			
			mesh.addBoundary();
		}
		for (int i=0; i<mesh.boundary.size(); ++i) {
			mesh.boundary[i].name = trim(boundary_name[i]);
			mesh.boundary[i].nFaces = boundary_nFaces[i];
			mesh.boundary[i].startFace = boundary_startFace[i];
			mesh.boundary[i].myProcNo = boundary_myProcNo[i];
			mesh.boundary[i].neighbProcNo = boundary_neighbProcNo[i];
			
			
			// cout << mesh.boundary.back().startFace << endl;
		}
	}
	
	if(rank==0) cout << "-> completed" << endl;
	if(rank==0) cout << "└────────────────────────────────────────────────────" << endl;
	

	mesh.check();
	
	mesh.buildCells();
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// mesh.checkLists();
	
	// cell's points, faces connection
	mesh.connectCelltoFaces();
	
	mesh.connectCelltoPoints();
	
	mesh.cellsGlobal();
		
	mesh.informations();
		
		
		
	// mesh.checkQualities();
	
	
	
}

