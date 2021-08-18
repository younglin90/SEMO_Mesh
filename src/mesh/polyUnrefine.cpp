#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>

#include "build.h"
#include "mpi.h"
#include "polyAMR.h"
#include "geometric.h"



void SEMO_Poly_AMR_Builder::sortCellCanUnrefine(
	vector<int>& vLevel, 
	int saveLevel, 
	int& num, 
	vector<vector<int>>& vGroupLevel, 
	vector<vector<int>>& vGroupNumber
	) {
		
		
		
		
	// cout << "AAAAAAA" << endl;
		
	vector<int> tmpLevels;
	vector<int> tmpNumbers;
	int limitNum = 7;
	if (saveLevel == 0) limitNum = 100;
	if (saveLevel == 1) limitNum = 100;
	while (1) {
		if (vLevel.size() - 1 >= num) {
			int level = vLevel[num];
			// cout << num << " " << tmpLevels.size()<< " " << saveLevel << " " << level << endl;

			if (
				saveLevel == level &&
				tmpLevels.size() < limitNum
				) {
				tmpLevels.push_back(saveLevel);
				tmpNumbers.push_back(num);
			}
			else if (
				saveLevel == level &&
				tmpLevels.size() == limitNum) {
				if (std::find(tmpNumbers.begin(), tmpNumbers.end(), -1) == tmpNumbers.end()) {
					tmpLevels.push_back(saveLevel);
					tmpNumbers.push_back(num);
					vGroupLevel.push_back(tmpLevels);
					vGroupNumber.push_back(tmpNumbers);
				}
				break;
			}
			else if (
				saveLevel < level
				) {
				// cout << "start" << endl;

				if (tmpLevels.size() == limitNum+1) {
					tmpLevels.clear();
					tmpNumbers.clear();
				}

				sortCellCanUnrefine(vLevel, saveLevel + 1, num, vGroupLevel, vGroupNumber);
				tmpLevels.push_back(-1);
				tmpNumbers.push_back(-1);
				// cout << tmpLevels.size() << endl;
			}
			else {
				--num;
				break;
			}
			++num;
		}
		else {
			break;
		}
	}
   
}




void SEMO_Poly_AMR_Builder::sortFaceCanUnrefine(
	vector<int>& vLevel, 
	int saveLevel, 
	int& num, 
	vector<vector<int>>& vGroupLevel, 
	vector<vector<int>>& vGroupNumber
	) {
		
	vector<int> tmpLevels;
	vector<int> tmpNumbers;
	int limitNum = 3;
	if (saveLevel == 0) limitNum = 100;
	if (saveLevel == 1) limitNum = 100;
	while (1) {
		if (vLevel.size() - 1 >= num) {
			int level = vLevel[num];
			// cout << num << " " << tmpLevels.size()<< " " << saveLevel << " " << level << endl;

			if (
				saveLevel == level &&
				tmpLevels.size() < limitNum
				) {
				tmpLevels.push_back(saveLevel);
				tmpNumbers.push_back(num);
			}
			else if (
				saveLevel == level &&
				tmpLevels.size() == limitNum) {
				if (std::find(tmpNumbers.begin(), tmpNumbers.end(), -1) == tmpNumbers.end()) {
					tmpLevels.push_back(saveLevel);
					tmpNumbers.push_back(num);
					vGroupLevel.push_back(tmpLevels);
					vGroupNumber.push_back(tmpNumbers);
				}
				break;
			}
			else if (
				saveLevel < level
				) {
				// cout << "start" << endl;

				if (tmpLevels.size() == limitNum+1) {
					tmpLevels.clear();
					tmpNumbers.clear();
				}

				sortFaceCanUnrefine(vLevel, saveLevel + 1, num, vGroupLevel, vGroupNumber);
				tmpLevels.push_back(-1);
				tmpNumbers.push_back(-1);
				// cout << tmpLevels.size() << endl;
			}
			else {
				--num;
				break;
			}
			++num;
		}
		else {
			break;
		}
	}
   
}









void SEMO_Poly_AMR_Builder::polyUnrefine(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	// cout << "000000" << endl;
	

	if(rank==0) cout << "----Unrefine start-----" << endl;
	
	int maxLevel = 10;
	
	
	
	
	SEMO_Mesh_Geometric geometric;
	

	
	// cout << "111111" << endl;

	//====================================================
	// 셀 Unrefine : Unrefine 되는 셀 & 면 조사
	
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	

	SEMO_Utility_Math math;
	vector<vector<double>> gradVF;
	// math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradVF);
	math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	
	vector<bool> boolCellUnrefine(mesh.cells.size(),false);
	for(int i=0; i<mesh.cells.size(); ++i){
		// if( distr(eng) < 0.9 ){
			// boolCellUnrefine[i] = true;
		// }
			// boolCellUnrefine[i] = true;
		
		if( mesh.cells[i].var[controls.indicatorAMR] < 0.5 ){
			boolCellUnrefine[i] = true;
		}
		
			// boolCellUnrefine[i] = true;
		
		
		// 만약 셀의 레벨이 0 이면, false
		if(mesh.cells[i].level == 0) boolCellUnrefine[i] = false;
	}
	
	
	//#################################################
	// prcoface 앞 셀은 Refine 안되게
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			boolCellUnrefine[face.owner] = false;
		}
	}
	
	
	// cout << "111111" << endl;
	
	//====================================================
	// 셀 Unrefine : MPI
	vector<int> cLevel_send, cLevel_recv;
	vector<int> cUnrefine_send, cUnrefine_recv;
	if(size>1){
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cLevel_send.push_back(mesh.cells[face.owner].level);
				
				if(boolCellUnrefine[face.owner]){
					cUnrefine_send.push_back(1);
				}
				else{
					cUnrefine_send.push_back(0);
				}
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					cLevel_send, cLevel_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
		
		mpi.setProcsFaceDatas(
					cUnrefine_send, cUnrefine_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					   
	}
	cLevel_send.clear();
	cUnrefine_send.clear();
	
	
	//====================================================
	// 셀 Unrefine : Unrefine 되는 셀 옆에 레벨보다 작으면 안됨
	int proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			if(mesh.cells[face.owner].level < mesh.cells[face.neighbour].level){
				boolCellUnrefine[face.owner] = false;
			}
			if(mesh.cells[face.owner].level > mesh.cells[face.neighbour].level){
				boolCellUnrefine[face.neighbour] = false;
			}
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			if(mesh.cells[face.owner].level < cLevel_recv[proc_num]){
				boolCellUnrefine[face.owner] = false;
			}
			++proc_num;
		}
	}
	
	
	//====================================================
	// 셀 Unrefine : 그룹 셀 일정량 넘으면 그룹 전체 다 언리파인 true, 아니면 false
	// 그룹 지정
	// 내부 면 찾기
	
	

	if(rank==0) cout << "cell groupping start" << endl;
	
	// 오리지널 그룹 
	vector<vector<int>> groupCellListsLevelZero;
	vector<vector<int>> groupCellLevelListsLevelZero;
	int saveGroup = mesh.cells[0].group;
	for(int i=0; i<mesh.cells.size(); ++i){
		
		vector<int> tmp(1,i);
		vector<int> tmp2(1,mesh.cells[i].level);
		groupCellListsLevelZero.push_back(tmp);
		groupCellLevelListsLevelZero.push_back(tmp2);
		
		if(mesh.cells[i].level==0) continue;
		
		auto& cellOrg = mesh.cells[i];
		int groupOrg = cellOrg.group;
		
		// cout << groupOrg<< endl;
		
		for(int j=i+1; j<mesh.cells.size(); ++j){
			auto& cell = mesh.cells[j];
			int group = cell.group;
			int level = cell.level;
			
			if(groupOrg != group) {
				i = j-1;
				break;
			}
			
			if(j+1==mesh.cells.size()) i = j;
			
			groupCellListsLevelZero.back().push_back(j);
			groupCellLevelListsLevelZero.back().push_back(level);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	// // for(auto& i : groupCellListsLevelZero){
	// for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		// // for(auto& j : i){
		// for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
			// cout << rank << " : " << groupCellListsLevelZero[i][0] << " " << groupCellListsLevelZero[i][j] << " " << mesh.cells.size() << " " << groupCellLevelListsLevelZero[i][j] << endl;
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// 오리지널 그룹 내에서 그룹핑
	// vector<vector<int>> groupCellListsAlongLevels;
	// vector<vector<int>> groupCellLevelListsAlongLevels;
	// vector<vector<int>> groupCellListsCanUnrefine;
	// for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		
		// bool boolAllLevelOne = true;
		// for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
			// if(groupCellLevelListsLevelZero[i][j] != 1) 
				// boolAllLevelOne = false;
		// }
		// if(boolAllLevelOne){
			// vector<int> tmp;
			// for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
				// tmp.push_back(groupCellListsLevelZero[i][j]);
			// }
			// groupCellListsCanUnrefine.push_back(tmp);
			// // continue;
		// }
		
		// int saveLevel = groupCellLevelListsLevelZero[i][0];
		
		// for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
			// int iCellOrg = groupCellListsLevelZero[i][j];
			// int levelOrg = groupCellLevelListsLevelZero[i][j];
			
			// vector<int> tmp(1,iCellOrg);
			// vector<int> tmp2(1,levelOrg);
			// groupCellListsAlongLevels.push_back(tmp);
			// groupCellLevelListsAlongLevels.push_back(tmp2);
			
			// int tmpNumb = 1;
			// for(int k=j+1; k<groupCellListsLevelZero[i].size(); ++k){
				// int iCell = groupCellListsLevelZero[i][k];
				// int level = groupCellLevelListsLevelZero[i][k];
				
				
				// if( levelOrg != level || (tmpNumb == 8 && level != 1) ) {
					
					// if(tmpNumb == 8) 
						// groupCellListsCanUnrefine.push_back(groupCellListsAlongLevels.back());
					
					
					// j = k-1;
					// break;
				// }
				
				// // if(k+1==groupCellListsLevelZero[i].size()) j = k;
				
				// groupCellListsAlongLevels.back().push_back(iCell);
				// groupCellLevelListsAlongLevels.back().push_back(level);
				
				
				// ++tmpNumb;
				
				// if(k+1==groupCellListsLevelZero[i].size()) {
					// if(tmpNumb == 8 && level != 1) 
						// groupCellListsCanUnrefine.push_back(groupCellListsAlongLevels.back());
					// j = k;
				// }
				
			// }
			
		// }
		
	// }
	
	
	
	// cout << endl;
	// cout << endl;
	
	// for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		// for(auto& j : groupCellListsLevelZero[i]){
			// cout << i << " " << j << " " << mesh.cells[j].level << endl;
		// }
	// }
	
	
	// cout << endl;
	// cout << endl;

	vector<vector<int>> groupCellListsCanUnrefine;
	for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		
		bool boolAllLevelOne = true;
		for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
			if(groupCellLevelListsLevelZero[i][j] != 1) {
				boolAllLevelOne = false;
				break;
			}
		}
		
		if(boolAllLevelOne==true){
			vector<int> tmpVal;
			for(int j=0; j<groupCellListsLevelZero[i].size(); ++j){
				int pushVal = groupCellListsLevelZero[i][j];
				tmpVal.push_back(pushVal);
			}
			groupCellListsCanUnrefine.push_back(tmpVal);
		}
		
		
		if(boolAllLevelOne==false){
			vector<vector<int>> tmpGroupLevel;
			vector<vector<int>> tmpGroupNumber;
			int tmpNum = 0;
			
			// cout << "AAA" << endl;
			
			sortCellCanUnrefine(
				groupCellLevelListsLevelZero[i],
				0,
				tmpNum,
				tmpGroupLevel,
				tmpGroupNumber
				);
			
			// cout << "BBB" << endl;
			
			int startCellNum = groupCellListsLevelZero[i][0];
			for(auto& j : tmpGroupNumber){
				for(auto& k : j){
					k += startCellNum;
				}
			}
			
			for(auto& j : tmpGroupNumber){
				groupCellListsCanUnrefine.push_back(j);
			}
			
		}
		
	}
	
	
	
	// for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		// for(auto& j : groupCellListsCanUnrefine[i]){
			// cout << i << " " << j << " " << mesh.cells[j].level << endl;
		// }
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	
	
	
	
	
	
	
	// // cout << " Setting Group " << endl;
	
	// for(int i=0; i<groupCellListsAlongLevels.size(); ++i){
		// // for(auto& j : i){
		// for(int j=0; j<groupCellListsAlongLevels[i].size(); ++j){
			// cout << rank << " : " << groupCellListsAlongLevels[i][0] << " " << groupCellListsAlongLevels[i][j] << " " << mesh.cells.size() << " " << groupCellLevelListsAlongLevels[i][j] << endl;
		// }
	// }
	
	// cout << " Setting Group " << endl;
	
	// for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		// // for(auto& j : i){
		// for(int j=0; j<groupCellListsCanUnrefine[i].size(); ++j){
			// cout << rank << " : " << groupCellListsCanUnrefine[i][0] << " " << groupCellListsCanUnrefine[i][j] << endl;
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	// cout << "111111" << endl;
	
	
	//====================================================
	// 셀 Unrefine : 그룹된 셀 기준으로 Unrefine 적용
	vector<vector<int>> groupCellNewUnrefine;
	for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		bool boolAllUnrefine = true;
		for(int j=0; j<groupCellListsCanUnrefine[i].size(); ++j){
			int iCell = groupCellListsCanUnrefine[i][j];
			if(boolCellUnrefine[iCell] == false) 
				boolAllUnrefine = false;
		}
		
		// cout << groupCellListsCanUnrefine[i].size() << endl;
		
		if(boolAllUnrefine){
			vector<int> tmpCells;
			
			for(int j=0; j<groupCellListsCanUnrefine[i].size(); ++j){
				int iCell = groupCellListsCanUnrefine[i][j];
				tmpCells.push_back(iCell);
			}
			groupCellNewUnrefine.push_back(tmpCells);
		}
	}
	
	std::fill(boolCellUnrefine.begin(), boolCellUnrefine.end(), false);
	
	for(auto& i : groupCellNewUnrefine){
		for(auto& j : i){
			boolCellUnrefine[j] = true;
		}
	}
	
	
	
	// for(int i=0; i<groupCellNewUnrefine.size(); ++i){
		// for(int j=0; j<groupCellNewUnrefine[i].size(); ++j){
			// cout << rank << " : " << groupCellNewUnrefine[i][0] << " " << groupCellNewUnrefine[i][j] << endl;
		// }
	// }
	
	// cout << "22222" << endl;

	
	//====================================================
	// 셀 Unrefine : MPI
	if(size>1){
		cUnrefine_send.clear();
	
		for(int i=0; i<mesh.faces.size(); ++i){
			auto& face = mesh.faces[i];
			
			if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				cLevel_send.push_back(mesh.cells[face.owner].level);
				
				if(boolCellUnrefine[face.owner]){
					cUnrefine_send.push_back(1);
				}
				else{
					cUnrefine_send.push_back(0);
				}
			}
		}
		
		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					cUnrefine_send, cUnrefine_recv,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
					
		cUnrefine_send.clear();
					   
	}
	
	
	
	
	//====================================================
	// 면 Unrefine : outer faces, internal faces searching
	// 포인트 Unrefine : cell center point searching
	vector<bool> boolFaceUnrefine(mesh.faces.size(),false);
	vector<bool> boolInternalFaces(mesh.faces.size(),false);
	vector<vector<int>> groupIntFacesUnrefine;
	vector<vector<int>> groupOutFacesUnrefine;
	for(auto& i : groupCellNewUnrefine){
		vector<int> internalFaces;
		vector<int> outerFaces;
		for(auto& cell0 : i){
			for(auto& cell1 : i){
				if(cell0==cell1) continue;
				for(auto& face0 : mesh.cells[cell0].faces){
					for(auto& face1 : mesh.cells[cell1].faces){
						if(face0==face1){
							if ( std::find( internalFaces.begin(), internalFaces.end(), face0 ) 
								== internalFaces.end() ) {
								internalFaces.push_back(face0);
							}
						}
					}
				}
			}
		}
		for(auto& cell0 : i){
			for(auto& face0 : mesh.cells[cell0].faces){
				if ( std::find( internalFaces.begin(), internalFaces.end(), face0 ) 
					== internalFaces.end() ) {
					outerFaces.push_back(face0);
				}
			}
		}
		groupIntFacesUnrefine.push_back(internalFaces);
		groupOutFacesUnrefine.push_back(outerFaces);
		
			// cout << internalFaces.size() << endl;
		
		for(auto& j : internalFaces){
			boolFaceUnrefine[j] = true;
			boolInternalFaces[j] = true;
		}
		
		for(auto& j : outerFaces){
			auto& face = mesh.faces[j];
			if(face.getType()==SEMO_Types::INTERNAL_FACE){
				if(
				mesh.cells[face.owner].level == mesh.cells[face.neighbour].level &&
				boolCellUnrefine[face.owner] == true &&
				boolCellUnrefine[face.neighbour] == true
				){
					boolFaceUnrefine[j] = true;
				}
				
				if(
				mesh.cells[face.owner].level > mesh.cells[face.neighbour].level &&
				boolCellUnrefine[face.owner] == true
				){
					boolFaceUnrefine[j] = true;
				}
				
				if(
				mesh.cells[face.owner].level < mesh.cells[face.neighbour].level &&
				boolCellUnrefine[face.neighbour] == true
				){
					boolFaceUnrefine[j] = true;
				}
			}
			else{
				boolFaceUnrefine[j] = true;
			}
		}
	}
	// proc faces
	proc_num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			if(
			mesh.cells[face.owner].level == cLevel_recv[proc_num] &&
			boolCellUnrefine[face.owner] == true &&
			cUnrefine_recv[proc_num] == true
			){
				boolFaceUnrefine[i] = true;
			}
			else if(
			mesh.cells[face.owner].level > cLevel_recv[proc_num] &&
			boolCellUnrefine[face.owner] == true
			){
				boolFaceUnrefine[i] = true;
			}
			else if(
			mesh.cells[face.owner].level < cLevel_recv[proc_num] &&
			cUnrefine_recv[proc_num] == true
			){
				boolFaceUnrefine[i] = true;
			}
			else{
				boolFaceUnrefine[i] = false;
			}
			
			++proc_num;
		}
	}
	

	// for(int i=0; i<groupOutFacesUnrefine.size(); ++i){
		// for(int j=0; j<groupOutFacesUnrefine[i].size(); ++j){
			// cout << rank << " : " << groupOutFacesUnrefine[i][0] << " " << groupOutFacesUnrefine[i][j] << endl;
		// }
	// }
	
	
	

	//====================================================
	// 면 Unrefine : 그룹핑
	
	
	if(rank==0) cout << "face groupping start" << endl;
	
	
	// 오리지널 그룹 
	vector<vector<int>> groupFaceListsLevelZero;
	vector<vector<int>> groupFaceLevelListsLevelZero;
	saveGroup = mesh.faces[0].group;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		vector<int> tmp(1,i);
		vector<int> tmp2(1,mesh.faces[i].level);
		groupFaceListsLevelZero.push_back(tmp);
		groupFaceLevelListsLevelZero.push_back(tmp2);
		
		if(mesh.faces[i].level==0) continue;
		
		auto& faceOrg = mesh.faces[i];
		int groupOrg = faceOrg.group;
		
		for(int j=i+1; j<mesh.faces.size(); ++j){
			auto& face = mesh.faces[j];
			int group = face.group;
			int level = face.level;
			
			if(groupOrg != group) {
				i = j-1;
				break;
			}
			
			if(j+1==mesh.faces.size()) i = j;
			
			groupFaceListsLevelZero.back().push_back(j);
			groupFaceLevelListsLevelZero.back().push_back(level);
			
		}
		
		// if(groupFaceListsLevelZero.back().size()==3) {
			// for(auto& j : groupFaceListsLevelZero.back()){
				// cout << i << " " << j << endl;
			// }
			
			// for(int j=i-5; j<i+5; ++j){
				// auto& face = mesh.faces[j];
				// int group = face.group;
				// int level = face.level;
				// cout << i << " " << j << " " << group << " " << level << endl;
				
			// }
			

			// auto fa = groupFaceListsLevelZero.back();
			// int face0 = fa[0];
			// int centerPoint = -1;
			// for(auto& point0 : mesh.faces[face0].points){
				// int nCentP=0;
				// for(int j=1; j<groupFaceListsLevelZero.back().size(); ++j){
					// int face1 = fa[j];
					// for(auto& point1 : mesh.faces[face1].points){
						// if(point0==point1) ++nCentP;
					// }
				// }
				// if(nCentP>=2){
					// centerPoint = point0;
					// break;
				// }
			// }
			
			// if(centerPoint==-1){
				// cout << endl;
				// cout << "| #Error 0 : No searching center point" << endl;
				// for(int j=0; j<fa.size(); ++j){
					// int face1 = fa[j];
					// for(auto& point1 : mesh.faces[face1].points){
						// cout << face1 << " " << point1 << endl;
					// }
				// }
				// cout << endl;
				// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			// }
		// }
	}
	
	
	
	
	
	
	
	
	
	
	
	// vector<vector<int>> groupFaceListsAlongLevels;
	// vector<vector<int>> groupFaceLevelListsAlongLevels;
	// vector<vector<int>> groupFaceListsCanUnrefine;
	// for(int i=0; i<groupFaceListsLevelZero.size(); ++i){
		
		// if(boolInternalFaces[groupFaceListsLevelZero[i][0]]==true){
			// i = i + groupFaceListsLevelZero[i].size()-1;
			// continue;
		// }
		
		// bool boolAllLevelOne = true;
		// for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
			// if(groupFaceLevelListsLevelZero[i][j] != 1) 
				// boolAllLevelOne = false;
		// }
		// if(boolAllLevelOne){
			// vector<int> tmp;
			// for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
				// tmp.push_back(groupFaceListsLevelZero[i][j]);
			// }
			// groupFaceListsCanUnrefine.push_back(tmp);
		// }
		
		// int saveLevel = groupFaceLevelListsLevelZero[i][0];
		
		// for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
			// int iFaceOrg = groupFaceListsLevelZero[i][j];
			// int levelOrg = groupFaceLevelListsLevelZero[i][j];
			
			// vector<int> tmp(1,iFaceOrg);
			// vector<int> tmp2(1,levelOrg);
			// groupFaceListsAlongLevels.push_back(tmp);
			// groupFaceLevelListsAlongLevels.push_back(tmp2);
			
			// int tmpNumb = 1;
			// for(int k=j+1; k<groupFaceListsLevelZero[i].size(); ++k){
				// int iCell = groupFaceListsLevelZero[i][k];
				// int level = groupFaceLevelListsLevelZero[i][k];
				
				
				// if( levelOrg != level || (tmpNumb == 4 && level != 1) ) {
					
					// if(tmpNumb == 4) 
						// groupFaceListsCanUnrefine.push_back(groupFaceListsAlongLevels.back());
					
					// j = k-1;
					// break;
				// }
				
				// groupFaceListsAlongLevels.back().push_back(iCell);
				// groupFaceLevelListsAlongLevels.back().push_back(level);
				
				
				// ++tmpNumb;
				
				// if(k+1==groupFaceListsLevelZero[i].size()) {
					// if(tmpNumb == 4 && level != 1) 
						// groupFaceListsCanUnrefine.push_back(groupFaceListsAlongLevels.back());
					// j = k;
				// }
			// }
		// }
	// }
	
	
	
	
	
	// cout << endl;
	// cout << endl;
	
	// for(int i=0; i<groupFaceListsLevelZero.size(); ++i){
		// for(auto& j : groupFaceListsLevelZero[i]){
			// cout << i << " " << j << " " << mesh.faces[j].level << endl;
		// }
	// }
	
	// cout << endl;
	// cout << endl; 
	
	
	

	vector<vector<int>> groupFaceListsCanUnrefine;
	for(int i=0; i<groupFaceListsLevelZero.size(); ++i){
		
		if(boolInternalFaces[groupFaceListsLevelZero[i][0]]==true){
			i = i + groupFaceListsLevelZero[i].size()-1;
			continue;
		}
		
		bool boolAllLevelOne = true;
		for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
			if(groupFaceLevelListsLevelZero[i][j] != 1) {
				boolAllLevelOne = false;
				break;
			}
		}

		if(boolAllLevelOne==true){
			vector<int> tmpVal;
			for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
				int pushVal = groupFaceListsLevelZero[i][j];
				tmpVal.push_back(pushVal);
			}
			groupFaceListsCanUnrefine.push_back(tmpVal);
		}
		
		
		if(boolAllLevelOne==false){
			vector<vector<int>> tmpGroupLevel;
			vector<vector<int>> tmpGroupNumber;
			int tmpNum = 0;
			sortFaceCanUnrefine(
				groupFaceLevelListsLevelZero[i],
				0,
				tmpNum,
				tmpGroupLevel,
				tmpGroupNumber
				);
			
			int startFaceNum = groupFaceListsLevelZero[i][0];
			for(auto& j : tmpGroupNumber){
				for(auto& k : j){
					k += startFaceNum;
				}
			}
			
			for(auto& j : tmpGroupNumber){
				groupFaceListsCanUnrefine.push_back(j);
			}
			
		}
		
	}




	
	// for(int i=0; i<groupFaceListsCanUnrefine.size(); ++i){
		// for(auto& j : groupFaceListsCanUnrefine[i]){
			// cout << i << " " << j << " " << mesh.faces[j].level << endl;
		// }
	// }


	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
	// // 내부 면인 것들 제거 
	// vector<bool> boolIntFaceOneFace(groupFaceListsCanUnrefine.size(),false);
	// for(int i=0; i<groupFaceListsCanUnrefine.size(); ++i){
		// if(groupFaceListsCanUnrefine[i].size() == 1){
			// int iFace = groupFaceListsCanUnrefine[i][0];
			// if(mesh.faces[iFace].level > 0){
				// boolIntFaceOneFace[i] = true;
			// }
		// }
	// }

	// int numN2 = 0;
	// groupFaceListsCanUnrefine.erase( 
		// std::remove_if( groupFaceListsCanUnrefine.begin(), groupFaceListsCanUnrefine.end(), 
		// [&boolIntFaceOneFace, &numN2](vector<int> const& v) { 
		// return boolIntFaceOneFace[numN2++]; 
		// }), groupFaceListsCanUnrefine.end());
		
	
	
	
	// for(int i=0; i<groupFaceListsCanUnrefine.size(); ++i){
		// for(auto& j : groupFaceListsCanUnrefine[i]){
			// cout << i << " " << j << " " << mesh.faces[j].group << " " << mesh.faces[j].level << endl;
		// }
	// }
	
	



	// numN2 = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
	// // for(int i=0; i<groupFaceListsCanUnrefine.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// if(boolInternalFaces[i]==true){
			// continue;
		// }
		
		// // cout << groupFaceListsCanUnrefine[numN2].size() << " " <<  groupFaceListsCanUnrefine[numN2][0] << endl;
		
		// // combine faces (outer faces)
		// // if(groupFaceListsCanUnrefine.size()>0){
			// int face0 = groupFaceListsCanUnrefine[numN2][0];
			// if(i==face0){
				// int centerPoint = -1;
				// for(auto& point0 : mesh.faces[face0].points){
					// int nCentP=0;
					// for(int j=1; j<groupFaceListsCanUnrefine[numN2].size(); ++j){
						// int face1 = groupFaceListsCanUnrefine[numN2][j];
						// for(auto& point1 : mesh.faces[face1].points){
							// if(point0==point1) ++nCentP;
						// }
					// }
					// if(nCentP>=3){
						// centerPoint = point0;
						// break;
					// }
				// }
				
				// if(centerPoint==-1){
					// cout << endl;
					// cout << "| #Error 2 : No searching center point" << endl;
					// for(int j=0; j<groupFaceListsCanUnrefine[numN2].size(); ++j){
						// int face1 = groupFaceListsCanUnrefine[numN2][j];
						// for(auto& point1 : mesh.faces[face1].points){
							// cout << face1 << " " << point1 << endl;
						// }
					// }
					// cout << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				// i = i + groupFaceListsCanUnrefine[numN2].size()-1;
				// ++numN2;
				// if(numN2 == groupFaceListsCanUnrefine.size()) numN2 = 0;
				
				
				// continue;
			// }
		// // }
	// }
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	
	//====================================================
	// 면 Unrefine : 그룹된 면 기준으로 Unrefine 적용
	vector<vector<int>> groupFaceNewUnrefine;
	for(int i=0; i<groupFaceListsCanUnrefine.size(); ++i){
		bool boolAllUnrefine = true;
		for(int j=0; j<groupFaceListsCanUnrefine[i].size(); ++j){
			int iFace = groupFaceListsCanUnrefine[i][j];
			if(boolFaceUnrefine[iFace] == false) 
				boolAllUnrefine = false;
		}
		
		if(boolAllUnrefine){
			vector<int> tmpFaces;
			for(int j=0; j<groupFaceListsCanUnrefine[i].size(); ++j){
				int iFace = groupFaceListsCanUnrefine[i][j];
				tmpFaces.push_back(iFace);
			}
			groupFaceNewUnrefine.push_back(tmpFaces);
		}
	}
	
	std::fill(boolFaceUnrefine.begin(), boolFaceUnrefine.end(), false);
	
	for(auto& i : groupFaceNewUnrefine){
		for(auto& j : i){
			boolFaceUnrefine[j] = true;
		}
	}
	
	for(auto& i : groupCellNewUnrefine){
		for(auto& j : i){
			if(boolCellUnrefine[j] == true){
				for(auto& k : mesh.cells[j].faces){
					if(boolInternalFaces[k] == true){
						boolFaceUnrefine[k] = true;
					}
				}
			}
		}
	}
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// for(int i=0; i<groupFaceListsLevelZero.size(); ++i){
		// for(int j=0; j<groupFaceListsLevelZero[i].size(); ++j){
			// int aaaa = groupFaceListsLevelZero[i][j];
			// cout << rank << " : " << groupFaceListsLevelZero[i][0] << " " << groupFaceListsLevelZero[i][j] << " " << 
			// mesh.faces[aaaa].level << " " << mesh.faces[aaaa].group << endl;
		// }
	// }
	
	
	// for(int i=0; i<groupFaceListsAlongLevels.size(); ++i){
		// for(int j=0; j<groupFaceListsAlongLevels[i].size(); ++j){
			// int aaaa = groupFaceListsAlongLevels[i][j];
			// cout << rank << " : " << groupFaceListsAlongLevels[i][0] << " " << groupFaceListsAlongLevels[i][j] << " " << 
			// mesh.faces[aaaa].level << " " << mesh.faces[aaaa].group << endl;
		// }
	// }
	
	
	
				
	// for(int i=0; i<groupFaceNewUnrefine.size(); ++i){
		// // for(int j=0; j<groupFaceNewUnrefine[i].size(); ++j){
			// // cout << rank << " : " << groupFaceNewUnrefine[i][0] << " " << groupFaceNewUnrefine[i][j] << endl;
		// // }
		// int face0 = groupFaceNewUnrefine[i][0];
		// int centerPoint = -1;
		// for(auto& point0 : mesh.faces[face0].points){
			// int nCentP=0;
			// for(int j=1; j<groupFaceNewUnrefine[i].size(); ++j){
				// int face1 = groupFaceNewUnrefine[i][j];
				// for(auto& point1 : mesh.faces[face1].points){
					// if(point0==point1) ++nCentP;
				// }
			// }
			// if(nCentP>=3){
				// centerPoint = point0;
				// break;
			// }
		// }
		
		// if(centerPoint==-1){
			// cout << endl;
			// cout << "| #Error 2 : No searching center point" << endl;
			// for(int j=0; j<groupFaceNewUnrefine[i].size(); ++j){
				// int face1 = groupFaceNewUnrefine[i][j];
				// for(auto& point1 : mesh.faces[face1].points){
					// cout << face1 << " " << point1 << endl;
				// }
			// }
			// cout << endl;
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
			
			
	// }
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	
	
	

	//====================================================
	// 엣지 생성
	
	if(rank==0) cout << "edge add start" << endl;
	
	
	vector<vector<int>> edgesPoints;
	vector<vector<int>> facesEdges(mesh.faces.size(),vector<int>(0,0));
	vector<vector<int>> pointsFaces(mesh.points.size(),vector<int>(0,0));
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		int pointSize = face.points.size();
		for(int j=0; j<pointSize; ++j){
			int ipoint0 = face.points[j];
			int ipoint1 = ( j+1 == pointSize ? face.points[0] : face.points[j+1] );
			vector<int> matchFaces;
			for(auto& k : pointsFaces[ipoint0]){
				for(auto& l : pointsFaces[ipoint1]){
					if( k == l ) {
						matchFaces.push_back(l);
					}
				}
			}
			
			if(matchFaces.size()==0){
				vector<int> tmpPoints;
				tmpPoints.push_back(ipoint0);
				tmpPoints.push_back(ipoint1);
				edgesPoints.push_back(tmpPoints);
				
				facesEdges[i].push_back(edgesPoints.size()-1);
				
			}
			else{
				int iFace = matchFaces[0];
				int iEdgeSave = -1;
				for(auto& iEdge : facesEdges[iFace]){
					if(
					(edgesPoints[iEdge][0]==ipoint0 && edgesPoints[iEdge][1]==ipoint1) ||
					(edgesPoints[iEdge][1]==ipoint0 && edgesPoints[iEdge][0]==ipoint1) 
					){
						iEdgeSave = iEdge;
						break;
					}
				}
				
				facesEdges[i].push_back(iEdgeSave);
			}
			pointsFaces[ipoint0].push_back(i);
			
		}
	}
	pointsFaces.clear();
	
	
	vector<vector<int>> edgesFaces(edgesPoints.size(),vector<int>(0,0));
	for(int i=0; i<mesh.faces.size(); ++i){
		for(auto& j : facesEdges[i]){
			edgesFaces[j].push_back(i);
		}
	}
	

	
	//====================================================
	// Unrefine : 엣지 레벨 (최대 포인트 레벨)
	vector<int> edgeLevel(edgesPoints.size(),0);
	for(int i=0; i<edgesPoints.size(); ++i){
		int point0 = edgesPoints[i][0];
		int point1 = edgesPoints[i][1];
		
		edgeLevel[i] = max(
			mesh.points[point0].level, 
			mesh.points[point1].level);
	}
	

	// cout << "444444" << endl;
	
	//====================================================
	// Unrefine : 엣지 Unrefine
	vector<bool> boolEdgeUnrefine(edgesPoints.size(),true);
	for(int i=0; i<edgesPoints.size(); ++i){
		int point0 = edgesPoints[i][0];
		int point1 = edgesPoints[i][1];
		
		bool tmpBool = true;
		for(auto& j : edgesFaces[i]){
			int level = mesh.faces[j].level;
			if(
			boolFaceUnrefine[j] == false &&
			edgeLevel[i] <= level 
			){
				tmpBool = false;
			}
		}
		
		if(tmpBool==false) boolEdgeUnrefine[i] = false;
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	//====================================================
	// Unrefine : 엣지 refine MPI
	
	// if(size>1){
		// vector<int> eUnrefine_send, eUnrefine_recv;
		// vector<int> countsEdge(size,0);
		// vector<int> displsEdge(size,0);
		
		// for(int ip=0; ip<size; ++ip){
			// for(auto& boundary : mesh.boundary){
				// if(ip == boundary.neighbProcNo){
					// int str = boundary.startFace;
					// int end = str + boundary.nFaces;
					// for(int i=str; i<end; ++i){
						// SEMO_Face& face = mesh.faces[i];
						// for(auto& j : facesEdges[i]){
							// if(boolEdgeUnrefine[j]){
								// eUnrefine_send.push_back(1);
							// }
							// else{
								// eUnrefine_send.push_back(0);
							// }
						// }
						// countsEdge[ip] += facesEdges[i].size();
					// }
					// break;
				// }
			// }
		// }
		
		// displsEdge[0] = 0;
		// for(int ip=1; ip<size; ++ip){
			// displsEdge[ip] = displsEdge[ip-1] + countsEdge[ip-1];
			
		// }
		
		// SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(eUnrefine_send, eUnrefine_recv,
                              // countsEdge, countsEdge, 
                              // displsEdge, displsEdge);

		// int tmp_proc_num = 0;
		// for(int ip=0; ip<size; ++ip){
			// for(auto& boundary : mesh.boundary){
				// if(ip == boundary.neighbProcNo){
					// int str = boundary.startFace;
					// int end = str + boundary.nFaces;
					// for(int i=str; i<end; ++i){
						// SEMO_Face& face = mesh.faces[i];
						// vector<int> tmpEdgeBool;
						// for(auto& j : facesEdges[i]){
							// tmpEdgeBool.push_back(eUnrefine_recv[tmp_proc_num]);
							// ++tmp_proc_num;
						// }
						// // std::reverse(tmpEdgeBool.begin(), tmpEdgeBool.end());
						// for(int j=0; j<facesEdges[i].size(); ++j){
							// if(tmpEdgeBool[j]==1){
								// boolEdgeUnrefine[facesEdges[i][j]] = true;
							// }
						// }
					// }
					// break;
				// }
			// }
		// }	
	// }
	
	
	
	
	
	//====================================================
	// Unrefine : 엣지 연결 포인트 제거 체크 (면 중심, 셀 중심 다 포함)
	// vector<bool> boolPointUnrefine(mesh.points.size(),false);
	// for(int i=0; i<mesh.faces.size(); ++i){
		// for(auto& j : facesEdges[i]){
			// if(boolEdgeUnrefine[j] == true){
				// int iPoint0 = edgesPoints[j][0];
				// int iPoint1 = edgesPoints[j][1];
				// int pLevel0 = mesh.points[iPoint0].level;
				// int pLevel1 = mesh.points[iPoint1].level;
				// if(pLevel0 > pLevel1){
					// boolPointUnrefine[iPoint0] = true;
				// }
				// else if(pLevel1 > pLevel0){
					// boolPointUnrefine[iPoint1] = true;
				// }
				// else if(pLevel0 == pLevel1){
					// boolPointUnrefine[iPoint0] = true;
					// boolPointUnrefine[iPoint1] = true;
				// }	
			// }
		// }
	// }
	vector<bool> boolPointUnrefine(mesh.points.size(),true);
	for(int i=0; i<edgesPoints.size(); ++i){
		int iPoint0 = edgesPoints[i][0];
		int iPoint1 = edgesPoints[i][1];
		int pLevel0 = mesh.points[iPoint0].level;
		int pLevel1 = mesh.points[iPoint1].level;
		if(boolEdgeUnrefine[i]==false){
			boolPointUnrefine[iPoint0] = false;
			boolPointUnrefine[iPoint1] = false;
		}
		else{
			bool tmpBool0 = true;
			bool tmpBool1 = true;
			for(auto& j : edgesFaces[i]){
				int level = mesh.faces[j].level;
				if(boolFaceUnrefine[j] == true){
					if(pLevel0 <= level-1) tmpBool0 = false;
					if(pLevel1 <= level-1) tmpBool1 = false;
				}
			}
			
			if(tmpBool0==false) boolPointUnrefine[iPoint0] = false;
			if(tmpBool1==false) boolPointUnrefine[iPoint1] = false;
			
		}
	}
	
	// int test22 = 0;
	// for(int i=0; i<edgesPoints.size(); ++i){
		// if(boolEdgeUnrefine[i]) ++test22;
	// }
	// cout << test22 << endl;
	// int test33 = 0;
	// for(int i=0; i<mesh.points.size(); ++i){
		// if(boolPointUnrefine[i]==true){
			// ++test33;
		// }
	// }
	
	// cout << test33 << endl;
	
	
	
	
	//====================================================
	// Unrefine : 새로운 포인트 넘버링
	
	if(rank==0) cout << "new point numbering start" << endl;
	
	
	vector<int> newPointNum(mesh.points.size(),-1);
	int newNumbering = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		auto& point = mesh.points[i];
		if(boolPointUnrefine[i]==false){
			newPointNum[i] = newNumbering;
			++newNumbering;
		}
	}
	
	
	//====================================================
	// Unrefine : 새로운 면 넘버링
	
	if(rank==0) cout << "new face numbering start" << endl;
	
	
	vector<int> newFaceNum(mesh.faces.size(),-1);
	newNumbering = 0;
	int numN = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		// delete faces (internal faces)
		if(boolInternalFaces[i]==true){
			continue;
		}
		
		// combine faces (outer faces)
		if(groupFaceNewUnrefine.size()>0){
			if(i==groupFaceNewUnrefine[numN][0]){
				
				for(int j=0; j<groupFaceNewUnrefine[numN].size(); ++j){
					newFaceNum[i] = newNumbering;
					++i;
				}
				--i;
				
				++numN;
				if(numN == groupFaceNewUnrefine.size()) numN = 0;
				++newNumbering;
				continue;
			}
		}
		
		
		// original faces
		newFaceNum[i] = newNumbering;
		
		++newNumbering;
	}
	
	
	
	//====================================================
	// Unrefine : 새로운 셀 넘버링
	// Unrefine : 셀 그룹핑
	
	if(rank==0) cout << "new cell numbering start" << endl;
	
	
	vector<int> cellGroupping;
	vector<int> cellLeveling;
	vector<vector<double>> cellVariables;
	vector<int> newCellNum(mesh.cells.size(),-1);
	newNumbering = 0;
	numN = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		int level = mesh.cells[i].level;
		int group = mesh.cells[i].group;
		
		// combine cells
		if(groupCellNewUnrefine.size()>0){
			if(i==groupCellNewUnrefine[numN][0]){
				
				// cout << "CELL VAR " << cell.var.size() << endl;
				
				vector<double> tmpVariables(cell.var.size(),0.0);
				double groupCellSize = (double)groupCellNewUnrefine[numN].size();
				for(int j=0; j<groupCellNewUnrefine[numN].size(); ++j){
					int iCell = groupCellNewUnrefine[numN][j];
					newCellNum[i] = newNumbering;
				// cout << "sub CELL VAR " << mesh.cells[iCell].var.size() << endl;
					
					for(int k=0; k<controls.nEq; ++k){
						double var = mesh.cells[iCell].var[k];
						tmpVariables[k] += var/groupCellSize;
					}
					
					++i;
				}
				
				i = i - 1;
				
				++numN;
				if(numN == groupCellNewUnrefine.size()) numN = 0;
				
				cellLeveling.push_back(level-1);
				cellGroupping.push_back(group);
				cellVariables.push_back(tmpVariables);
				++newNumbering;
				continue;
			}
		}
		
		// original cells
		newCellNum[i] = newNumbering;
		cellLeveling.push_back(level);
		cellGroupping.push_back(group);
		cellVariables.push_back(cell.var);
		
		++newNumbering;
		
	}
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// cout << newCellNum[i] << endl;
	// }
	
		
	
	//====================================================
	// Unrefine : 면 합침 , 면 - 포인트 연결
	
	

	if(rank==0) cout << "face combined start" << endl;
	
	vector<bool> deleteFaces(mesh.faces.size(),false);
	numN = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(boolInternalFaces[i]==true){
			if(boolFaceUnrefine[i]==true){
				deleteFaces[i] = true;
			}
			continue;
		}
		
		// combine faces (outer faces)
		if(groupFaceNewUnrefine.size()>0){
			int face0 = groupFaceNewUnrefine[numN][0];
			// cout << face0 << endl;
			if(i==face0){
				
				
				// cout << "AAAAAAAAA" << endl;

				for(int j=1; j<groupFaceNewUnrefine[numN].size(); ++j){
					int face1 = groupFaceNewUnrefine[numN][j];
					deleteFaces[face1] = true;
				}
				
				int centerPoint = -1;
				for(auto& point0 : mesh.faces[face0].points){
					int nCentP=0;
					for(int j=1; j<groupFaceNewUnrefine[numN].size(); ++j){
						int face1 = groupFaceNewUnrefine[numN][j];
						for(auto& point1 : mesh.faces[face1].points){
							if(point0==point1) ++nCentP;
						}
					}
					if(nCentP>=2){
						centerPoint = point0;
						break;
					}
				}
				
				if(centerPoint==-1){
					cout << endl;
					cout << "| #Error : No searching center point" << endl;
					for(int j=0; j<groupFaceNewUnrefine[numN].size(); ++j){
						int face1 = groupFaceNewUnrefine[numN][j];
						int level1 = mesh.faces[face1].level;
						int group1 = mesh.faces[face1].group;
						for(auto& point1 : mesh.faces[face1].points){
							int levelp = mesh.points[point1].level;
							cout << face1 << " " << point1 << " " << level1 << " " << group1 << " " << levelp << endl;
						}
					}
					cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				
				// cout << "BBBBBB" << endl;
				// cout << centerPoint << endl;
				// for(int j=0; j<groupFaceNewUnrefine[numN].size(); ++j){
					// int face1 = groupFaceNewUnrefine[numN][j];
					// for(auto& k : mesh.faces[face1].points){
						// cout << face1 << " " << k << endl;
					// }
					// // cout << mesh.faces[face1].points.size() << " " << centerPoint << " " << mesh.points.size() << endl;
				// }
				
				vector<int> tmpPoint;
				for(int j=0; j<1; ++j){
					int face1 = groupFaceNewUnrefine[numN][j];
					auto it0 = mesh.faces[face1].points.begin();
					auto it4 = mesh.faces[face1].points.end();
					auto it = std::find( it0, it4, centerPoint );
					auto it1 = it-1;
					tmpPoint.insert(tmpPoint.end(), it0, it1);
					// cout << std::distance(it0, it4) << endl;
					// cout << std::distance(it0, it) << endl;
				}
				
				// cout << "CCCCC" << endl;
				
					// cout << tmpPoint.size() << " " << centerPoint << endl;
				for(int j=1; j<groupFaceNewUnrefine[numN].size(); ++j){
					int face1 = groupFaceNewUnrefine[numN][j];
					// cout << "11" << endl;
					auto it0 = mesh.faces[face1].points.begin();
					// cout << "22" << endl;
					auto it4 = mesh.faces[face1].points.end();
					// cout << "33" << endl;
					auto it = std::find( it0, it4, centerPoint );
					// cout << "44" << endl;
					auto it1 = it-1;
					// cout << "55" << endl;
					auto it3 = it+1;
					// cout << "66" << endl;
					// cout << std::distance(it0,it) << " " << mesh.faces[face1].points.size() << endl;
					tmpPoint.insert(tmpPoint.end(), it3, it4);
					// cout << "77" << endl;
					tmpPoint.insert(tmpPoint.end(), it0, it1);
					// cout << "88" << endl;
				}
				
				// cout << "DDDD" << endl;
				
				
				for(int j=0; j<1; ++j){
					int face1 = groupFaceNewUnrefine[numN][j];
					auto it0 = mesh.faces[face1].points.begin();
					auto it4 = mesh.faces[face1].points.end();
					auto it = std::find( it0, it4, centerPoint );
					auto it3 = it+1;
					tmpPoint.insert(tmpPoint.end(), it3, it4);
				}
				// for(auto& j : tmpPoint){
					// cout << face0 << " " << j << endl;
				// }
				
				// cout << "EEEEE" << endl;
				
				
				
				tmpPoint.erase( std::remove_if( tmpPoint.begin(), tmpPoint.end(), 
					[&boolPointUnrefine](int const& n) { 
					return boolPointUnrefine[n]; 
					}), tmpPoint.end());
				
				mesh.faces[face0].points.clear();
				for(auto& j : tmpPoint){
					int newP = newPointNum[j];
					mesh.faces[face0].points.push_back(newP);
				}
				
				--mesh.faces[face0].level;
				
				if(mesh.faces[face0].points.size()<3){
					cout << endl;
					cout << "| #Error 2 : face's points size < 3" << endl;
					cout << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				// for(auto& j : mesh.faces[face0].points){
					// cout << face0 << " " << j << endl;
				// }
				
				// cout << "DDDDDDDD" << endl;
				
				i = i + groupFaceNewUnrefine[numN].size()-1;
				++numN;
				if(numN == groupFaceNewUnrefine.size()) numN = 0;
				
				// cout << "EEEEEEE" << endl;
				
				
				continue;
			}
		}
		
		// original faces
		vector<int> tmpPointsSave = face.points;
		face.points.erase( std::remove_if( face.points.begin(), face.points.end(), 
			[&boolPointUnrefine](int const& n) { 
			return boolPointUnrefine[n]; 
			}), face.points.end());
		if(face.points.size()<3){
			cout << endl;
			cout << "| #Error 3 : face's points size < 3" << endl;
			cout << endl;
			cout << i << " " << boolFaceUnrefine[i] << " " << boolInternalFaces[i] << " " << mesh.faces[i].level << endl;
			for(auto& j : tmpPointsSave){
				cout << j << " " << mesh.points[j].x << " " << mesh.points[j].y << " " << mesh.points[j].z << " " << boolPointUnrefine[j] << endl;
				cout << mesh.points[j].level << endl;
			}
			MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		for(auto& j : face.points){
			j = newPointNum[j];
		}
		
	}
	
		
	
	// cout << "666666" << endl;
	//====================================================
	// Unrefine : 포인트 삭제
	
	
	// SEMO_Mesh_Builder newMesh;
	
	// for(int i=0; i<mesh.points.size(); ++i){
		// auto& point = mesh.points[i];
		// if(boolPointUnrefine[i]==false){
			// newMesh.addPoint();
			// newMesh.points.back().x = point.x;
			// newMesh.points.back().y = point.y;
			// newMesh.points.back().z = point.z;
			// newMesh.points.back().level = point.level;
		// }
	// }
	
	
	if(rank==0) cout << "point erase start" << endl;
	
	
	
	numN = 0;
	mesh.points.erase( std::remove_if( mesh.points.begin(), mesh.points.end(), 
		[&boolPointUnrefine, &numN](SEMO_Point const& v) { 
		return boolPointUnrefine[numN++]; 
		}), mesh.points.end());
		
	
	
	
		
	//====================================================
	// Unrefine : 바운더리 셋팅
	
	
	if(rank==0) cout << "boundary setting start" << endl;
	
	vector<int> BCFacestart;
	vector<int> BCFaceNfaces;
	int BCnum = 0;
	numN = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		if(BCnum < mesh.boundary.size()){
			if( mesh.boundary[BCnum].startFace == i ){
				BCFacestart.push_back(numN);
				if(BCnum != 0) {
					BCFaceNfaces.push_back(numN-BCFacestart[BCnum-1]);
				}
				++BCnum;
			}
		}
		
		if(deleteFaces[i] == false){
			++numN;
		}
		
	}
	BCFaceNfaces.push_back(numN-BCFacestart[BCnum-1]);
	
	
	numN = 0;
	for(auto& boundary : mesh.boundary){
		boundary.startFace = BCFacestart[numN];
		boundary.nFaces = BCFaceNfaces[numN];
		++numN;
	}
	
	// cout << "777777" << endl;

	//====================================================
	// Unrefine : 면 삭제 , 면 - 셀
	
	if(rank==0) cout << "face erase, face-cell connect start" << endl;
	
	
	numN = 0;
	mesh.faces.erase( std::remove_if( mesh.faces.begin(), mesh.faces.end(), 
		[&deleteFaces, &numN](SEMO_Face const& v) { 
		return deleteFaces[numN++]; 
		}), mesh.faces.end());
		
		
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		face.var.resize(controls.nTotalFaceVar,0.0);
		
		face.owner = newCellNum[face.owner];
		if(face.getType()==SEMO_Types::INTERNAL_FACE){
			face.neighbour = newCellNum[face.neighbour];
		}
		else{
			face.neighbour = -1;
		}
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	

	// cout << "888888" << endl;
	
	
	//====================================================
	
	if(rank==0) cout << "mesh clear, setting start" << endl;
	
	
	// 원래 메쉬 클리어 및 복사
	mesh.cells.clear();
	

	mesh.check();
	
	
	
	// mesh.checkMatchingProcessorFace();
	
	
	
	mesh.buildCells();
	
	// mesh.setFaceTypes();

	mesh.buildLists();
	
	mesh.connectCelltoFaces();
	
	
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		// int minLevel = 999999999;
		// for(int j=0; j<cell.faces.size(); ++j){
			// auto& face = mesh.faces[cell.faces[j]];
			// minLevel = min(minLevel,face.level);
		// }
		// cell.level = minLevel;
		
		cell.level = cellLeveling[i];
		cell.group = cellGroupping[i];
		
		cell.var.resize(controls.nTotalCellVar,0.0); 
		int tmpNum = 0;
		for(auto& k : cellVariables[i]){
			cell.var[tmpNum] = k;
			++tmpNum;
		}
	}
	
	
	
	mesh.connectCelltoPoints();
	
	// set processor face counts
	mesh.setCountsProcFaces();
	
	// set processor face displacements
	mesh.setDisplsProcFaces(); 
		
	// mesh.informations();
	
	
	// cout << "99999" << endl;
	
	// SEMO_Mesh_Save save;
	// string tmpFile = "./Uf" + to_string(iter);
	// // string tmpFile = "./";
	// save.vtu(tmpFile, rank, mesh);
	
	
	
	if(rank==0) cout << "----Unrefine finished-----" << endl;

	// cout << "SIZE : " << cellGroupping.size() << " " << cellLeveling.size() << " " << mesh.cells.size() << endl;
		
	// cout << "CELL" << endl;
	// for(int i=0; i<mesh.cells.size(); ++i){
		// auto& cell = mesh.cells[i];
		// cout << i << " " << cell.level << " " << cell.group << endl;
	// }
	
	// cout << "FACE" << endl;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		// cout << i << " " << face.level << " " << face.group << endl;
	// }
		

	// for(int i=0; i<mesh.cells.size(); ++i){
		// cout << mesh.cells[i].var[controls.VF[0]] << endl;
		// // mesh.cells[i].var[controls.U] = 0.0;
		// // mesh.cells[i].var[controls.V] = 0.0;
		// // mesh.cells[i].var[controls.W] = 0.0;
		// // mesh.cells[i].var[controls.VF[0]] = 1.0;
		
	// }
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	

	
	
}