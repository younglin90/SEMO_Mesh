#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <random>
#include <set>

#include "build.h"
#include "mpi.h"
#include "polyAMR.h"
#include "geometric.h"


class faces_Unrefind{

private:

public:
	int id;
	bool boolCombine;
	vector<int> child;
	vector<int> points;
	int owner;
	int neighbour;
};



class groupCells_Unrefine{

private:

public:
	vector<int> child;
	// vector<int> outChildFaces;
	vector<int> intChildFaces;
	vector<int> groupOutFaces_id;
	vector<int> faces;
	vector<int> points;
	vector<int> owner;
	int neighbour;
	int cellCenterPoint;
	bool boolCombine;
	vector<int> faceCenterPoints;
};




// void SEMO_Poly_AMR_Builder::sortFaceCanUnrefine(
	// vector<int>& vLevel, 
	// int saveLevel, 
	// int& num, 
	// vector<vector<int>>& vGroupLevel, 
	// vector<vector<int>>& vGroupNumber
	// ) {
		
	// // vector<int> tmpLevels;
	// // vector<int> tmpNumbers;
	// // int limitNum = 3;
	// // if (saveLevel == 0) limitNum = 100;
	// // if (saveLevel == 1) limitNum = 100;
	// // while (1) {
		// // if (vLevel.size() - 1 >= num) {
			// // int level = vLevel[num];
			// // // cout << num << " " << tmpLevels.size()<< " " << saveLevel << " " << level << endl;

			// // if (
				// // saveLevel == level &&
				// // tmpLevels.size() < limitNum
				// // ) {
				// // tmpLevels.push_back(saveLevel);
				// // tmpNumbers.push_back(num);
			// // }
			// // else if (
				// // saveLevel == level &&
				// // tmpLevels.size() == limitNum) {
				// // if (std::find(tmpNumbers.begin(), tmpNumbers.end(), -1) == tmpNumbers.end()) {
					// // tmpLevels.push_back(saveLevel);
					// // tmpNumbers.push_back(num);
					// // vGroupLevel.push_back(tmpLevels);
					// // vGroupNumber.push_back(tmpNumbers);
				// // }
				// // break;
			// // }
			// // else if (
				// // saveLevel < level
				// // ) {
				// // // cout << "start" << endl;

				// // if (tmpLevels.size() == limitNum+1) {
					// // tmpLevels.clear();
					// // tmpNumbers.clear();
				// // }

				// // sortFaceCanUnrefine(vLevel, saveLevel + 1, num, vGroupLevel, vGroupNumber);
				// // tmpLevels.push_back(-1);
				// // tmpNumbers.push_back(-1);
				// // // cout << tmpLevels.size() << endl;
			// // }
			// // else {
				// // --num;
				// // break;
			// // }
			// // ++num;
		// // }
		// // else {
			// // break;
		// // }
	// // }
   
// }







void sortCellCanUnrefine(
	vector<int>& vLevel, 
	int saveLevel, 
	int& num, 
	vector<vector<int>>& vGroupLevel, 
	vector<vector<int>>& vGroupNumber
	) {
		
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
		
	vector<int> tmpLevels;
	vector<int> tmpNumbers;
	int limitNum = 7;
	if (saveLevel <= 1) limitNum = 2147483647;
	while (1) {
		if (vLevel.size() - 1 >= num) {
			int level = vLevel.at(num);
			// cout << num << " " << tmpLevels.size()<< " " << saveLevel << " " << level << endl;

			if(
			saveLevel == level &&
			tmpLevels.size() < limitNum
			){
				tmpLevels.push_back(saveLevel);
				tmpNumbers.push_back(num);
			}
			else if(
			saveLevel == level &&
			tmpLevels.size() == limitNum
			){
				if (std::find(tmpNumbers.begin(), tmpNumbers.end(), -1) == tmpNumbers.end()) {
					tmpLevels.push_back(saveLevel);
					tmpNumbers.push_back(num);
					vGroupLevel.push_back(tmpLevels);
					vGroupNumber.push_back(tmpNumbers);
				}
				break;
			}
			else if(
			saveLevel < level
			){
				// cout << "start" << endl;

				if (tmpLevels.size() == limitNum+1) {
					tmpLevels.clear();
					tmpNumbers.clear();
				}

				sortCellCanUnrefine(vLevel, saveLevel + 1, num, vGroupLevel, vGroupNumber);
				tmpLevels.push_back(-1);
				tmpNumbers.push_back(-1);
				
				if (tmpLevels.size() == limitNum + 1) {
					break;
				}
				// cout << tmpLevels.size() << endl;
			}
			else {
				--num;
				break;
			}
			++num;
		}
		else {
			
			if(saveLevel>1 && tmpLevels.size()!=8){
				cout << "WWWWWWWWWWWWWWW" << " " << rank << " " << saveLevel << " " << tmpLevels.size() << " " << vLevel.size() << " " << num << endl;
				
				cout << rank << " ";
				for(auto& k : vLevel){
					cout << k << " ";
				}
				cout << endl;
				
				MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
			}
			
			break;
		}
	}
   
}






void restrictCellUnrefine(
	SEMO_Mesh_Builder& mesh, 
	vector<bool>& boolCellUnrefine,
	vector<int>& cLevel_recv, 
	vector<int>& cUnrefine_recv){

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
}



void extractGroupCellListsCanUnrefine(
	SEMO_Mesh_Builder& mesh, 
	vector<vector<int>>& groupCellListsCanUnrefine
){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	vector<vector<int>> groupCellListsLevelZero;
	vector<vector<int>> groupCellLevelListsLevelZero;
	for(int i=0; i<mesh.cells.size(); ++i){
		
		groupCellListsLevelZero.push_back(vector<int>(1,i));
		groupCellLevelListsLevelZero.push_back(vector<int>(1,mesh.cells[i].level));
		
		if(mesh.cells[i].level<=0) continue;
		
		int groupOrg = mesh.cells[i].group;
		
		++i;
		while(groupOrg==mesh.cells[i].group && i<mesh.cells.size()){
			groupCellListsLevelZero.back().push_back(i);
			groupCellLevelListsLevelZero.back().push_back(mesh.cells[i].level);
			++i;
		}
		--i;
	}
	
	
	int testNum = 0;
	// for(auto& j : groupCellListsLevelZero[0]){
		// // if(rank==0) cout << j << " ";
		// ++testNum;
	// }
			
	
	
	
	vector<int> test;
	for(int i=0; i<groupCellListsLevelZero.size(); ++i){
		// cout << i << endl;
		
		
		
		
		if(groupCellLevelListsLevelZero[i].size()==1) continue;
		
		
		
		
		bool boolAllLevelOne = true;
		for(auto& j : groupCellLevelListsLevelZero[i]){
			if(j != 1) {
				boolAllLevelOne = false;
				break;
			}
		}
		
		
		
		if(boolAllLevelOne==true){
			test.push_back(i);
			groupCellListsCanUnrefine.push_back(groupCellListsLevelZero[i]);
		}
		else{
			vector<vector<int>> tmpGroupLevel;
			vector<vector<int>> tmpGroupNumber;
			int tmpNum2 = 0;
			
				
		// cout << " 1 : " << groupCellListsLevelZero[0].size() << endl;
		
	// testNum = 0;
	// for(auto& j : groupCellListsLevelZero[0]){
		// if(rank==0) cout << j << " ";
		// ++testNum;
	// }	
	// cout << endl;
			
			// cout << rank << " " << groupCellListsLevelZero[i][0] << " " << i << endl;
			
			sortCellCanUnrefine(
				groupCellLevelListsLevelZero[i],
				0,
				tmpNum2,
				tmpGroupLevel,
				tmpGroupNumber
				);
		
		// cout << " 2 : " << groupCellListsLevelZero[0].size() << endl;
			
	// testNum = 0;
	// for(auto& j : groupCellListsLevelZero[0]){
		// if(rank==0) cout << j << " ";
		// ++testNum;
	// }	
	// cout << endl;
			
			int startCellNum = groupCellListsLevelZero[i][0];
			for(auto& j : tmpGroupNumber){
				for(auto& k : j){
					k += startCellNum;
				}
			}
			
			for(auto& j : tmpGroupNumber){
				test.push_back(i);
				groupCellListsCanUnrefine.push_back(j);
			}
		}
		
	}
	
	
	
	
	for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){

		auto& group = groupCellListsCanUnrefine[i];
		
	// if(rank==0) cout << "11111111111" << endl;
		// search cell center point
		auto& targetCell = mesh.cells[group[0]];
		int cellCenterPoint = -1;
		for(auto& targetPoint : targetCell.points){
			bool boolCellCenterPoint = true;
			for(int j=1; j<group.size(); ++j){
				auto& cell = mesh.cells[group[j]];
				if(
				std::find(cell.points.begin(),cell.points.end(),targetPoint)
				== cell.points.end()
				){
					boolCellCenterPoint = false;
					break;
				}
			}
			if(boolCellCenterPoint==true){
				cellCenterPoint = targetPoint;
				break;
			}
		}
		
		
		if(cellCenterPoint==-1){
			
			
			
			
			for(int j=0; j<group.size(); ++j){
				auto& cell = mesh.cells[group[j]];
				if(rank==0) cout << cell.level << endl;
				for(auto& k : cell.points){
					if(rank==0) cout << rank << " " << group[j] << " " << k << endl;
				}
				
			}
			cout << endl;
			
			if(rank==0) cout << "cellLevel" << endl;
			
			
			testNum = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << groupCellListsLevelZero[test[i]][testNum] << " " << groupCellLevelListsLevelZero[test[i]][testNum] << " ";
				++testNum;
			}
			
			cout << endl;
			
			testNum = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << j << " ";
				++testNum;
			}
			cout << endl;
			
			testNum = 0;
			for(auto& j : groupCellListsLevelZero[0]){
				if(rank==0) cout << j << " ";
				++testNum;
			}	
			cout << endl;
			
			testNum = 0;
			int saveFirst = groupCellListsLevelZero[test[i]][0];
			int saveLast = 0;
			if(rank==0) cout << groupCellListsLevelZero[test[i]][0] << " " << groupCellListsLevelZero[test[i]].size() << endl;
			for(auto& j : groupCellListsLevelZero[test[i]]){
				if(rank==0) cout << mesh.cells[j].group << " ";
				saveLast = j;
				++testNum;
			}
			
			
			cout << endl;
			if(rank==0) cout << mesh.cells[saveFirst-1].group << endl;
			if(rank==0) cout << mesh.cells[saveLast+1].group << endl;
			if(rank==0) cout << mesh.cells.size() << endl;
			
			for(int j=0; j<mesh.cells.size(); ++j){
				if(rank==0) cout << mesh.cells[j].level << " ";
			}
			
			if(rank==0) MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		}
		
		
		
	}
	
	
	
	
	
	
	
	
	
	
}





void extractGroupUnrefineCells(
	SEMO_Mesh_Builder& mesh, 
	vector<bool>& boolCellUnrefine,
	vector<vector<int>>& groupCellListsCanUnrefine,
	vector<groupCells_Unrefine>& groupCellsUnrefine
	){

	for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		bool boolAllUnrefine = true;
		for(auto& j : groupCellListsCanUnrefine[i]){
			if(boolCellUnrefine[j] == false) {
				boolAllUnrefine = false;
				break;
			}
		}
		
		if(boolAllUnrefine){
			groupCellsUnrefine.push_back(groupCells_Unrefine());
			for(auto& j : groupCellListsCanUnrefine[i]){
				groupCellsUnrefine.back().child.push_back(j);
			}
		}
	}
}











void SEMO_Poly_AMR_Builder::polyUnrefine(
	SEMO_Mesh_Builder& mesh, 
	SEMO_Controls_Builder& controls,
	int iter){

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int proc_num = 0;
	
	// cout << rank << " CCCCC " << mesh.cells.size() << endl;
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// for(auto& j : mesh.faces[i].points){
			// if(rank==0) cout << i << " " << mesh.faces[i].owner << " " << mesh.faces[i].neighbour << " " << j << endl;
		// }
	// }
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// ++proc_num;
		// }
	// }
	// cout << "11111 : " << rank << " " << proc_num << " " << mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1] << endl;
	
	
	// for(auto& boundary : mesh.boundary){
		// if(rank==0) cout << boundary.startFace << " " << boundary.nFaces << endl;
	// }
	
	

	if(rank==0) cout << "┌────────────────────────────────────────────────────" << endl;
	if(rank==0) cout << "| execute AMR - Unrefine " << endl;
	
	
	SEMO_Mesh_Geometric geometric;
	
	//====================================================
	// 셀 Unrefine : Unrefine 되는 셀 & 면 조사
	
	// random
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0.0, 1.0);
	
	

	// SEMO_Utility_Math math;
	// vector<vector<double>> gradVF;
	// // math.calcLeastSquare2nd(mesh, controls.VF[0], controls.fVF[0], gradVF);
	// math.calcGaussGreen(mesh, controls.VF[0], controls.fVF[0], gradVF);
	
	// vector<bool> boolCanNotUnrefineCells(mesh.cells.size(),false);
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(mesh.cells[i].level < 0) boolCanNotUnrefineCells[i] = true;
	// }
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].level < 0) mesh.faces[i].level = 0;
	// }
	
	
	double indicatorUnrefine_AMR = controls.indicatorUnrefine;
	
	vector<bool> boolCellUnrefine(mesh.cells.size(),false);
	for(int i=0; i<mesh.cells.size(); ++i){
		// if( distr(eng) < 0.9 ){
			// boolCellUnrefine[i] = true;
		// }
		// boolCellUnrefine[i] = true;
		// boolCellUnrefine[i] = false;
		
		// if( mesh.cells[i].var[controls.indicatorAMR] < 0.5 ){
			// boolCellUnrefine[i] = true;
		// }
		
			// boolCellUnrefine[i] = true;
		

		if(mesh.cells[i].var[controls.indicatorAMR[0]] < indicatorUnrefine_AMR) 
			boolCellUnrefine[i] = true;
		
		
		// 만약 셀의 레벨이 0 이면, false
		if(mesh.cells[i].level == 0) boolCellUnrefine[i] = false;
		
		
		if(mesh.cells[i].level < 0) boolCellUnrefine[i] = false;
		// if(boolCanNotUnrefineCells[i]==true) {
			// boolCellUnrefine[i] = false;
			// mesh.cells[i].level = 0;
		// }
		
		
	}
	
	
		// boolCellUnrefine[0] = true;
		// boolCellUnrefine[1] = true;
		// boolCellUnrefine[2] = true;
		// boolCellUnrefine[3] = true;
		// boolCellUnrefine[4] = true;
		// boolCellUnrefine[5] = true;
		// boolCellUnrefine[6] = true;
		// boolCellUnrefine[7] = true;
	
	vector<int> cLevel_recv;
	vector<int> cUnrefine_recv;
	this->mpiLevelRefine(mesh, boolCellUnrefine, cLevel_recv, cUnrefine_recv);
	
	restrictCellUnrefine(mesh, boolCellUnrefine, cLevel_recv, cUnrefine_recv);
	
	this->mpiRefines(mesh, boolCellUnrefine, cUnrefine_recv);
	
	
	
	//====================================================
	// cout << rank << " : cell groupping start" << endl;
	
	// 오리지널 그룹 
	// int saveGroup = mesh.cells[0].group;
	
	vector<vector<int>> groupCellListsCanUnrefine;
	
	extractGroupCellListsCanUnrefine(mesh, groupCellListsCanUnrefine);
	
	
	//====================================================
	
	vector<groupCells_Unrefine> groupCellsUnrefine;
	// vector<vector<int>> groupCellsUnrefine;
	extractGroupUnrefineCells(mesh, boolCellUnrefine, groupCellListsCanUnrefine, 
		groupCellsUnrefine);
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(rank==0) cout << rank << " " << i << " " << mesh.cells[i].level << " " << mesh.cells[i].group << endl;
	// }
	// if(rank==0) cout << "AA : " << groupCellListsCanUnrefine.size() << endl;
	// for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){
		// for(auto& j : groupCellListsCanUnrefine[i]){
			// if(rank==0) cout << rank << " " << i << " " << j << endl;
		// }
	// }
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(rank==1) cout << rank << " " << mesh.cells[i].level << " " << boolCellUnrefine[i] << endl;
	// }
	
	
	std::fill(boolCellUnrefine.begin(), boolCellUnrefine.end(), false);
	vector<int> groupCells_id(mesh.cells.size(),-1);
	for(int i=0; i<groupCellsUnrefine.size(); ++i){
		for(auto& j : groupCellsUnrefine[i].child){
			boolCellUnrefine[j] = true;
			groupCells_id[j] = i;
		}
	}
	
	
	// for(int i=0; i<mesh.cells.size(); ++i){
		// if(rank==0) cout << rank << " " << mesh.cells[i].level << " " << boolCellUnrefine[i] << endl;
	// }
	
	
	this->mpiRefines(mesh, boolCellUnrefine, cUnrefine_recv);
	
	
	
	// for(int i=0; i<groupCellsUnrefine.size(); ++i){
		// // int sumGroup = 0;
		// // for(auto& j : groupCellsUnrefine[i].child){
			// // sumGroup += mesh.cells[j].group;
		// // }
		
		// // if(
		// // (groupCellsUnrefine[i].child.size() != 8) ||
		// // (mesh.cells[groupCellsUnrefine[i].child[0]].group*8 != sumGroup)
		// // ){
			// // for(auto& j : groupCellsUnrefine[i].child){
				// // if(rank==0) cout << i << " " << j << " " << mesh.cells[j].level << " " << mesh.cells[j].group << " " << endl;
			// // }
		// // }
		// auto& group = groupCellsUnrefine[i];
		
	// // if(rank==0) cout << "11111111111" << endl;
		// // search cell center point
		// auto& targetCell = mesh.cells[group.child[0]];
		// int cellCenterPoint = -1;
		// for(auto& targetPoint : targetCell.points){
			// bool boolCellCenterPoint = true;
			// for(int j=1; j<group.child.size(); ++j){
				// auto& cell = mesh.cells[group.child[j]];
				// if(
				// std::find(cell.points.begin(),cell.points.end(),targetPoint)
				// == cell.points.end()
				// ){
					// boolCellCenterPoint = false;
					// break;
				// }
			// }
			// if(boolCellCenterPoint==true){
				// cellCenterPoint = targetPoint;
				// break;
			// }
		// }
		
		
		// if(cellCenterPoint==-1){
			
			// for(int j=0; j<group.child.size(); ++j){
				// auto& cell = mesh.cells[group.child[j]];
				// cout << cell.level << endl;
				// for(auto& k : cell.points){
					// cout << j << " " << k << endl;
				// }
				
			// }
			
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		
		
	// }
	
	
	
	// for(int i=0; i<groupCellListsCanUnrefine.size(); ++i){

		// auto& group = groupCellListsCanUnrefine[i];
		
	// // if(rank==0) cout << "11111111111" << endl;
		// // search cell center point
		// auto& targetCell = mesh.cells[group[0]];
		// int cellCenterPoint = -1;
		// for(auto& targetPoint : targetCell.points){
			// bool boolCellCenterPoint = true;
			// for(int j=1; j<group.size(); ++j){
				// auto& cell = mesh.cells[group[j]];
				// if(
				// std::find(cell.points.begin(),cell.points.end(),targetPoint)
				// == cell.points.end()
				// ){
					// boolCellCenterPoint = false;
					// break;
				// }
			// }
			// if(boolCellCenterPoint==true){
				// cellCenterPoint = targetPoint;
				// break;
			// }
		// }
		
		
		// if(cellCenterPoint==-1){
			
			// for(int j=0; j<group.size(); ++j){
				// auto& cell = mesh.cells[group[j]];
				// if(rank==0) cout << cell.level << endl;
				// for(auto& k : cell.points){
					// if(rank==0) cout << rank << " " << j << " " << k << endl;
				// }
				
			// }
			
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		// }
		
		
		
	// }
	
			// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// if( boolCellUnrefine[i] ){
			// cout << rank << " CELL UNREFINE : " << i << endl;
		// }
	// } 
	
	
	
	
	
	//====================================================
	// 셀 넘버링
	
	vector<int> newCellsNumber(mesh.cells.size(),-1);
	vector<int> cellsLevel(mesh.cells.size(),-1);
	vector<int> cellsGroup(mesh.cells.size(),-1);
	int newCellNum = 0;
	for(int i=0; i<mesh.cells.size(); ++i){
		if(groupCells_id[i] == -1){
			cellsLevel[newCellNum] = mesh.cells[i].level;
			// if(boolCanNotUnrefineCells[i]==true) cellsLevel[newCellNum] = -1;
			cellsGroup[newCellNum] = mesh.cells[i].group;
			
			newCellsNumber[i] = newCellNum;
			
			++newCellNum;
		}
		else{
			cellsLevel[newCellNum] = mesh.cells[i].level-1;
			// if(boolCanNotUnrefineCells[i]==true) cellsLevel[newCellNum] = -1;
			cellsGroup[newCellNum] = mesh.cells[i].group;
			
			for(auto& j : groupCellsUnrefine[groupCells_id[i]].child){
				newCellsNumber[j] = newCellNum;
				++i;
			}
			--i;
			
			++newCellNum;
		}
	}
	
	int totalCellNum = newCellNum;
	
	
	//====================================================
	// 포인트 삭제 및 넘버링

	vector<bool> boolDeletePoints(mesh.points.size(),true);
	for(int i=0; i<mesh.cells.size(); ++i){
		auto& cell = mesh.cells[i];
		int cellLevel = cell.level;
		
		
		if(cellLevel==-1) cellLevel=0;
		
		
		if(boolCellUnrefine[i] == true){
			for(auto& j : cell.points){
				int pointLevel = mesh.points[j].level;
				if(pointLevel <= cellLevel-1){
					boolDeletePoints[j] = false;
				}
			}
		}
		else{
			for(auto& j : cell.points){
				int pointLevel = mesh.points[j].level;
				if(pointLevel <= cellLevel)
					boolDeletePoints[j] = false;
			}
		}
	}
	
	// proc
	proc_num = 0;
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			auto& face = mesh.faces[i];
			int faceLevel = face.level;
			int ownLevel = mesh.cells[face.owner].level;
			int ngbLevel = cLevel_recv[proc_num];
			bool ownUnrefine = boolCellUnrefine[face.owner];
			bool ngbUnrefine = cUnrefine_recv[proc_num];
			
			
			if(faceLevel==-1) faceLevel=0;
			if(ownLevel==-1) ownLevel=0;
			if(ngbLevel==-1) ngbLevel=0;
			
			
			if(
			(ownLevel<ngbLevel && ngbUnrefine==false) ||
			(ownLevel==ngbLevel && ownUnrefine==true && ngbUnrefine==false) 
			){
				for(auto& j : face.points){
					int pointLevel = mesh.points[j].level;
					if(pointLevel <= faceLevel){
						boolDeletePoints[j] = false;
					}
					// else{
						// cout << "CCC " << faceLevel << " " << pointLevel << endl;
						// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
					// }
				}
			}
			++proc_num;
		}
	}
	
	
	
	vector<int> newPointsNumber(mesh.points.size(),-1);
	int newPointNum = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		if(boolDeletePoints[i] == false){
			newPointsNumber[i] = newPointNum;
			++newPointNum;
		}
	}
	
	
	
	//====================================================
	// Unrefine : 면 합침 , 면 - 포인트 연결
	// if(rank==0) cout << "face combined start" << endl;
	
	// cout << rank << " : face combined start" << endl;
	
	int proc_total_num = 0;
	
	vector<bool> boolIntFaces_Deleted(mesh.faces.size(),false);
	vector<int> groupFaces_id(mesh.faces.size(),-1);
	vector<faces_Unrefind> groupOutFaces;
	for(int i=0; i<groupCellsUnrefine.size(); ++i){
		auto& group = groupCellsUnrefine[i];
		
	// if(rank==0) cout << "11111111111" << endl;
		// search cell center point
		auto& targetCell = mesh.cells[group.child[0]];
		group.cellCenterPoint = -1;
		for(auto& targetPoint : targetCell.points){
			bool boolCellCenterPoint = true;
			for(int j=1; j<group.child.size(); ++j){
				auto& cell = mesh.cells[group.child[j]];
				if(
				std::find(cell.points.begin(),cell.points.end(),targetPoint)
				== cell.points.end()
				){
					boolCellCenterPoint = false;
					break;
				}
			}
			if(boolCellCenterPoint==true){
				group.cellCenterPoint = targetPoint;
				break;
			}
		}
		
		
	// if(rank==0) cout << "222222222" << endl;
		// search cell internal & outer faces
		int cellCenterPoint = group.cellCenterPoint;
		vector<int> outChildFaces;
		set<int> intChildFaces;
		for(auto& j : group.child){
			auto& cell = mesh.cells[j];
			for(auto& k : cell.faces){
				auto& face = mesh.faces[k];
				if(
				std::find(face.points.begin(),face.points.end(),cellCenterPoint)
				== face.points.end()
				){
					outChildFaces.push_back(k);
				}
				else{
					boolIntFaces_Deleted[k] = true;
					groupFaces_id[k] = i;
					intChildFaces.insert(k);
				}
			}
		}
		
		group.intChildFaces.assign(intChildFaces.begin(),intChildFaces.end());
		
		// group.intChildFaces.erase(
			// std::unique(group.intChildFaces.begin(),group.intChildFaces.end()),
			// group.intChildFaces.end());
		
		// for(auto& j : group.intChildFaces){
			// if(rank==0) cout << "151515 " << j << endl;
		// }
		
	// if(rank==0) cout << "3333333 " << group.child.size() << " " << 
				// group.intChildFaces.size() << " " << outChildFaces.size() << endl;
				
		// search face center points
		set<int> tmpFaceCenterPoints;
		for(int j=0; j<group.intChildFaces.size(); ++j){
			auto& face0 = mesh.faces[group.intChildFaces[j]];
			for(int k=j+1; k<group.intChildFaces.size(); ++k){
				auto& face1 = mesh.faces[group.intChildFaces[k]];
				for(auto& p0 : face0.points){
					for(auto& p1 : face1.points){
						if(p0==p1){
							tmpFaceCenterPoints.insert(p0);
						}
					}
				}
			}
		}
		
	// if(rank==0) cout << "4444444 " << group.intChildFaces.size() << endl;
		// tmpFaceCenterPoints.erase(
			// std::find(tmpFaceCenterPoints.begin(),
				// tmpFaceCenterPoints.end(),cellCenterPoint));
		
		// tmpFaceCenterPoints.erase(
			// std::unique(tmpFaceCenterPoints.begin(),
				// tmpFaceCenterPoints.end()),tmpFaceCenterPoints.end());
		
	// if(rank==0) cout << cellCenterPoint << endl;
		if(tmpFaceCenterPoints.find(cellCenterPoint)==tmpFaceCenterPoints.end()){
			cout << "| #WARNING : NO searching cellCenterPoint " << rank << " " << group.cellCenterPoint << endl;
		}
		
		tmpFaceCenterPoints.erase(tmpFaceCenterPoints.find(cellCenterPoint));
		group.faceCenterPoints.assign(tmpFaceCenterPoints.begin(),tmpFaceCenterPoints.end());
		
	// if(rank==0) cout << "5555555555 " << group.faceCenterPoints.size() << endl;
		
		// groupping outer faces, save points, own, ngb
		for(auto& centP : group.faceCenterPoints){
	// if(rank==0) cout << "1-1" << endl;
			vector<int> tmpChildFaces;
			for(auto& j : outChildFaces){
				auto& face = mesh.faces[j];
				if(
				std::find(face.points.begin(),face.points.end(),centP)
				!= face.points.end()
				){
					tmpChildFaces.push_back(j);
				}
			}
			
			
			
			// if(rank==0) cout << "1-2" << " " << tmpChildFaces.size() << " " << endl;
			// if(mesh.faces[ tmpChildFaces[0] ].getType() == SEMO_Types::PROCESSOR_FACE){
				// if(rank==0) cout << "PROC" << endl;
			// }
			
			
	// if(rank==0) cout << "1-2" << " " << tmpChildFaces.size() << " " << tmpFaceCenterPoints.size()<<endl;
			int face0 = tmpChildFaces[0];
			if(groupFaces_id[face0]==-1){
	// if(rank==0) cout << "1-3" << endl;
	
				// for(auto& j : tmpChildFaces){
					// auto& face = mesh.faces[j];
					// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
						// ++proc_total_num;
					// }
				// }
				
				
				
					
				for(auto& j : tmpChildFaces){
					groupFaces_id[j] = groupOutFaces.size();
				}
				group.groupOutFaces_id.push_back(groupOutFaces.size());
				groupOutFaces.push_back(faces_Unrefind());
				groupOutFaces.back().child = tmpChildFaces;
				
				// points
				std::sort(groupOutFaces.back().child.begin(), 
						  groupOutFaces.back().child.end());
						  
				vector<vector<int>> vertexLists;
				for(auto& j : groupOutFaces.back().child){
					auto& face = mesh.faces[j];
					int faceLevel = face.level;
					vertexLists.push_back(vector<int>());
					for(auto& p0 : face.points){
						auto& point = mesh.points[p0];
						int pointLevel = point.level;
						if(pointLevel<=faceLevel){
							vertexLists.back().push_back(p0);
						}
					}
				}
				
				
				// cout << "00000000000000" << endl;
				vector<int> groupFacePoints0;
				int tmptmp = 0;
				for(auto& j : groupOutFaces.back().child){
					auto& face = mesh.faces[j];
					auto& points = face.points;
					
					if(tmptmp!=0){
						
						auto iStr = std::find(points.begin(),points.end(),vertexLists[tmptmp][3]);
						
						// std::copy( iStr, points.end(), groupFacePoints0.end()-1 );
						groupFacePoints0.insert( groupFacePoints0.end(), iStr, points.end() );
						
					}
					
					auto iEnd = std::find(points.begin(),points.end(),vertexLists[tmptmp][1]);
					// cout << "1-1" << endl;
					groupFacePoints0.insert( groupFacePoints0.end(), points.begin(), iEnd );
					// cout << "1-2" << endl;
					
					
					++tmptmp;
				}
				// cout << "111111111111111111" << endl;
				auto& faceZero = mesh.faces[groupOutFaces.back().child[0]];
				auto& pointsZero = faceZero.points;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists[0][3]);
				// std::copy( iStr0, pointsZero.end(), groupFacePoints0.end()-1 );
				groupFacePoints0.insert( groupFacePoints0.end(), iStr0, pointsZero.end() );
				// cout << "222222222222222222" << endl;
				
				vector<int> groupFacePoints;
				for(auto& j : groupFacePoints0){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						groupFacePoints.push_back(newPN);
					}
				}
							
						  
				// vector<int> groupFacePoints;
				// for(auto& j : groupOutFaces.back().child){
					// auto& face = mesh.faces[j];
					// int faceLevel = face.level;
					// for(auto& p0 : face.points){
						// auto& point = mesh.points[p0];
						// int pointLevel = point.level;
						// if(pointLevel<=faceLevel-1){
							// int newPN = newPointsNumber[p0];
							// groupFacePoints.push_back(newPN);
						// }
					// }
				// }
				
				
	// if(rank==0) cout << "1-4" << endl;
				groupOutFaces.back().points = groupFacePoints;
				
				// own, ngb
				groupOutFaces.back().owner = newCellsNumber[mesh.faces[face0].owner];
				if(mesh.faces[face0].neighbour==-1){
					groupOutFaces.back().neighbour = -1;
				}
				else{
					groupOutFaces.back().neighbour = newCellsNumber[mesh.faces[face0].neighbour];
				}
				
	// if(rank==0) cout << "1-5" << endl;
				
			}
			else{
	// if(rank==0) cout << "1-6" << endl;
				int group_id = groupFaces_id[face0];
				group.groupOutFaces_id.push_back(group_id);
	// if(rank==0) cout << "1-7" << endl;
			}
			
		}
	// if(rank==0) cout << "6666666666" << endl;
	}
	
	
	
	
	// cout << rank << " : 222222222" << endl;
	
	

	// proc
	vector<bool> boolProcFaces_Combine(proc_num,false);
	vector<int> groupProcFaces_id(proc_num,-1);
	vector<faces_Unrefind> groupProcFaces;
	
	
	proc_num = 0;
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			auto& face = mesh.faces[i];
			int faceLevel = face.level;
			int ownLevel = mesh.cells[face.owner].level;
			int ngbLevel = cLevel_recv[proc_num];
			bool ownUnrefine = boolCellUnrefine[face.owner];
			bool ngbUnrefine = cUnrefine_recv[proc_num];
			if(
			(ownLevel<ngbLevel && ngbUnrefine==true)
			){
				
				if(boolCellUnrefine[face.owner]==true){
					cout << "NONONONONONONONOON" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				// face center point
				auto& face1 = mesh.faces[i+1];
				auto& face2 = mesh.faces[i+2];
				int centerP = -1;
				for(auto& p0 : face.points){
					if(std::find(face1.points.begin(),face1.points.end(),p0)==face1.points.end()){
						continue;
					}
					if(std::find(face2.points.begin(),face2.points.end(),p0)==face2.points.end()){
						continue;
					}
					centerP = p0;
					break;
				}
				
				if(centerP==-1){
					cout << "NONONONONONONONOON" << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				
				
				vector<int> tmpProcFace;
				int ownTarget = face.owner;
				for(int j=i; j<end; ++j){
					auto& face_tmp = mesh.faces[j];
					if(std::find(face_tmp.points.begin(),face_tmp.points.end(),centerP)!=face_tmp.points.end()){
						tmpProcFace.push_back(j);
					}
					else{
						break;
					}
				}
				
				
				
				// proc_total_num += tmpProcFace.size();
				
				
				int groupNum = groupProcFaces.size();
				groupProcFaces.push_back(faces_Unrefind());
				for(auto& j : tmpProcFace){
					boolProcFaces_Combine[proc_num] = true;
					groupProcFaces_id[proc_num] = groupNum;
					groupProcFaces.back().child.push_back(j);
					++i;
					++proc_num;
				}
				--i;
				--proc_num;
				
				// combined face points point
				// vector<int> tmpPoints;
				// for(auto& j : tmpProcFace){
					// for(auto& k : mesh.faces[j].points){
						// if(mesh.points[k].level <= faceLevel-1){
							// int newPN = newPointsNumber[k];
							// tmpPoints.push_back(newPN);
						// }
					// }
				// }
				// groupProcFaces.back().points = tmpPoints;
				
				
				// proc 면 포인트 순서 바꾸기
				vector<int> tmpPoints;
				if(boundary.neighbProcNo > rank){
					
					std::reverse(tmpProcFace.begin()+1, tmpProcFace.end());
						
				}
				

				vector<int> vertexLists0;
				int tmptmp = 0;
				for(auto& j : tmpProcFace){
					auto& face_tmp = mesh.faces[j];
					int faceLevel_tmp = face_tmp.level;
					auto& points_tmp = face_tmp.points;
					
					vector<int> vertexLists;
					for(auto& p0 : points_tmp){
						if(mesh.points[p0].level<=faceLevel_tmp){
							vertexLists.push_back(p0);
						}
					}
					if(tmptmp==0){
						vertexLists0 = vertexLists;
					}
			
					if(tmptmp!=0){
						auto iStr = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[3]);
						tmpPoints.insert( tmpPoints.end(), iStr, points_tmp.end() );
					}
					auto iEnd = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[1]);
					tmpPoints.insert( tmpPoints.end(), points_tmp.begin(), iEnd );
					++tmptmp;
				}
				auto& faceZero = mesh.faces[tmpProcFace[0]];
				auto& pointsZero = faceZero.points;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists0[3]);
				tmpPoints.insert( tmpPoints.end(), iStr0, pointsZero.end() );
						
				vector<int> groupFacePoints;
				for(auto& j : tmpPoints){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						groupFacePoints.push_back(newPN);
					}
				}
				
				groupProcFaces.back().points = groupFacePoints;
				
				
			}
			++proc_num;
		}
	}
	
	
	
	
	
	
	
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(groupFaces_id[i]==-1){
			// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
				// if(rank==0) cout << "PROC" << endl;
			// }
		// }
	// }
	
	
	
	
	
	
	// cout << rank << " : 3333333" << endl;
	
	
	
	//====================================================
	// proc 면 포인트 순서 바꾸기 & 재부여
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo == -1) continue;
		if(boundary.neighbProcNo <= rank) continue;
		int str = boundary.startFace;
		int end = str + boundary.nFaces;
		for(int i=str; i<end; ++i){
			if(groupFaces_id[i] != -1){
				
				auto& groupFace = groupOutFaces[groupFaces_id[i]];
				
				
				// proc_total_num += groupFace.child.size();
				
				
				
				vector<int> tmpProcFace = groupFace.child;
				
				std::reverse(tmpProcFace.begin()+1, tmpProcFace.end());
				
				
				vector<int> tmpPoints;
				vector<int> vertexLists0;
				int tmptmp = 0;
				for(auto& j : tmpProcFace){
					auto& face_tmp = mesh.faces[j];
					int faceLevel_tmp = face_tmp.level;
					auto& points_tmp = face_tmp.points;
					
					vector<int> vertexLists;
					for(auto& p0 : points_tmp){
						if(mesh.points[p0].level<=faceLevel_tmp){
							vertexLists.push_back(p0);
						}
					}
					if(tmptmp==0){
						vertexLists0 = vertexLists;
					}
			
					if(tmptmp!=0){
						auto iStr = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[3]);
						tmpPoints.insert( tmpPoints.end(), iStr, points_tmp.end() );
					}
					auto iEnd = std::find(points_tmp.begin(),points_tmp.end(),vertexLists[1]);
					tmpPoints.insert( tmpPoints.end(), points_tmp.begin(), iEnd );
					++tmptmp;
				}
				auto& faceZero = mesh.faces[tmpProcFace[0]];
				auto& pointsZero = faceZero.points;
				auto iStr0 = std::find(pointsZero.begin(),pointsZero.end(),vertexLists0[3]);
				tmpPoints.insert( tmpPoints.end(), iStr0, pointsZero.end() );
						
						
				auto& points = groupFace.points;
				points.clear();
				for(auto& j : tmpPoints){
					int newPN = newPointsNumber[j];
					if(newPN != -1){
						points.push_back(newPN);
					}
				}
				
				
				i += groupFace.child.size();
				--i;
			}
		}
	}
	
	

	
	// proc_total_num = 0;
	// proc_num = 0;
	// for(auto& boundary : mesh.boundary){
		// if(boundary.neighbProcNo == -1) continue;
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		// for(int i=str; i<end; ++i){
			// auto& face = mesh.faces[i];
			// int faceLevel = face.level;
			// int ownLevel = mesh.cells[face.owner].level;
			// int ngbLevel = cLevel_recv[proc_num];
			// bool ownUnrefine = boolCellUnrefine[face.owner];
			// bool ngbUnrefine = cUnrefine_recv[proc_num];
			// if(ownLevel<ngbLevel && ngbUnrefine==true){
				// ++proc_total_num;
			// }
			// if(ownLevel==ngbLevel && ownUnrefine==true && ngbUnrefine==true){
				// ++proc_total_num;
			// }
		// }
	// }
	
	// proc_total_num = 0;
	// proc_num = 0;
	// for(auto& boundary : mesh.boundary){
		// if(boundary.neighbProcNo == -1) continue;
		// int str = boundary.startFace;
		// int end = str + boundary.nFaces;
		// for(int i=str; i<end; ++i){
			// auto& face = mesh.faces[i];
			// int faceLevel = face.level;
			// int ownLevel = mesh.cells[face.owner].level;
			// int ngbLevel = cLevel_recv[proc_num];
			// bool ownUnrefine = boolCellUnrefine[face.owner];
			// bool ngbUnrefine = cUnrefine_recv[proc_num];
			// if(ownLevel<ngbLevel && ngbUnrefine==true){
				// // if(groupFaces_id[i]!=-1){
					// // cout << "AAAAAAAAAAAAAAAAAAAAA" << endl;
				// // }
				// ++proc_total_num;
			// }
			// // else if(ownLevel>ngbLevel && ownUnrefine==true){
				// // if(groupFaces_id[i]==-1){
					// // cout << "AAAAAAAAAAAAAAAAAAAAA" << endl;
				// // }
				// // ++proc_total_num;
			// // }
			// // else if(ownLevel==ngbLevel && ownUnrefine==true && ngbUnrefine==true){
				// // if(groupFaces_id[i]==-1){
					// // cout << "AAAAAAAAAAAAAAAAAAAAA" << endl;
				// // }
				// // ++proc_total_num;
			// // }
			// else{
				// if(groupFaces_id[i]!=-1){
					// ++proc_total_num;
				// }
				
				// // if(groupFaces_id[i]!=-1){
					// // cout << "AAAAAAAAAAAAAAAAAAAAA " << ownUnrefine << " " << ngbUnrefine << " " << ownLevel << " " << ngbLevel << endl;
				// // }
			// }
			// // else if(groupFaces_id[i] != -1){
				// // ++proc_total_num;
			// // }
			// ++proc_num;
		// }
	// }
	
	proc_total_num = 0;
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// if(boolProcFaces_Combine[proc_num]==true){
				// ++proc_total_num;
			// }
			// else{
				// if(groupFaces_id[i]!=-1){
					int faceLevel = face.level;
					
					bool ownUnrefine = boolCellUnrefine[face.owner];
					bool ngbUnrefine = false;
					if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
					int ownLevel = mesh.cells[face.owner].level;
					int ngbLevel = cLevel_recv[proc_num];
					
					if(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel){
						++proc_total_num;
					}
					else if(ownUnrefine==true && ownLevel>ngbLevel){
						++proc_total_num;
					}
					else if(ngbUnrefine==true && ownLevel<ngbLevel){
						
						if(boolProcFaces_Combine[proc_num]==false){
							cout << "AAAAAAAAAAAAAA" << endl;
						}
						++proc_total_num;
					}
					else{
						
					}
				// }
			// }
			++proc_num;
			
		}
		
	}
	
	// cout << rank << " prco_total_num = " << proc_total_num << endl;
	
	
	//====================================================
	// 포인트 삭제
	int numN = 0;
	mesh.points.erase( std::remove_if( mesh.points.begin(), mesh.points.end(), 
		[&boolDeletePoints, &numN](SEMO_Point const& v) { 
		return boolDeletePoints[numN++]; 
		}), mesh.points.end());
		
		
	mesh.points.shrink_to_fit();
	
	
	//====================================================
	// 면 재정립
	
	
	// cout << rank << " : 444444" << endl;
	
	proc_total_num=0;
	int nBC = 0;
	int saveI;
	numN = 0;
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		int orgnBC = nBC;
		for(int j=orgnBC; j<mesh.boundary.size(); ++j){
			if(
			numN==mesh.boundary[j].startFace ||
			mesh.boundary[j].startFace == 0 ||
			mesh.boundary[j].nFaces == 0){
				mesh.boundary[j].startFace = i;
				++nBC;
			}
			else{
				break;
			}
		}
		
		// internal faces
		if(boolIntFaces_Deleted[numN]==true){
			auto& group = groupCellsUnrefine[ groupFaces_id[numN] ];
			// cout << "ADD : " << rank << " " << group.intChildFaces.size() << endl;
			numN += group.intChildFaces.size();
			--i;
			
		}
		// outer faces
		else{
			
			
			
			
			if(mesh.faces[numN].getType() != SEMO_Types::PROCESSOR_FACE){
			
			
			
				
				// let it be
				if(groupFaces_id[numN]==-1){
					auto& face_copy = mesh.faces[numN];
					
					// proc faces
					if(face_copy.getType() == SEMO_Types::PROCESSOR_FACE){
						
						if(boolProcFaces_Combine[proc_num]==true){
							
							++proc_total_num;
							
							auto& groupPrcocFace = groupProcFaces[ groupProcFaces_id[proc_num] ];
							
							face.points = groupPrcocFace.points;
							face.owner = newCellsNumber[ face_copy.owner ];
							face.neighbour = -1;
							face.setType(face_copy.getType());
							
							for(auto& j : groupPrcocFace.child){
								++numN;
								++proc_num;
							}
							--numN;
							--proc_num;
						}
						else{
							vector<int> tmpPoints;
							for(auto& j : face_copy.points){
								if(boolDeletePoints[j]==false){
									tmpPoints.push_back( newPointsNumber[j] );
								}
							}
							face.points.clear();
							face.points = tmpPoints;
							face.owner = newCellsNumber[ face_copy.owner ];
							face.neighbour = -1;
							face.setType(face_copy.getType());
							
						}
						
						
						++proc_num;
						
					}
					else{
						
						vector<int> tmpPoints;
						for(auto& j : face_copy.points){
							if(boolDeletePoints[j]==false){
								tmpPoints.push_back( newPointsNumber[j] );
							}
						}
						face.points.clear();
						face.points = tmpPoints;
						face.owner = newCellsNumber[ face_copy.owner ];
						if(face_copy.neighbour != -1){
							face.neighbour = newCellsNumber[ face_copy.neighbour ];
						}
						else{
							face.neighbour = -1;
						}
						
						// if(face.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
						
						face.setType(face_copy.getType());
					
					}
					
					++numN;
					
					
				}
				// Unrefine
				else{
					
					auto& groupFace = groupOutFaces[ groupFaces_id[numN] ];
					
					// cout << "TQA : " << numN << " " << groupFace.child[0] << endl;
					
					auto& face_copy = mesh.faces[numN];
					
					int faceLevel = face_copy.level;
					
					bool ownUnrefine = boolCellUnrefine[face_copy.owner];
					bool ngbUnrefine;
					int ownLevel = mesh.cells[face_copy.owner].level;
					int ngbLevel;
					if(face_copy.getType() == SEMO_Types::INTERNAL_FACE){
						ngbUnrefine = boolCellUnrefine[face_copy.neighbour];
						ngbLevel = mesh.cells[face_copy.neighbour].level;
					}
					else if(face_copy.getType() == SEMO_Types::PROCESSOR_FACE){
						ngbUnrefine = false;
						if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
						ngbLevel = cLevel_recv[proc_num];
						proc_num += groupFace.child.size();
					}
					else if(face_copy.getType() == SEMO_Types::BOUNDARY_FACE){
						ngbUnrefine = true;
						ngbLevel = 0;
					}
					
					
					if( 
					(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel) ||
					(ownUnrefine==true && ownLevel>ngbLevel) ||
					(ngbUnrefine==true && ownLevel<ngbLevel)
					){
						
						if(face_copy.getType() == SEMO_Types::PROCESSOR_FACE){
							++proc_total_num;
						}
						
						// if(rank==0) cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
						
						face.points.clear();
						for(auto& j : groupFace.points){
							face.points.push_back( j );
						}
						face.owner = groupFace.owner;
						face.neighbour = groupFace.neighbour;
						face.setType(face_copy.getType());
						
						// if(face.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
						// if(face.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
						
						numN += groupFace.child.size();
						
					}
					else{
						// face.points.clear();
						// for(auto& j : groupFace.points){
							// face.points.push_back( j );
						// }
						// face.owner = groupFace.owner;
						// face.neighbour = groupFace.neighbour;
						// face.setType(face_copy.getType());
						
						// cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAA" << endl;
						
						for(auto& j : groupFace.child){
							auto& face_target = mesh.faces[i];
							auto& face_child = mesh.faces[j];
							
							vector<int> tmpPoints;
							for(auto& k : face_child.points){
								if(boolDeletePoints[k]==false){
									tmpPoints.push_back( newPointsNumber[k] );
								}
							}
							
							
							// if(tmpPoints.size()<4){
								// cout << rank << " GGG : " << tmpPoints.size() << " " << face_child.points.size() << endl;
								// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
							// }
							
							face_target.points.clear();
							face_target.points = tmpPoints;
							face_target.owner = newCellsNumber[ face_child.owner ];
							if(face_child.neighbour != -1){
								face_target.neighbour = newCellsNumber[ face_child.neighbour ];
							}
							else{
								face_target.neighbour = -1;
							}
							
						// if(face_target.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
							
							
							face_target.setType(face_child.getType());
					
							++numN;
							++i;
						}
						--i;
						
						// numN += groupFace.child.size();
						
					}
					
				}
				
				
				
				
						
			}
			else{
				
				auto& face_copy = mesh.faces[numN];
				
			
				int faceLevel = face_copy.level;
				
				bool ownUnrefine = boolCellUnrefine[face_copy.owner];
				bool ngbUnrefine = false;
				if(cUnrefine_recv[proc_num]==1) ngbUnrefine = true;
				int ownLevel = mesh.cells[face_copy.owner].level;
				int ngbLevel = cLevel_recv[proc_num];
				
				if(
				(ownUnrefine==true && ngbUnrefine==true && ownLevel==ngbLevel) ||
				(ownUnrefine==true && ownLevel>ngbLevel)
				){
					
					
					// cout << "1-1" << endl;
					
					auto& groupFace = groupOutFaces[ groupFaces_id[numN] ];
					
					face.points.clear();
					for(auto& j : groupFace.points){
						face.points.push_back( j );
					}
					face.owner = groupFace.owner;
					face.neighbour = groupFace.neighbour;
					face.setType(face_copy.getType());
					
					// if(face.neighbour==9){
						// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
					// }
					// if(face.neighbour==9){
						// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
					// }
					
					numN += groupFace.child.size();
					proc_num += groupFace.child.size();
					
					// if(groupFace.child.size()!=4){
						// cout << rank << " AAAAAAAAAAAAAAA " << groupFace.child.size() << endl;
					// }
					proc_total_num += groupFace.child.size();
					
					
					// cout << "1-2" << endl;
					
				}
				else if(ngbUnrefine==true && ownLevel<ngbLevel){
					
					// if(boolProcFaces_Combine[proc_num]==false){
						// cout << "AAAAAAAAAAAAAA" << endl;
					// }
					
					// cout << "2-1 " << groupProcFaces_id[proc_num] << " " << boolProcFaces_Combine[proc_num] << endl;
					
					auto& groupProcFace = groupProcFaces[ groupProcFaces_id[proc_num] ];
					
					face.points = groupProcFace.points;
					face.owner = newCellsNumber[ face_copy.owner ];
					face.neighbour = -1;
					face.setType(face_copy.getType());
					
					
					numN += groupProcFace.child.size();
					proc_num += groupProcFace.child.size();
						
					// if(groupProcFace.child.size()!=4){
						// cout << rank << " BBBBBBBBBBBBBB " << groupProcFace.child.size() << " " << ownUnrefine << " " << groupProcFaces_id[proc_num] << " " << boolProcFaces_Combine[proc_num] << endl;
					// }
					proc_total_num += groupProcFace.child.size();
					
					// cout << "2-2" << endl;
					
				}
				else{
					
					
					
					// cout << "3-1" << endl;
					
					if(groupFaces_id[numN]==-1){
					
						vector<int> tmpPoints;
						for(auto& j : face_copy.points){
							if(boolDeletePoints[j]==false){
								tmpPoints.push_back( newPointsNumber[j] );
							}
						}
						face.points.clear();
						face.points = tmpPoints;
						face.owner = newCellsNumber[ face_copy.owner ];
						if(face_copy.neighbour != -1){
							face.neighbour = newCellsNumber[ face_copy.neighbour ];
						}
						else{
							face.neighbour = -1;
						}
						
						// if(face.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
						
						face.setType(face_copy.getType());
					
						++numN;
						++proc_num;
						
					
					}
					else{
						
						auto& groupFace = groupOutFaces[ groupFaces_id[numN] ];
						
						for(auto& j : groupFace.child){
							auto& face_target = mesh.faces[i];
							auto& face_child = mesh.faces[j];
							
							vector<int> tmpPoints;
							for(auto& k : face_child.points){
								if(boolDeletePoints[k]==false){
									tmpPoints.push_back( newPointsNumber[k] );
								}
							}
							
							
							// if(tmpPoints.size()<4){
								// cout << rank << " GGG : " << tmpPoints.size() << " " << face_child.points.size() << endl;
								// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
							// }
							
							face_target.points.clear();
							face_target.points = tmpPoints;
							face_target.owner = newCellsNumber[ face_child.owner ];
							if(face_child.neighbour != -1){
								face_target.neighbour = newCellsNumber[ face_child.neighbour ];
							}
							else{
								face_target.neighbour = -1;
							}
							
						// if(face_target.neighbour==9){
							// cout << rank << "GGGGGGGGGGGGGGGGGGGGGGGG" << endl;
						// }
							
							
							face_target.setType(face_child.getType());
					
							++numN;
							++proc_num;
							
							++i;
							
						}
						--i;
						
						
					}
					
					
					// cout << "3-2" << endl;
					
				}
				
				
				// ++proc_num;
				
				
			}
			
			
		}
			
		if(numN > mesh.faces.size()-1){
			// cout << "NUMN : " << numN << " " << mesh.faces.size() << endl;
			saveI = i;
			break;
		}
	}
	
	
	
	
	
	if(nBC!=mesh.boundary.size()){
		cout << rank << " NO BOUNDARY MATCHING" << endl;
		for(auto& boundary : mesh.boundary){
			cout << rank << " : " << boundary.neighbProcNo << " " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
		}
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	// cout << rank << " prco_total_num = " << proc_total_num << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 55555555" << endl;
	
	
	
	
	int orgFacesSize = mesh.faces.size();
	for(int i=orgFacesSize-1; i>saveI; --i){
		mesh.faces.pop_back();
	}
	
	// mesh.faces.erase(mesh.faces.begin()+saveI, mesh.faces.end());
	
	
	mesh.faces.shrink_to_fit();
	
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// ++proc_num;
		// }
	// }
	// cout << "11111 : " << rank << " " << proc_num << endl;
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(rank==0) cout << rank << " " << i << " " << mesh.faces[i].points.size() << endl;
	// }
	
	
	// if(rank==0) cout << "2222222222222" << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 66666666" << endl;

	//====================================================
	// boundary setting
	// for (int i=0; i<mesh.boundary.size(); ++i) {
		// mesh.boundary[i].startFace = startFaces[ mesh.boundary[i].startFace ];
	// }
	
	

	int maxBCnum = mesh.boundary.size()-1;
	if(mesh.boundary[maxBCnum].nFaces == 0){
		mesh.boundary[maxBCnum].startFace = mesh.faces.size();
	}
	for(int i=maxBCnum-1; i>=0; --i){
		if(mesh.boundary[i].nFaces == 0){
			mesh.boundary[i].startFace = mesh.boundary[i+1].startFace;
		}
	} 
	
	for (int i=0; i<mesh.boundary.size()-1; ++i) {
		mesh.boundary[i].nFaces = mesh.boundary[i+1].startFace-mesh.boundary[i].startFace;
	}
	int maxBDsize = mesh.boundary.size()-1;
	mesh.boundary[maxBDsize].nFaces = mesh.faces.size()-mesh.boundary[maxBDsize].startFace;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(auto& boundary : mesh.boundary){
		// if(rank==0) cout << rank << " : " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(auto& boundary : mesh.boundary){
		// if(rank==1) cout << rank << " : " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(auto& boundary : mesh.boundary){
		// if(rank==2) cout << rank << " : " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(auto& boundary : mesh.boundary){
		// if(rank==3) cout << rank << " : " << boundary.startFace << " " << boundary.nFaces << " " << mesh.faces.size() << endl;
	// }
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 77777777" << endl;
	//====================================================
	// Unrefine : 면 삭제 , 면 - 셀
	
	// if(rank==0) cout << "face erase, face-cell connect start" << endl;
	
	// cout << rank << " : face erase, face-cell connect start" << endl;
		
	// for(int i=0; i<mesh.faces.size(); ++i){
		// auto& face = mesh.faces[i];
		
		// face.varL.resize(controls.nTotalFaceLRVar,0.0);
		// face.varR.resize(controls.nTotalFaceLRVar,0.0);
		// face.var.resize(controls.nTotalFaceVar,0.0);
		
		// face.owner = newCellNum[face.owner];
		// if(face.getType()==SEMO_Types::INTERNAL_FACE){
			// face.neighbour = newCellNum[face.neighbour];
		// }
		// else{
			// face.neighbour = -1;
		// }
	// }
	

	
	//====================================================
	
	// if(rank==0) cout << "mesh clear, setting start" << endl;
	
	
	// cout << rank << " : mesh clear, setting start" << endl;

	
	
	
	
	// for(int i=0; i<mesh.faces.size(); ++i){
		// for(auto& j : mesh.faces[i].points){
			// if(rank==0) cout << i << " " << mesh.faces[i].owner << " " << mesh.faces[i].neighbour << " " << j << endl;
		// }
	// }
	
	
	
	// mesh.setCountsProcFaces();
	
	// mesh.setDisplsProcFaces(); 
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// ++proc_num;
		// }
	// }
	// if(rank==0) cout << "11111 : " << rank << " " << proc_num << " " << mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1] << endl;
	
	// for(auto& boundary : mesh.boundary){
		// if(rank==0) cout << boundary.startFace << " " << boundary.nFaces << endl;
	// }
	
	
	
	
	
	

	// mesh.buildCells();
	int orgCellSize = mesh.cells.size();
	int tmpCellNum = 0;
	for(int i=0; i<orgCellSize; ++i){
		if(groupCells_id[i] == -1){
			
			mesh.cells[tmpCellNum].var.assign(
				mesh.cells[i].var.begin(),mesh.cells[i].var.end());
			
		}
		else{
			
			int varSize = mesh.cells[tmpCellNum].var.size();
			int subCellSize = groupCellsUnrefine[groupCells_id[i]].child.size();
			double dSubCellSize = (double)subCellSize;
			vector<double> tmpVars(varSize,0.0);
			for(auto& j : groupCellsUnrefine[groupCells_id[i]].child){
				for(int k=0; k<varSize; ++k){
					tmpVars[k] += mesh.cells[i].var[k]/dSubCellSize;
				}
				++i;
			}
			--i;
			
			mesh.cells[tmpCellNum].var.assign(tmpVars.begin(),tmpVars.end());
			
		}
		
		mesh.cells[tmpCellNum].points.clear();
		mesh.cells[tmpCellNum].faces.clear();
		
		++tmpCellNum;
		
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 888888888" << endl;
	
	for(int i=totalCellNum; i<orgCellSize; ++i){
		mesh.cells.pop_back();
	}
	mesh.cells.shrink_to_fit();
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 999999999" << endl;
	
	// cout << rank << " : SETTING 0" << endl;
	
	
	mesh.buildLists();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 10" << endl;
	
	mesh.connectCelltoFaces();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 11" << endl;
	
	mesh.connectCelltoPoints();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 12" << endl;
	
	mesh.setCountsProcFaces();
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 13" << endl;
	
	mesh.setDisplsProcFaces(); 
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 14" << endl;
	
	// proc_num = 0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			// ++proc_num;
		// }
	// }
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "15 : " << rank << " " << proc_num << " " << mesh.displsProcFaces[size-1] + mesh.countsProcFaces[size-1] << " " << mesh.displsProcFaces[size-1] << " " <<  mesh.countsProcFaces[size-1] << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(auto& boundary : mesh.boundary){
		// cout << boundary.startFace << " " << boundary.nFaces << endl;
	// }

	
	// if(rank==0) cout << rank << " : SETTING 3" << endl;
	
	
	
	
	// level setting
	// if(rank==0) cout << rank << " CCCCC " << mesh.cells.size() << " " << cellsLevel.size() << " " << cellsGroup.size() << endl;
	
	for(int i=0; i<mesh.cells.size(); ++i){
		mesh.cells[i].level = cellsLevel[i];
		mesh.cells[i].group = cellsGroup[i];
	}
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 16" << endl;
	
	
	this->mpiLevels(mesh, cLevel_recv);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << rank << " : 17" << endl;
	
	proc_num = 0;
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		
		face.var.resize(controls.nTotalFaceVar,0.0);
		face.varL.resize(controls.nTotalFaceLRVar,0.0);
		face.varR.resize(controls.nTotalFaceLRVar,0.0);
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			
			int maxLevel = 
				max(mesh.cells[face.owner].level,
					mesh.cells[face.neighbour].level);
			face.level = maxLevel;
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			
			int maxLevel = 
				max(mesh.cells[face.owner].level,
					cLevel_recv[proc_num]);
			face.level = maxLevel;
			
			++proc_num;
			
		}
		else if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			face.level = mesh.cells[face.owner].level;
		}
	}
	
	// cout << rank << " : SETTING 6" << endl;
	
		
	// mesh.informations();
	
	// SEMO_Mesh_Save save;
	// string tmpFile = "./Urf" + to_string(iter);
	// // string tmpFile = "./";
	// save.vtu(tmpFile, rank, mesh);
	
	
	if(rank==0){
		cout << "| AMR - Unrefine completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
}