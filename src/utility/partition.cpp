#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <numeric>
#include <string>
#include <map>
using namespace std;

#include "parmetis.h" 
#include "scotch.h" 
// #include "metis.h" 
// #include "scotchf.h" 

#include "../mesh/build.h" 
#include "../mesh/geometric.h" 

void parMETIS_Graph_Partition(int nBlocks, vector<int>& idBlockCell);
void parMETIS_Mesh_Partition(int nBlocks, vector<int>& idBlockCell);
void partitionFromSerial(int nBlocks, vector<int>& idBlockCell, vector<SEMO_Mesh_Builder>& newMesh);

SEMO_Mesh_Builder mesh;

int main(int argc, char* argv[]) {
	

	MPI::Init(); 
    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

	// for(int i=0; i<argc; ++i){
		// cout << argv[i] << endl;
	// }
	

	// SEMO_Controls_Builder controls;
	// controls.readConfigures();
	
	
	map<string,string> mapArgv;
	
	for(int i=1; i<argc; i+=2){
		string first = argv[i];
		string second;
		if(i+1==argc){
			second = "nan";
		}
		else{
			second = argv[i+1];
		}
		mapArgv.insert(make_pair(first, second));
		// cout <<  first<< " " << second << endl;
	}
	
	if( 
	mapArgv.find("-help") != mapArgv.end() ||
	mapArgv.find("-h") != mapArgv.end()
	){
		if(rank==0){
			cout << endl;
			cout << endl;
			cout << "┌─────── Partitioning helper ─────────────────────────────── " << endl;
			cout << "| -n \"int num.\"   : # of partition" << endl;
			cout << "| -s \"real num.\"  : scale of mesh" << endl;
			cout << "└───────────────────────────────────────────────────────────────── " << endl;
			cout << endl;
			cout << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		return 1;
	}
	
	
	
	int nBlocks = stoi(mapArgv["-n"]);
	double scaleMesh = stod(mapArgv["-s"]);
	
	
	
	SEMO_Mesh_Load load;
	load.OpenFoam(mesh, "./grid/");
	
	// bool boolGeo = true;
	// if(boolGeo){
		// SEMO_Utility_Math math;
		// SEMO_Mesh_Geometric geometric;
		// geometric.init(mesh);
	// }
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	SEMO_Mesh_Save save;

	// save.vtu("./grid/0/", 0, mesh);
	// return 0;
	
	
	vector<int> idBlockCell(mesh.cells.size(),0);

	// cout << "PARMETIS" << endl;
	
	parMETIS_Graph_Partition(nBlocks, idBlockCell);
	// parMETIS_Mesh_Partition(nBlocks, idBlockCell);

	
	// METIS_PartGraphKway(
		// NULL);



	// cout << "PARMETIS END" << endl;

	vector<SEMO_Mesh_Builder> newMesh;
	for(int i=0; i<nBlocks; ++i) {
		SEMO_Mesh_Builder tmpMesh;
		newMesh.push_back(tmpMesh);
	}

	partitionFromSerial(nBlocks, idBlockCell, newMesh);
	
	
	for(int ip=0; ip<nBlocks; ++ip) {
		for(auto& point : newMesh[ip].points){
			point.x *= scaleMesh;
			point.y *= scaleMesh;
			point.z *= scaleMesh;
		}
	}

	cout.precision(20);
	for(int ip=0; ip<nBlocks; ++ip){

		SEMO_Utility_Math math;
		SEMO_Mesh_Geometric geometric;
		geometric.init(newMesh[ip]);
		
		
		save.vtu("./grid/0/", ip, newMesh[ip]);
	}

	return 0;
}



void libSCOTCH_PartGraph (
	vector<int>& xadj,
	vector<int>& adjncy,
	int* nparts,
	vector<int>& part) {
	/* Scotch graph object to interface with libScotch */
	SCOTCH_Graph grafdat;
	SCOTCH_Strat stradat;
	int vertnbr, edgenbr;

	SCOTCH_graphInit (&grafdat);

	vertnbr = xadj.size()-1;
	edgenbr = xadj[vertnbr];
	
	
	SCOTCH_graphBuild (&grafdat,
					   0,
					   vertnbr, // 6, 
					   xadj.data(), // tmpxadj, 
					   NULL, 
					   NULL, 
					   NULL,
					   edgenbr, // 14, 
					   adjncy.data(), // tmpadjncy, 
					   NULL);

	SCOTCH_graphCheck(&grafdat);
	SCOTCH_stratInit(&stradat);
	SCOTCH_graphPart(&grafdat, *nparts, &stradat, part.data());
	SCOTCH_stratExit(&stradat);
	SCOTCH_graphExit(&grafdat);

}




void parMETIS_Graph_Partition(int nBlocks, vector<int>& idBlockCell){

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute ParMETIS ... ";
	}
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	// options[METIS_OPTION_DBGLVL]=1;
	// options[0] = 0;
	options[0] = 0;
	options[1] = 0;
	options[2] = 0;
	

	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	// real_t ubvec[ncon];
	// std::fill_n(ubvec, ncon, 1.02);
	
	real_t ubvec = 1.02;
	
	
	int temp_vtxdist[nBlocks];
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	int vtxdist[nBlocks+1];
	vtxdist[0] = 0;
	for(int i=1; i<nBlocks+1; ++i){
		vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	}
	
	
	
	
	vector<int> xadj(ncells+1,0);
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		for(auto& face : cell.faces){
			if(
			mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE ||
			mesh.faces[face].getType() == SEMO_Types::PROCESSOR_FACE) ++numt2;
		}
		++numt;
		xadj[numt] = xadj[numt-1] + numt2;
	}
	
	
	
	
	// int nSize = mesh.displsProcFaces[nBlocks-1] + mesh.countsProcFaces[nBlocks-1];
	// vector<int> idCell_Send;
	// for(auto& face : mesh.faces){
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			// idCell_Send.push_back(face.owner);
		// }
	// }
	
	
	// vector<int> idCell_Recv(nSize,0);
	// MPI_Alltoallv( idCell_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // idCell_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // MPI_COMM_WORLD);
				   
				   
	
	// vector<int> iRank_Send;
	// for(auto& face : mesh.faces){
		
		// if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			// iRank_Send.push_back(rank);
		// }
	// }
	// // if(nSize==0) nSize = 1;
	// vector<int> iRank_Recv(nSize,0);;
	// MPI_Alltoallv( iRank_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // iRank_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   // MPI_COMM_WORLD);
	
	
	vector<int> adjncy(xadj[ncells],0);
	
	
	vector<int> num_ncell(ncells,0);
	int proc_num = 0;
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			int tnum = xadj[face.owner] + num_ncell[face.owner];
			adjncy[tnum] = vtxdist[rank] + face.neighbour;
			++num_ncell[face.owner];
			
			// cout << face.neighbour << endl;
			
			tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			adjncy[tnum] = vtxdist[rank] + face.owner;
			++num_ncell[face.neighbour];
			
		}
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// int tnum = xadj[face.owner] + num_ncell[face.owner];
			// int rank_neighbour = iRank_Recv[proc_num];
			// // cout << proc_num << " " << iRank_Recv.size() << endl;
			// adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// ++num_ncell[face.owner];
			
			// ++proc_num;
		// }
	}
	

	ParMETIS_V3_PartKway(
		vtxdist, xadj.data(), adjncy.data(), NULL, NULL, &wgtflag, &numflag,
		&ncon, &nBlocks, tpwgts, &ubvec,
		options, &objval, idBlockCell.data(), &comm);

		
	// libSCOTCH_PartGraph(xadj, adjncy, &nBlocks, idBlockCell);


	// for(int i=0; i<ncells; ++i){
		// idBlockCell[i] = 1;
	// }
	// for(int i=0; i<ncells/2; ++i){
		// idBlockCell[i] = 0;
	// }

	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}





void parMETIS_Mesh_Partition(int nBlocks, vector<int>& idBlockCell){

    int rank = MPI::COMM_WORLD.Get_rank(); 
    int size = MPI::COMM_WORLD.Get_size();

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute ParMETIS ... ";
	}
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_DBGLVL]=1;
	options[0] = 0;
	

	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	real_t ubvec[ncon];
	std::fill_n(ubvec, ncon, 1.05);
	
	
	
	int temp_vtxdist[nBlocks];
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	int vtxdist[nBlocks+1];
	vtxdist[0] = 0;
	for(int i=1; i<nBlocks+1; ++i){
		vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	}
	
	
	
	
	vector<int> xadj(ncells+1,0);
	xadj[0] = 0;
	for(int i=1; i<ncells+1; ++i){
		xadj[i] = xadj[i-1] + mesh.cells[i-1].points.size();
	}
	
	
	
	// vector<int> num_ncell(ncells,0);
	// vector<int> adjncy(xadj[ncells],0);
	vector<int> adjncy;
	int proc_num = 0;
	for(int i=0; i<ncells; ++i){
		for(auto& point : mesh.cells[i].points){
			adjncy.push_back(point);
		}
	}
		
	// for(auto& face : mesh.points){
		
		// if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			// int tnum = xadj[face.owner] + num_ncell[face.owner];
			// adjncy[tnum] = vtxdist[rank] + face.neighbour;
			// ++num_ncell[face.owner];
			
			// // cout << face.neighbour << endl;
			
			// tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			// adjncy[tnum] = vtxdist[rank] + face.owner;
			// ++num_ncell[face.neighbour];
			
		// }
		// // else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			// // int tnum = xadj[face.owner] + num_ncell[face.owner];
			// // int rank_neighbour = iRank_Recv[proc_num];
			// // // cout << proc_num << " " << iRank_Recv.size() << endl;
			// // adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			// // ++num_ncell[face.owner];
			
			// // ++proc_num;
		// // }
	// }
	

	int ncommonnodes = 3;

	// ParMETIS_V3_PartMeshKway(
		// vtxdist, xadj.data(), adjncy.data(), NULL, &wgtflag, &numflag,
		// &ncon, &ncommonnodes, &nBlocks, tpwgts, ubvec,
		// options, &objval, idBlockCell.data(), &comm);

	// objval = adjncy.size()+1;

	vector<int> npart(npoints,0);

	METIS_PartMeshDual(
		&ncells, &npoints, xadj.data(), adjncy.data(), NULL, NULL,
		&ncommonnodes, &nBlocks, NULL, NULL,
		&objval, idBlockCell.data(), npart.data());

	
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	
	
	
	
}







void partitionFromSerial(int nBlocks, vector<int>& idBlockCell, vector<SEMO_Mesh_Builder>& newMesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int ncells = mesh.cells.size();
	
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0)); // point block id (copies)
	vector<int> nCellsLocal(nBlocks,0); // local total cells
	vector<int> idCellLocal(mesh.cells.size(),0); // local cell id
	// std::fill_n(nCellsLocal, nBlocks, 0);
	for(int i=0; i<mesh.cells.size(); ++i) {
		for(auto& j : mesh.cells[i].points){
			int l=0;
			for(int k=0; k<idBlockPoint[j].size(); ++k){
				if(idBlockPoint[j][k] == idBlockCell[i]) ++l;
			}
			if(l==0) idBlockPoint[j].push_back( idBlockCell[i] );
		}
		idCellLocal[i] = nCellsLocal[ idBlockCell[i] ];
		++nCellsLocal[ idBlockCell[i] ];
	}
	
	int nTotalLocalPointSize = 0; // total all # of local point
	for(auto& i : idBlockPoint) 
		nTotalLocalPointSize += i.size();

	// point distribution (CSR format)
	//
	// / gP : global points / lP : local points / BL : block id /
	//
	//     gP[0]     gP[1]      gP[2]      gP[3]   gP[4] ...
	// - - - - - - - - - - - - - - - - - - - - - - - - - ....
	// |  lP[0,5]  |lP[6,9]| lP[10,15] |  ....   |     |
	// |  BL[0,5]  |BL[6,9]| BL[10,15] |  ....   |     |
	// |           |       |           |         |     |
	// strPoints[0]|  strPoints[2]  strPoints[3] |  strPoints[5]  ...
	//         strPoints[1]                  strPoints[4]
	//
	vector<int> nPointsLocal(nBlocks,0); // local total points
	vector<int> idPointLocal(nTotalLocalPointSize,0); // local point id
	vector<int> idBlockPointLocal(nTotalLocalPointSize,0); // local point block id
	vector<int> strPoints(mesh.points.size()+1,0); // start of each global point
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(int j=0; j<idBlockPoint[i].size(); ++j){
			int idBlock = idBlockPoint[i][j];
			idPointLocal[nIndex] = nPointsLocal[ idBlock ];
			++nPointsLocal[ idBlock ];
			idBlockPointLocal[nIndex] = idBlock;
			++nIndex;
		}
	}
	strPoints[mesh.points.size()] = nIndex;
	
	nPointsLocal.clear();
	nPointsLocal.resize(nBlocks,0);
	for(int i=0; i<nTotalLocalPointSize; ++i)
		++nPointsLocal[ idBlockPointLocal[i] ];


	vector<int> nDisplPoint(nBlocks,0);
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			// int displPoint = sendDispls[idBlock];
			// displPoint += nDisplPoint[idBlock]*3;
			// xyz[displPoint] = mesh.points[i].x;
			// xyz[displPoint+1] = mesh.points[i].y;
			// xyz[displPoint+2] = mesh.points[i].z;
			
			newMesh[idBlock].addPoint();
			newMesh[idBlock].points.back().x = mesh.points[i].x;
			newMesh[idBlock].points.back().y = mesh.points[i].y;
			newMesh[idBlock].points.back().z = mesh.points[i].z;
			
			++nDisplPoint[idBlock];
		}
	}
	nDisplPoint.clear();
	
	
	// face setting
	vector<int> nFacesLocal(nBlocks,0); // total local faces
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal[idBlockOwner];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				face.setType(SEMO_Types::PROCESSOR_FACE); // set processor face
				++nFacesLocal[idBlockNeighbour];
			}
		}
	}
	
	
	vector<int> idFaceLocal(nBlocks,0);  // temporary local face id
	
	// internal faces
	vector<vector<int>> idFacePoint_Send(nBlocks,vector<int>(0,0));
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlock = idBlockCell[face.owner];
			++idFaceLocal[idBlock];
			
			
			vector<int> idFacePoint; // face's points id
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			// save face's points size , local points id
			// idFacePoint_Send[idBlock].push_back(face.points.size());
			// newMesh[idBlock].faces.back().points.push_back(face.points.size());
			newMesh[idBlock].addFace();
			for(int i=0; i<idFacePoint.size(); ++i){
				newMesh[idBlock].faces.back().points.push_back(idFacePoint[i]);
				// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
		}
	}
	
	
	// boundary faces
	// local boundary face size
	vector<int> nFacesBoundaryLocal(nBlocks,0); 
	// local each boundary face size
	vector<vector<int>> nFacesEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0)); 
	// local each boundary face start
	vector<vector<int>> nStartFaceEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0)); 
	vector<vector<bool>> boolFaceEachBoundaryLocal(nBlocks,vector<bool>(mesh.boundary.size(),false)); 
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int idBlock = idBlockCell[face.owner];
			++nFacesBoundaryLocal[ idBlock ];
			
			++nFacesEachBoundaryLocal[ idBlock ][ face.getTypeBC() ];
			if(boolFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] == false){
				nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = idFaceLocal[ idBlock ];
				boolFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = true;
			}
			++idFaceLocal[ idBlock ];

			// wirte bc
			vector<int> idFacePoint;
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			// save face's points size , local points id
			// idFacePoint_Send[idBlock].push_back(face.points.size());
			newMesh[idBlock].addFace();
			for(int i=0; i<idFacePoint.size(); ++i){
				newMesh[idBlock].faces.back().points.push_back(idFacePoint[i]);
				// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
		
	}
	
	
	// processor faces
	// local processor face size
	vector<int> nFacesProcessorLocal(nBlocks,0); 
	 // sending size from each processor to each processor 
	vector<vector<int>> sendCountsProcs(nBlocks,vector<int>(nBlocks,0));
	vector<int> idFacesProcessor; // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			++sendCountsProcs[idBlockOwner][idBlockNeighbour];
			++sendCountsProcs[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			int k = idFacesProcessor[jp];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			if(n==ip) {
				int idBlock = m;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].points){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				// idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				newMesh[idBlock].addFace();
				for(int i=0; i<idFacePoint.size(); ++i){
					// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					newMesh[idBlock].faces.back().points.push_back(idFacePoint[i]);
					
				}
				
			}
			else if(m==ip) {
				int idBlock = n;
				vector<int> idFacePoint;
				for(auto& i : mesh.faces[k].points){
					for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
						if(idBlockPointLocal[j] == idBlock){
							idFacePoint.push_back( idPointLocal[j] );
							break;
						}
					}
				}
			
				// idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				newMesh[idBlock].addFace();
				std::reverse(idFacePoint.begin()+1,idFacePoint.end());
				for(int i=0; i<idFacePoint.size(); ++i){
					// idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
					newMesh[idBlock].faces.back().points.push_back(idFacePoint[i]);
				}
					
			}
		}
	}
	
	
	
	
	// write owner of internal face
	vector<vector<int>> idOwnerLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			idOwnerLocal[i].push_back(idCellLocal[j]);
		}
	}
	for(int i=0; i<nBlocks; ++i){
		for(int j=0; j<idFacesProcessor.size(); ++j){
			int k = idFacesProcessor[j];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			
			if(n==i) {
				idOwnerLocal[m].push_back( idCellLocal[ mesh.faces[k].owner ] );
			}
			else if(m==i) {
				idOwnerLocal[n].push_back( idCellLocal[ mesh.faces[k].neighbour ] );
			}
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		int i=0;
		for(auto& owner : idOwnerLocal[ip]){
			newMesh[ip].faces[i].owner = owner;
			++i;
		}
	}
	
	// write neighbour of internal face
	vector<vector<int>> idNeighbourLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(
		face.getType() == SEMO_Types::INTERNAL_FACE
		){
			int j = face.neighbour;
			int i = idBlockCell[j];
			idNeighbourLocal[i].push_back(idCellLocal[j]);
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		int i=0;
		for(auto& neighbour : idNeighbourLocal[ip]){
			newMesh[ip].faces[i].neighbour = neighbour;
			++i;
		}
	}
	

	// write of boundary faces
	// vector<int> nFaceBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
	// vector<int> nStartBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
	// int temp_bound=0;
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
			newMesh[ip].addBoundary();
			
			string bcnames = mesh.boundary[ibcs].name;
			
			bcnames.erase(
			std::find_if(bcnames.rbegin(), 
			bcnames.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), 
			bcnames.end());
			
			bcnames.erase(
			bcnames.begin(), 
			std::find_if(bcnames.begin(), bcnames.end(), 
			std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			newMesh[ip].boundary.back().name = bcnames;
			
			if( nFacesEachBoundaryLocal[ip][ibcs] > 0) {
				
				newMesh[ip].boundary.back().nFaces = nFacesEachBoundaryLocal[ip][ibcs];
				newMesh[ip].boundary.back().startFace = nStartFaceEachBoundaryLocal[ip][ibcs];
			}
			else{
				newMesh[ip].boundary.back().nFaces = 0;
				newMesh[ip].boundary.back().startFace = 0;
			}
		}
	}
	
	
	
	
	// write of processor faces
	for(int ip=0; ip<nBlocks; ++ip){
		int n = nFacesLocal[ip];
		for(int jp=0; jp<nBlocks; ++jp){
			n -= sendCountsProcs[ip][jp];
		}
		
		for(int jp=0; jp<nBlocks; ++jp){
			if( sendCountsProcs[ip][jp] > 0) {
				newMesh[ip].addBoundary();
				
				string bcnames = "procBoundary" + to_string(ip) + "to" + to_string(jp);
				
				newMesh[ip].boundary.back().name = bcnames;
				
				newMesh[ip].boundary.back().nFaces = sendCountsProcs[ip][jp];
				newMesh[ip].boundary.back().startFace = n;
				newMesh[ip].boundary.back().myProcNo = ip;
				newMesh[ip].boundary.back().neighbProcNo = jp;
				
				n += sendCountsProcs[ip][jp];
			}
		}
	}
	
	// //==========================================
	
	vector<int> nPointsInf;
	vector<int> nFacesInf;
	vector<int> nCellsInf;
	vector<int> nFacesIntInf;
	vector<int> nFacesBCInf;
	vector<int> nFacesProcInf;
	vector<int> nTriangleInf;
	vector<int> nQuadrangleInf;
	vector<int> nPolygonInf;
	vector<int> nTetrahedronInf;
	vector<int> nHexahedronInf;
	vector<int> nPrismInf;
	vector<int> nPyramidInf;
	vector<int> nPolyhedronInf;
	
	int nbcs=0;
	
	for(int ip=0; ip<nBlocks; ++ip){

		newMesh[ip].check();
		
		newMesh[ip].buildCells();
		
		newMesh[ip].setFaceTypes();
		
		// create list
		newMesh[ip].buildLists();
		
		// check list
		// newMesh[ip].checkLists();
		
		// cell's faces connection
		newMesh[ip].connectCelltoFaces();
		
		// cell's points connection
		newMesh[ip].connectCelltoPoints();
		
		// set processor face counts
		newMesh[ip].setCountsProcFaces();
		
		// set processor face displacements
		newMesh[ip].setDisplsProcFaces(); 
		
		// newMesh[ip].informations();
		
		nPointsInf.push_back(newMesh[ip].points.size());
		nFacesInf.push_back(newMesh[ip].faces.size());
		nCellsInf.push_back(newMesh[ip].cells.size());
	
		nbcs=0;
		int nFacesBC = 0;
		int nFacesProc = 0;
		for(auto& boundary : newMesh[ip].boundary){
			
			if(boundary.neighbProcNo == -1){
				nFacesBC += boundary.nFaces;
				++nbcs;
			}
			else{
				nFacesProc += boundary.nFaces;
			}
		}
		
		
		// faces type
		int nFacesInt = 0;
		int nTriangle = 0;
		int nQuadrangle = 0;
		int nPolygon = 0;
		for(auto& face : newMesh[ip].faces){
			
			if(face.neighbour != -1){
				++nFacesInt;
			}
			
			if(face.points.size() == 3){
				++nTriangle;
			}
			else if(face.points.size() == 4){
				++nQuadrangle;
			}
			else{
				++nPolygon;
			}
		}
		
		
		
		// cells type
		int nTetrahedron = 0;
		int nHexahedron = 0;
		int nPrism = 0;
		int nPyramid = 0;
		int nPolyhedron = 0;
		for(auto& cell : newMesh[ip].cells){
			if( cell.points.size() == 4 && cell.faces.size() == 4 ){
				++nTetrahedron;
			}
			else if( cell.points.size() == 5 && cell.faces.size() == 5 ){
				++nPyramid;
			}
			else if( cell.points.size() == 6 && cell.faces.size() == 5 ){
				++nPrism;
			}
			else if( cell.points.size() == 8 && cell.faces.size() == 6 ){
				++nHexahedron;
			}
			else {
				++nPolyhedron;
			}
		}
		
		nFacesIntInf.push_back(nFacesInt);
		nFacesBCInf.push_back(nFacesBC);
		nFacesProcInf.push_back(nFacesProc);
		nTriangleInf.push_back(nTriangle);
		nQuadrangleInf.push_back(nQuadrangle);
		nPolygonInf.push_back(nPolygon);
		nTetrahedronInf.push_back(nTetrahedron);
		nHexahedronInf.push_back(nHexahedron);
		nPrismInf.push_back(nPrism);
		nPyramidInf.push_back(nPyramid);
		nPolyhedronInf.push_back(nPolyhedron);
	}
	

    if(rank == 0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| partitioned MPI size : " << nBlocks << endl;
		cout << "| points size : ";
		for(auto& i : nPointsInf) cout << i << " | ";
		cout << endl;
		// gatherValue = mesh.faces.size();
		cout << "| faces size : ";
		for(auto& i : nFacesInf) cout << i << " | ";
		cout << endl;
		// gatherValue = mesh.cells.size();
		cout << "| cells size : ";
		for(auto& i : nCellsInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nFacesInt;
		cout << "| internal faces size : ";
		for(auto& i : nFacesIntInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nFacesBC;
		cout << "| boundary faces size : ";
		for(auto& i : nFacesBCInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nFacesProc;
		cout << "| processor faces size : ";
		for(auto& i : nFacesProcInf) cout << i << " | ";
		cout << endl;
		cout << "| boundary types : " << nbcs << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nTriangle;
		cout << "| Triangle faces : ";
		for(auto& i : nTriangleInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nQuadrangle;
		cout << "| Quadrangle faces : ";
		for(auto& i : nQuadrangleInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPolygon;
		cout << "| Polygon faces : ";
		for(auto& i : nPolygonInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
		cout << "┌────────────────────────────────────────────────────" << endl;
		// gatherValue = nTetrahedron;
		cout << "| Tetrahedron cells : ";
		for(auto& i : nTetrahedronInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nHexahedron;
		cout << "| Hexahedron cells : ";
		for(auto& i : nHexahedronInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPrism;
		cout << "| Prism cells : ";
		for(auto& i : nPrismInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPyramid;
		cout << "| Pyramid cells : ";
		for(auto& i : nPyramidInf) cout << i << " | ";
		cout << endl;
		// gatherValue = nPolyhedron;
		cout << "| Polyhedron cells : ";
		for(auto& i : nPolyhedronInf) cout << i << " | ";
		cout << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
    }
	
	
}








