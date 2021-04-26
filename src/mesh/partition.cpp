#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "build.h"
#include <cmath>
#include <array>
#include "mpi.h"

#include "parmetis.h" 
#include "scotch.h" 

#include "../mpi/build.h"

void SEMO_Mesh_Builder::distributeOneToAll(string type){

		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	SEMO_Mesh_Builder& mesh = *this;
	
	
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	int nBlocks = size;
		
		
		
		
		
		
	mesh.distributeEvenlyOneToAll();
		
		
		
		
		
		
		
	mesh.check();
	
	mesh.buildCells();
	
	mesh.setFaceTypes();

	// create list
	mesh.buildLists();
	
	// check list
	mesh.checkLists();
	
	// cell's faces connection
	mesh.connectCelltoFaces();
	
	// cell's points connection
	mesh.connectCelltoPoints();
	
	// set processor face counts
	mesh.setCountsProcFaces();
	
	// set processor face displacements
	mesh.setDisplsProcFaces(); 
	
	// mesh.saveFile("vtu");
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	int idBlockCell[mesh.cells.size()];
	
	
	mesh.parMETIS(nBlocks, idBlockCell);
	
	
	//======================================================================
	

	vector<int> sendValues(0,0);
	vector<int> recvValues(0,0);
	vector<int> sendCounts(nBlocks,0);
	vector<int> recvCounts(nBlocks,0);
	vector<int> sendDisps(nBlocks,0);
	vector<int> recvDisps(nBlocks,0);
	sendDisps[0] = 0;
	recvDisps[0] = 0;
	int sendSize = 0;
	int recvSize = 0;
	
	SEMO_Mesh_Builder newMesh[nBlocks];
	

	vector<int> rank_Send;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			rank_Send.push_back(rank);
		}
	}
	vector<int> rank_Recv(rank_Send.size(),0);
	MPI_Alltoallv( rank_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   rank_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	// each points block id , local cell id , local N cells
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0)); // point block id (copies)
	vector<int> nCellsLocal(nBlocks,0); // local total cells
	vector<int> idCellLocal(mesh.cells.size(),0); // local cell id
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
	for(auto& i : idBlockPoint) nTotalLocalPointSize += i.size();

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
	vector<int> nPointsLocal_Send(nBlocks,0); // local total points
	vector<int> idPointLocal(nTotalLocalPointSize,0); // local point id
	vector<int> idBlockPointLocal(nTotalLocalPointSize,0); // local point block id
	vector<int> strPoints(mesh.points.size()+1,0); // start of each global point
	int nIndex = 0;
	for(int i=0; i<mesh.points.size(); ++i){
		strPoints[i] = nIndex;
		for(int j=0; j<idBlockPoint[i].size(); ++j){
			int idBlock = idBlockPoint[i][j];
			idPointLocal[nIndex] = nPointsLocal_Send[ idBlock ];
			++nPointsLocal_Send[ idBlock ];
			idBlockPointLocal[nIndex] = idBlock;
			++nIndex;
		}
	}
	strPoints[mesh.points.size()] = nIndex;
	
	
	// MPI local npoints
	//     local point's x y z
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> nPointsLocal_Recv(nBlocks,0);
	
	MPI_Alltoallv( nPointsLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nPointsLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	

	for(int i=0; i<nBlocks; ++i) sendCounts[i] = nPointsLocal_Send[i]*3;
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendSize = 0;
	for(auto& i : sendCounts) sendSize += i;
	double* xyz_Send = new double[sendSize]; 
	int nDisplPoint[nBlocks];
	std::fill_n(nDisplPoint, nBlocks, 0);
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			int displPoint = sendDisps[idBlock];
			displPoint += nDisplPoint[idBlock]*3;
			if(abs(mesh.points[i].x) < 1.e-300) mesh.points[i].x = 0.0;
			if(abs(mesh.points[i].y) < 1.e-300) mesh.points[i].y = 0.0;
			if(abs(mesh.points[i].z) < 1.e-300) mesh.points[i].z = 0.0;
			xyz_Send[displPoint] = mesh.points[i].x;
			xyz_Send[displPoint+1] = mesh.points[i].y;
			xyz_Send[displPoint+2] = mesh.points[i].z;
			++nDisplPoint[idBlock];
		}
	}
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nPointsLocal_Recv[i]*3;
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	vector<double> xyz_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0.0);
	
	
	MPI_Alltoallv( xyz_Send, sendCounts.data(), sendDisps.data(), MPI_DOUBLE, 
				   xyz_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_DOUBLE, 
				   MPI_COMM_WORLD);
			
	for(int ip=0; ip<nBlocks; ++ip){
		int i=recvDisps[ip];
		int num = 0;
		while(i<recvDisps[ip] + recvCounts[ip]){
			// cout << recvDisps[ip] + recvCounts[ip] << " " << i << endl;
			newMesh[ip].addPoint();
			newMesh[ip].points[num].x = xyz_Recv[i]; ++i;
			newMesh[ip].points[num].y = xyz_Recv[i]; ++i;
			newMesh[ip].points[num].z = xyz_Recv[i]; ++i;
			++num;
		}
	}
	delete[] xyz_Send; 
	xyz_Recv.clear();
	
	
	
	// local processor face id order
	
	vector<int> idBlockCell_Send;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idBlockCell_Send.push_back(idBlockCell[face.owner]);
		}
	}
	vector<int> idBlockCell_Recv(idBlockCell_Send.size(),0);
	MPI_Alltoallv( idBlockCell_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   idBlockCell_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	int proc_num = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			face.neighbour = proc_num;
			++proc_num;
		}
	}

	// ====================   Internal face 
	// face setting
	proc_num = 0;
	vector<int> nFacesLocal_Send(nBlocks,0); // total local faces
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal_Send[idBlockOwner];
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				// face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				++nFacesLocal_Send[idBlockNeighbour];
			}
		}
	}
	
	
	
	
	
	vector<int> nFacesLocal_Recv(nBlocks,0);

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	MPI_Alltoallv( nFacesLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFacesLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	vector<int> idFaceLocal_Send(nBlocks,0);  // temporary local face id
	vector<vector<int>> idFacePoint_Send(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::INTERNAL_FACE || 
		   // face.getType() == SEMO_Types::TO_BE_INTERNAL_FACE){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int idBlock = idBlockCell[face.owner];
			++idFaceLocal_Send[idBlock];
			
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
			idFacePoint_Send[idBlock].push_back(face.points.size());
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
		}
	}
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		sendValues[i] = idFacePoint_Send[i].size();
	}
	recvValues.resize(nBlocks,0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	// sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);

	MPI_Alltoallv( sendValues.data(),       sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idFacePoint_Recv(nBlocks+1,0);
	str_idFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// ====================   boundary face 
	
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
	
	// boundary faces
	vector<int> nFacesBoundaryLocal(nBlocks,0); // local boundary face size
	vector<vector<int>> nFacesEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0));
	vector<vector<int>> nStartFaceEachBoundaryLocal(nBlocks,vector<int>(mesh.boundary.size(),0));
	vector<vector<bool>> boolStartFaceEachBoundaryLocal(nBlocks,vector<bool>(mesh.boundary.size(),false));
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int idBlock = idBlockCell[face.owner];
			++nFacesBoundaryLocal[ idBlock ];
			
			// mesh.getBoundaryNumber();
			++nFacesEachBoundaryLocal[ idBlock ][ face.getTypeBC() ];
			
			// if(idBlock==1) cout << face.getTypeBC() << endl;
			
			if(boolStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] == false){
				nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = idFaceLocal_Send[ idBlock ];
				boolStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = true;
			}
			++idFaceLocal_Send[ idBlock ];
			
			// cout << face.getTypeBC() << endl;
			
			// wirte bc
			vector<int> idFacePoint;
			idFacePoint.clear();
			for(auto& i : face.points){
				for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
					if(idBlockPointLocal[j] == idBlock){
						idFacePoint.push_back( idPointLocal[j] );
						break;
					}
				}
			}
			
			// save face's points size , local points id
			idFacePoint_Send[idBlock].push_back(face.points.size());
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
		
	}
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		sendValues[i] = idFacePoint_Send[i].size();
	}
	recvValues.resize(nBlocks,0);

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];

	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idBoundaryFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idBoundaryFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idBoundaryFacePoint_Recv(nBlocks+1,0);
	str_idBoundaryFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idBoundaryFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];


	// //===============
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
	
	// PROCESSOR_FACE
	vector<int> nFacesOriginProcessorLocal(nBlocks,0); 
	int sendCounts_temp[nBlocks][nBlocks]; // sending size from each processor to each processor , i<->j
	for(int i=0; i<nBlocks; ++i){ for(int j=0; j<nBlocks; ++j) { sendCounts_temp[i][j]=0; } }
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::PROCESSOR_FACE ){
			int idBlock = idBlockCell[face.owner];
			++nFacesOriginProcessorLocal[idBlock];
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell_Recv[face.neighbour];
			
			++sendCounts_temp[idBlockOwner][idBlockNeighbour];
			// ++sendCounts_temp[idBlockNeighbour][idBlockOwner];
			
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
			idFacePoint_Send[idBlock].push_back(face.points.size());
			
			//######################################
			if( 
			idBlock == idBlockCell_Recv[face.neighbour] &&
			rank > rank_Send[face.neighbour]
			){
				std::reverse(idFacePoint.begin(),idFacePoint.end());
			}
			
			//######################################
			
			for(int i=0; i<idFacePoint.size(); ++i){
				idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
			}
			
		}
	}
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
	recvValues.resize(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idOriginProcessorFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),                sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idOriginProcessorFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idOriginProcessorFacePoint_Recv(nBlocks+1,0);
	str_idOriginProcessorFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idOriginProcessorFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];


	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	// //====================
	
	
	
	// ====================   processor face 
	// TO_BE_PROCESSOR_FACE
	vector<int> nFacesProcessorLocal(nBlocks,0); // local processor face size
	int sendCounts_temp2[nBlocks][nBlocks]; // sending size from each processor to each processor , i<->j
	for(int i=0; i<nBlocks; ++i){ for(int j=0; j<nBlocks; ++j) { sendCounts_temp2[i][j]=0; } }
	vector<int> idFacesProcessor(0,0); // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::TO_BE_PROCESSOR_FACE ){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			++sendCounts_temp2[idBlockOwner][idBlockNeighbour];
			++sendCounts_temp2[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	idFacePoint_Send.clear();
	idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
	for(int ip=0; ip<nBlocks; ++ip){
		for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			int k = idFacesProcessor[jp];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			// cout<<k<<endl;
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
			
				idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				for(int i=0; i<idFacePoint.size(); ++i){
					// cout << idFacePoint[i] << endl;
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
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
			
				idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
				std::reverse(idFacePoint.begin(),idFacePoint.end());
				for(int i=0; i<idFacePoint.size(); ++i){
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
				}
					
			}
		}
	}
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
	recvValues.resize(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	
	sendValues.resize(sendDisps[nBlocks-1] + sendCounts[nBlocks-1],0);
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : idFacePoint_Send){ for(auto& j : i){ sendValues.push_back(j); } }
	
	vector<int> idProcessorFacePoint_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),                sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idProcessorFacePoint_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);

	vector<int> str_idProcessorFacePoint_Recv(nBlocks+1,0);
	str_idProcessorFacePoint_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idProcessorFacePoint_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// //====================
	
	
	
	
	// write owner of internal face
	vector<vector<int>> idOwnerLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::INTERNAL_FACE ){
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
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
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
	
	
	 
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idOwnerLocal[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	sendValues.clear();
	for(auto& i : idOwnerLocal){ for(auto& j : i){ sendValues.push_back(j); } }
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nFacesLocal_Recv[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idOwnerLocal_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idOwnerLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	vector<int> str_idOwnerLocal_Recv(nBlocks+1,0);
	str_idOwnerLocal_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idOwnerLocal_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// //====================

	// write # of neighbour
	vector<vector<int>> idNeighbourLocal(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int j = face.neighbour;
			int i = idBlockCell[j];
			idNeighbourLocal[i].push_back(idCellLocal[j]);
		}
	}
	
	vector<int> nNeighbourLocal(nBlocks,0);
	for(int i=0; i<nBlocks; ++i){
		nNeighbourLocal[i] = nFacesLocal_Send[i] 
		- nFacesBoundaryLocal[i] - nFacesOriginProcessorLocal[i] - nFacesProcessorLocal[i];
	}
	// cout<< idNeighbourLocal[0].size() << " " << nNeighbourLocal[0] << endl;
	

	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : nNeighbourLocal){ sendValues.push_back(i); }
	
	vector<int> nNeighbourLocal_Recv(nBlocks,0);
	

	MPI_Alltoallv( sendValues.data(),           sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nNeighbourLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = idNeighbourLocal[i].size();
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	sendValues.clear();
	for(auto& i : idNeighbourLocal){ for(auto& j : i){ sendValues.push_back(j); } }
	
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = nNeighbourLocal_Recv[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	vector<int> idNeighbourLocal_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);
	
	MPI_Alltoallv( sendValues.data(),            sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idNeighbourLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	

	vector<int> str_idNeighbourLocal_Recv(nBlocks+1,0);
	str_idNeighbourLocal_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idNeighbourLocal_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	// ==================
	
	
	
	// write of boundary faces
	int bound_num=0;
	for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
		if(mesh.boundary[ibcs].myProcNo == -1){
			++bound_num;
		}
	}
	
	int nFaceBoundaryFaceLocal_Send[nBlocks*bound_num];
	int nStartBoundaryFaceLocal_Send[nBlocks*bound_num];
	int temp_bound=0;
	for(int i=0; i<nBlocks; ++i){
		for(int ibcs=0; ibcs<bound_num; ++ibcs){
			if( nFacesEachBoundaryLocal[i][ibcs] > 0) {
				nFaceBoundaryFaceLocal_Send[temp_bound] = nFacesEachBoundaryLocal[i][ibcs];
				nStartBoundaryFaceLocal_Send[temp_bound] = nStartFaceEachBoundaryLocal[i][ibcs];
			}
			else{
				nFaceBoundaryFaceLocal_Send[temp_bound] = 0;
				nStartBoundaryFaceLocal_Send[temp_bound] = 0;
			}
			++temp_bound;
		}
	}
		
	
	std::fill(sendCounts.begin(),sendCounts.end(),bound_num);
	std::fill(recvCounts.begin(),recvCounts.end(),bound_num);
	
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	for(auto& i : nFaceBoundaryFaceLocal_Send){ 
		sendValues.push_back(i); 
		// cout << i << endl;
	}
	
	vector<int> nFaceBoundaryFaceLocal_Recv(nBlocks*bound_num,0);

	MPI_Alltoallv( sendValues.data(),                  sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceBoundaryFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	sendValues.clear();
	for(auto& i : nStartBoundaryFaceLocal_Send){ 
		sendValues.push_back(i);
		
		// cout << i << endl;
	}
	
	vector<int> startFaceBoundaryFaceLocal_Recv(nBlocks*bound_num,0);
	

	MPI_Alltoallv( sendValues.data(),                      sendCounts.data(), sendDisps.data(), MPI_INT, 
				   startFaceBoundaryFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int ibcs=0; ibcs<bound_num; ++ibcs){
			int numm = ip*bound_num + ibcs;
			// if(rank==2) cout << rank << " " << ip << " " << ibcs << " " << nBlocks << endl;
			newMesh[ip].addBoundary();
			newMesh[ip].boundary.back().name = mesh.boundary[ibcs].name;
			newMesh[ip].boundary.back().nFaces = nFaceBoundaryFaceLocal_Recv[numm];
			newMesh[ip].boundary.back().startFace = startFaceBoundaryFaceLocal_Recv[numm];
			newMesh[ip].boundary.back().myProcNo = -1;
			newMesh[ip].boundary.back().neighbProcNo = -1;
			// // cout << nFaceBoundaryFaceLocal_Recv[numm] << endl;
			// if(rank==2) cout << rank << " " << numm << " " << newMesh[ip].points.size() << " " << bound_num << " " << startFaceBoundaryFaceLocal_Recv[numm] << " " << nFaceBoundaryFaceLocal_Recv[numm] << endl;
		}
	}
	
	
	//=======================
	
	
	// write of processor faces, PROCESSOR_FACE
	vector<int> nFaceProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> nStartProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> myProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> neighbProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	int temp_proc=0;
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal_Send[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp[i][j];
		}
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp2[i][j];
		}
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts_temp[i][j] > 0) {
				nFaceProcessorFaceLocal_Send[temp_proc] = sendCounts_temp[i][j];
				nStartProcessorFaceLocal_Send[temp_proc] = n;
				myProcNoProcessorFaceLocal_Send[temp_proc] = i;
				neighbProcNoProcessorFaceLocal_Send[temp_proc] = j;
				
				// cout << "========= " << sendCounts_temp[i][j] << " " << n << endl;
				
				n += sendCounts_temp[i][j];
			}
			else {
				nFaceProcessorFaceLocal_Send[temp_proc]=0;
				nStartProcessorFaceLocal_Send[temp_proc]=0;
				myProcNoProcessorFaceLocal_Send[temp_proc]=i;
				neighbProcNoProcessorFaceLocal_Send[temp_proc]=j;
			}
			++temp_proc;
		}
	}
	
	

	// write of processor faces, TO_BE_PROCESSOR_FACE
	vector<int> nFaceToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> nStartToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> myProcNoToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	vector<int> neighbProcNoToBeProcessorFaceLocal_Send(nBlocks*nBlocks,0);
	temp_proc=0;
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal_Send[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts_temp2[i][j];
		}
		// cout << n << endl;
		
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts_temp2[i][j] > 0) {
				nFaceToBeProcessorFaceLocal_Send[temp_proc] = sendCounts_temp2[i][j];
				nStartToBeProcessorFaceLocal_Send[temp_proc] = n;
				myProcNoToBeProcessorFaceLocal_Send[temp_proc] = i;
				neighbProcNoToBeProcessorFaceLocal_Send[temp_proc] = j;
				
				
				n += sendCounts_temp2[i][j];
			}
			else {
				nFaceToBeProcessorFaceLocal_Send[temp_proc]=0;
				nStartToBeProcessorFaceLocal_Send[temp_proc]=0;
				myProcNoToBeProcessorFaceLocal_Send[temp_proc]=i;
				neighbProcNoToBeProcessorFaceLocal_Send[temp_proc]=j;
			}
			++temp_proc;
		}
	}
	
	
	
	std::fill(sendCounts.begin(),sendCounts.end(),nBlocks);
	std::fill(recvCounts.begin(),recvCounts.end(),nBlocks);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];

	vector<int> nFaceProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nFaceProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> nStartProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nStartProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nStartProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> myProcNoProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( myProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   myProcNoProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> neighbProcNoProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( neighbProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   neighbProcNoProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);


	vector<vector<int>> ibcProcessorFace(nBlocks,vector<int>(0,0));
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=0; i<nBlocks; ++i){
			int numm = ip*nBlocks + i;
			if(nFaceProcessorFaceLocal_Recv[numm] > 0){
				
				ibcProcessorFace[ip].push_back( newMesh[ip].boundary.size() );
				
				newMesh[ip].addBoundary();
				string bcnames = "procBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[numm]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[numm]);
				newMesh[ip].boundary.back().name = bcnames;
				newMesh[ip].boundary.back().nFaces = nFaceProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().startFace = nStartProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().myProcNo = myProcNoProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().neighbProcNo = neighbProcNoProcessorFaceLocal_Recv[numm];
			}
		}
	}
	

	
	vector<int> nFaceToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nFaceToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nFaceToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> nStartToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( nStartToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   nStartToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> myProcNoToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( myProcNoToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   myProcNoToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	vector<int> neighbProcNoToBeProcessorFaceLocal_Recv(nBlocks*nBlocks,0);
	MPI_Alltoallv( neighbProcNoToBeProcessorFaceLocal_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   neighbProcNoToBeProcessorFaceLocal_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);


	vector<vector<int>> ibcToBeProcessorFace(nBlocks,vector<int>(0,0));

	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=0; i<nBlocks; ++i){
			int numm = ip*nBlocks + i;
			if(nFaceToBeProcessorFaceLocal_Recv[numm] > 0){
				
				ibcToBeProcessorFace[ip].push_back( newMesh[ip].boundary.size() );
				
				newMesh[ip].addBoundary();
				string bcnames = "procBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[numm]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[numm]);
				newMesh[ip].boundary.back().name = bcnames;
				newMesh[ip].boundary.back().nFaces = nFaceToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().startFace = nStartToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().myProcNo = myProcNoToBeProcessorFaceLocal_Recv[numm];
				newMesh[ip].boundary.back().neighbProcNo = neighbProcNoToBeProcessorFaceLocal_Recv[numm];
			}
		}
	}
	
	
	// startFace re calc.
	
	for(int ip=0; ip<nBlocks; ++ip){
		int nbc=0;
		for(auto& ibc : newMesh[ip].boundary){
			if(ibc.neighbProcNo == -1) ++nbc;
		}
		
		for(int ibc=nbc; ibc<newMesh[ip].boundary.size(); ++ibc){
			int ostart = newMesh[ip].boundary[ibc-1].startFace;
			int onface = newMesh[ip].boundary[ibc-1].nFaces;
			int nstart = ostart + onface;
			newMesh[ip].boundary[ibc].startFace = nstart;
		}
	}
	
	
	
	
	
	
	
	
	//==========================================
	
	
	// push face points
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idFacePoint_Recv[ip]; i<str_idFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idFacePoint_Recv[i] );
			}
		}
	}
	

	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idBoundaryFacePoint_Recv[ip]; i<str_idBoundaryFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idBoundaryFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idBoundaryFacePoint_Recv[i] );
				// cout << num1 << " " << idBoundaryFacePoint_Recv[i] << endl;
			}
		}
	}


	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idOriginProcessorFacePoint_Recv[ip]; i<str_idOriginProcessorFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idOriginProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idOriginProcessorFacePoint_Recv[i] );
			}
		}
	}
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int i=str_idProcessorFacePoint_Recv[ip]; i<str_idProcessorFacePoint_Recv[ip+1]; ++i){
			newMesh[ip].addFace();
			int num1 = idProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				newMesh[ip].faces.back().points.push_back( idProcessorFacePoint_Recv[i] );
				// cout << num1 << " " << idProcessorFacePoint_Recv[i] << endl;
			}
		}
	}

// if(rank==1){
	// cout << rank << " " << newMesh[0].faces.size() << endl;
	// cout << rank << " " << newMesh[1].faces.size() << endl;
	// cout << rank << " " << newMesh[0].boundary[0].nFaces << " " <<  newMesh[1].boundary[0].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[0].startFace << " " <<  newMesh[1].boundary[0].startFace << endl;
	
	// cout << rank << " " << newMesh[0].boundary[1].nFaces << " " <<  newMesh[1].boundary[1].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[1].startFace << " " <<  newMesh[1].boundary[1].startFace << endl;
	
	// cout << rank << " " << newMesh[0].boundary[2].nFaces << " " <<  newMesh[1].boundary[2].nFaces << endl;
	// cout << rank << " " << newMesh[0].boundary[2].startFace << " " <<  newMesh[1].boundary[2].startFace << endl;

// }
	//==========================================
	
	
	// push owner 
	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idOwnerLocal_Recv[ip]; i<str_idOwnerLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].owner = idOwnerLocal_Recv[i];
			
			// if(idOwnerLocal_Recv[i]>100000000)  cout<< idOwnerLocal_Recv[i] << endl;
			// cout << ip << " " << i << " " << idOwnerLocal_Recv[i] << endl;
			++numm;
		}
	}
	
	
	// push neighbour 
	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idNeighbourLocal_Recv[ip]; i<str_idNeighbourLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].neighbour = idNeighbourLocal_Recv[i];
			// cout << ip << " " << i << " " << idNeighbourLocal_Recv[i] << endl;
			++numm;
		}
	}
	
	
	// cout << newMesh[1].points.size() << 
	// " " << newMesh[1].faces.size() <<
	// " " <<	newMesh[1].boundary.size() <<
	// " " <<	newMesh[1].boundary[0].myProcNo <<
	// " " <<	newMesh[1].boundary[1].myProcNo <<
	// " " <<	newMesh[1].boundary[2].myProcNo  << endl;
	// // cout << newMesh[1].boundary[2].name  << endl;
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	for(int ip=0; ip<nBlocks; ++ip){
		int numm=0;
		for(int i=str_idOwnerLocal_Recv[ip]; i<str_idOwnerLocal_Recv[ip+1]; ++i){
			newMesh[ip].faces[numm].owner = idOwnerLocal_Recv[i];
			if( newMesh[ip].faces[numm].neighbour == -1 ){
				newMesh[ip].faces[numm].setType(SEMO_Types::BOUNDARY_FACE);
			}
			else{
				newMesh[ip].faces[numm].setType(SEMO_Types::INTERNAL_FACE);
			}
			++numm;
		}
		
		// newMesh[ip].check();
		
		// newMesh[ip].buildCells();
		
		// newMesh[ip].setFaceTypes();
	
	}
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// if(rank==1){
		// int ip=3;
		
		// // for(auto& i : newMesh[ip].faces){
			// // for(auto& j : i.points){
				// // cout << " face's points : " << j << endl;
			// // }
		// // }
		// // for(auto& i : newMesh[ip].faces){
			// // cout << " face's owner : " << i.owner << endl;
		// // }
		// // for(auto& i : newMesh[ip].faces){
			// // cout << " face's neighbour : " << i.neighbour << endl;
		// // }
		// // for(auto& i : newMesh[ip].boundary){
			// // cout << " boundary name : " << i.name << endl;
			// // cout << " boundary startFace : " << i.startFace << endl;
			// // cout << " boundary nFaces : " << i.nFaces << endl;
			// // cout << " boundary neighbProcNo : " << i.neighbProcNo << endl;
		// // }
		
		// newMesh[ip].check();
		
		// newMesh[ip].buildCells();
		
		// newMesh[ip].setFaceTypes();
		
		// // create list
		// newMesh[ip].buildLists();
		
		// // check list
		// newMesh[ip].checkLists();
		
		// // cell's faces connection
		// newMesh[ip].connectCelltoFaces();
		
		// // cell's points connection
		// newMesh[ip].connectCelltoPoints();
		
		// // set processor face counts
		// newMesh[ip].setCountsProcFaces();
		
		// // set processor face displacements
		// newMesh[ip].setDisplsProcFaces(); 

		// newMesh[ip].saveFile("vtu");
	
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	//==========================================
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	
	vector<int> nToBeInternalFace_Send(nBlocks,0);
	vector<int> startToBeInternalFace_Send(nBlocks,0);
	vector<vector<int>> idToBeInternalFace(nBlocks,vector<int>(0,0));
	vector<vector<int>> neighbMeshNoToBeInternalFace(nBlocks,vector<int>(0,0));
	for(auto& face : mesh.faces){
		if( face.getType() == SEMO_Types::INTERNAL_FACE ){
			int j = face.owner;
			int i = idBlockCell[j];
			++startToBeInternalFace_Send[i];
		}
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			++startToBeInternalFace_Send[i];
		}
	}
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			int k = idBlockCell_Recv[face.neighbour];
			if(i == k){
				idToBeInternalFace[i].push_back( startToBeInternalFace_Send[i] );
				neighbMeshNoToBeInternalFace[i].push_back( rank_Recv[face.neighbour] );
			}
			++startToBeInternalFace_Send[i];
		}
	}
	
	
	
	std::fill(sendCounts.begin(),sendCounts.end(),1);
	std::fill(recvCounts.begin(),recvCounts.end(),1);
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	sendValues.resize(nBlocks,0);
	for(int i=0; i<nBlocks; ++i) sendValues[i] = idToBeInternalFace[i].size();
	
	recvValues.clear();
	recvValues.resize(nBlocks,0);
	
	vector<int> nToBeInternalFace_Recv(nBlocks,0);
	MPI_Alltoallv( sendValues.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   recvValues.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	
	for(int i=0; i<nBlocks; ++i) sendCounts[i] = sendValues[i];
	for(int i=1; i<nBlocks; ++i) sendDisps[i] = sendDisps[i-1] + sendCounts[i-1];
	for(int i=0; i<nBlocks; ++i) recvCounts[i] = recvValues[i];
	for(int i=1; i<nBlocks; ++i) recvDisps[i] = recvDisps[i-1] + recvCounts[i-1];
	
	sendValues.clear();
	// sendValues.resize(nBlocks,0);
	recvValues.clear();
	// recvValues.resize(nBlocks,0);
	
	vector<int> idToBeInternalFace_Send;
	for(int ip=0; ip<nBlocks; ++ip) {for(auto& i : idToBeInternalFace[ip]) { idToBeInternalFace_Send.push_back(i); } } 
	vector<int> idToBeInternalFace_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);

	MPI_Alltoallv( idToBeInternalFace_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   idToBeInternalFace_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
			
			
	vector<int> neighbMeshNoToBeInternalFace_Send;
	for(int ip=0; ip<nBlocks; ++ip) {for(auto& i : neighbMeshNoToBeInternalFace[ip]) { neighbMeshNoToBeInternalFace_Send.push_back(i); } } 
	vector<int> neighbMeshNoToBeInternalFace_Recv(recvDisps[nBlocks-1] + recvCounts[nBlocks-1],0);	  
	
	MPI_Alltoallv( neighbMeshNoToBeInternalFace_Send.data(), sendCounts.data(), sendDisps.data(), MPI_INT, 
				   neighbMeshNoToBeInternalFace_Recv.data(), recvCounts.data(), recvDisps.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	vector<int> str_idToBeInternalFace_Recv(nBlocks+1,0);
	str_idToBeInternalFace_Recv[0] = 0;
	for(int i=0; i<nBlocks+1; ++i) str_idToBeInternalFace_Recv[i] = recvDisps[i-1] + recvCounts[i-1];
	
	
	
	vector<int> nCellsEach(nBlocks,0);
	vector<int> startCells(nBlocks+1,0);
	
	for(int ip=0; ip<nBlocks; ++ip){
		int maxCells=-1;
		for(auto& face : newMesh[ip].faces){
			maxCells = max(maxCells, face.owner);
		}
		if(maxCells==-1){
			nCellsEach[ip] = 0; 
		}
		else{
			nCellsEach[ip] = maxCells+1; 
		}
	}
	
	startCells[0] = 0;
	for(int ip=1; ip<nBlocks+1; ++ip){
		startCells[ip] = startCells[ip-1] + nCellsEach[ip-1];
	}
	
	
	vector<vector<int>> erasePoints(nBlocks,vector<int>(0,0));
	vector<vector<int>> eraseFaces(nBlocks,vector<int>(0,0));
	vector<vector<int>> replacePoints(nBlocks,vector<int>(0,0));
	vector<vector<int>> replacePointsRank(nBlocks,vector<int>(0,0));
	
	int nnbcs = 0;
	for(int ibcs=0; ibcs<newMesh[0].boundary.size(); ++ibcs){
		if(newMesh[0].boundary[ibcs].myProcNo == -1){
			// mesh.addBoundary();
			++nnbcs;
		}
	}
	vector<vector<int>> temp_bcfacepoint(nnbcs,vector<int>(0,0));
	vector<vector<int>> temp_procfacepoint(nBlocks,vector<int>(0,0));
	vector<vector<int>> temp_Toprocfacepoint(nBlocks,vector<int>(0,0));
	
	vector<vector<int>> temp_bcfacepointOwner(nnbcs,vector<int>(0,0));
	vector<vector<int>> temp_procfacepointOwner(nBlocks,vector<int>(0,0));
	vector<vector<int>> temp_ToprocfacepointOwner(nBlocks,vector<int>(0,0));
	
	
	

	mesh.points.clear();
	mesh.faces.clear();
	mesh.cells.clear();
	mesh.boundary.clear();
	
	
	
	// point id reset
	vector<int> startpointsize(nBlocks,0);
	vector<int> iii_temp(nBlocks,0);
	for(int ip=0; ip<nBlocks; ++ip){
		if(ip>0){
			startpointsize[ip] = startpointsize[ip-1] 
				+ newMesh[ip-1].points.size() 
				- erasePoints[ip-1].size() ;
		}
		

		// points
		vector<int> replacePointsId(newMesh[ip].points.size(),0);
		int iipoint = 0;
		for(int i=0; i<newMesh[ip].points.size(); ++i){
			
			int l=0;
			bool eraseOn = true;
			for(auto& j : erasePoints[ip]){
				if( i == j ) {
					eraseOn = false;
					break;
				}
				++l;
			}
			if(eraseOn) {
				mesh.addPoint();
				mesh.points.back().x = newMesh[ip].points[i].x;
				mesh.points.back().y = newMesh[ip].points[i].y;
				mesh.points.back().z = newMesh[ip].points[i].z;
				replacePointsId[i] = startpointsize[ip] + iipoint;
				++iipoint;
			}
			else{
				int repP = replacePoints[ip][l];
				int repR = replacePointsRank[ip][l];
				if(repR >= ip){
					cout << " ERROR " << endl;
					MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				}
				// replacePointsId[i] = startpointsize[repR] + repP;
				replacePointsId[i] = repP;
				
			// if(rank==1) cout << "aaaa====== " << ip << " " << repR << " " << replacePointsId[i] << endl;
				// if(ip==1 && i==0) {
				// }
			}
		}
		
		
		for(int i=str_idToBeInternalFace_Recv[ip]; i<str_idToBeInternalFace_Recv[ip+1]; ++i){
			int j = idToBeInternalFace_Recv[i]; // id
			int ipNgb = neighbMeshNoToBeInternalFace_Recv[i]; // neighbour face's Processor
			if( ip < ipNgb ) {
				int ii = str_idToBeInternalFace_Recv[ipNgb] + iii_temp[ipNgb];
				int jj = idToBeInternalFace_Recv[ii];
				
				int owner = newMesh[ip].faces[j].owner;
				int ownerNgb = newMesh[ipNgb].faces[jj].owner;
				
				newMesh[ip].faces[j].neighbour = startCells[ipNgb] - startCells[ip] + ownerNgb;
				
				int opoint = 0;
				
				// if( jj >= newMesh[ipNgb].faces.size() ) {
					// cout << endl;
					// cout << "#ERROR : " << ipNgb << " " << jj << " " << 
					// newMesh[ipNgb].faces.size() << " " << newMesh[ipNgb].faces[jj].points.size() << endl;
					// cout << ii << " " <<  iii[ipNgb] << " " << idToBeInternalFace_Recv.size() << endl;
					// cout << str_idToBeInternalFace_Recv[ipNgb+1] << endl;
					// cout << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				for(auto& point : newMesh[ipNgb].faces[jj].points){
					bool eraseOn = true;
					for(auto& jjj : erasePoints[ipNgb]){
						if(jjj == point) {
							eraseOn = false;
							break;
						}
					}
					if(eraseOn){
						erasePoints[ipNgb].push_back(point);
						// reverse connection for each other
						int reverse_opoint = newMesh[ip].faces[j].points.size() - opoint -1;
						int addIPoint = newMesh[ip].faces[j].points[ reverse_opoint ];
						int addIPointReplace = replacePointsId[ addIPoint ];
						// replacePoints[ipNgb].push_back( addIPoint );
						replacePoints[ipNgb].push_back( addIPointReplace );
						replacePointsRank[ipNgb].push_back(ip);
						
						// if(rank==1) cout << addIPoint << " " << replacePointsId[ addIPoint ] << endl;
						// if(rank==1){
							// cout << endl;
							// cout << newMesh[ip].points[ addIPoint ].x << " " << newMesh[ip].points[ addIPoint ].y << " " << newMesh[ip].points[ addIPoint ].z << endl;
							
							// cout << newMesh[ipNgb].points[point].x << " " << newMesh[ipNgb].points[point].y << " " << newMesh[ipNgb].points[point].z << endl;
						// }
						// cout << newMesh[ip].faces[j].points[ reverse_opoint ] << endl;
					}
					++opoint;
				}
				eraseFaces[ipNgb].push_back(jj);
				
				newMesh[ip].faces[j].setType(SEMO_Types::INTERNAL_FACE);
				newMesh[ipNgb].faces[jj].setType(SEMO_Types::TO_BE_DELETE_FACE);
				++iii_temp[ipNgb];
			}
		}
		
			
		
		for(auto& i : replacePointsId){
			if(rank==1) cout << i << endl;
		}

		// faces , faces points
		for(auto& face : newMesh[ip].faces){
			if(face.getType() == SEMO_Types::INTERNAL_FACE){
				mesh.addFace();
				for(auto& point : face.points){
					mesh.faces.back().points.push_back( replacePointsId[point] );
				}
				
				// owner
				mesh.faces.back().owner = startCells[ip] + face.owner;
				
				// neighbour
				mesh.faces.back().neighbour = startCells[ip] + face.neighbour;
			}
		}
		
		
		//boundary
		int nbcs = 0;
		for(int ibcs=0; ibcs<newMesh[ip].boundary.size(); ++ibcs){
			if(newMesh[ip].boundary[ibcs].myProcNo == -1){
				// mesh.addBoundary();
				++nbcs;
			}
		}
		// cout << "aaaa" << nbcs << endl;
		for(int ibcs=0; ibcs<nbcs; ++ibcs){
			int strF = newMesh[ip].boundary[ibcs].startFace;
			int endF = strF + newMesh[ip].boundary[ibcs].nFaces;
			for(int i=strF; i<endF; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
					temp_bcfacepointOwner[ibcs].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				temp_bcfacepoint[ibcs].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_bcfacepoint[ibcs].push_back( replacePointsId[point] );
					
					// if(rank==1) cout << " " << point << " " << replacePointsId[point] << endl;
				}
				
				// // owner
				// mesh.faces.back().owner = startCells[ip] + face.owner;
			}
		}
		
	

		//procs
		for(auto& ibcs : ibcProcessorFace[ip]){
			
			int strf = newMesh[ip].boundary[ibcs].startFace;
			int endf = strf + newMesh[ip].boundary[ibcs].nFaces;
			int negh = newMesh[ip].boundary[ibcs].neighbProcNo;
			for(int i=strf; i<endf; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
				temp_procfacepointOwner[negh].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				temp_procfacepoint[negh].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_procfacepoint[negh].push_back( replacePointsId[point] );
					// if(rank==1) cout << replacePointsId[point] << endl;
				}
			}
		}
		
		for(auto& ibcs : ibcToBeProcessorFace[ip]){
			int strf = newMesh[ip].boundary[ibcs].startFace;
			int endf = strf + newMesh[ip].boundary[ibcs].nFaces;
			int negh = newMesh[ip].boundary[ibcs].neighbProcNo;
	// if(rank==0) cout << "aaa " << strf << " " << endf <<  endl;
			for(int i=strf; i<endf; ++i){
				if(newMesh[ip].faces[i].getType() == SEMO_Types::INTERNAL_FACE) continue;
				if(newMesh[ip].faces[i].getType() == SEMO_Types::TO_BE_DELETE_FACE) continue;
				
				temp_ToprocfacepointOwner[negh].push_back( startCells[ip] + newMesh[ip].faces[i].owner );
				
				// if(rank==0) cout << startCells[ip] << endl;
				
				temp_Toprocfacepoint[negh].push_back( newMesh[ip].faces[i].points.size() );
				// if(rank==1) cout << endl;
				for(auto& point : newMesh[ip].faces[i].points){
					temp_Toprocfacepoint[negh].push_back( replacePointsId[point] );
					// if(rank==1) cout << replacePointsId[point] << endl;
				}
			}
		}
		
	}
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	

	// boundary
	int nbcs = 0;
	for(int ibcs=0; ibcs<newMesh[0].boundary.size(); ++ibcs){
		if(newMesh[0].boundary[ibcs].myProcNo == -1){
			mesh.addBoundary();
			mesh.boundary[ibcs].name = newMesh[0].boundary[ibcs].name;
			++nbcs;
		}
	}
	
	
	// if(rank==0){
		
		// for(int ip=0; ip<nbcs; ++ip){
			// for(auto& i : temp_bcfacepoint[ip]){
				
				// // cout << i << endl;
				
			// }
		// }
	// }
	
	
		
	for(int ibcs=0; ibcs<nbcs; ++ibcs){
		
		mesh.boundary[ibcs].startFace = mesh.faces.size();
		int nFaces = 0;
		for(int i=0; i<temp_bcfacepoint[ibcs].size(); ++i){
			mesh.addFace();
			
			
			int nnn = temp_bcfacepoint[ibcs][i];
			for(int j=0; j<nnn; ++j){
				++i;
				mesh.faces.back().points.push_back( temp_bcfacepoint[ibcs][i] );
			}
			++nFaces;
			
			
		}
		mesh.boundary[ibcs].nFaces = nFaces;
		mesh.boundary[ibcs].myProcNo = -1;
		mesh.boundary[ibcs].neighbProcNo = -1;
		
		
		// owner
		int iii=0;
		for(int i=mesh.faces.size()-temp_bcfacepointOwner[ibcs].size(); 
		i<mesh.faces.size(); ++i){
			mesh.faces[i].owner = temp_bcfacepointOwner[ibcs][iii];
			++iii;
		}
		
	}
	

	// MPI_Barrier(MPI_COMM_WORLD);
	// if(rank==1){
		// for(int ip=0; ip<nBlocks; ++ip){
			// cout << "ip = " << ip << endl;
			// for(auto& i : temp_Toprocfacepoint[ip]){
				// cout << i << endl;
			// }
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	
	// processor face
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(auto& i : temp_Toprocfacepoint[ip]){
			temp_procfacepoint[ip].push_back( i );
		}
		temp_Toprocfacepoint[ip].clear();
		
		for(auto& i : temp_ToprocfacepointOwner[ip]){
			temp_procfacepointOwner[ip].push_back( i );
		}
		temp_ToprocfacepointOwner[ip].clear();
	}
	
	
	// if(rank==1){
		// cout << endl;
		// for(auto& i : temp_procfacepoint[0]){
			// cout << i << endl;
		// }
	// }
	
	
	
	int iii_proc=0;
	for(int ip=0; ip<nBlocks; ++ip){
		
		// if(ip==rank) continue;
		
		if(temp_procfacepoint[ip].size()>0){
			
			mesh.addBoundary();
		
			mesh.boundary[nbcs+iii_proc].startFace = mesh.faces.size();
			int nFaces = 0;
			for(int i=0; i<temp_procfacepoint[ip].size(); ++i){
				mesh.addFace();
				int nnn = temp_procfacepoint[ip][i];
				for(int j=0; j<nnn; ++j){
					++i;
					mesh.faces.back().points.push_back( temp_procfacepoint[ip][i] );
				}
				++nFaces;
			}
			mesh.boundary[nbcs+iii_proc].nFaces = nFaces;
			mesh.boundary[nbcs+iii_proc].myProcNo = rank;
			mesh.boundary[nbcs+iii_proc].neighbProcNo = ip;
			
			// // owner
			// mesh.faces.back().owner = 
			
			++iii_proc;
			
			// owner
			int iii=0;
			for(int i=mesh.faces.size()-temp_procfacepointOwner[ip].size(); 
			i<mesh.faces.size(); ++i){
				mesh.faces[i].owner = temp_procfacepointOwner[ip][iii];
				++iii;
			}
			
		}
	}
	
	//==========================================
	

	
	// startFace re calc.
	
	for(int ip=0; ip<nBlocks; ++ip){
		int nbc=0;
		for(auto& ibc : mesh.boundary){
			if(ibc.neighbProcNo == -1) ++nbc;
		}
		
		for(int ibc=nbc; ibc<mesh.boundary.size(); ++ibc){
			int ostart = mesh.boundary[ibc-1].startFace;
			int onface = mesh.boundary[ibc-1].nFaces;
			int nstart = ostart + onface;
			mesh.boundary[ibc].startFace = nstart;
		}
	}
	
	
	
	
	
	
	

	// if(rank==1){
		// int ip=0;
		
		// for(auto& i : temp_procfacepoint[ip]){
			// cout << i << endl;
		// }
		
		
		// int face = 0;
		// for(auto& i : mesh.faces){
			// cout << " face : " << face << endl;
			// for(auto& j : i.points){
				// cout << " face's points : " << j << endl;
			// }
			// ++face;
		// }
		// for(auto& i : mesh.faces){
			// cout << " face's owner : " << i.owner << endl;
		// }
		// for(auto& i : mesh.faces){
			// cout << " face's neighbour : " << i.neighbour << endl;
		// }
		// for(auto& i : mesh.boundary){
			// cout << " boundary name : " << i.name << endl;
			// cout << " boundary startFace : " << i.startFace << endl;
			// cout << " boundary nFaces : " << i.nFaces << endl;
			// cout << " boundary neighbProcNo : " << i.neighbProcNo << endl;
		// }
	// }
	
	
	

		mesh.check();
		
		mesh.buildCells();
		
		mesh.setFaceTypes();
		
		// create list
		mesh.buildLists();
		
		// check list
		mesh.checkLists();
		
		// cell's faces connection
		mesh.connectCelltoFaces();
		
		// cell's points connection
		mesh.connectCelltoPoints();
		
		// set processor face counts
		mesh.setCountsProcFaces();
		
		// set processor face displacements
		mesh.setDisplsProcFaces(); 

		mesh.saveFile("vtu");
	
	
	
	
	
	// if(rank==3){
		// for(auto& i : newMesh[3].points){
			// cout << i.x << " " << i.y << " " << i.z << endl;
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// if(rank==0){
	// // for(int ip=3; ip<4; ++ip){
		// int ip = 1;
		
		// newMesh[ip].check();
		
		// newMesh[ip].buildCells();
		
		// newMesh[ip].setFaceTypes();
		
		// // create list
		// newMesh[ip].buildLists();
		
		// // check list
		// newMesh[ip].checkLists();
		
		// // cell's faces connection
		// newMesh[ip].connectCelltoFaces();
		
		// // cell's points connection
		// newMesh[ip].connectCelltoPoints();
		
		// // set processor face counts
		// newMesh[ip].setCountsProcFaces();
		
		// // set processor face displacements
		// newMesh[ip].setDisplsProcFaces(); 

		// newMesh[ip].saveFile("vtu");
	
	// // }
	// }
		

		
		
		
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	// //======================================================================
	
	
	// // mesh.saveFile("vtu");
	
	
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}


// void SEMO_Mesh_Builder::partitionOneToAll(){
	
	
	
	
	
	
	
	
	
	
	
	
	
// }



void SEMO_Mesh_Builder::partitionInit(string type){
	
	int ncells = (*this).cells.size();
	int idBlockCell[ncells];
	
	
	int nBlocks = (*this).mpi.getSize();
	
	// if(type=="METIS"){
		SEMO_Mesh_Builder::parMETIS(nBlocks, idBlockCell);
	// }


	SEMO_Mesh_Builder::distribute(nBlocks, idBlockCell);
	
}




void SEMO_Mesh_Builder::partition(string type){
	
	int ncells = (*this).cells.size();
	int idBlockCell[ncells];
	
	
	int nBlocks = 2;
	
	// if(type=="METIS"){
		SEMO_Mesh_Builder::parMETIS(nBlocks, idBlockCell);
	// }


	SEMO_Mesh_Builder::distribute(nBlocks, idBlockCell);
	
}




void SEMO_Mesh_Builder::parMETIS(int nBlocks, int idBlockCell[]){
		
		
	SEMO_Mesh_Builder& mesh = *this;
		
	int ncells = mesh.cells.size();
	int npoints = mesh.points.size();
	int nfaces = mesh.faces.size();
	int ncon=1;
	int ncommon=3;
	int objval;
	
	int options[METIS_NOPTIONS];
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_DBGLVL]=1;
	
	// int idBlockCell[ncells];
	// int dummy1[npoints];

	// cout << ncells << endl;
	

	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
	
	// cout << "111" << endl;
	
	int wgtflag=0;
	int numflag=0;
	real_t tpwgts[nBlocks*ncon];
	for(int i=0;i<nBlocks*ncon;++i) 
		tpwgts[i]=1.0/nBlocks;

	// tpwgts[0]=0.4;
	// tpwgts[1]=0.9;
	
	real_t ubvec[ncon];
	std::fill_n(ubvec, ncon, 1.05);
	
	
	
	int temp_vtxdist[nBlocks];
    MPI_Allgather(&ncells, 1, MPI_INT, temp_vtxdist, 1, MPI_INT, MPI_COMM_WORLD);
	
	int vtxdist[nBlocks+1];
	vtxdist[0] = 0;
	for(int i=1; i<nBlocks+1; ++i){
		vtxdist[i] = vtxdist[i-1] + temp_vtxdist[i-1];
	}
	
	// int str_Cell = 0;
	// for(int i=0; i<rank; ++i){
		// str_Cell += 
	// }
	
	int xadj[ncells+1];
	xadj[0] = 0;
	int numt = 0;
	for(auto& cell : mesh.cells){
		int numt2 = 0;
		for(auto& face : cell.faces){
			if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE){
				++numt2;
			}
			
			
			else if(mesh.faces[face].getType() == SEMO_Types::PROCESSOR_FACE) {
				++numt2;
			}
		}
		++numt;
		xadj[numt] = xadj[numt-1] + numt2;
	}
	
	
	
	
	
	
	int nSize = mesh.displsProcFaces[nBlocks-1] + mesh.countsProcFaces[nBlocks-1];
	vector<int> idCell_Send;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			idCell_Send.push_back(face.owner);
		}
	}
	// cout << nSize << endl;
	// int nnn = 0;
	// for(auto& i : mesh.boundary){
		// ++nnn;
		// cout << nnn << " " << i.neighbProcNo << endl;
	// }
	
	// cout << nSize << " " << idCell_Send.size() << 
	// " " << mesh.displsProcFaces[0] << 
	// " " << mesh.displsProcFaces[1] << 
	// " " << mesh.displsProcFaces[2] << 
	// " " << mesh.displsProcFaces[3] << 
	// " " << mesh.countsProcFaces[0] << 
	// " " << mesh.countsProcFaces[1] << 
	// " " << mesh.countsProcFaces[2] << 
	// " " << mesh.countsProcFaces[3] << 
	// endl;
	
	
	// if(nSize==0) nSize = 1;
	vector<int> idCell_Recv(nSize,0);
	MPI_Alltoallv( idCell_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   idCell_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
				   
				   
				   
	int rank = mesh.mpi.getRank();
	
	vector<int> iRank_Send;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			iRank_Send.push_back(rank);
		}
	}
	// if(nSize==0) nSize = 1;
	vector<int> iRank_Recv(nSize,0);;
	MPI_Alltoallv( iRank_Send.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   iRank_Recv.data(), mesh.countsProcFaces.data(), mesh.displsProcFaces.data(), MPI_INT, 
				   MPI_COMM_WORLD);
	
	
	int adjncy[ xadj[ncells] ];
	
	
	numt = 0;
	
	
	// cout << xadj[ncells] << endl;
	// cout << idCell_Recv.size() << endl;
	// cout << iRank_Recv.size() << endl;
	
	int num_ncell[ncells];
	std::fill_n(num_ncell, ncells, 0);
	int proc_num = 0;
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::INTERNAL_FACE) {
			
			int tnum = xadj[face.owner] + num_ncell[face.owner];
			adjncy[tnum] = vtxdist[rank] + face.neighbour;
			++num_ncell[face.owner];
			
			tnum = xadj[face.neighbour] + num_ncell[face.neighbour];
			adjncy[tnum] = vtxdist[rank] + face.owner;
			++num_ncell[face.neighbour];
			
		}
		else if(face.getType() == SEMO_Types::PROCESSOR_FACE) {
			
			int tnum = xadj[face.owner] + num_ncell[face.owner];
			int rank_neighbour = iRank_Recv[proc_num];
			// cout << proc_num << " " << iRank_Recv.size() << endl;
			adjncy[tnum] = vtxdist[rank_neighbour] + idCell_Recv[proc_num];
			++num_ncell[face.owner];
			
			++proc_num;
		}
	}
	
	
	// for(auto& cell : mesh.cells){
		// for(auto& face : cell.faces){
			// if(mesh.faces[face].getType() == SEMO_Types::INTERNAL_FACE) {
				// if( &cell == &mesh.cells[ mesh.faces[face].owner ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].neighbour;
				// }
				// else if( &cell == &mesh.cells[ mesh.faces[face].neighbour ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].owner;
				// }
				// else {
					// cerr << "#error, not matching cell != face owner or neighbour" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				// ++numt;
			// }
			
			
			// else if(mesh.faces[face].getType() == SEMO_Types::PROCESSOR_FACE) {
				
				// if( &cell == &mesh.cells[ mesh.faces[face].owner ] ){
					// // int idNeighbour = mesh.faces[face].neighbour;
					// int rank_neighbour = idCell_Recv[i];
					
					// adjncy[numt] = vtxdist[rank_neighbour] + idNeighbour;
					
					
				// }
				// else if( &cell == &mesh.cells[ mesh.faces[face].neighbour ] ){
					// adjncy[numt] = vtxdist[rank] + mesh.faces[face].owner;
				// }
				// else {
					// cerr << "#error, not matching cell != face owner or neighbour" << endl;
					// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
				// }
				
				
				// ++numt;
			// }
			
			
		// }
	// }
	


	ParMETIS_V3_PartKway(
		vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag,
		&ncon, &nBlocks, tpwgts, ubvec,
		options, &objval, idBlockCell, &comm);

	
}







void SEMO_Mesh_Builder::distribute(int nBlocks, int idBlockCell[]){
	
	SEMO_MPI_Builder mpi;
	mpi.setRank();
	mpi.setSize();
	
	SEMO_Mesh_Builder& mesh = *this;
	
	vector<vector<int>> idBlockPoint(mesh.points.size(),vector<int>(0,0)); // point block id (copies)
	int nCellsLocal[nBlocks]; // local total cells
	int idCellLocal[mesh.cells.size()]; // local cell id
	std::fill_n(nCellsLocal, nBlocks, 0);
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
	int nPointsLocal[nBlocks]; // local total points
	std::fill_n(nPointsLocal, nBlocks, 0);
	int idPointLocal[nTotalLocalPointSize]; // local point id
	int idBlockPointLocal[nTotalLocalPointSize]; // local point block id
	int strPoints[mesh.points.size()+1]; // start of each global point
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
	
	
	std::fill_n(nPointsLocal, nBlocks, 0);
	for(int i=0; i<nTotalLocalPointSize; ++i)
		++nPointsLocal[ idBlockPointLocal[i] ];
	
	// open & write local npoints
	ofstream decomPointsName[nBlocks];
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/points." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nPointsLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// wirte local point's x y z
	for(int i=0; i<mesh.points.size(); ++i){
		for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
			int idBlock = idBlockPointLocal[j];
			decomPointsName[idBlock] << "(" << mesh.points[i].x;
			decomPointsName[idBlock] << " " << mesh.points[i].y;
			decomPointsName[idBlock] << " " << mesh.points[i].z << ")" << endl;
		}
	}
	
	// close points files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	// //==========================================
	// // @@@ at MPI
	// int* idBlockCell_Send;
	// int* idBlockCell_Recv;
	// if(mpi.getSize() > 0){
		// // MPI_Allgather(&ncells, 1, MPI_INT, myvtxdist, 1, MPI_INT, MPI_COMM_WORLD);
		// // int numb_RecvNode[size];
		// // MPI_Alltoall(&numb_SendNode, 1, MPI_INT, numb_RecvNode, 1, MPI_INT, MPI_COMM_WORLD);
		// idBlockCell_Send = new int[mesh.mpi.nSend];
		// idBlockCell_Recv = new int[mesh.mpi.nRecv];
		// MPI_Alltoallv( idBlockCell_Send, mesh.mpi.sendCounts, mesh.mpi.sendDispls, MPI_INT, 
					   // idBlockCell_Recv, mesh.mpi.recvCounts, mesh.mpi.recvDispls, MPI_INT, 
					   // MPI_COMM_WORLD);
	// }
	// //==========================================
				   
	
	
	// face setting
	int nFacesLocal[nBlocks]; // total local faces
	std::fill_n(nFacesLocal, nBlocks, 0);
	for(auto& face : mesh.faces){
		int idBlockOwner = idBlockCell[face.owner];
		
		++nFacesLocal[idBlockOwner];
		
		int idBlockNeighbour;
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			idBlockNeighbour = idBlockCell[face.neighbour];
			if(idBlockOwner != idBlockNeighbour){
				face.setType(SEMO_Types::PROCESSOR_FACE); // set processor face
				++nFacesLocal[idBlockNeighbour];
			}
		}
		// //==========================================
		// // @@@ at MPI 
		// else if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// idBlockNeighbour = idBlockCell_Recv[face.neighbour]; // MPI neighbour
			// if(idBlockOwner != idBlockNeighbour){
				// face.setType(SEMO_Types::TO_BE_PROCESSOR_FACE); // set processor face
				// ++nFacesLocal[idBlockNeighbour];
			// }
			// else {
				// face.setType(SEMO_Types::TO_BE_INTERNAL_FACE); // set unknown face
				// ++nFacesLocal[idBlockNeighbour];
			// }
		// }
		// //==========================================
		
	}
	
	// write # of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/faces." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	
	int idFaceLocal[nBlocks];  // temporary local face id
	std::fill_n(idFaceLocal, nBlocks, 0);
	
	// internal faces
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
			decomPointsName[idBlock] << face.points.size() << "(";
			decomPointsName[idBlock] << idFacePoint[0];
			for(int i=1; i<idFacePoint.size(); ++i){
				decomPointsName[idBlock] << " " << idFacePoint[i];
			}
			decomPointsName[idBlock] << ")" << endl;
		}
	}
	
	
	// //==========================================
	// // @@@ at MPI
	// // TO_BE_INTERNAL_FACE 
	
	// vector<int> idBlockCell_Recv;
	// if(mpi.getSize() > 0){
		
		// int nBlockCell_Send[nBlocks];
		// int nBlockCell_Recv[nBlocks];
		// std::fill_n(nBlockCell_Send, nBlocks, 0);
		// std::fill_n(nBlockCell_Recv, nBlocks, 0);
		
		// vector<int> idBlockCell_vec;
		// for(auto& face : mesh.faces){
			// if(face.getType() == SEMO_Types::PROCESSOR_FACE){
				// int idBlock = idBlockCell[face.owner];
				// idBlockCell_vec.push_back(idBlock);
				// ++nBlockCell_Send[idBlock];
				// ++nBlockCell_Recv[idBlock];
			// }
		// }
		
		// int disBlockCell_Send[nBlocks];
		// int disBlockCell_Recv[nBlocks];
		// disBlockCell_Send[0]=0;
		// disBlockCell_Recv[0]=0;
		// for(int i=1; i<nBlocks; ++i){
			// disBlockCell_Send[i] = disBlockCell_Send[i-1] + nBlockCell_Send[i-1];
			// disBlockCell_Recv[i] = disBlockCell_Recv[i-1] + nBlockCell_Recv[i-1];
		// }
		// int idBlockCell_Send[idBlockCell_vec.size()];
		// // int idBlockCell_Recv[idBlockCell_vec.size()];
		
		// for(int i=0; i<idBlockCell_vec.size(); ++i){
			// idBlockCell_Send[i] = idBlockCell_vec[i];
		// }
		
		// MPI_Alltoallv( idBlockCell_Send, nBlockCell_Send, disBlockCell_Send, MPI_INT, 
					   // idBlockCell_Recv.data(), nBlockCell_Recv, disBlockCell_Recv, MPI_INT, 
					   // MPI_COMM_WORLD);
		
		
	// }
	
	
	
	
	// vector<int> idFaceToBeInt_Send;
	// int nFaceToBeInt_Send[nBlocks];
	// int nFaceToBeInt_Recv[nBlocks];
	// std::fill_n(nFaceToBeInt_Send, nBlocks, 0);
	// std::fill_n(nFaceToBeInt_Recv, nBlocks, 0);
	// int temp_num_to_be_internal = 0;
	// for(auto& face : mesh.faces){
		// if(face.getType() == SEMO_Types::TO_BE_INTERNAL_FACE){
			
			// int idBlockOwner = idBlockCell[face.owner];
			// int idBlockNeighbour = idBlockCell_Recv[face.neighbour];
			
			// if( idBlockOwner != mpi.getRank() ){
				// idFaceToBeInt_Send.push_back(temp_num_to_be_internal);
				// // local send face size [from block]
				// ++nFaceToBeInt_Send[idBlockNeighbour]; 
			// }
			// else {
				// // local recv face size [from block]
				// ++nFaceToBeInt_Recv[idBlockNeighbour]; 
			// }
			
		// }
		// ++temp_num_to_be_internal;
	// }
	
	// vector<int> idCellLocal_Recv;
	// if(mpi.getSize() > 0){
		// // int sendCounts[nBlocks];
		// // int recvCounts[nBlocks];
		// // std::fill_n(sendCounts, nBlocks, 0);
		// // std::fill_n(recvCounts, nBlocks, 0);
		// // MPI_Allgather(&ncells, 1, MPI_INT, myvtxdist, 1, MPI_INT, MPI_COMM_WORLD);
		// // MPI_Alltoall(&numb_SendNode, 1, MPI_INT, numb_RecvNode, 1, MPI_INT, MPI_COMM_WORLD);
		// // for(int i=0; i<mpi.getSize(); ++i){
			// // sendCounts[i] = nFaceNeighbourLocalToBeInternal[i];
		// // }
		// // int nSend = nFaceNeighbourLocalToBeInternal;
		// // int* idBlockCell_Send = new int[mesh.mpi.nSend];
		// // int* idBlockCell_Recv = new int[mesh.mpi.nRecv];
		// int sendDispls[nBlocks];
		// int recvDispls[nBlocks];
		// sendDispls[0]=0;
		// recvDispls[0]=0;
		// for(int i=1; i<nBlocks; ++i){
			// sendDispls[i] = sendDispls[i-1] + nFaceToBeInt_Send[i-1];
			// recvDispls[i] = recvDispls[i-1] + nFaceToBeInt_Recv[i-1];
		// }
		// int nSizeSend = sendDispls[nBlocks-1]+nFaceToBeInt_Send[nBlocks-1];
		// int nSizeRecv = recvDispls[nBlocks-1]+nFaceToBeInt_Recv[nBlocks-1];
		// int idCellLocal_Send[nSizeSend];
		
		// for(int i=0; i<nBlocks; ++i){
			// int strB=sendDispls[i];
			// int endB=strB+nFaceToBeInt_Send[i];
			// for(int j=strB; j<endB; ++j){
				// int idFace = idFaceToBeInt_Send[j];
				// int idOwner = mesh.faces[idFace].owner;
				// int idOwnerLocal = idCellLocal[idOwner];
				// idCellLocal_Send[j] = idOwnerLocal;
			// }
		// }
		
		// MPI_Alltoallv( idCellLocal_Send, nFaceToBeInt_Send, sendDispls, MPI_INT, 
					   // idCellLocal_Recv.data(), nFaceToBeInt_Recv, recvDispls, MPI_INT, 
					   // MPI_COMM_WORLD);
					   
		// // for(int i=1; i<nBlocks; ++i){
			// // int strB=recvDispls[i];
			// // int endB=strB+nFaceToBeInt_Recv[i];
			// // for(int j=strB; j<endB; ++j){
				// // int idOwnerLocal = idCellLocal_Recv[j];
			// // }
		// // }
		
	// }
	
	
	
	// // idCellLocal_Recv : TO_BE_INTERNAL_FACE connection neighbour local cell
	
	
	// //==========================================
	
	
	
	// boundary faces
	int nFacesBoundaryLocal[nBlocks]; // local boundary face size
	int nFacesEachBoundaryLocal[nBlocks][mesh.boundary.size()]; // local each boundary face size
	int nStartFaceEachBoundaryLocal[nBlocks][mesh.boundary.size()]; // local each boundary face start
	for(int i=0; i<nBlocks; ++i){
		nFacesBoundaryLocal[i] = 0;
		for(int j=0; j<mesh.boundary.size(); ++j){
			nFacesEachBoundaryLocal[i][j]=0;
			nStartFaceEachBoundaryLocal[i][j]=0;
		}
	}
	for(auto& face : mesh.faces){
		
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int idBlock = idBlockCell[face.owner];
			++nFacesBoundaryLocal[ idBlock ];
			
			// mesh.getBoundaryNumber();
			++nFacesEachBoundaryLocal[ idBlock ][ face.getTypeBC() ];
			if(nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] == 0)
				nStartFaceEachBoundaryLocal[ idBlock ][ face.getTypeBC() ] = idFaceLocal[ idBlock ];
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
		
			decomPointsName[idBlock] << face.points.size() << "(";
			decomPointsName[idBlock] << idFacePoint[0];
			for(int i=1; i<idFacePoint.size(); ++i){
				decomPointsName[idBlock] << " " << idFacePoint[i];
			}
			decomPointsName[idBlock] << ")" << endl;
		}
		
	}
	
	
	
	// processor faces
	int nFacesProcessorLocal[nBlocks]; // local processor face size
	int sendCounts[nBlocks][nBlocks]; // sending size from each processor to each processor 
	for(int i=0; i<nBlocks; ++i){
		nFacesProcessorLocal[i] = 0;
		for(int j=0; j<nBlocks; ++j)
			sendCounts[i][j]=0;
	}
	vector<int> idFacesProcessor; // local processor face id
	int temp_num_proc_face = 0;
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			idFacesProcessor.push_back(temp_num_proc_face);
			
			int idBlockOwner = idBlockCell[face.owner];
			int idBlockNeighbour = idBlockCell[face.neighbour];
			
			++nFacesProcessorLocal[idBlockOwner];
			++nFacesProcessorLocal[idBlockNeighbour];
			++sendCounts[idBlockOwner][idBlockNeighbour];
			++sendCounts[idBlockNeighbour][idBlockOwner];
			
		}
		++temp_num_proc_face;
	}
	
	
	for(int ip=0; ip<nBlocks; ++ip){
		for(int jp=0; jp<idFacesProcessor.size(); ++jp){
			int k = idFacesProcessor[jp];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			// cout<<k<<endl;
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
			
				decomPointsName[idBlock] << mesh.faces[k].points.size() << "(";
				decomPointsName[idBlock] << idFacePoint[0];
				for(int i=1; i<idFacePoint.size(); ++i){
					decomPointsName[idBlock] << " " << idFacePoint[i];
				}
				decomPointsName[idBlock] << ")" << endl;
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
			
				decomPointsName[idBlock] << mesh.faces[k].points.size() << "(";
				
                // inverse face vector direction, because ngbr cell index convert to owner cell index
				std::reverse(idFacePoint.begin(),idFacePoint.end());
				decomPointsName[idBlock] << idFacePoint[0];
				for(int i=1; i<idFacePoint.size(); ++i){
					decomPointsName[idBlock] << " " << idFacePoint[i];
				}
				
				decomPointsName[idBlock] << ")" << endl;
			}
		}
	}
	
	// close faces files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	
	//==========================================
	// @@@ at MPI
	// TO_BE_PROCESSOR_FACE
	
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::TO_BE_PROCESSOR_FACE){
			
		}
	}
	
	
	
	//==========================================
	
	
	// write # of owner, same as number of faces
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/owner." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// write owner of internal face
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// write owner of boundary face
	for(auto& face : mesh.faces){
		if(face.getType() == SEMO_Types::BOUNDARY_FACE){
			int j = face.owner;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// write owner of processor face
	for(int i=0; i<nBlocks; ++i){
		for(int j=0; j<idFacesProcessor.size(); ++j){
			int k = idFacesProcessor[j];
			int m = idBlockCell[ mesh.faces[k].owner ];
			int n = idBlockCell[ mesh.faces[k].neighbour ];
			
			if(n==i) {
				decomPointsName[m] << idCellLocal[ mesh.faces[k].owner ] << endl;
			}
			else if(m==i) {
				decomPointsName[n] << idCellLocal[ mesh.faces[k].neighbour ] << endl;
			}
		}
	}
	
	
	// close owner files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	// write # of neighbour
	for(int i=0; i<nBlocks; ++i){
		string sFilename = "./grid/neighbour." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << nFacesLocal[i] - nFacesBoundaryLocal[i] - nFacesProcessorLocal[i] << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// write neighbour of internal face
	for(auto& face : mesh.faces){
		if(
		face.getType() == SEMO_Types::INTERNAL_FACE
		){
			int j = face.neighbour;
			int i = idBlockCell[j];
			decomPointsName[i] << idCellLocal[j] << endl;
		}
	}
	
	// close neighbour files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	// write # of boundary and processor faces
	for(int i=0; i<nBlocks; ++i){
		int n=0;
		for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
			if( nStartFaceEachBoundaryLocal[i][ibcs] > 0) ++n;
		}
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts[i][j] > 0) ++n;
		}
		
		string sFilename = "./grid/boundary." + to_string(i);
		decomPointsName[i].open(sFilename);
		decomPointsName[i] << n << endl;
		decomPointsName[i] << "(" << endl;
	}
	
	// write of boundary faces
	for(int i=0; i<nBlocks; ++i){
		for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
			if( nStartFaceEachBoundaryLocal[i][ibcs] > 0) {
				decomPointsName[i] << "   " << mesh.boundary[ibcs].name << endl;
				decomPointsName[i] << "   {" << endl;
				decomPointsName[i] << "      type Unspecified;" << endl;
				decomPointsName[i] << "      nFaces " << nFacesEachBoundaryLocal[i][ibcs] << endl;
				decomPointsName[i] << "      startFace " << nStartFaceEachBoundaryLocal[i][ibcs] << endl;
				decomPointsName[i] << "   }" << endl;
				decomPointsName[i] << endl;
				
			}
		}
	}
	
	
	// write of processor faces
	for(int i=0; i<nBlocks; ++i){
		int n = nFacesLocal[i];
		for(int j=0; j<nBlocks; ++j){
			n -= sendCounts[i][j];
		}
		
		for(int j=0; j<nBlocks; ++j){
			if( sendCounts[i][j] > 0) {
				string bcnames = "   procBoundary" + to_string(i) + "to" + to_string(j);
				decomPointsName[i] << bcnames << endl;
				decomPointsName[i] << "   {" << endl;
				decomPointsName[i] << "      type processor;" << endl;
				decomPointsName[i] << "      nFaces " << sendCounts[i][j] << endl;
				decomPointsName[i] << "      startFace " << n << endl;
				decomPointsName[i] << "      myProcNo " << i << endl;
				decomPointsName[i] << "      neighbProcNo " << j << endl;
				decomPointsName[i] << "   }" << endl;
				decomPointsName[i] << endl;
				
				n += sendCounts[i][j];
			}
		}
	}
	
	
	// close boundary & processor face files
	for(int i=0; i<nBlocks; ++i) {
		decomPointsName[i] << ")";
		decomPointsName[i].close();
	}
	
	
	
	
}