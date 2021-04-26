#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "build.h"
// #include <cmath>
// #include <array>
#include "mpi.h"

#include "../mpi/build.h"

void SEMO_Mesh_Builder::distributeEvenlyOneToAll(){
	
	SEMO_Mesh_Builder& mesh = *this;
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	int nBlocks = size;
	
	
	if(size == 1){
		cerr << "#ERROR : mpi size = 1, use more > 1" << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	if(rank == 0){
		int ncells = mesh.cells.size();
		int idBlockCell[ncells];
		int nDistribute[size+1];
		nDistribute[0] = 0;
		for(int i=1; i<size; ++i){
			nDistribute[i] = nDistribute[i-1] + ncells/size;
		}
		nDistribute[size] = ncells;
		
		for(int i=0; i<size; ++i){
			for(int j=nDistribute[i]; j<nDistribute[i+1]; ++j){
				idBlockCell[j] = i;
			}
		}
		
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
		
		vector<int> sendCounts(nBlocks,1);
		vector<int> sendDispls(nBlocks,0);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		int nPointsLocal_Recv;
		
		MPI_Scatterv(nPointsLocal.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
		             &nPointsLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		sendCounts[0] = nPointsLocal[0]*3;
		for(int i=1; i<nBlocks; ++i){
			sendCounts[i] = nPointsLocal[i]*3;
			sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		}
		int sendSize = 0;
		for(auto& i : sendCounts) sendSize += i;
		vector<double> xyz(sendSize,0.0); 
		vector<int> nDisplPoint(nBlocks,0);
		for(int i=0; i<mesh.points.size(); ++i){
			for(int j=strPoints[i]; j<strPoints[i+1]; ++j){
				int idBlock = idBlockPointLocal[j];
				int displPoint = sendDispls[idBlock];
				displPoint += nDisplPoint[idBlock]*3;
				xyz[displPoint] = mesh.points[i].x;
				xyz[displPoint+1] = mesh.points[i].y;
				xyz[displPoint+2] = mesh.points[i].z;
				++nDisplPoint[idBlock];
			}
		}
		nDisplPoint.clear();
		
		vector<double> xyz_Recv(nPointsLocal_Recv*3,0.0);
		
		MPI_Scatterv(xyz.data(), sendCounts.data(), sendDispls.data(), MPI_DOUBLE, 
					 xyz_Recv.data(), nPointsLocal_Recv*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					 
					
		// delete & create points
		mesh.points.clear();
		// idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
		
		for(int i=0; i<nPointsLocal_Recv; ++i){
			mesh.addPoint();
			mesh.points[i].x = xyz_Recv[i*3];
			mesh.points[i].y = xyz_Recv[i*3+1];
			mesh.points[i].z = xyz_Recv[i*3+2];
		}
		
		xyz.clear(); 
		xyz_Recv.clear();
		
		
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
		
		
		int nFacesLocal_Recv;
		std::fill(sendCounts.begin(),sendCounts.end(),1);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		MPI_Scatterv(nFacesLocal.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
		             &nFacesLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
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
				idFacePoint_Send[idBlock].push_back(face.points.size());
				for(int i=0; i<idFacePoint.size(); ++i){
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
				}
			}
		}
		
		// cout << "aaaaaa" << endl;
		// for(auto& i : idFacePoint_Send[2]){
			// cout << i << endl;
		// }
		
		std::fill(sendCounts.begin(),sendCounts.end(),1);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		vector<int> displIdFacePoint_Send(nBlocks,0);
		for(int i=0; i<nBlocks; ++i) displIdFacePoint_Send[i] = idFacePoint_Send[i].size();
		int displIdFacePoint_Recv;
		MPI_Scatterv(displIdFacePoint_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 &displIdFacePoint_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		vector<int> sendValues(sendDispls[nBlocks-1] + sendCounts[nBlocks-1],0);
		int num=0;
		for(auto& i : idFacePoint_Send){
			for(auto& j : i){
				sendValues[num] = j;
				++num;
			}
		}
		idFacePoint_Send.clear();
		idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
		
		vector<int> idFacePoint_Recv(displIdFacePoint_Recv,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 idFacePoint_Recv.data(), displIdFacePoint_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		sendValues.clear();
		
		
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
				idFacePoint_Send[idBlock].push_back(face.points.size());
				for(int i=0; i<idFacePoint.size(); ++i){
					idFacePoint_Send[idBlock].push_back(idFacePoint[i]);
				}
				
			}
			
		}
		
		std::fill(sendCounts.begin(),sendCounts.end(),1);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(nBlocks,0);
		for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
		vector<int> recvValues(1,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 recvValues.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(sendDispls[nBlocks-1] + sendCounts[nBlocks-1],0);
		num=0;
		for(auto& i : idFacePoint_Send){
			for(auto& j : i){
				sendValues[num] = j;
				++num;
			}
		}
		idFacePoint_Send.clear();
		idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
		
		int nBoundaryFacePoint_Recv = recvValues[0];
		vector<int> idBoundaryFacePoint_Recv(nBoundaryFacePoint_Recv,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 idBoundaryFacePoint_Recv.data(), nBoundaryFacePoint_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		sendValues.clear();
		recvValues.clear();
		
		
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
				
					idFacePoint_Send[idBlock].push_back(mesh.faces[k].points.size());
					for(int i=0; i<idFacePoint.size(); ++i){
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
		
		
		std::fill(sendCounts.begin(),sendCounts.end(),1);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(nBlocks,0);
		for(int i=0; i<nBlocks; ++i) sendValues[i] = idFacePoint_Send[i].size();
		recvValues.resize(1,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 recvValues.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		sendValues.clear();

		for(int i=0; i<nBlocks; ++i) sendCounts[i] = idFacePoint_Send[i].size();
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(sendDispls[nBlocks-1] + sendCounts[nBlocks-1],0);
		num=0;
		for(auto& i : idFacePoint_Send){
			for(auto& j : i){
				sendValues[num] = j;
				++num;
			}
		}
		idFacePoint_Send.clear();
		idFacePoint_Send.resize(nBlocks,vector<int>(0,0));
		
		int nProcessorFacePoint_Recv = recvValues[0];
		vector<int> idProcessorFacePoint_Recv(nProcessorFacePoint_Recv,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 idProcessorFacePoint_Recv.data(), nProcessorFacePoint_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		sendValues.clear();
		recvValues.clear();
		
		
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
		
		
		for(int i=0; i<nBlocks; ++i) sendCounts[i] = idOwnerLocal[i].size();
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(sendDispls[nBlocks-1] + sendCounts[nBlocks-1],0);
		num=0;
		for(auto& i : idOwnerLocal){
			for(auto& j : i){
				sendValues[num] = j;
				++num;
			}
		}
		idOwnerLocal.clear();
		idOwnerLocal.resize(nBlocks,vector<int>(0,0));
		
		int nOwnerLocal_Recv = nFacesLocal_Recv;
		vector<int> idOwnerLocal_Recv(nOwnerLocal_Recv,0);
		
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 idOwnerLocal_Recv.data(), nOwnerLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		sendValues.clear();
		
		
		

		// write # of neighbour
		vector<int> nNeighbourLocal_Send(nBlocks,0);
		for(int i=0; i<nBlocks; ++i){
			nNeighbourLocal_Send[i] = nFacesLocal[i] - nFacesBoundaryLocal[i] - nFacesProcessorLocal[i];
		}
		std::fill(sendCounts.begin(),sendCounts.end(),1);
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		sendValues.resize(nBlocks,0);
		for(int i=0; i<nBlocks; ++i) sendValues[i] = nNeighbourLocal_Send[i];
		recvValues.resize(1,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 recvValues.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		int nNeighbourLocal_Recv = recvValues[0];
		
		sendValues.clear();
		recvValues.clear();
		
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
		for(int i=0; i<nBlocks; ++i) sendCounts[i] = idNeighbourLocal[i].size();
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		// cout << displIdFacePoint_Recv << endl;
		sendValues.resize(sendDispls[nBlocks-1] + sendCounts[nBlocks-1],0);
		num=0;
		for(auto& i : idNeighbourLocal){
			for(auto& j : i){
				sendValues[num] = j;
				++num;
			}
		}
		idNeighbourLocal.clear();
		idNeighbourLocal.resize(nBlocks,vector<int>(0,0));
		
		vector<int> idNeighbourLocal_Recv(nNeighbourLocal_Recv,0);
		MPI_Scatterv(sendValues.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 idNeighbourLocal_Recv.data(), nNeighbourLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		sendValues.clear();
		
		

		// write of boundary faces
		vector<int> nFaceBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
		vector<int> nStartBoundaryFaceLocal_Send(nBlocks*mesh.boundary.size(),0);
		int temp_bound=0;
		for(int i=0; i<nBlocks; ++i){
			for(int ibcs=0; ibcs<mesh.boundary.size(); ++ibcs){
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
		
		
		int nBoundaryFaceLocal_Recv = mesh.boundary.size();
		MPI_Bcast(&nBoundaryFaceLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		string namesBoundaryFaceLocal_Recv[nBoundaryFaceLocal_Recv];
		for(int i=0; i<nBoundaryFaceLocal_Recv; ++i){
			int temp_size = mesh.boundary[i].name.length();
			MPI_Bcast(&temp_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			char tmp_ptr[50];
			memmove(tmp_ptr, mesh.boundary[i].name.c_str(), mesh.boundary[i].name.length());
			MPI_Bcast(tmp_ptr, temp_size, MPI_BYTE, 0, MPI_COMM_WORLD);
			
			string strrrr(tmp_ptr);
			namesBoundaryFaceLocal_Recv[i] = strrrr;
		}
		

		for(int i=0; i<nBlocks; ++i) sendCounts[i] = nBoundaryFaceLocal_Recv;
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		vector<int> nFaceBoundaryFaceLocal_Recv(nBoundaryFaceLocal_Recv,0);
		MPI_Scatterv(nFaceBoundaryFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 nFaceBoundaryFaceLocal_Recv.data(), nBlocks*nBoundaryFaceLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> nStartBoundaryFaceLocal_Recv(nBoundaryFaceLocal_Recv,0);
		MPI_Scatterv(nStartBoundaryFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 nStartBoundaryFaceLocal_Recv.data(), nBoundaryFaceLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
		
		// write of processor faces
		vector<int> nFaceProcessorFaceLocal_Send(nBlocks*nBlocks,0);
		vector<int> nStartProcessorFaceLocal_Send(nBlocks*nBlocks,0);
		vector<int> myProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
		vector<int> neighbProcNoProcessorFaceLocal_Send(nBlocks*nBlocks,0);
		int temp_proc=0;
		for(int i=0; i<nBlocks; ++i){
			int n = nFacesLocal[i];
			for(int j=0; j<nBlocks; ++j){
				n -= sendCountsProcs[i][j];
			}
			
			for(int j=0; j<nBlocks; ++j){
				if( sendCountsProcs[i][j] > 0) {
					nFaceProcessorFaceLocal_Send[temp_proc] = sendCountsProcs[i][j];
					nStartProcessorFaceLocal_Send[temp_proc] = n;
					myProcNoProcessorFaceLocal_Send[temp_proc] = i;
					neighbProcNoProcessorFaceLocal_Send[temp_proc] = j;
					
					n += sendCountsProcs[i][j];
				}
				else {
					nFaceProcessorFaceLocal_Send[temp_proc]=0;
					nStartProcessorFaceLocal_Send[temp_proc]=0;
					myProcNoProcessorFaceLocal_Send[temp_proc]=0;
					neighbProcNoProcessorFaceLocal_Send[temp_proc]=0;
				}
				++temp_proc;
			}
		}
		
		for(int i=0; i<nBlocks; ++i) sendCounts[i] = nBlocks;
		for(int i=1; i<nBlocks; ++i) sendDispls[i] = sendDispls[i-1] + sendCounts[i-1];
		
		vector<int> nFaceProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(nFaceProcessorFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 nFaceProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> nStartProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(nStartProcessorFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 nStartProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> myProcNoProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(myProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 myProcNoProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> neighbProcNoProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(neighbProcNoProcessorFaceLocal_Send.data(), sendCounts.data(), sendDispls.data(), MPI_INT, 
					 neighbProcNoProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
	
	

		
		//==========================================
		
		
		
		
		// push face points
		for(auto& i : mesh.faces){
			i.points.clear();
		}
		mesh.faces.clear();
		for(int i=0; i<displIdFacePoint_Recv; ++i){
			mesh.addFace();
			int num1 = idFacePoint_Recv[i];
			
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idFacePoint_Recv[i] );
			}
		}
		
		for(int i=0; i<nBoundaryFacePoint_Recv; ++i){
			mesh.addFace();
			int num1 = idBoundaryFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idBoundaryFacePoint_Recv[i] );
			}
		}
		
		for(int i=0; i<nProcessorFacePoint_Recv; ++i){
			mesh.addFace();
			int num1 = idProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idProcessorFacePoint_Recv[i] );
			}
		}
		
		
		// push owner 
		for(int i=0; i<nOwnerLocal_Recv; ++i){
			mesh.faces[i].owner = idOwnerLocal_Recv[i];
		}
		// push neighbour 
		for(int i=0; i<nNeighbourLocal_Recv; ++i){
			mesh.faces[i].neighbour = idNeighbourLocal_Recv[i];
		}
		
		mesh.boundary.clear();
		for(int ibcs=0; ibcs<nBoundaryFaceLocal_Recv; ++ibcs){
			mesh.addBoundary();
			mesh.boundary.back().name = namesBoundaryFaceLocal_Recv[ibcs];
			mesh.boundary.back().nFaces = nFaceBoundaryFaceLocal_Recv[ibcs];
			mesh.boundary.back().startFace = nStartBoundaryFaceLocal_Recv[ibcs];
		}
		
		
		for(int i=0; i<nBlocks; ++i){
			if(nFaceProcessorFaceLocal_Recv[i] > 0){
				mesh.addBoundary();
				string bcnames = "procBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[i]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[i]);
				mesh.boundary.back().name = bcnames;
				mesh.boundary.back().nFaces = nFaceProcessorFaceLocal_Recv[i];
				mesh.boundary.back().startFace = nStartProcessorFaceLocal_Recv[i];
				mesh.boundary.back().myProcNo = myProcNoProcessorFaceLocal_Recv[i];
				mesh.boundary.back().neighbProcNo = neighbProcNoProcessorFaceLocal_Recv[i];
			}
		}
		
		
		//==========================================
		
		
		
           
	}
	else{
		
		int nPointsLocal_Recv;            
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
		             &nPointsLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<double> xyz_Recv(nPointsLocal_Recv*3,0.0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, 
					 xyz_Recv.data(), nPointsLocal_Recv*3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
					 

		// create points
		mesh.points.clear();
		
		for(int i=0; i<nPointsLocal_Recv; ++i){
			mesh.addPoint();
			mesh.points[i].x = xyz_Recv[i*3];
			mesh.points[i].y = xyz_Recv[i*3+1];
			mesh.points[i].z = xyz_Recv[i*3+2];
		}
		
		
		
		xyz_Recv.clear();

		int nFacesLocal_Recv;            
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
		             &nFacesLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		int displIdFacePoint_Recv;
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 &displIdFacePoint_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> idFacePoint_Recv(displIdFacePoint_Recv,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 idFacePoint_Recv.data(), displIdFacePoint_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
					 

		vector<int> recvValues(1,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 recvValues.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> idBoundaryFacePoint_Recv(recvValues[0],0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 idBoundaryFacePoint_Recv.data(), recvValues[0], MPI_INT, 0, MPI_COMM_WORLD);
					 
		
		
		
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 recvValues.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> idProcessorFacePoint_Recv(recvValues[0],0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 idProcessorFacePoint_Recv.data(), recvValues[0], MPI_INT, 0, MPI_COMM_WORLD);
					 
		recvValues.clear();




		int nOwnerLocal_Recv = nFacesLocal_Recv;
		vector<int> idOwnerLocal_Recv(nOwnerLocal_Recv,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 idOwnerLocal_Recv.data(), nOwnerLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		int nNeighbourLocal_Recv;
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 &nNeighbourLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> idNeighbourLocal_Recv(nNeighbourLocal_Recv,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 idNeighbourLocal_Recv.data(), nNeighbourLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
		
		
		
		
		
		// boundary face
		int nBoundaryFaceLocal_Recv;
		MPI_Bcast(&nBoundaryFaceLocal_Recv, 1, MPI_INT, 0, MPI_COMM_WORLD);
					 
		string namesBoundaryFaceLocal_Recv[nBoundaryFaceLocal_Recv];
		
		char* temp_char;
		for(int i=0; i<nBoundaryFaceLocal_Recv; ++i){
			int temp_size;
			MPI_Bcast(&temp_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
			char tmp_ptr[50];
			MPI_Bcast(tmp_ptr, temp_size, MPI_BYTE, 0, MPI_COMM_WORLD);
			// mesh.addBoundary();
			string strrrr(tmp_ptr);
			namesBoundaryFaceLocal_Recv[i] = strrrr;
		}
		
		vector<int> nFaceBoundaryFaceLocal_Recv(nBoundaryFaceLocal_Recv,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT,  
					 nFaceBoundaryFaceLocal_Recv.data(), nBoundaryFaceLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> nStartBoundaryFaceLocal_Recv(nBoundaryFaceLocal_Recv,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT,  
					 nStartBoundaryFaceLocal_Recv.data(), nBoundaryFaceLocal_Recv, MPI_INT, 0, MPI_COMM_WORLD);
		
		

		
		
		
		
		// processor face
		vector<int> nFaceProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT,  
					 nFaceProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> nStartProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT,  
					 nStartProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
		
		vector<int> myProcNoProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 myProcNoProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
					 
		vector<int> neighbProcNoProcessorFaceLocal_Recv(nBlocks,0);
		MPI_Scatterv(NULL, NULL, NULL, MPI_INT, 
					 neighbProcNoProcessorFaceLocal_Recv.data(), nBlocks, MPI_INT, 0, MPI_COMM_WORLD);
		
		

		//==========================================
		
		// push face points
		for(auto& i : mesh.faces) i.points.clear();
		mesh.faces.clear();
		for(int i=0; i<displIdFacePoint_Recv; ++i){
			mesh.addFace();
			int num1 = idFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idFacePoint_Recv[i] );
			// if(rank==2) cout << idFacePoint_Recv[i] << endl;
			}
		}
		
			// if(rank==2) cout << idBoundaryFacePoint_Recv.size() << endl;
		for(int i=0; i<idBoundaryFacePoint_Recv.size(); ++i){
			mesh.addFace();
			int num1 = idBoundaryFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idBoundaryFacePoint_Recv[i] );
			// if(rank==2) cout << idBoundaryFacePoint_Recv[i] << endl;
			}
		}
		
		for(int i=0; i<idProcessorFacePoint_Recv.size(); ++i){
			mesh.addFace();
			int num1 = idProcessorFacePoint_Recv[i];
			for(int j=0; j<num1; ++j){
				++i;
				mesh.faces.back().points.push_back( idProcessorFacePoint_Recv[i] );
			}
		}
		
		
		// push owner 
		for(int i=0; i<idOwnerLocal_Recv.size(); ++i){
			mesh.faces[i].owner = idOwnerLocal_Recv[i];
		}
		// push neighbour 
		for(int i=0; i<idNeighbourLocal_Recv.size(); ++i){
			mesh.faces[i].neighbour = idNeighbourLocal_Recv[i];
		}
		
		mesh.boundary.clear();
		for(int ibcs=0; ibcs<nBoundaryFaceLocal_Recv; ++ibcs){
			mesh.addBoundary();
			mesh.boundary.back().name = namesBoundaryFaceLocal_Recv[ibcs];
			mesh.boundary.back().nFaces = nFaceBoundaryFaceLocal_Recv[ibcs];
			mesh.boundary.back().startFace = nStartBoundaryFaceLocal_Recv[ibcs];
		}
		


		for(int i=0; i<nBlocks; ++i){
			if(nFaceProcessorFaceLocal_Recv[i] > 0){
				mesh.addBoundary();
				string bcnames = "procBoundary" + to_string(myProcNoProcessorFaceLocal_Recv[i]) + "to" + to_string(neighbProcNoProcessorFaceLocal_Recv[i]);
				mesh.boundary.back().name = bcnames;
				mesh.boundary.back().nFaces = nFaceProcessorFaceLocal_Recv[i];
				mesh.boundary.back().startFace = nStartProcessorFaceLocal_Recv[i];
				mesh.boundary.back().myProcNo = myProcNoProcessorFaceLocal_Recv[i];
				mesh.boundary.back().neighbProcNo = neighbProcNoProcessorFaceLocal_Recv[i];
			}
		}
		//==========================================
		
		
		
    }
	
	
	
	
	
	// // check
	// // if(rank==2){
		// // for(auto& i : mesh.faces){
			// // cout << i.neighbour << endl;
		// // }
		// mesh.check();
		
		// mesh.buildCells();
		
		// mesh.setFaceTypes();
		// // for(auto& i : mesh.faces){
			// // cout << i.owner << endl;
		// // }

		// // create list
		// mesh.buildLists();
		
		// // check list
		// mesh.checkLists();
		
		// // cell's faces connection
		// mesh.connectCelltoFaces();
		
		// // cell's points connection
		// mesh.connectCelltoPoints();
		
		
		// mesh.saveFile("vtu");
	// // }
	
	
	
	
	
	
	
	
	
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}