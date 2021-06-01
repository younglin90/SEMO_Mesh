#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include "mpi.h"

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_krylov.h"

void SEMO_Solvers_Builder::solveHYPRE(
	SEMO_Mesh_Builder& mesh,
	vector<double>& resiVar, 
	vector<double>& linA, vector<double>& linAL, vector<double>& linAR, 
	vector<double>& linB,
	string solver, double tolerance, double relTol, string preconditioner,
	int maxIter
	){
		
	int ierr;
	
	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();
	
	int ncells = linA.size();
	vector<int> procNcells(size,0);
	MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	int myStart=0;
	for(int i=0; i<rank; ++i){
		myStart += procNcells[i];
	}
	
	int ncellTot;
	MPI_Allreduce(&ncells, &ncellTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	vector<int> procStart(size,0);
	MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	
	
	
	HYPRE_IJMatrix A;
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_IJVector b;
	HYPRE_ParVector par_b;
	HYPRE_IJVector x;
	HYPRE_ParVector par_x;
	
	HYPRE_Solver sol, pre;
	
	// if(HYPRE_Init()) cout << "# Error 1" << endl;
	
	HYPRE_Init();
	
	int ilower = myStart;
	int iupper = ilower + procNcells[rank] - 1;
	
	// if(HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A)) cout << "# Error 1" << endl;
	// if(HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR)) cout << "# Error 1" << endl;
	// if(HYPRE_IJMatrixInitialize(A)) cout << "# Error 1" << endl;
	
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);
	
	int cols;
	int rows;
	double values;
	
	int i1 = 1;

	for(int i=0; i<mesh.cells.size(); ++i){
		cols = myStart + i;
		values = linA[i];
		HYPRE_IJMatrixSetValues(A, 1, &i1, &cols, &cols, &values);
		
	}
	
	
	
	
	vector<int> neighbProcNo;
	vector<int> ownerNo;
	vector<int> neighbNo;
	
	for(auto& boundary : mesh.boundary){
		if(boundary.neighbProcNo != -1){
			int str=boundary.startFace;
			int end=str+boundary.nFaces;
			for(int i=str; i<end; ++i){
				neighbProcNo.push_back(boundary.neighbProcNo);
				ownerNo.push_back(mesh.faces[i].owner);
				// cout << mesh.faces[i].owner << endl;
			}
		}
	}
	
	if(size>1){

		SEMO_MPI_Builder mpi;
		
		mpi.setProcsFaceDatas(
					ownerNo, neighbNo,
					mesh.countsProcFaces, mesh.countsProcFaces, 
					mesh.displsProcFaces, mesh.displsProcFaces);
	
	}
	
	int num=0;
	for(int i=0; i<mesh.faces.size(); ++i){
		
		if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			
			int own = mesh.faces[i].owner;
			int ngb = mesh.faces[i].neighbour;
			
			cols = myStart + ngb;
			rows = myStart + own;
			values = linAR[i];
			
			HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			cols = myStart + own;
			rows = myStart + ngb;
			values = linAL[i];
			
			HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
		}
		if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			
			int own = myStart + mesh.faces[i].owner;
			int ngb = procStart[neighbProcNo[num]] + neighbNo[num];
			
			// cout << ngb << " " << own << " " << ncellTot << endl;
			
			cols = ngb;
			rows = own;
			values = linAR[i];
			
			HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			cols = own;
			rows = ngb;
			values = linAL[i];
			
			HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			
			++num;
		}
		
		
	}
	// cout << ilower << " " << iupper << endl;
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b );
	HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(b);
	
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x );
	HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(x);
	
	int *save_rows = new int[mesh.cells.size()];
	double *BB = new double[mesh.cells.size()];
	double *RR = new double[mesh.cells.size()];
	for(int i=0; i<mesh.cells.size(); ++i){
		save_rows[i] = myStart + i;
		BB[i] = linB[i];
		RR[i] = resiVar[i];
		// rows = myStart + i;
		// HYPRE_IJVectorSetValues(b, 1, &rows, &linB[i] );
		// HYPRE_IJVectorSetValues(x, 1, &rows, &resiVar[i] );
	}
	HYPRE_IJVectorSetValues(b, ncells, save_rows, BB );
	HYPRE_IJVectorSetValues(x, ncells, save_rows, RR );
	
	delete[] save_rows;
	delete[] BB;
	delete[] RR;
	
	HYPRE_IJVectorAssemble( b);
	HYPRE_IJVectorAssemble( x);
	
	HYPRE_IJVectorGetObject( b, (void **) &par_b);
	HYPRE_IJVectorGetObject( x, (void **) &par_x);


	ierr = 0;
	if(solver == "cg"){
		// PCG + AMG
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &sol);
		HYPRE_ParCSRPCGSetTol(sol, tolerance);
		HYPRE_ParCSRPCGSetMaxIter(sol, maxIter);

		HYPRE_BoomerAMGCreate(&pre);

		HYPRE_PCGSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
							
		HYPRE_ParCSRPCGSetup(sol, parcsr_A, par_b, par_x);
		ierr = HYPRE_ParCSRPCGSolve(sol, parcsr_A, par_b, par_x);
		 
		
		HYPRE_ParCSRPCGDestroy( sol );
		HYPRE_BoomerAMGDestroy( pre );

	}
	else if(solver == "gmres"){
		// GMRES + AMG
		HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol);
		HYPRE_ParCSRFlexGMRESSetTol(sol, tolerance);
		HYPRE_ParCSRFlexGMRESSetMaxIter(sol, maxIter);
		
		HYPRE_BoomerAMGCreate(&pre);
		HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							(HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
		
		HYPRE_ParCSRFlexGMRESSetup(sol, parcsr_A, par_b, par_x);
		ierr = HYPRE_ParCSRFlexGMRESSolve(sol, parcsr_A, par_b, par_x);
		
		HYPRE_ParCSRFlexGMRESDestroy(sol);
		HYPRE_BoomerAMGDestroy( pre );
	
	}
	else{
		// AMG
		HYPRE_BoomerAMGCreate(&sol);
		HYPRE_BoomerAMGSetTol(sol, tolerance);
		HYPRE_BoomerAMGSetMaxIter(sol, maxIter);
		// HYPRE_BoomerAMGSetMinCoarseSize(sol, 2);
		// HYPRE_BoomerAMGSetMaxLevels(sol, 25);
		// HYPRE_BoomerAMGSetNumSweeps(sol, 1);
		// HYPRE_BoomerAMGSetRelaxType(sol, 18);
		
		HYPRE_BoomerAMGSetup(sol, parcsr_A, par_b, par_x);
		ierr = HYPRE_BoomerAMGSolve(sol, parcsr_A, par_b, par_x);
		
		HYPRE_BoomerAMGDestroy( sol );
	}
	





	// if(ierr){
		// if(rank==0) cerr << "| #Error : HYPRE error" << endl;
		// MPI_Barrier(MPI_COMM_WORLD);
		// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// }
	

	for(int i=0; i<mesh.cells.size(); ++i){
		rows = myStart + i;
		HYPRE_IJVectorGetValues(x, 1, &rows, &resiVar[i] );
	}


	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(b);
	HYPRE_IJVectorDestroy(x);
	
	HYPRE_Finalize();
	
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
}