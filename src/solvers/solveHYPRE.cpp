#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include "mpi.h"

// #include "HYPRE.h"
// #include "HYPRE_parcsr_ls.h"
// #include "HYPRE_krylov.h"


// #include <petscvec.h>
// #include <petscmat.h>
// #include <petscksp.h>



void SEMO_Solvers_Builder::solveHypre_Test_Coupled(
	SEMO_Mesh_Builder& mesh,
	vector<double>& resiVar, 
	vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	vector<double>& B_vals,
	int ncells, int ncellTot,
	string solver, double tolerance, double relTol, string preconditioner,
	int maxIter
	){
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 

	// // ofstream outputFile;
	// // string filenamePlot = "Amat." + to_string(rank) + ".txt";
	
	// // outputFile.open(filenamePlot);
	
	// // for(int i=0; i<A_rows.size(); ++i){
		// // outputFile << A_rows[i] << " " << A_cols[i] << " " << A_vals[i] << " " << endl;
	// // }
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
								 
	// int ierr;
	
	// int myStart = mesh.startCellGlobal;
	// int myEnd   = mesh.startCellGlobal + ncells - 1;
	

	// HYPRE_IJMatrix A_hypre;
	// HYPRE_ParCSRMatrix parcsr_A;
	// HYPRE_IJVector b_hypre;
	// HYPRE_ParVector par_b;
	// HYPRE_IJVector x_hypre;
	// HYPRE_ParVector par_x;
	
	// HYPRE_Solver sol, pre;
	
	// HYPRE_Init();
	
	// HYPRE_IJMatrixCreate(MPI_COMM_WORLD, myStart, myEnd, myStart, myEnd, &A_hypre);
	// HYPRE_IJMatrixSetObjectType(A_hypre, HYPRE_PARCSR);
	// HYPRE_IJMatrixInitialize(A_hypre);
	
	// for(int i=0; i<A_vals.size(); ++i){
		// int il = 1;
		// HYPRE_IJMatrixSetValues(A_hypre,1, &il, &A_rows[i], &A_cols[i], &A_vals[i]);
	// }
	
	
	// HYPRE_IJMatrixAssemble(A_hypre);
	// HYPRE_IJMatrixGetObject(A_hypre, (void**) &parcsr_A);
	

	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, myStart, myEnd, &b_hypre );
	// HYPRE_IJVectorSetObjectType(b_hypre, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(b_hypre);
	
	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, myStart, myEnd, &x_hypre );
	// HYPRE_IJVectorSetObjectType(x_hypre, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(x_hypre);
	
	// // // cout << "11" << endl;
	// // int *save_rows = new int[n];
	// // double *BB = new double[n];
	// // double *RR = new double[n];
	// // for(int i=0; i<n; ++i){
		// // save_rows[i] = i;
		
		// // double Bvalues = B_vals[i];
		// // double Xvalues = 0.0;
		// // BB[i] = 1.0;
		
		// // RR[i] = 0.0;
		// // // rows = myStart + i;
		// // // HYPRE_IJVectorSetValues(b, 1, &rows, &linB[i] );
		// // // HYPRE_IJVectorSetValues(x, 1, &rows, &resiVar[i] );
	// // }
	// // HYPRE_IJVectorSetValues(b_hypre, n, save_rows, BB );
	// // HYPRE_IJVectorSetValues(x_hypre, n, save_rows, RR );
	
	// // delete[] save_rows;
	// // delete[] BB;
	// // delete[] RR;
	
	// for(int i=0; i<ncells; ++i){
		// int cols = myStart + i;
		// double Bvalues = B_vals[i];
		// double Xvalues = 0.0;
		// HYPRE_IJVectorSetValues(b_hypre,1,&cols,&Bvalues);
		// HYPRE_IJVectorSetValues(x_hypre,1,&cols,&Xvalues);
	// }
	
	
	
	// HYPRE_IJVectorAssemble( b_hypre);
	// HYPRE_IJVectorAssemble( x_hypre);
	
	// HYPRE_IJVectorGetObject( b_hypre, (void **) &par_b);
	// HYPRE_IJVectorGetObject( x_hypre, (void **) &par_x);


	// // // AMG
		// HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol);
		// HYPRE_ParCSRFlexGMRESSetTol(sol, tolerance);
		// HYPRE_ParCSRFlexGMRESSetMaxIter(sol, maxIter);
		
		// // HYPRE_BoomerAMGCreate(&pre);
		// // HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
		
		// HYPRE_ParCSRFlexGMRESSetup(sol, parcsr_A, par_b, par_x);
		// ierr = HYPRE_ParCSRFlexGMRESSolve(sol, parcsr_A, par_b, par_x);
		
		// HYPRE_ParCSRFlexGMRESDestroy(sol);
		// // HYPRE_BoomerAMGDestroy( pre );
		

		// // // AMG
		// // HYPRE_BoomerAMGCreate(&sol);
		// // HYPRE_BoomerAMGSetTol(sol, tolerance);
		// // HYPRE_BoomerAMGSetMaxIter(sol, maxIter);
		// // // HYPRE_BoomerAMGSetMinCoarseSize(sol, 2);
		// // // HYPRE_BoomerAMGSetMaxLevels(sol, 25);
		// // // HYPRE_BoomerAMGSetNumSweeps(sol, 1);
		// // // HYPRE_BoomerAMGSetRelaxType(sol, 18);
		
		// // HYPRE_BoomerAMGSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_BoomerAMGSolve(sol, parcsr_A, par_b, par_x);
		
		// // HYPRE_BoomerAMGDestroy( sol );
		
	

	// for(int i=0; i<ncells; ++i){
		// int cols = myStart + i;
		// double values;
		// HYPRE_IJVectorGetValues(x_hypre,1,&cols,&values);
		// resiVar[i] = values;
	// }


	// HYPRE_IJMatrixDestroy(A_hypre);
	// HYPRE_IJVectorDestroy(b_hypre);
	// HYPRE_IJVectorDestroy(x_hypre);
	
	// HYPRE_Finalize();
	
	
	
	
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
	
	
	
	
}






// void SEMO_Solvers_Builder::solveHYPRE_Coupled(
	// SEMO_Mesh_Builder& mesh,
	// vector<double>& resiVar, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B_vals,
	// int ncells, int ncellTot,
	// string solver, double tolerance, double relTol, string preconditioner,
	// int maxIter
	// ){
		
		
								 
	// int ierr;
	
	// ierr = PetscInitialize(NULL,NULL,(char*)0,NULL);




	// Vec              x,b,u;
	// Mat              A;
	// KSP              ksp;
	// PC               pc;
	  
	  
	  
	// MatCreate(PETSC_COMM_WORLD,&A);
	// MatSetSizes(A,ncells,ncells,ncellTot,ncellTot);
	// MatSetFromOptions(A);
	
	
    // MatSeqAIJSetPreallocation(A, 60, NULL);
    // MatMPIAIJSetPreallocation(A, 60, NULL, 59, NULL);

	// int myStart, myEnd;
    // MatGetOwnershipRange(A, &myStart, &myEnd);
	
	// for(int i=0; i<A_vals.size(); ++i){
		// MatSetValues(A,1,&A_rows[i],1,&A_cols[i],&A_vals[i],INSERT_VALUES);
	// }
	
	// VecCreate(PETSC_COMM_WORLD,&x);
	// VecSetSizes(x,ncells,ncellTot);
	// VecSetFromOptions(x);
	// VecDuplicate(x,&b);
	
	// for(int i=0; i<ncells; ++i){
		// int cols = myStart + i;
		// double Bvalues = B_vals[i];
		// double Xvalues = 0.0;
		// VecSetValues(b,1,&cols,&Bvalues,ADD_VALUES);
		// VecSetValues(x,1,&cols,&Xvalues,ADD_VALUES);
	// }
	
 
	// MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	// MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	// VecAssemblyBegin(x);
	// VecAssemblyEnd(x);
	// VecAssemblyBegin(b);
	// VecAssemblyEnd(b);
	
	
	// KSPCreate(PETSC_COMM_WORLD,&ksp);

	// KSPSetOperators(ksp,A,A);
	
		// // KSPSetType(ksp,KSPGMRES);
		// // KSPSetType(ksp,KSPCGNE);
		// // KSPSetType(ksp,KSPPIPEFGMRES);
		// KSPSetType(ksp,KSPFGMRES);
		// // KSPSetType(ksp,KSPLGMRES);
		// // KSPSetType(ksp,KSPDGMRES);
		// // KSPSetType(ksp,KSPPGMRES);
		// // KSPSetType(ksp,KSPTCQMR);
		// // KSPSetType(ksp,KSPBCGS); // BICGSTAB
		// // KSPSetType(ksp,KSPIBCGS);
		// // KSPSetType(ksp,KSPFBCGS);
		// // KSPSetType(ksp,KSPFBCGSR);
		// // KSPSetType(ksp,KSPBCGSL);
		// // KSPSetType(ksp,KSPPIPEBCGS);
		// // KSPSetType(ksp,KSPCGS);
		// // KSPSetType(ksp,KSPTFQMR);
		// // KSPSetType(ksp,KSPLSQR);
		// // KSPSetType(ksp,KSPBICG);
		// // KSPSetType(ksp,KSPLCD );
		// // KSPSetType(ksp,KSPGCR);
		// // KSPSetType(ksp,KSPPIPEGCR);
		// // KSPSetType(ksp,KSPCGLS);
				 

	// KSPGetPC(ksp,&pc);
	
		// // PCSetType(pc,PCNONE);
		// PCSetType(pc,PCJACOBI);
		// // PCSetType(pc,PCSOR);
		// // PCSetType(pc,PCMG);
		// // PCSetType(pc,PCSVD);
		// // PCSetType(pc,PCGAMG);
		// // PCSetType(pc,PCKACZMARZ);
		// // PCSetType(pc,PCHMG);
	
	
		
	// KSPSetTolerances(ksp, 1.e-200, 1.e-4, 1.e+5, 50);
	

	// KSPSetFromOptions(ksp);
	// KSPSetUp(ksp);

	
	// ierr = KSPSolve(ksp,b,x);

	// // if(permute) {
		// // VecPermute(x,rowperm,PETSC_TRUE);
	// // }
	
	// for(int i=0; i<ncells; ++i){
		// int cols = myStart + i;
		// double values;
		// VecGetValues(x,1,&cols,&values);
		
		// resiVar[i] = values;
		
		// // cout << i << " " << values << endl;
	// }
	// // cout << endl;
	
	

	// KSPDestroy(&ksp);
	// VecDestroy(&x);
	// VecDestroy(&b);
	// MatDestroy(&A);
	
	// PetscFinalize();
	
// }



// void SEMO_Solvers_Builder::solveHYPRE_Coupled(
	// SEMO_Mesh_Builder& mesh,
	// vector<double>& resiVar, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B_vals,
	// int ncells, int ncellTot,
	// string solver, double tolerance, double relTol, string preconditioner,
	// int maxIter
	// ){
		
		
		
    // int    n = 8;
    // int    ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
    // int    ja[20] = { 0,    2,       5, 6, 
                         // 1, 2,    4,
                            // 2,             7,
                               // 3,       6,
                         // 1,
                            // 2,       5,    7,
                         // 1,             6,
                            // 2,          6, 7 };
    // double  a[20] = { 7.0,      1.0,           2.0, 7.0, 
                          // -4.0, 8.0,      2.0,
                                // 1.0,                     5.0,
                                     // 7.0,           9.0,
                          // -4.0,
                                // 7.0,           3.0,      8.0,
                           // 1.0,                    11.0,
                               // -3.0,                2.0, 5.0 };
								 
	// int ierr;
	
	// ierr = PetscInitialize(NULL,NULL,(char*)0,NULL);




	// Vec              x,b,u;
	// Mat              A;
	// KSP              ksp;
	// PC               pc;
	  
	  
	  
	// MatCreate(PETSC_COMM_WORLD,&A);
	// MatSetSizes(A,n,n,n,n);
	// MatSetFromOptions(A);
	
	
    // MatSeqAIJSetPreallocation(A, 60, NULL);
    // MatMPIAIJSetPreallocation(A, 60, NULL, 59, NULL);

	// int myStart, myEnd;
    // MatGetOwnershipRange(A, &myStart, &myEnd);
	
	// int i1 = 1;

	// for(int i=0; i<n; ++i){
		// for(int j=ia[i]; j<ia[i+1]; ++j){
			// MatSetValues(A,i1,&i,i1,&ja[j],&a[j],INSERT_VALUES);
		// }
	// }
	
	// VecCreate(PETSC_COMM_WORLD,&x);
	// VecSetSizes(x,n,n);
	// VecSetFromOptions(x);
	// VecDuplicate(x,&b);
	
	// for(int i=0; i<n; ++i){
		// int cols = myStart + i;
		// double Bvalues = 1.0;
		// double Xvalues = 0.0;
		// VecSetValues(b,i1,&cols,&Bvalues,ADD_VALUES);
		// VecSetValues(x,i1,&cols,&Xvalues,ADD_VALUES);
	// }
	
 
	// MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	// MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	// VecAssemblyBegin(x);
	// VecAssemblyEnd(x);
	// VecAssemblyBegin(b);
	// VecAssemblyEnd(b);
	

	// // // reordering (ex, the reverse cuthill-mckee)
	// // bool permute = true;
	// // if(permute) {
		// // char ordering[256] = MATORDERINGRCM;
		// // Mat Aperm;
		// // MatGetOrdering(A,ordering,&rowperm,&colperm);
		// // MatPermute(A,rowperm,colperm,&Aperm);
		// // VecPermute(b,colperm,PETSC_FALSE);
		// // MatDestroy(&A);
		// // A = Aperm;               /* Replace original operator with permuted version */
	// // }
	
	
	
	// KSPCreate(PETSC_COMM_WORLD,&ksp);

	// KSPSetOperators(ksp,A,A);
	
		// // KSPSetType(ksp,KSPGMRES);
		// // KSPSetType(ksp,KSPCGNE);
		// // KSPSetType(ksp,KSPPIPEFGMRES);
		// // KSPSetType(ksp,KSPFGMRES);
		// // KSPSetType(ksp,KSPLGMRES);
		// // KSPSetType(ksp,KSPDGMRES);
		// // KSPSetType(ksp,KSPPGMRES);
		// // KSPSetType(ksp,KSPTCQMR);
		// KSPSetType(ksp,KSPBCGS); // BICGSTAB
		// // KSPSetType(ksp,KSPIBCGS);
		// // KSPSetType(ksp,KSPFBCGS);
		// // KSPSetType(ksp,KSPFBCGSR);
		// // KSPSetType(ksp,KSPBCGSL);
		// // KSPSetType(ksp,KSPPIPEBCGS);
		// // KSPSetType(ksp,KSPCGS);
		// // KSPSetType(ksp,KSPTFQMR);
		// // KSPSetType(ksp,KSPLSQR);
		// // KSPSetType(ksp,KSPBICG);
		// // KSPSetType(ksp,KSPLCD );
		// // KSPSetType(ksp,KSPGCR);
		// // KSPSetType(ksp,KSPPIPEGCR);
		// // KSPSetType(ksp,KSPCGLS);
				 

	// KSPGetPC(ksp,&pc);
	
		// // PCSetType(pc,PCNONE);
		// // PCSetType(pc,PCJACOBI);
		// // PCSetType(pc,PCSOR);
		// PCSetType(pc,PCMG);
		// // PCSetType(pc,PCSVD);
		// // PCSetType(pc,PCGAMG);
		// // PCSetType(pc,PCKACZMARZ);
		// // PCSetType(pc,PCHMG);
	
	
		
	// // KSPSetTolerances(ksp, 1.e-6, 1.e-6, 1.e+5, 6);
	

	// KSPSetFromOptions(ksp);
	// KSPSetUp(ksp);

	
	// ierr = KSPSolve(ksp,b,x);

	// // if(permute) {
		// // VecPermute(x,rowperm,PETSC_TRUE);
	// // }
	
	// for(int i=0; i<n; ++i){
		// int cols = myStart + i;
		// double values;
		// VecGetValues(x,i1,&cols,&values);
		
		// cout << i << " " << values << endl;
	// }
	// cout << endl;


	// HYPRE_IJMatrix A_hypre;
	// HYPRE_ParCSRMatrix parcsr_A;
	// HYPRE_IJVector b_hypre;
	// HYPRE_ParVector par_b;
	// HYPRE_IJVector x_hypre;
	// HYPRE_ParVector par_x;
	
	// HYPRE_Solver sol_hypre, pre_hypre;
	
	// HYPRE_Init();
	
	// HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, n, 0, n, &A_hypre);
	// HYPRE_IJMatrixSetObjectType(A_hypre, HYPRE_PARCSR);
	// HYPRE_IJMatrixInitialize(A_hypre);
	
	// // int i1 = 1;

	// for(int i=0; i<n; ++i){
		// for(int j=ia[i]; j<ia[i+1]; ++j){
			// HYPRE_IJMatrixSetValues(A_hypre, 1, &i1, &i, &ja[j], &a[j]);
		// }
	// }
	
	
	// // for(int i=0; i<n; ++i){
		// // double values = 1.0;
		// // HYPRE_IJMatrixSetValues(A, 1, &i1, &i, &i, &values);
	// // }
	
	
	// HYPRE_IJMatrixAssemble(A_hypre);
	// HYPRE_IJMatrixGetObject(A_hypre, (void**) &parcsr_A);
	

	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n, &b_hypre );
	// HYPRE_IJVectorSetObjectType(b_hypre, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(b_hypre);
	
	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n, &x_hypre );
	// HYPRE_IJVectorSetObjectType(x_hypre, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(x_hypre);
	
	// // cout << "11" << endl;
	// int *save_rows = new int[n];
	// double *BB = new double[n];
	// double *RR = new double[n];
	// for(int i=0; i<n; ++i){
		// save_rows[i] = i;
		
		// BB[i] = 1.0;
		
		// RR[i] = 0.0;
		// // rows = myStart + i;
		// // HYPRE_IJVectorSetValues(b, 1, &rows, &linB[i] );
		// // HYPRE_IJVectorSetValues(x, 1, &rows, &resiVar[i] );
	// }
	// HYPRE_IJVectorSetValues(b_hypre, n, save_rows, BB );
	// HYPRE_IJVectorSetValues(x_hypre, n, save_rows, RR );
	
	
	// delete[] save_rows;
	// delete[] BB;
	// delete[] RR;
	
	// HYPRE_IJVectorAssemble( b_hypre);
	// HYPRE_IJVectorAssemble( x_hypre);
	
	// HYPRE_IJVectorGetObject( b_hypre, (void **) &par_b);
	// HYPRE_IJVectorGetObject( x_hypre, (void **) &par_x);


	// // AMG
	// ierr=0;
		// HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol_hypre);
		// HYPRE_ParCSRFlexGMRESSetTol(sol_hypre, tolerance);
		// HYPRE_ParCSRFlexGMRESSetMaxIter(sol_hypre, maxIter);
		
		// // HYPRE_BoomerAMGCreate(&pre);
		// // HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
		
		// HYPRE_ParCSRFlexGMRESSetup(sol_hypre, parcsr_A, par_b, par_x);
		// ierr = HYPRE_ParCSRFlexGMRESSolve(sol_hypre, parcsr_A, par_b, par_x);
		
		// HYPRE_ParCSRFlexGMRESDestroy(sol_hypre);
		// // HYPRE_BoomerAMGDestroy( pre );
	

	// for(int i=0; i<n; ++i){
		// double valuess = 0.0;
		// HYPRE_IJVectorGetValues(x_hypre, 1, &i, &valuess );
		
		// cout << i << " " << valuess << endl;
	// }


	// HYPRE_IJMatrixDestroy(A_hypre);
	// HYPRE_IJVectorDestroy(b_hypre);
	// HYPRE_IJVectorDestroy(x_hypre);
	
	// HYPRE_Finalize();
	
	
	
	
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
		
		
		
		
		
		
		
		
		
	// // int ierr;
	
	// // int rank = MPI::COMM_WORLD.Get_rank();
	// // int size = MPI::COMM_WORLD.Get_size();
	
	// // // int ncells = linA.size();
	// // vector<int> procNcells(size,0);
	// // MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	// // int myStart=0;
	// // for(int i=0; i<rank; ++i){
		// // myStart += procNcells[i];
	// // }
	
	// // // int ncellTot;
	// // MPI_Allreduce(&ncells, &ncellTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	// // vector<int> procStart(size,0);
	// // MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	
	
	
	// // HYPRE_IJMatrix A;
	// // HYPRE_ParCSRMatrix parcsr_A;
	// // HYPRE_IJVector b;
	// // HYPRE_ParVector par_b;
	// // HYPRE_IJVector x;
	// // HYPRE_ParVector par_x;
	
	// // HYPRE_Solver sol, pre;
	
	// // // if(HYPRE_Init()) cout << "# Error 1" << endl;
	
	// // HYPRE_Init();
	
	// // int ilower = myStart;
	// // int iupper = ilower + procNcells[rank] - 1;
	
	// // // if(HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A)) cout << "# Error 1" << endl;
	// // // if(HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR)) cout << "# Error 1" << endl;
	// // // if(HYPRE_IJMatrixInitialize(A)) cout << "# Error 1" << endl;
	
	// // HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	// // HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	// // HYPRE_IJMatrixInitialize(A);
	
	// // int cols;
	// // int rows;
	// // double values;
	
	// // int i1 = 1;

	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // // cols = myStart + i;
		// // // values = linA[i];
		// // // HYPRE_IJMatrixSetValues(A, 1, &i1, &cols, &cols, &values);
		// // HYPRE_IJMatrixSetValues(A, 1, &i1, &A_rows[i], &A_cols[i], &A_vals[i]);
		
	// // }
	
	
	
	
	// // vector<int> neighbProcNo;
	// // vector<int> ownerNo;
	// // vector<int> neighbNo;
	
	// // for(auto& boundary : mesh.boundary){
		// // if(boundary.neighbProcNo != -1){
			// // int str=boundary.startFace;
			// // int end=str+boundary.nFaces;
			// // for(int i=str; i<end; ++i){
				// // neighbProcNo.push_back(boundary.neighbProcNo);
				// // ownerNo.push_back(mesh.faces[i].owner);
				// // // cout << mesh.faces[i].owner << endl;
			// // }
		// // }
	// // }
	
	// // if(size>1){

		// // SEMO_MPI_Builder mpi;
		
		// // mpi.setProcsFaceDatas(
					// // ownerNo, neighbNo,
					// // mesh.countsProcFaces, mesh.countsProcFaces, 
					// // mesh.displsProcFaces, mesh.displsProcFaces);
	
	// // }
	
	// // int num=0;
	// // // for(int i=0; i<mesh.faces.size(); ++i){
		
		// // // if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			
			// // // int own = mesh.faces[i].owner;
			// // // int ngb = mesh.faces[i].neighbour;
			
			// // // cols = myStart + ngb;
			// // // rows = myStart + own;
			// // // values = linAR[i];
			
			// // // HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			// // // cols = myStart + own;
			// // // rows = myStart + ngb;
			// // // values = linAL[i];
			
			// // // HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
		// // // }
		// // // if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			
			// // // int own = myStart + mesh.faces[i].owner;
			// // // int ngb = procStart[neighbProcNo[num]] + neighbNo[num];
			
			// // // // cout << ngb << " " << own << " " << ncellTot << endl;
			
			// // // cols = ngb;
			// // // rows = own;
			// // // values = linAR[i];
			
			// // // HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			// // // cols = own;
			// // // rows = ngb;
			// // // values = linAL[i];
			
			// // // HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			
			// // // ++num;
		// // // }
		
		
	// // // }
	// // // cout << ilower << " " << iupper << endl;
	
	// // // MPI_Barrier(MPI_COMM_WORLD);
	// // // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);



	// // HYPRE_IJMatrixAssemble(A);
	// // HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	

	// // HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b );
	// // HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	// // HYPRE_IJVectorInitialize(b);
	
	// // HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x );
	// // HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	// // HYPRE_IJVectorInitialize(x);
	
	// // // cout << "11" << endl;
	// // int *save_rows = new int[ncells];
	// // double *BB = new double[ncells];
	// // double *RR = new double[ncells];
	// // for(int i=0; i<ncells; ++i){
		// // save_rows[i] = myStart + i;
		
		// // BB[i] = B_vals[i];
		
		// // RR[i] = resiVar[i];
		// // // rows = myStart + i;
		// // // HYPRE_IJVectorSetValues(b, 1, &rows, &linB[i] );
		// // // HYPRE_IJVectorSetValues(x, 1, &rows, &resiVar[i] );
	// // }
	// // HYPRE_IJVectorSetValues(b, ncells, save_rows, BB );
	// // HYPRE_IJVectorSetValues(x, ncells, save_rows, RR );
	
	
	// // delete[] save_rows;
	// // delete[] BB;
	// // delete[] RR;
	
	// // HYPRE_IJVectorAssemble( b);
	// // HYPRE_IJVectorAssemble( x);
	
	// // HYPRE_IJVectorGetObject( b, (void **) &par_b);
	// // HYPRE_IJVectorGetObject( x, (void **) &par_x);


	// // ierr = 0;
	// // if(solver == "cg"){
		// // // PCG + AMG
		// // HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &sol);
		// // HYPRE_ParCSRPCGSetTol(sol, tolerance);
		// // HYPRE_ParCSRPCGSetMaxIter(sol, maxIter);

		// // HYPRE_BoomerAMGCreate(&pre);

		// // HYPRE_PCGSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
							
		// // HYPRE_ParCSRPCGSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_ParCSRPCGSolve(sol, parcsr_A, par_b, par_x);
		 
		
		// // HYPRE_ParCSRPCGDestroy( sol );
		// // HYPRE_BoomerAMGDestroy( pre );

	// // }
	// // else if(solver == "gmres"){
		// // // GMRES + AMG
		// // HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol);
		// // HYPRE_ParCSRFlexGMRESSetTol(sol, tolerance);
		// // HYPRE_ParCSRFlexGMRESSetMaxIter(sol, maxIter);
		
		// // HYPRE_BoomerAMGCreate(&pre);
		// // HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
		
		// // HYPRE_ParCSRFlexGMRESSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_ParCSRFlexGMRESSolve(sol, parcsr_A, par_b, par_x);
		
		// // HYPRE_ParCSRFlexGMRESDestroy(sol);
		// // HYPRE_BoomerAMGDestroy( pre );
	
	// // }
	// // else{
		// // // AMG
		// // HYPRE_BoomerAMGCreate(&sol);
		// // HYPRE_BoomerAMGSetTol(sol, tolerance);
		// // HYPRE_BoomerAMGSetMaxIter(sol, maxIter);
		// // // HYPRE_BoomerAMGSetMinCoarseSize(sol, 2);
		// // // HYPRE_BoomerAMGSetMaxLevels(sol, 25);
		// // // HYPRE_BoomerAMGSetNumSweeps(sol, 1);
		// // // HYPRE_BoomerAMGSetRelaxType(sol, 18);
		
		// // HYPRE_BoomerAMGSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_BoomerAMGSolve(sol, parcsr_A, par_b, par_x);
		
		// // HYPRE_BoomerAMGDestroy( sol );
	// // }
	





	// // // if(ierr){
		// // // if(rank==0) cerr << "| #Error : HYPRE error" << endl;
		// // // MPI_Barrier(MPI_COMM_WORLD);
		// // // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // // }
	

	// // for(int i=0; i<4*mesh.cells.size(); ++i){
		// // rows = myStart + i;
		// // HYPRE_IJVectorGetValues(x, 1, &rows, &resiVar[i] );
	// // }


	// // HYPRE_IJMatrixDestroy(A);
	// // HYPRE_IJVectorDestroy(b);
	// // HYPRE_IJVectorDestroy(x);
	
	// // HYPRE_Finalize();
	
	
	// // // MPI_Barrier(MPI_COMM_WORLD);
	// // // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
// }






void SEMO_Solvers_Builder::solveHYPRE(
	SEMO_Mesh_Builder& mesh,
	vector<double>& resiVar, 
	vector<double>& linA, vector<double>& linAL, vector<double>& linAR, 
	vector<double>& linB,
	string solver, double tolerance, double relTol, string preconditioner,
	int maxIter
	){
		
	// int ierr;
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	// int ncells = linA.size();
	// vector<int> procNcells(size,0);
	// MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	// int myStart=0;
	// for(int i=0; i<rank; ++i){
		// myStart += procNcells[i];
	// }
	
	// int ncellTot;
	// MPI_Allreduce(&ncells, &ncellTot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
	// vector<int> procStart(size,0);
	// MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,MPI_COMM_WORLD);
	
	
	
	
	// HYPRE_IJMatrix A;
	// HYPRE_ParCSRMatrix parcsr_A;
	// HYPRE_IJVector b;
	// HYPRE_ParVector par_b;
	// HYPRE_IJVector x;
	// HYPRE_ParVector par_x;
	
	// HYPRE_Solver sol, pre;
	
	// // if(HYPRE_Init()) cout << "# Error 1" << endl;
	
	// HYPRE_Init();
	
	// int ilower = myStart;
	// int iupper = ilower + procNcells[rank] - 1;
	
	// // if(HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A)) cout << "# Error 1" << endl;
	// // if(HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR)) cout << "# Error 1" << endl;
	// // if(HYPRE_IJMatrixInitialize(A)) cout << "# Error 1" << endl;
	
	// HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	// HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	// HYPRE_IJMatrixInitialize(A);
	
	// int cols;
	// int rows;
	// double values;
	
	// int i1 = 1;

	// for(int i=0; i<mesh.cells.size(); ++i){
		// cols = myStart + i;
		// values = linA[i];
		// HYPRE_IJMatrixSetValues(A, 1, &i1, &cols, &cols, &values);
		
	// }
	
	
	
	
	// vector<int> neighbProcNo;
	// vector<int> ownerNo;
	// vector<int> neighbNo;
	
	// for(auto& boundary : mesh.boundary){
		// if(boundary.neighbProcNo != -1){
			// int str=boundary.startFace;
			// int end=str+boundary.nFaces;
			// for(int i=str; i<end; ++i){
				// neighbProcNo.push_back(boundary.neighbProcNo);
				// ownerNo.push_back(mesh.faces[i].owner);
				// // cout << mesh.faces[i].owner << endl;
			// }
		// }
	// }
	
	// if(size>1){

		// SEMO_MPI_Builder mpi;
		
		// mpi.setProcsFaceDatas(
					// ownerNo, neighbNo,
					// mesh.countsProcFaces, mesh.countsProcFaces, 
					// mesh.displsProcFaces, mesh.displsProcFaces);
	
	// }
	
	// int num=0;
	// for(int i=0; i<mesh.faces.size(); ++i){
		
		// if(mesh.faces[i].getType() == SEMO_Types::INTERNAL_FACE){
			
			// int own = mesh.faces[i].owner;
			// int ngb = mesh.faces[i].neighbour;
			
			// cols = myStart + ngb;
			// rows = myStart + own;
			// values = linAR[i];
			
			// HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			// cols = myStart + own;
			// rows = myStart + ngb;
			// values = linAL[i];
			
			// HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
		// }
		// if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			
			// int own = myStart + mesh.faces[i].owner;
			// int ngb = procStart[neighbProcNo[num]] + neighbNo[num];
			
			// // cout << ngb << " " << own << " " << ncellTot << endl;
			
			// cols = ngb;
			// rows = own;
			// values = linAR[i];
			
			// HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			// cols = own;
			// rows = ngb;
			// values = linAL[i];
			
			// HYPRE_IJMatrixSetValues(A, 1, &i1, &rows, &cols, &values);
			
			
			// ++num;
		// }
		
		
	// }
	// // cout << ilower << " " << iupper << endl;
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	// HYPRE_IJMatrixAssemble(A);
	// HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	

	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b );
	// HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(b);
	
	// HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x );
	// HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	// HYPRE_IJVectorInitialize(x);
	
	// int *save_rows = new int[mesh.cells.size()];
	// double *BB = new double[mesh.cells.size()];
	// double *RR = new double[mesh.cells.size()];
	// for(int i=0; i<mesh.cells.size(); ++i){
		// save_rows[i] = myStart + i;
		// BB[i] = linB[i];
		// RR[i] = resiVar[i];
		// // rows = myStart + i;
		// // HYPRE_IJVectorSetValues(b, 1, &rows, &linB[i] );
		// // HYPRE_IJVectorSetValues(x, 1, &rows, &resiVar[i] );
	// }
	// HYPRE_IJVectorSetValues(b, ncells, save_rows, BB );
	// HYPRE_IJVectorSetValues(x, ncells, save_rows, RR );
	
	// delete[] save_rows;
	// delete[] BB;
	// delete[] RR;
	
	// HYPRE_IJVectorAssemble( b);
	// HYPRE_IJVectorAssemble( x);
	
	// HYPRE_IJVectorGetObject( b, (void **) &par_b);
	// HYPRE_IJVectorGetObject( x, (void **) &par_x);


	// ierr = 0;
	// // if(solver == "cg"){
		// // PCG + AMG
		
		
		
		// // HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &sol);
		// // HYPRE_ParCSRPCGSetTol(sol, tolerance);
		// // HYPRE_ParCSRPCGSetMaxIter(sol, maxIter);
        // // HYPRE_ParCSRPCGSetTwoNorm(sol, 1);         /* use the two norm as the stopping criteria */
        // // HYPRE_ParCSRPCGSetLogging(sol, 1); /* needed to get run info later */


		// HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol);
		// HYPRE_ParCSRFlexGMRESSetTol(sol, tolerance);
		// HYPRE_ParCSRFlexGMRESSetMaxIter(sol, maxIter);
        // HYPRE_ParCSRFlexGMRESSetLogging(sol, 1);





		// HYPRE_BoomerAMGCreate(&pre);  











		// // AMG coarsening options:
		// int coarsen_type = 10; // 10 = HMIS, 8 = PMIS, 6 = Falgout, 0 = CLJP
		// int agg_levels = 1;    // number of aggressive coarsening levels
		// double theta = 0.25;   // strength threshold: 0.25, 0.5, 0.8

		// // AMG interpolation options:
		// int interp_type = 6; // 6 = extended+i, 0 = classical
		// int Pmax = 4;        // max number of elements per row in P

		// // AMG relaxation options:
		// int relax_type = 8;   // 8 = l1-GS, 6 = symm. GS, 3 = GS, 18 = l1-Jacobi
		// int relax_sweeps = 1; // relaxation sweeps on each level

		// // Additional options:
		// int print_level = 0; // print AMG iterations? 1 = no, 2 = yes
		// int max_levels = 25; // max number of levels in AMG hierarchy

		// HYPRE_BoomerAMGSetCoarsenType(pre, coarsen_type);
		// HYPRE_BoomerAMGSetAggNumLevels(pre, agg_levels);
		// HYPRE_BoomerAMGSetRelaxType(pre, relax_type);
		// HYPRE_BoomerAMGSetNumSweeps(pre, relax_sweeps);
		// HYPRE_BoomerAMGSetStrongThreshold(pre, theta);
		// HYPRE_BoomerAMGSetInterpType(pre, interp_type);
		// HYPRE_BoomerAMGSetPMaxElmts(pre, Pmax);
		// HYPRE_BoomerAMGSetPrintLevel(pre, print_level);
		// HYPRE_BoomerAMGSetMaxLevels(pre, max_levels);

		// // Use as a preconditioner (one V-cycle, zero tolerance)
		// HYPRE_BoomerAMGSetMaxIter(pre, 1);
		// HYPRE_BoomerAMGSetTol(pre, 0.0);




		// // // Make sure the systems AMG options are set1
		// // int dim = 3;
		// // HYPRE_BoomerAMGSetNumFunctions(pre, dim);

		// // // More robust options with respect to convergence
		// // HYPRE_BoomerAMGSetAggNumLevels(pre, 0);
		// // HYPRE_BoomerAMGSetStrongThreshold(pre, 0.5);

		// // // Nodal coarsening options (nodal coarsening is required for this solver)
		// // // See hypre's new_ij driver and the paper for descriptions.
		// // int nodal = 4;        // strength reduction norm: 1, 3 or 4
		// // int nodal_diag = 1;   // diagonal in strength matrix: 0, 1 or 2
		// // int relax_coarse = 8; // smoother on the coarsest grid: 8, 99 or 29

		// // // Elasticity interpolation options
		// // int interp_vec_variant = 2;    // 1 = GM-1, 2 = GM-2, 3 = LN
		// // int q_max = 4;                 // max elements per row for each Q
		// // int smooth_interp_vectors = 1; // smooth the rigid-body modes?

		// // // Optionally pre-process the interpolation matrix through iterative weight
		// // // refinement (this is generally applicable for any system)
		// // int interp_refine = 1;

		// // // HYPRE_BoomerAMGSetNodal(pre, nodal);
		// // // HYPRE_BoomerAMGSetNodalDiag(pre, nodal_diag);
		// // // HYPRE_BoomerAMGSetCycleRelaxType(pre, relax_coarse, 3);
		// // HYPRE_BoomerAMGSetInterpVecVariant(pre, interp_vec_variant);
		// // HYPRE_BoomerAMGSetInterpVecQMax(pre, q_max);












		// // HYPRE_PCGSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
							
		// // HYPRE_ParCSRPCGSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_ParCSRPCGSolve(sol, parcsr_A, par_b, par_x);
		 
		
		// // HYPRE_ParCSRPCGDestroy( sol );
		// // HYPRE_BoomerAMGDestroy( pre );
		
		
		
							
		// HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
							
		// HYPRE_ParCSRFlexGMRESSetup(sol, parcsr_A, par_b, par_x);
		// ierr = HYPRE_ParCSRFlexGMRESSolve(sol, parcsr_A, par_b, par_x);
		
		// HYPRE_ParCSRFlexGMRESDestroy(sol);
		// HYPRE_BoomerAMGDestroy( pre );
		
		
		
		
		
		
		
		
		

	// // }
	// // else if(solver == "gmres"){
		// // // GMRES + AMG
		// // HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD,&sol);
		// // HYPRE_ParCSRFlexGMRESSetTol(sol, tolerance);
		// // HYPRE_ParCSRFlexGMRESSetMaxIter(sol, maxIter);
		
		// // HYPRE_BoomerAMGCreate(&pre);
		// // HYPRE_FlexGMRESSetPrecond(sol, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							// // (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, pre);
		
		// // HYPRE_ParCSRFlexGMRESSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_ParCSRFlexGMRESSolve(sol, parcsr_A, par_b, par_x);
		
		// // HYPRE_ParCSRFlexGMRESDestroy(sol);
		// // HYPRE_BoomerAMGDestroy( pre );
	
	// // }
	// // else{
		// // // AMG
		// // HYPRE_BoomerAMGCreate(&sol);
		// // HYPRE_BoomerAMGSetTol(sol, tolerance);
		// // HYPRE_BoomerAMGSetMaxIter(sol, maxIter);
		// // // HYPRE_BoomerAMGSetMinCoarseSize(sol, 2);
		// // // HYPRE_BoomerAMGSetMaxLevels(sol, 25);
		// // // HYPRE_BoomerAMGSetNumSweeps(sol, 1);
		// // // HYPRE_BoomerAMGSetRelaxType(sol, 18);
		
		// // HYPRE_BoomerAMGSetup(sol, parcsr_A, par_b, par_x);
		// // ierr = HYPRE_BoomerAMGSolve(sol, parcsr_A, par_b, par_x);
		
		// // HYPRE_BoomerAMGDestroy( sol );
	// // }
	





	// // if(ierr){
		// // if(rank==0) cerr << "| #Error : HYPRE error" << endl;
		// // MPI_Barrier(MPI_COMM_WORLD);
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // }
	

	// for(int i=0; i<mesh.cells.size(); ++i){
		// rows = myStart + i;
		// HYPRE_IJVectorGetValues(x, 1, &rows, &resiVar[i] );
	// }


	// HYPRE_IJMatrixDestroy(A);
	// HYPRE_IJVectorDestroy(b);
	// HYPRE_IJVectorDestroy(x);
	
	// HYPRE_Finalize();
	
	
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);

	
}