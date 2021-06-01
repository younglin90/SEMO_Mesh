#include "build.h"
#include <cmath>
#include <array>
#include <numeric>
#include "mpi.h"

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

void SEMO_Solvers_Builder::solvePETSc(
	SEMO_Mesh_Builder& mesh,
	vector<double>& resiVar, 
	vector<double>& linA, vector<double>& linAL, vector<double>& linAR, 
	vector<double>& linB,
	string solver, double tolerance, double relTol, string preconditioner,
	int maxIter
	){
		
	
	int rank = MPI::COMM_WORLD.Get_rank();
	// int size = MPI::COMM_WORLD.Get_size();
	
	Vec              x,b,u;
	Mat              A;
	KSP              ksp;
	PC               pc;
	PetscReal        norm,tol;
	PetscReal        rtol,abstol,dtol;
	PetscErrorCode   ierr;
	PetscInt i,ncells,col[3],its,i1,i2,i3,cols,rows;
	PetscInt own,ngb;
	PetscInt myStart,myEnd;
	PetscBool  flg;
	PetscMPIInt size;
	PetscScalar      none,one,value[3],values,zero;
	PetscLogStage    stages[2];
	IS rowperm = NULL, colperm = NULL;
	  
	// PetscInitialize(PETSC_NULL_CHARACTER,PETSC_NULL_CHARACTER,(char*)0,help);
	// PetscInitialize(PETSC_NULL_CHARACTER,ierr);
	
	

	vector<clock_t> startTime;
	bool checkTime = false;
	if(rank==0 && checkTime) cout << endl;
	
	if(rank==0 && checkTime) startTime.push_back(clock());
	
	
	ierr = PetscInitialize(NULL,NULL,(char*)0,NULL);
	

	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	
	none = -1.0;
	one  = 1.0;
	zero = 0.0;
	ncells    = linA.size();
	i1 = 1;
	i2 = 2;
	i3 = 3;
	
	int ncellTot;
	
	// vector<int> eachNCells(size,0);
	// MPI_Allgather(&ncells,1,MPI_INT,eachNCells.data(),1,MPI_INT,PETSC_COMM_WORLD);
	
	MPI_Allreduce(&ncells, &ncellTot,1,MPI_INTEGER,MPI_SUM,PETSC_COMM_WORLD);
	
	MatCreate(PETSC_COMM_WORLD,&A);
	MatSetSizes(A,ncells,ncells,ncellTot,ncellTot);
	MatSetFromOptions(A);
	
    MatSeqAIJSetPreallocation(A, 24, NULL);
    MatMPIAIJSetPreallocation(A, 24, NULL, 23, NULL);

    MatGetOwnershipRange(A, &myStart, &myEnd);
	
	vector<int> procStart(size,0);
	vector<int> procNcells(size,0);
	MPI_Allgather(&myStart,1,MPI_INT,procStart.data(),1,MPI_INT,PETSC_COMM_WORLD);
	MPI_Allgather(&ncells,1,MPI_INT,procNcells.data(),1,MPI_INT,PETSC_COMM_WORLD);
	
	
	// MatSetUp(A);

	// resiVar.resize(linB.size(),0.0);
	
	// cout << myStart << " " << ncells << " " << mesh.cells.size() << endl;
	for(int i=0; i<mesh.cells.size(); ++i){
		cols = myStart + i;
		values = linA[i];
		
		MatSetValues(A,i1,&cols,i1,&cols,&values,INSERT_VALUES);
		
		// resiVar[i] = 0.0;
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
			
			MatSetValues(A,i1,&rows,i1,&cols,&values,INSERT_VALUES);
			
			cols = myStart + own;
			rows = myStart + ngb;
			values = linAL[i];
			
			MatSetValues(A,i1,&rows,i1,&cols,&values,INSERT_VALUES);
		}
		if(mesh.faces[i].getType() == SEMO_Types::PROCESSOR_FACE){
			
			int own = myStart + mesh.faces[i].owner;
			int ngb = procStart[neighbProcNo[num]] + neighbNo[num];
			// cout << neighbNo[num] << endl;
			cols = ngb;
			rows = own;
			values = linAR[i];
			
			MatSetValues(A,i1,&rows,i1,&cols,&values,INSERT_VALUES);
			
			cols = own;
			rows = ngb;
			values = linAL[i];
			
			MatSetValues(A,i1,&rows,i1,&cols,&values,INSERT_VALUES);
			
			++num;
		}
		
		
	}



	VecCreate(PETSC_COMM_WORLD,&x);
	VecSetSizes(x,ncells,ncellTot);
	VecSetFromOptions(x);
	VecDuplicate(x,&b);
	


	for(int i=0; i<mesh.cells.size(); ++i){
		cols = myStart + i;
		VecSetValues(b,i1,&cols,&linB[i],ADD_VALUES);
		VecSetValues(x,i1,&cols,&resiVar[i],ADD_VALUES);
	}
	
 
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
	

	if(rank==0 && checkTime) cout<< "1st : " << (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	if(rank==0 && checkTime) startTime.pop_back();
	

	if(rank==0 && checkTime) startTime.push_back(clock());
	
	
	
	
	// reordering (ex, the reverse cuthill-mckee)
	bool permute = true;
	if(permute) {
		char ordering[256] = MATORDERINGRCM;
		Mat Aperm;
		MatGetOrdering(A,ordering,&rowperm,&colperm);
		MatPermute(A,rowperm,colperm,&Aperm);
		VecPermute(b,colperm,PETSC_FALSE);
		MatDestroy(&A);
		A = Aperm;               /* Replace original operator with permuted version */
	}
	
	
	
	KSPCreate(PETSC_COMM_WORLD,&ksp);

	KSPSetOperators(ksp,A,A);
	
	if(solver == "gmres"){
		KSPSetType(ksp,KSPGMRES);
	} 
	else if(solver == "richardson"){
		KSPSetType(ksp,KSPRICHARDSON);
	} 
	else if(solver == "chebyshev"){
		KSPSetType(ksp,KSPCHEBYSHEV);
	} 
	else if(solver == "cg"){
		KSPSetType(ksp,KSPCG);
	} 
	else if(solver == "groppcg"){
		KSPSetType(ksp,KSPGROPPCG);
	} 
	else if(solver == "pipecg"){
		KSPSetType(ksp,KSPPIPECG);
	} 
	else if(solver == "pipecgrr"){
		KSPSetType(ksp,KSPPIPECGRR);
	} 
	else if(solver == "pipelcg"){
		KSPSetType(ksp,KSPPIPELCG);
	} 
	else if(solver == "pipeprcg"){
		KSPSetType(ksp,KSPPIPEPRCG);
	} 
	else if(solver == "pipecg2"){
		KSPSetType(ksp,KSPPIPECG2);
	} 
	else if(solver == "cgne"){
		KSPSetType(ksp,KSPCGNE);
	} 
	else if(solver == "nash"){
		KSPSetType(ksp,KSPNASH);
	} 
	else if(solver == "stcg"){
		KSPSetType(ksp,KSPSTCG);
	} 
	else if(solver == "gltr"){
		KSPSetType(ksp,KSPGLTR );
	} 
	else if(solver == "fcg"){
		KSPSetType(ksp,KSPFCG);
	} 
	else if(solver == "pipefcg"){
		KSPSetType(ksp,KSPPIPEFCG);
	} 
	else if(solver == "pipefgmres"){
		KSPSetType(ksp,KSPPIPEFGMRES);
	} 
	else if(solver == "fgmres"){
		KSPSetType(ksp,KSPFGMRES);
	} 
	else if(solver == "lgmres"){
		KSPSetType(ksp,KSPLGMRES);
	} 
	else if(solver == "dgmres"){
		KSPSetType(ksp,KSPDGMRES);
	} 
	else if(solver == "pgmres"){
		KSPSetType(ksp,KSPPGMRES);
	} 
	else if(solver == "tcqmr"){
		KSPSetType(ksp,KSPTCQMR);
	} 
	else if(solver == "bcgs"){
		KSPSetType(ksp,KSPBCGS);
	} 
	else if(solver == "ibcgs"){
		KSPSetType(ksp,KSPIBCGS);
	} 
	else if(solver == "fbcgs"){
		KSPSetType(ksp,KSPFBCGS);
	} 
	else if(solver == "fbcgsr"){
		KSPSetType(ksp,KSPFBCGSR);
	} 
	else if(solver == "bcgsl"){
		KSPSetType(ksp,KSPBCGSL);
	} 
	else if(solver == "pipebcgs"){
		KSPSetType(ksp,KSPPIPEBCGS);
	} 
	else if(solver == "cgs"){
		KSPSetType(ksp,KSPCGS);
	} 
	else if(solver == "tfqmr"){
		KSPSetType(ksp,KSPTFQMR);
	} 
	else if(solver == "cr"){
		KSPSetType(ksp,KSPCR);
	} 
	else if(solver == "pipecr"){
		KSPSetType(ksp,KSPPIPECR);
	} 
	else if(solver == "lsqr"){
		KSPSetType(ksp,KSPLSQR);
	} 
	else if(solver == "preonly"){
		KSPSetType(ksp,KSPPREONLY);
	} 
	else if(solver == "qcg"){
		KSPSetType(ksp,KSPQCG);
	} 
	else if(solver == "bicg"){
		KSPSetType(ksp,KSPBICG);
	} 
	else if(solver == "minres"){
		KSPSetType(ksp,KSPMINRES);
	} 
	else if(solver == "symmlq"){
		KSPSetType(ksp,KSPSYMMLQ);
	} 
	else if(solver == "lcd"){
		KSPSetType(ksp,KSPLCD );
	} 
	else if(solver == "python"){
		KSPSetType(ksp,KSPPYTHON);
	} 
	else if(solver == "gcr"){
		KSPSetType(ksp,KSPGCR);
	} 
	else if(solver == "pipegcr"){
		KSPSetType(ksp,KSPPIPEGCR);
	} 
	else if(solver == "tsirm"){
		KSPSetType(ksp,KSPTSIRM);
	} 
	else if(solver == "cgls"){
		KSPSetType(ksp,KSPCGLS);
	} 
	else if(solver == "fetidp"){
		KSPSetType(ksp,KSPFETIDP);
	} 
	// else if(solver == "hpddm"){
		// KSPSetType(ksp,SPHPDDM);
	// } 
	else{
		cerr << "| #Error : no defined PETSc solver name -> " << solver << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	

	KSPGetPC(ksp,&pc);
	
	
	if(preconditioner == "none"){
		PCSetType(pc,PCNONE);
	} 
	else if(preconditioner == "jacobi"){
		PCSetType(pc,PCJACOBI);
	} 
	else if(preconditioner == "sor"){
		PCSetType(pc,PCSOR);
	} 
	else if(preconditioner == "lu"){
		PCSetType(pc,PCLU );
	} 
	else if(preconditioner == "shell"){
		PCSetType(pc,PCSHELL);
	} 
	else if(preconditioner == "bjacobi"){
		PCSetType(pc,PCBJACOBI);
	} 
	else if(preconditioner == "mg"){
		PCSetType(pc,PCMG);
	} 
	else if(preconditioner == "eisenstat"){
		PCSetType(pc,PCEISENSTAT);
	} 
	else if(preconditioner == "ilu"){
		PCSetType(pc,PCILU);
	} 
	else if(preconditioner == "icc"){
		PCSetType(pc,PCICC );
	} 
	else if(preconditioner == "asm"){
		PCSetType(pc,PCASM);
	} 
	else if(preconditioner == "gasm"){
		PCSetType(pc,PCGASM);
	} 
	else if(preconditioner == "ksp"){
		PCSetType(pc,PCKSP);
	} 
	else if(preconditioner == "composite"){
		PCSetType(pc,PCCOMPOSITE);
	} 
	else if(preconditioner == "redundant"){
		PCSetType(pc,PCREDUNDANT);
	} 
	else if(preconditioner == "spai"){
		PCSetType(pc,PCSPAI);
	} 
	else if(preconditioner == "nn"){
		PCSetType(pc,PCNN);
	} 
	else if(preconditioner == "cholesky"){
		PCSetType(pc,PCCHOLESKY);
	} 
	else if(preconditioner == "pbjacobi"){
		PCSetType(pc,PCPBJACOBI);
	} 
	else if(preconditioner == "vpbjacobi"){
		PCSetType(pc,PCVPBJACOBI);
	} 
	else if(preconditioner == "mat"){
		PCSetType(pc,PCMAT);
	} 
	else if(preconditioner == "hypre"){
		PCSetType(pc,PCHYPRE);
	} 
	else if(preconditioner == "parms"){
		PCSetType(pc,PCPARMS);
	} 
	else if(preconditioner == "fieldsplit"){
		PCSetType(pc,PCFIELDSPLIT);
	} 
	else if(preconditioner == "tfs"){
		PCSetType(pc,PCTFS);
	} 
	else if(preconditioner == "ml"){
		PCSetType(pc,PCML);
	} 
	else if(preconditioner == "galerkin"){
		PCSetType(pc,PCGALERKIN );
	} 
	else if(preconditioner == "exotic"){
		PCSetType(pc,PCEXOTIC);
	} 
	else if(preconditioner == "cp"){
		PCSetType(pc,PCCP);
	} 
	else if(preconditioner == "bfbt"){
		PCSetType(pc, PCBFBT);
	} 
	else if(preconditioner == "lsc"){
		PCSetType(pc,PCLSC );
	} 
	else if(preconditioner == "python"){
		PCSetType(pc,PCPYTHON);
	} 
	else if(preconditioner == "pfmg"){
		PCSetType(pc,PCPFMG );
	} 
	else if(preconditioner == "syspfmg"){
		PCSetType(pc,PCSYSPFMG);
	} 
	else if(preconditioner == "redistribute"){
		PCSetType(pc,PCREDISTRIBUTE);
	} 
	else if(preconditioner == "svd"){
		PCSetType(pc,PCSVD);
	} 
	else if(preconditioner == "gamg"){
		PCSetType(pc,PCGAMG);
	} 
	else if(preconditioner == "chowiluviennacl"){
		PCSetType(pc,PCCHOWILUVIENNACL);
	} 
	else if(preconditioner == "rowscalingviennacl"){
		PCSetType(pc,PCROWSCALINGVIENNACL);
	} 
	else if(preconditioner == "saviennacl"){
		PCSetType(pc,PCSAVIENNACL);
	} 
	else if(preconditioner == "bddc"){
		PCSetType(pc,PCBDDC);
	} 
	else if(preconditioner == "kaczmarz"){
		PCSetType(pc,PCKACZMARZ);
	} 
	else if(preconditioner == "telescope"){
		PCSetType(pc,PCTELESCOPE );
	} 
	else if(preconditioner == "patch"){
		PCSetType(pc,PCPATCH);
	} 
	else if(preconditioner == "lmvm"){
		PCSetType(pc,PCLMVM );
	} 
	else if(preconditioner == "hmg"){
		PCSetType(pc,PCHMG);
	} 
	else if(preconditioner == "deflation"){
		PCSetType(pc,PCDEFLATION);
	} 
	else if(preconditioner == "hpddm"){
		PCSetType(pc,PCHPDDM);
	} 
	else if(preconditioner == "hara"){
		PCSetType(pc,PCHARA);
	} 
	else{
		cerr << "| #Error : no defined PETSc preconditioner name -> " << preconditioner << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	
	
	// // KSPSetType(ksp,KSPCG);
	// // KSPSetType(ksp,KSPGMRES);
	// // KSPSetType(ksp,KSPBCGS);
	// KSPSetType(ksp,KSPFGMRES);
	
	// // KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	
	// KSPGetPC(ksp,&pc);
	
	// PCSetType(pc,PCMG);
	// // PCSetType(pc,PCPBJACOBI);
	
	
	// // PCSetType(pc,PCGAMG);
	// // PCHMGSetInnerPCType(pc,PCGAMG);
	// // PCHMGSetReuseInterpolation(pc,PETSC_TRUE);
	// // PCHMGSetUseSubspaceCoarsening(pc,PETSC_TRUE);
	// // PCHMGUseMatMAIJ(pc,PETSC_FALSE);
	// // PCHMGSetCoarseningComponent(pc,0);
	
	

	if(rank==0 && checkTime) cout<< "2nd : " << (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	if(rank==0 && checkTime) startTime.pop_back();
	

	if(rank==0 && checkTime) startTime.push_back(clock());
	
	// KSPSetTolerances(ksp, relTol, tolerance, PETSC_DEFAULT, maxIter);
	KSPSetTolerances(ksp, relTol, tolerance, 1.e+5, maxIter);
	

	KSPSetFromOptions(ksp);
	KSPSetUp(ksp);

	
	ierr = KSPSolve(ksp,b,x);
	
	if(ierr){
		if(rank==0) cerr << "| #Error : PETSc error" << endl;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	
	if(permute) {
		VecPermute(x,rowperm,PETSC_TRUE);
	}
	

	if(rank==0 && checkTime) cout<< "3rd : " << (double)(clock()-startTime.back())/1000.0 << " sec" << endl;
	if(rank==0 && checkTime) startTime.pop_back();
	
	for(int i=0; i<mesh.cells.size(); ++i){
		cols = myStart + i;
		VecGetValues(x,i1,&cols,&resiVar[i]);
	}


	// KSPGetTotalIterations(ksp,&its);
	// KSPGetResidualNorm(ksp,&norm);



	

	// if(rank==0 && checkTime) startTime.push_back(clock());

	KSPDestroy(&ksp);
	VecDestroy(&x);
	VecDestroy(&b);
	MatDestroy(&A);
	
	PetscFinalize();
	 
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
}