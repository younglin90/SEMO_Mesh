#include "build.h"
#include "mpi.h"

// #define AMGCL_NO_BOOST


// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>
// #include <boost/foreach.hpp>
// #include <boost/range/iterator_range.hpp>
// #include <boost/scope_exit.hpp>


// #include <amgcl/io/binary.hpp>
// #include <amgcl/io/mm.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
// #include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/schur_pressure_correction.hpp>
#include <amgcl/mpi/block_preconditioner.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/mpi/cpr.hpp>


// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/solver/runtime.hpp>
// #include <amgcl/coarsening/runtime.hpp>
// #include <amgcl/relaxation/runtime.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>

#include <amgcl/mpi/coarsening/runtime.hpp>
#include <amgcl/mpi/relaxation/runtime.hpp>
#include <amgcl/mpi/partition/runtime.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>
#include <amgcl/mpi/solver/runtime.hpp>


#include <amgcl/mpi/relaxation/as_preconditioner.hpp>

#include <amgcl/mpi/distributed_matrix.hpp>
#include <amgcl/mpi/amg.hpp>
#include <amgcl/mpi/relaxation/spai0.hpp>
#include <amgcl/mpi/relaxation/spai1.hpp>
#include <amgcl/mpi/coarsening/aggregation.hpp>
#include <amgcl/mpi/coarsening/smoothed_aggregation.hpp>
#include <amgcl/mpi/solver/bicgstab.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/runtime.hpp>


#include <amgcl/mpi/relaxation/ilu0.hpp>
#include <amgcl/mpi/relaxation/gauss_seidel.hpp>
#include <amgcl/mpi/solver/preonly.hpp>
#include <amgcl/mpi/solver/idrs.hpp>
// #include <amgcl/mpi/solver/cg.hpp>
// #include <amgcl/mpi/relaxation/ilut.hpp>
#include <amgcl/mpi/relaxation/iluk.hpp>
// #include <amgcl/mpi/relaxation/ilup.hpp>
#include <amgcl/mpi/solver/fgmres.hpp>
#include <amgcl/mpi/solver/lgmres.hpp>
#include <amgcl/mpi/solver/gmres.hpp>
#include <amgcl/mpi/solver/richardson.hpp>

#include <amgcl/mpi/direct_solver/skyline_lu.hpp>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/scope_exit.hpp>

// #include <amgcl/io/binary.hpp>
// #include <amgcl/io/mm.hpp>
// #include <amgcl/adapter/crs_tuple.hpp>
// #include <amgcl/backend/builtin.hpp>
// #include <amgcl/mpi/make_solver.hpp>
// #include <amgcl/mpi/cpr.hpp>
// #include <amgcl/mpi/amg.hpp>
// #include <amgcl/mpi/coarsening/runtime.hpp>
// #include <amgcl/mpi/relaxation/runtime.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>
// #include <amgcl/mpi/relaxation/as_preconditioner.hpp>
// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/mpi/partition/runtime.hpp>
// #include <amgcl/profiler.hpp>


// // mpi
// #include <amgcl/mpi/distributed_matrix.hpp>
// #include <amgcl/mpi/make_solver.hpp>

// #include <amgcl/mpi/coarsening/runtime.hpp>
// #include <amgcl/mpi/relaxation/runtime.hpp>
// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/mpi/partition/runtime.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>
		
// // preconditioners
// #include <amgcl/amg.hpp>
// #include <amgcl/mpi/amg.hpp>
// #include <amgcl/mpi/relaxation/as_preconditioner.hpp>
// #include <amgcl/mpi/cpr.hpp>
// #include <amgcl/mpi/schur_pressure_correction.hpp>
// #include <amgcl/mpi/block_preconditioner.hpp>

// // relaxation
// #include <amgcl/relaxation/runtime.hpp>
// #include <amgcl/mpi/relaxation/spai0.hpp>
// #include <amgcl/mpi/relaxation/ilu0.hpp>
// #include <amgcl/relaxation/as_preconditioner.hpp>
// // #include <amgcl/relaxation/damped_jacobi.hpp>
// // #include <amgcl/relaxation/gauss_seidel.hpp>
// // #include <amgcl/relaxation/chebyshev.hpp>
// // #include <amgcl/relaxation/ilu0.hpp>
// // #include <amgcl/relaxation/iluk.hpp>
// // #include <amgcl/relaxation/ilup.hpp>
// // #include <amgcl/relaxation/ilut.hpp>
// // #include <amgcl/relaxation/spai0.hpp>
// // #include <amgcl/relaxation/spai1.hpp>

// // coarsening strategies
// #include <amgcl/mpi/coarsening/aggregation.hpp>
// #include <amgcl/mpi/coarsening/smoothed_aggregation.hpp>
// #include <amgcl/coarsening/runtime.hpp>

// // solver
// #include <amgcl/mpi/subdomain_deflation.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>
// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/mpi/solver/idrs.hpp>
// #include <amgcl/mpi/solver/bicgstab.hpp>
// #include <amgcl/mpi/solver/preonly.hpp>
// // #include <amgcl/solver/cg.hpp>
// // #include <amgcl/solver/bicgstab.hpp>
// // #include <amgcl/solver/bicgstabl.hpp>
// // #include <amgcl/solver/gmres.hpp>
// // #include <amgcl/solver/lgmres.hpp>
// // #include <amgcl/solver/fgmres.hpp>
// // #include <amgcl/solver/preonly.hpp>


// // #include <amgcl/mpi/distributed_matrix.hpp>


// // #include <boost/program_options.hpp>
// // #include <boost/property_tree/ptree.hpp>
// // #include <boost/property_tree/json_parser.hpp>
// // #include <boost/range/iterator_range.hpp>
// // #include <boost/scope_exit.hpp>

// // #include <amgcl/io/binary.hpp>
// // #include <amgcl/io/mm.hpp>
// // #include <amgcl/adapter/crs_tuple.hpp>
// // #include <amgcl/amg.hpp>
// // #include <amgcl/coarsening/runtime.hpp>
// // #include <amgcl/relaxation/runtime.hpp>
// // #include <amgcl/relaxation/as_preconditioner.hpp>
// // #include <amgcl/mpi/make_solver.hpp>
// // #include <amgcl/mpi/schur_pressure_correction.hpp>
// // #include <amgcl/mpi/block_preconditioner.hpp>
// // #include <amgcl/mpi/subdomain_deflation.hpp>
// // #include <amgcl/mpi/solver/runtime.hpp>
// // #include <amgcl/mpi/direct_solver/runtime.hpp>
// // #include <amgcl/profiler.hpp>






void SEMO_Solvers_Builder::solveAMGCL(
	string equation,
	SEMO_Mesh_Builder& mesh,
	int B_n, 
	vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	vector<double>& B_vals,
	vector<double>& resiVar){
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	

	// ofstream outputFile;
	// string filenamePlot = "Amat." + to_string(rank) + ".txt";
	
	// outputFile.open(filenamePlot);
	// outputFile << "%%MatrixMarket matrix coordinate real general" << endl;
	// outputFile << B_vals.size() << " " << B_vals.size() << " " << A_rows.size() << endl;
	// for(int i=0; i<A_rows.size(); ++i){
		// outputFile << A_rows[i]+1 << " " << A_cols[i]+1 << " " << A_vals[i] << " " << endl;
	// }
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	
	
	// COO to CSR
	int strRow = mesh.startCellGlobal * B_n;
	
    //compute number of non-zero entries per row of A 
	int n_row = B_vals.size();
	int nnz = A_rows.size();
	
	vector<double> val(nnz,0.0);
	vector<int> col(nnz,0);
	vector<int> ptr(n_row+1,0);

    for (int n = 0; n < nnz; n++){            
        ptr[A_rows[n]-strRow]++;
    }

	// cout << "AAAAAAAAAAAAAAA" << endl;
    //cumsum the nnz per row to get ptr[]
    for(int i = 0, cumsum = 0; i < n_row; i++){     
        int temp = ptr[i];
        ptr[i] = cumsum;
        cumsum += temp;
    }
    ptr[n_row] = nnz; 

	// cout << "BBBBBBBB" << endl;
    //write col,val into col,val
    for(int n = 0; n < nnz; n++){
        int row  = A_rows[n] - strRow;
        int dest = ptr[row];
	// cout << row << " " << dest << endl;

        col[dest] = A_cols[n];
        val[dest] = A_vals[n];
		
		if(A_cols[n]>= mesh.ncellsTotal*B_n){
			cout << "ERROR " << A_cols[n] << " " << mesh.ncellsTotal*B_n <<  endl;
		}

        ptr[row]++;
    }
	// MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0, last = 0; i <= n_row; i++){
        int temp = ptr[i];
        ptr[i]  = last;
        last   = temp;
    }

	// cout << "1" << endl;
    //now ptr,col,val form a CSR representation (with possible duplicates)
	
	amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // The profiler:
    amgcl::profiler<> prof("Flows Solve");
	
	
	ptrdiff_t chunk = B_vals.size();

    // Compose the solver type
    typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend

	// cout << ptr.size() << " " << col.size() << " " << val.size() << endl;

	auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     world, std::tie(chunk, ptr, col, val));


    // typedef amgcl::mpi::make_solver<
        // amgcl::mpi::amg<
            // PBackend,
            // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
            // amgcl::mpi::relaxation::ilu0<PBackend>
            // >,
        // amgcl::mpi::solver::idrs<SBackend>
        // > Solver;

			
    // typedef
        // amgcl::mpi::make_solver<
            // amgcl::mpi::schur_pressure_correction<
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::block_preconditioner<
                        // amgcl::relaxation::as_preconditioner<SBackend, amgcl::runtime::relaxation::wrapper>
                        // >,
                    // amgcl::runtime::mpi::solver::wrapper<SBackend>
                    // >,
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::amg<
                        // SBackend,
                        // amgcl::runtime::mpi::coarsening::wrapper<SBackend>,
                        // amgcl::runtime::mpi::relaxation::wrapper<SBackend>,
                        // amgcl::runtime::mpi::direct::solver<double>,
                        // amgcl::runtime::mpi::partition::wrapper<SBackend>
                        // >,
                        // amgcl::runtime::mpi::solver::wrapper<SBackend>
                    // >
                // >,
                // amgcl::runtime::mpi::solver::wrapper<SBackend>
            // > Solver;




				
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// UBackend,
						// amgcl::mpi::coarsening::smoothed_aggregation<UBackend>,
						// amgcl::mpi::relaxation::ilu0<UBackend>
					// >,
					// amgcl::mpi::solver::bicgstab<UBackend>
				// > U1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// PBackend,
						// amgcl::mpi::coarsening::aggregation<PBackend>,
						// amgcl::mpi::relaxation::ilu0<PBackend>
					// >,
					// amgcl::mpi::solver::preonly<PBackend>
				// > U1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<PBackend, amgcl::relaxation::ilu0>
					// >,
					// amgcl::mpi::solver::preonly<PBackend>
				// > U1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// UBackend,
						// amgcl::mpi::coarsening::aggregation<UBackend>,
						// amgcl::mpi::relaxation::ilut<UBackend>
					// >,
					// amgcl::mpi::solver::preonly<UBackend>
				// > U1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// SBackend,
						// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
						// amgcl::mpi::relaxation::ilu0<SBackend>
					// >,
					// amgcl::mpi::solver::preonly<SBackend>
				// > U1;
					
	// typedef amgcl::mpi::make_solver<
				// amgcl::mpi::relaxation::as_preconditioner<
					// PBackend,
					// amgcl::mpi::relaxation::spai0<PBackend>
					// >,
				// amgcl::mpi::solver::preonly<PBackend>
				// > P1;
				
	// typedef amgcl::mpi::make_solver<
				// amgcl::mpi::amg<
					// PBackend,
					// amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
					// amgcl::mpi::relaxation::spai0<PBackend>
					// >,
					// amgcl::mpi::solver::preonly<PBackend>
				// > P1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<
							// SBackend, amgcl::relaxation::spai0
						// >
					// >,
					// amgcl::mpi::solver::preonly<SBackend>
				// > P1;
	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// UBackend,
						// amgcl::mpi::coarsening::aggregation<UBackend>,
						// amgcl::mpi::relaxation::spai0<UBackend>
					// >,
					// amgcl::mpi::solver::preonly<UBackend>
				// > P1;
				
	// typedef amgcl::mpi::schur_pressure_correction<P1,P1> SCHUR;
	// typedef amgcl::mpi::schur_pressure_correction<U1,U1> SCHUR;
	
	// typedef amgcl::mpi::solver::idrs<SBackend> S1;
	// typedef amgcl::mpi::solver::fgmres<SBackend> S1;
	// typedef amgcl::mpi::solver::gmres<SBackend> S1;
	
    // typedef amgcl::mpi::make_solver<SCHUR,S1> Solver;
	

	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// SBackend,
						// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
						// amgcl::mpi::relaxation::iluk<SBackend>
					// >,
					// amgcl::mpi::solver::bicgstab<SBackend>
				// > Solver;
	
	
	
    // typedef
        // amgcl::mpi::make_solver<
            // amgcl::mpi::schur_pressure_correction<
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::block_preconditioner<
                        // amgcl::relaxation::as_preconditioner<PBackend, amgcl::relaxation::ilu0<PBackend>>
                        // >,
                    // amgcl::mpi::solver::preonly<PBackend>
                    // >,
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::amg<
                        // PBackend,
                        // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
                        // amgcl::mpi::relaxation::ilu0<PBackend>
                        // >,
                        // amgcl::mpi::solver::bicgstab<PBackend>
                    // >
                // >,
                // amgcl::mpi::solver::bicgstab<SBackend>
            // > Solver;

			
					
		// typedef amgcl::mpi::make_solver<
						// amgcl::mpi::amg<
							// PBackend,
							// amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
							// amgcl::mpi::relaxation::ilu0<PBackend>
						// >,
						// amgcl::mpi::solver::preonly<PBackend>
					// > U1;
				
		// typedef amgcl::mpi::make_solver<
						// amgcl::mpi::amg<
							// PBackend,
							// amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
							// amgcl::mpi::relaxation::spai0<PBackend>
						// >,
						// amgcl::mpi::solver::preonly<PBackend>
					// > P1;

		// typedef amgcl::mpi::schur_pressure_correction<U1,P1> SCHUR;
		
		// typedef amgcl::mpi::solver::idrs<SBackend> S1;
		
		// typedef amgcl::mpi::make_solver<SCHUR,S1> Solver;
		
		
		// Solver::params prm;

		// prm.precond.usolver.precond.npre = 0;
		// prm.precond.usolver.precond.npost = 3;
		// prm.precond.usolver.precond.ncycle = 1;
		// // prm.precond.usolver.precond.relax.k = 2;
		
		// prm.precond.psolver.precond.npre = 0;
		// prm.precond.psolver.precond.npost = 3;
		// prm.precond.psolver.precond.ncycle = 1;
		// // prm.precond.psolver.precond.relax.k = 2;

		// // prm.precond.simplec_dia = true;
		// // prm.precond.adjust_p = 1;
		// // prm.precond.adjust_p = 2;
		
		// // prm.precond.coarsening.estimate_spectral_radius = true;

		// // prm.solver.ns_search = false;
		// // prm.solver.verbose = false;
		// prm.solver.maxiter = 5;
		// prm.solver.tol = 1e-6;
		// // prm.solver.replacement = true;
		// // prm.solver.smoothing = true;
		// // prm.solver.s = 5;
		// prm.precond.pmask.resize(n_row,0);
		// for(int i=0; i<n_row/B_n; ++i) prm.precond.pmask[i] = 1;
		// // prm.precond.pmask_size = n_row/B_n;

		

	if(equation == "coupled"){
		
		
		// runtime
		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::schur_pressure_correction<
				// amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<SBackend, amgcl::relaxation::iluk>
						// >,
					// amgcl::mpi::solver::preonly<SBackend>
					// >,
				// amgcl::mpi::subdomain_deflation<
					// amgcl::amg<SBackend, amgcl::coarsening::aggregation, amgcl::relaxation::iluk>,
					// amgcl::mpi::solver::preonly<SBackend>,
					// amgcl::mpi::direct::skyline_lu<double>
					// >
				// >,
			// amgcl::mpi::solver::gmres<SBackend>
			// >
		// Solver;
		
		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::schur_pressure_correction<
				// amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<SBackend, amgcl::relaxation::iluk>
						// // amgcl::amg<SBackend, amgcl::coarsening::aggregation, amgcl::relaxation::ilu0>
						// >,
					// amgcl::mpi::solver::preonly<SBackend>
					// >,
				// amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<SBackend, amgcl::relaxation::iluk>
						// >,
					// amgcl::mpi::solver::preonly<SBackend>
					// >
				// >,
			// amgcl::mpi::solver::gmres<SBackend>
			// >
		// Solver;
		

		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::schur_pressure_correction<
				// amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// // amgcl::relaxation::as_preconditioner<SBackend, amgcl::relaxation::iluk>
						// amgcl::amg<SBackend, amgcl::coarsening::aggregation, amgcl::relaxation::ilu0>
						// >,
					// amgcl::mpi::solver::idrs<SBackend>
					// >,
				// amgcl::mpi::make_solver<
					// amgcl::mpi::block_preconditioner<
						// amgcl::relaxation::as_preconditioner<SBackend, amgcl::relaxation::iluk>
						// >,
					// amgcl::mpi::solver::gmres<SBackend>
					// >
				// >,
			// amgcl::mpi::solver::gmres<SBackend>
			// >
		// Solver;
		
		
		
		// Solver::params prm;
		
		// prm.precond.usolver.precond.coarsening.over_interp = 1.0;
		// prm.precond.usolver.solver.maxiter = 10;
		// prm.precond.usolver.solver.tol = 1.e-6;
		// prm.precond.usolver.solver.replacement = true;
		// prm.precond.usolver.solver.smoothing = true;
		
		// prm.precond.psolver.solver.maxiter = 10;
		// prm.precond.psolver.solver.tol = 1.e-9;
		
		// // prm.precond.block_size = 6;
		// // prm.precond.active_rows = n_row/B_n;
		
		// prm.solver.maxiter = 1000;
		// prm.solver.tol = 1.e-9;
		
		// prm.precond.pmask.clear();
		// prm.precond.pmask.resize(n_row,1);
		// // std::fill(prm.precond.pmask.begin(), prm.precond.pmask.begin() + n_row/B_n, 0);
		// for(int i=0; i<n_row; ++i) prm.precond.pmask[i] = (i < n_row/B_n ? 0 : 1);
		

		// Solver solve(world, A, prm);
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		

		// if (world.rank == 0) {
			// // cout << solve << endl;
			// // cout << prof << endl;
			// cout << "|--- coupled. : " << iters << " " << error << endl;
		// }
		
		
		


		// iterative solver
		typedef amgcl::mpi::make_solver<
			amgcl::mpi::relaxation::as_preconditioner<
			amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			>,
			// amgcl::mpi::amg<
			// SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// amgcl::mpi::solver::bicgstab<SBackend>
			amgcl::mpi::solver::fgmres<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// amgcl::mpi::solver::preonly<SBackend>
			// amgcl::mpi::solver::lgmres<SBackend>
			// amgcl::mpi::solver::fgmres<SBackend>
		> Solver; 
		
		
		Solver::params prm;
		
		// prm.precond.direct_coarse = true;
		// prm.precond.npre = 1;
		// prm.precond.npost = 5;
		// prm.precond.pre_cycles = 5;
		
		// prm.precond.coarsening.over_interp = 1.0;
		
		prm.precond.k = 0;
		prm.precond.damping = 3.0;
		// prm.precond.relax.k = 2;
		// prm.solver.K = 5;
		// prm.solver.M = 50;
		// prm.solver.s = 8;
		// prm.solver.omega = 0.7;
		prm.solver.maxiter = 100;//1000;
		prm.solver.tol = 1.e-9;
		// prm.solver.replacement = true;
		// prm.solver.smoothing = true;
		
		
		prof.tic("setup");
		Solver solve(world, A, prm);
		prof.toc("setup");
		prof.tic("solve");
		int iters; double error;
		std::tie(iters, error) = solve(*A, B_vals, resiVar);
		prof.toc("solve");
		
		if (world.rank == 0) {
			cout << solve << endl;
			cout << prof << endl;
			cout << "|--- coupled. : " << iters << " " << error << endl;
		}
		
		
		

		// typedef amgcl::mpi::relaxation::as_preconditioner<
					// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
				// > PPrecond;

		// typedef amgcl::mpi::relaxation::as_preconditioner<
					// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
				// > SPrecond;

		// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::cpr<PPrecond, SPrecond>,
					// amgcl::mpi::solver::bicgstab<SBackend>
				// > Solver;

		// Solver::params prm;
		
		
		// prm.precond.block_size = B_n;
		// // prm.precond.active_rows = n_row/B_n;
		
		
		// // prm.precond.k = 2;
		// // prm.precond.relax.k = 2;
		// // prm.solver.K = 3;
		// // prm.solver.M = 30;
		// // prm.solver.maxiter = 500;
		
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		
		// if (world.rank == 0) {
			// cout << solve << endl;
			// cout << prof << endl;
			// cout << "|--- coupled. : " << iters << " " << error << endl;
		// }
		
		
		
		
		
		

		// // direct solver
		// typedef amgcl::mpi::direct::skyline_lu<double> Solver;
		// Solver::params prm;
		// Solver solve(world, *A);
		// solve(B_vals, resiVar);


		// // iterative solver
		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::amg<
			// SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// amgcl::mpi::relaxation::spai0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// // amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
			// // amgcl::mpi::solver::lgmres<SBackend>
			// amgcl::mpi::solver::fgmres<SBackend>
		// > Solver; 
		// Solver::params prm;
		
		
		// prm.precond.npre = 15;
		// prm.precond.npost = 15;
		// // prm.precond.ncycle = 2;
		// // prm.precond.pre_cycles = 2;
		// // prm.precond.coarsening.over_interp = 1.0;
		// // prm.precond.coarsening.relax = 0.7;

		// prm.solver.maxiter = 500;
		// // prm.solver.K = 5;
		// prm.solver.M = 50;
		// prm.solver.ns_search = true;
		// prm.solver.tol = 1e-9;
		// // prm.solver.replacement = true;
		// // prm.solver.smoothing = true;
	

		// // prm.solver.verbose = false;
		// // if(world.rank == 0) prm.solver.verbose = true;
		// // prm.solver.maxiter = 50;
		// // prm.solver.replacement = true;
		// // prm.solver.smoothing = true;
		
		
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		
		// if (world.rank == 0) {
			// cout << solve << endl;
			// cout << prof << endl;
			// cout << "|--- coupled. : " << iters << " " << error << endl;
		// }
		
		
		
		
		
		
		
		

		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::relaxation::as_preconditioner<
			// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::amg<
			// // SBackend,
			// // // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// // amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// amgcl::mpi::solver::preonly<SBackend>
		// > U1; 
		

		// typedef amgcl::mpi::make_solver<
			// amgcl::mpi::relaxation::as_preconditioner<
			// amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::amg<
			// // SBackend,
			// // // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// // amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// amgcl::mpi::solver::preonly<SBackend>
		// > P1; 
		
		// // typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::block_preconditioner<
				// // amgcl::mpi::relaxation::as_preconditioner<
				// // amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
				// // >
			// // >,
			// // amgcl::mpi::solver::preonly<SBackend>
		// // > U1; 
		
		// // typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::amg<
			// // SBackend,
			// // // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// // amgcl::mpi::solver::idrs<SBackend>
		// // > P1; 
		
		// // typedef amgcl::mpi::subdomain_deflation<
			// // amgcl::mpi::amg<
			// // SBackend,
			// // // amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// // amgcl::mpi::coarsening::aggregation<SBackend>,
			// // amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// // >,
			// // amgcl::mpi::solver::preonly<SBackend>
		// // > P1; 
		
		// // typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::block_preconditioner<
				// // amgcl::mpi::relaxation::as_preconditioner<
				// // amgcl::mpi::relaxation::iluk<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
				// // >
			// // >,
			// // amgcl::mpi::solver::preonly<SBackend>
		// // > P1; 
		
		// typedef amgcl::mpi::schur_pressure_correction<U1,P1> SCHUR;
		
		// typedef amgcl::mpi::solver::idrs<SBackend> S1;
		
		// typedef amgcl::mpi::make_solver<SCHUR,S1> Solver;
		
		
	
		// Solver::params prm;

		// // prm.precond.psolver.precond.npre = 5;
		// // prm.precond.psolver.precond.npost = 5;
		// // prm.precond.usolver.precond.npre = 5;
		// // prm.precond.usolver.precond.npost = 5;
		// // // prm.solver.maxiter = 100;
		// // // prm.solver.tol = 1e-8;
		// // // prm.solver.M = 50;
		// // prm.precond.usolver.solver.maxiter = 15;
		// // prm.precond.psolver.solver.maxiter = 15;
		// // // prm.precond.usolver.solver.tol = 1e-6;
		// // // prm.precond.psolver.solver.tol = 1e-6;
		// // prm.precond.psolver.solver.smoothing = true;
		
		// // prm.solver.smoothing = true;

		// // prm.precond.type = 2;
		// // prm.precond.adjust_p = 0;
		// // prm.precond.approx_schur = false;
		// // prm.precond.simplec_dia = false;
		// prm.precond.pmask.resize(n_row,0);
		// for(int i=0; i<n_row/B_n; ++i) prm.precond.pmask[i] = 1;
		

		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		
		// if (world.rank == 0) {
			// // cout << solve << endl;
			// // cout << prof << endl;
			// cout << "|--- coupled. : " << iters << " " << error << endl;
		// }
		
		
		
		
		
	
	}
	else if(equation == "volume_fraction"){
		
		typedef amgcl::mpi::make_solver<
			// amgcl::mpi::relaxation::as_preconditioner<
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			amgcl::mpi::amg<
			SBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			amgcl::mpi::coarsening::aggregation<SBackend>,
			amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			>,
			// amgcl::mpi::solver::bicgstab<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// amgcl::mpi::solver::preonly<SBackend>
			amgcl::mpi::solver::fgmres<SBackend>
			// amgcl::mpi::solver::lgmres<SBackend>
		> Solver; 
		Solver::params prm;

		// prm.solver.tol = 1e-10;
		// prm.solver.tol = 1e-12;
		prm.precond.npre = 1;
		prm.precond.npost = 5;
		// prm.precond.ncycle = 1;
		// prm.precond.relax.k = 2;
		// // prm.precond.max_levels = 40;
		prm.precond.coarsening.over_interp = 1.8;
		// prm.solver.K = 4;
		// prm.solver.M = 35;
		prm.solver.tol = 1e-9;
		// prm.solver.relTol = 1.e-200;
		prm.solver.maxiter = 100;//1000;
		
		// prm.solver.ns_search = true;
		// prm.solver.verbose = false;
		// if(world.rank == 0) prm.solver.verbose = true;
		// prm.solver.maxiter = 50;
		// prm.solver.replacement = true;
		// prm.solver.smoothing = true;
		
		
		prof.tic("setup");
		Solver solve(world, A, prm);
		prof.toc("setup");
		prof.tic("solve");
		int iters; double error;
		std::tie(iters, error) = solve(*A, B_vals, resiVar);
		prof.toc("solve");
		
		
		// if (world.rank == 0) {
			// cout << solve << endl;
			// cout << prof << endl;
			// cout << "|--- volume frac. : " << iters << " " << error << endl;
		// }
		

	}
	else if(equation == "pressure"){
		
		typedef amgcl::mpi::make_solver<
			// amgcl::mpi::relaxation::as_preconditioner<
			// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			amgcl::mpi::amg<
			SBackend,
			amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
			// amgcl::mpi::coarsening::aggregation<SBackend>,
			amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			>,
			// amgcl::mpi::solver::bicgstab<SBackend>
			amgcl::mpi::solver::idrs<SBackend>
			// amgcl::mpi::solver::preonly<SBackend>
			// amgcl::mpi::solver::lgmres<SBackend>
			// amgcl::mpi::solver::fgmres<SBackend>
		> Solver; 
		Solver::params prm;

		// prm.precond.relax.k = 2;
		prm.precond.npre = 1;
		prm.precond.npost = 5;
		// prm.precond.coarsening.over_interp = 1.2;
		// prm.precond.coarsening.relax = 0.7;

		// prm.solver.ns_search = true;
		// prm.solver.verbose = false;
		// if(world.rank == 0) prm.solver.verbose = true;
		prm.solver.tol = 1e-9;
		// prm.solver.replacement = true;
		prm.solver.maxiter = 1000;
		// prm.solver.smoothing = true;
		

		
		prof.tic("setup");
		Solver solve(world, A, prm);
		prof.toc("setup");
		prof.tic("solve");
		int iters; double error;
		std::tie(iters, error) = solve(*A, B_vals, resiVar);
		prof.toc("solve");
		

		// if (world.rank == 0) {
			// cout << solve << endl;
			// cout << prof << endl;
			// cout << "|--- pressure : " << iters << " " << error << endl;
		// }

	}
	

	
	
	// else if(equation == "momentum"){
		
		// typedef amgcl::mpi::make_solver<
			// // amgcl::mpi::block_preconditioner<
				// // amgcl::relaxation::as_preconditioner<
					// // PBackend, 
					// // amgcl::relaxation::spai0
				// // >
			// // >,
			// amgcl::mpi::amg<
			// PBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
			// amgcl::mpi::relaxation::ilu0<PBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// // amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
			// amgcl::mpi::solver::lgmres<SBackend>
			// // amgcl::mpi::solver::fgmres<SBackend>
		// > Solver; 
		// Solver::params prm;

		// // prm.solver.maxiter = 50;
		// // prm.solver.tol = 1e-6;
		
		
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");
		

	// }



		// typedef amgcl::mpi::make_solver<
			// /*
			// amgcl::relaxation::as_preconditioner<
			// PBackend,
			// amgcl::relaxation::ilu0
			// >,
			// */
			// amgcl::mpi::amg<
			// PBackend,
			// amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
			// amgcl::mpi::relaxation::iluk<PBackend>   // work : gauss_seidel, iluk, ilup, spai1 /// not work : ilu0 
			// >,
			// // amgcl::mpi::solver::bicgstab<SBackend>
			// amgcl::mpi::solver::idrs<SBackend>
			// // amgcl::mpi::solver::preonly<SBackend>
		// > Solver; 
			

		// Solver::params prm;

		// prm.precond.npre = 1;
		// prm.precond.npost = 3;
		// prm.precond.ncycle = 1;
		// prm.precond.relax.k = 2;
		// prm.precond.relax.damping = 0.7;
		// prm.precond.coarsening.estimate_spectral_radius = true;

		// // prm.solver.ns_search = false;
		// // prm.solver.verbose = false;
		// prm.solver.replacement = true;
		// prm.solver.smoothing = true;
		// prm.solver.s = 5;
		// prm.solver.maxiter = 15;
		// prm.solver.tol = 1e-6;
		// // prm.precond.simplec_dia = true;
		// // // prm.precond.adjust_p = 1;
		// prof.tic("setup");
		// Solver solve(world, A, prm);
		// prof.toc("setup");
		// prof.tic("solve");
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, B_vals, resiVar);
		// prof.toc("solve");

	
	
	

    // typedef amgcl::mpi::solver::skyline_lu<double> Solver;
	
	// prof.tic("setup");
	// Solver solve(world, A);
	// prof.toc("setup");
	// int iters;
	// double error;
	// std::vector<double> x(rows, 0.0);
	// prof.tic("solve");
	// solve(B_vals, resiVar);
	// prof.toc("solve");

	// // std::vector<double> r(rows);
	// // amgcl::backend::residual(B_vals, A, x, r);

	// // double rrrr = sqrt(amgcl::backend::inner_product(r, r));
	// // cout << rrrr << endl;

        // // Output the number of iterations, the relative error,
        // // and the profiling data:
        // //std::cout << "Iters: " << iters << std::endl
        // //    << "Error: " << error << std::endl
        // //    << prof << std::endl;
	
	
	
	
}




//+++++++++++++++++++++++++++
// Momentum
//+++++++++++++++++++++++++++
void SEMO_Solvers_Builder::solveAMGCL(
	string equation,
	SEMO_Mesh_Builder& mesh,
	int B_n, 
	vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	vector<double>& B0_vals, vector<double>& B1_vals, vector<double>& B2_vals,
	vector<double>& resiVar0, vector<double>& resiVar1, vector<double>& resiVar2){
		
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	
	// COO to CSR
	int strRow = mesh.startCellGlobal * B_n;
	
    //compute number of non-zero entries per row of A 
	int n_row = B0_vals.size();
	int nnz = A_rows.size();
	
	vector<double> val(nnz,0.0);
	vector<int> col(nnz,0);
	vector<int> ptr(n_row+1,0);

    for (int n = 0; n < nnz; n++){            
        ptr[A_rows[n]-strRow]++;
    }

	// cout << "AAAAAAAAAAAAAAA" << endl;
    //cumsum the nnz per row to get ptr[]
    for(int i = 0, cumsum = 0; i < n_row; i++){     
        int temp = ptr[i];
        ptr[i] = cumsum;
        cumsum += temp;
    }
    ptr[n_row] = nnz; 

	// cout << "BBBBBBBB" << endl;
    //write col,val into col,val
    for(int n = 0; n < nnz; n++){
        int row  = A_rows[n] - strRow;
        int dest = ptr[row];
	// cout << row << " " << dest << endl;

        col[dest] = A_cols[n];
        val[dest] = A_vals[n];
		
		if(A_cols[n]>= mesh.ncellsTotal*B_n){
			cout << "ERROR " << A_cols[n] << " " << mesh.ncellsTotal*B_n <<  endl;
		}

        ptr[row]++;
    }
	// MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0, last = 0; i <= n_row; i++){
        int temp = ptr[i];
        ptr[i]  = last;
        last   = temp;
    }

	// cout << "1" << endl;
    //now ptr,col,val form a CSR representation (with possible duplicates)
	
	amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // The profiler:
    amgcl::profiler<> prof("Flows Solve");
	
	
	ptrdiff_t chunk = B0_vals.size();

    // Compose the solver type
    typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend

	// cout << ptr.size() << " " << col.size() << " " << val.size() << endl;

	auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     world, std::tie(chunk, ptr, col, val));

	typedef amgcl::mpi::make_solver<
		// amgcl::mpi::relaxation::as_preconditioner<
		// amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
		// >,
		amgcl::mpi::amg<
		SBackend,
		// amgcl::mpi::coarsening::smoothed_aggregation<SBackend>,
		amgcl::mpi::coarsening::aggregation<SBackend>,
		amgcl::mpi::relaxation::ilu0<SBackend>   //gauss_seidel, ilu0, iluk, ilup, spai1
		>,
		// amgcl::mpi::solver::bicgstab<SBackend>
		// amgcl::mpi::solver::idrs<SBackend>
		// amgcl::mpi::solver::preonly<SBackend>
		// amgcl::mpi::solver::lgmres<SBackend>
		amgcl::mpi::solver::fgmres<SBackend>
	> Solver; 
	Solver::params prm;
	
	
	prm.precond.npre = 1;
	prm.precond.npost = 5;
	// prm.precond.ncycle = 2;
	// prm.precond.pre_cycles = 2;
	prm.precond.coarsening.over_interp = 1.0;
	// prm.precond.coarsening.relax = 0.7;

	// prm.solver.K = 4;
	// prm.solver.M = 35;
	// prm.solver.ns_search = true;
	prm.solver.tol = 1e-9;
	prm.solver.maxiter = 100;//1000;
	// prm.solver.replacement = true;
	// prm.solver.smoothing = true;
	
	
	// bool reorder = true;
	
    // if (reorder) {
        // amgcl::adapter::reorder<> perm(A);

        // Solver solve(world, perm(A), prm);
		// int iters; double error;
		// std::tie(iters, error) = solve(*A, perm(B0_vals), resiVar0);
		// std::tie(iters, error) = solve(*A, perm(B1_vals), resiVar1);
		// std::tie(iters, error) = solve(*A, perm(B2_vals), resiVar2);

        // perm.inverse(resiVar0, resiVar0);
        // perm.inverse(resiVar1, resiVar1);
        // perm.inverse(resiVar2, resiVar2);
	// }
	// else{
		prof.tic("setup");
		Solver solve(world, A, prm);
		prof.toc("setup");
		prof.tic("solve");
		int iters0; double error0;
		int iters1; double error1;
		int iters2; double error2;
		std::tie(iters0, error0) = solve(*A, B0_vals, resiVar0);
		std::tie(iters1, error1) = solve(*A, B1_vals, resiVar1);
		std::tie(iters2, error2) = solve(*A, B2_vals, resiVar2);
		prof.toc("solve");
	// }
	
	

	// if (world.rank == 0) {
		// cout << solve << endl;
		// cout << prof << endl;
		// cout << "|--- momentum : " << iters0+iters1+iters2 << " " << error0+error1+error2 << endl;
	// }
	
	
}



















void SEMO_Solvers_Builder::solveAMGCL_VF(
	SEMO_Mesh_Builder& mesh,
	int B_n, 
	vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	vector<double>& B_vals,
	vector<double>& resiVar){
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 
	
	// // COO to CSR
	// int strRow = mesh.startCellGlobal * B_n;
	
    // //compute number of non-zero entries per row of A 
	// int n_row = B_vals.size();
	// int nnz = A_rows.size();
	
	// vector<double> val(nnz,0.0);
	// vector<int> col(nnz,0);
	// vector<int> ptr(n_row,0);

    // for (int n = 0; n < nnz; n++){            
        // ptr[A_rows[n]-strRow]++;
    // }

    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = ptr[i];
        // ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // ptr[n_row] = nnz; 

    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = A_rows[n] - strRow;
        // int dest = ptr[row];

        // col[dest] = A_cols[n];
        // val[dest] = A_vals[n];

        // ptr[row]++;
    // }
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // cout << "AAAAAAAAAAAAAAA" << endl;

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = ptr[i];
        // ptr[i]  = last;
        // last   = temp;
    // }

    // //now ptr,col,val form a CSR representation (with possible duplicates)
	
	// amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // // The profiler:
    // amgcl::profiler<> prof("Volume Fraction Solve");
	
	
	// ptrdiff_t chunk = B_vals.size();
		

    // // Compose the solver type
    // typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    // typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    // //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    // typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend



    // // typedef amgcl::mpi::make_solver<
        // // amgcl::mpi::amg<
            // // PBackend,
            // // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
            // // amgcl::mpi::relaxation::spai0<PBackend>
            // // >,
        // // amgcl::mpi::solver::bicgstab<SBackend>
        // // > Solver;

	// typedef amgcl::mpi::make_solver<
					// amgcl::mpi::amg<
						// UBackend,
						// amgcl::mpi::coarsening::aggregation<UBackend>,
						// amgcl::mpi::relaxation::ilut<UBackend>
					// >,
					// amgcl::mpi::solver::preonly<SBackend>
				// > Solver;
			
			
			
    // Solver::params prm;
    // // prm.solver.tol = 1.e-12;
	// // prm.solver.maxiter = 400;
	
    // // prm.solver.replacement = true;
    // // prm.solver.smoothing = true;
    // // prm.solver.s = 5;

	
    // prof.tic("setup");
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, ptr, col, val));
    // Solver solve(world, A, prm);
    // prof.toc("setup");
    // prof.tic("solve");
    // int iters; double error;
    // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // prof.toc("solve");

    // if (world.rank == 0) {
        // std::cout
            // << "Iterations: " << iters << std::endl
            // << "Error:      " << error << std::endl
            // << std::endl
            // << prof << std::endl;
	// }
	


	
}












// void SEMO_Solvers_Builder::solveAMGCL(
	// SEMO_Mesh_Builder& mesh,
	// int B_n, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B_vals,
	// vector<double>& resiVar){
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 
	
	// // COO to CSR
	// int strRow = mesh.startCellGlobal * B_n;
	
    // //compute number of non-zero entries per row of A 
	// int n_row = B_vals.size();
	// int nnz = A_rows.size();
	
	// vector<double> val(nnz,0.0);
	// vector<int> col(nnz,0);
	// vector<int> ptr(n_row,0);

    // for (int n = 0; n < nnz; n++){            
        // ptr[A_rows[n]-strRow]++;
    // }

    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = ptr[i];
        // ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // ptr[n_row] = nnz; 

    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = A_rows[n] - strRow;
        // int dest = ptr[row];

        // col[dest] = A_cols[n];
        // val[dest] = A_vals[n];

        // ptr[row]++;
    // }
	// // MPI_Barrier(MPI_COMM_WORLD);
	// // cout << "AAAAAAAAAAAAAAA" << endl;

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = ptr[i];
        // ptr[i]  = last;
        // last   = temp;
    // }

    // //now ptr,col,val form a CSR representation (with possible duplicates)
	
	// amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // // The profiler:
    // amgcl::profiler<> prof("AMGCL MPI");
	
	
	// ptrdiff_t chunk = B_vals.size();
		
	
    // // Compose the solver type
    // typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    // typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend


    // // typedef amgcl::mpi::make_solver<
        // // amgcl::mpi::amg<
            // // PBackend,
            // // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
            // // amgcl::mpi::relaxation::ilu0<PBackend>
            // // >,
        // // amgcl::mpi::solver::idrs<SBackend>
        // // > Solver;


	
    // // SBackend::params bprm;
	
    // MPI_Barrier(world);
	

    // prof.tic("setup");
    // // typedef
        // // amgcl::mpi::make_solver<
            // // amgcl::mpi::schur_pressure_correction<
                // // amgcl::mpi::make_solver<
                    // // amgcl::mpi::block_preconditioner<
                        // // amgcl::relaxation::as_preconditioner<PBackend, amgcl::runtime::relaxation::wrapper>
                        // // >,
                    // // amgcl::runtime::mpi::solver::wrapper<PBackend>
                    // // >,
                // // amgcl::mpi::subdomain_deflation<
                    // // amgcl::amg<PBackend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>,
                    // // amgcl::runtime::mpi::solver::wrapper<PBackend>,
                    // // amgcl::runtime::mpi::direct::solver<double>
                    // // >
                // // >,
            // // amgcl::runtime::mpi::solver::wrapper<SBackend>
            // // > Solver;
			

    // // typedef
        // // amgcl::mpi::make_solver<
            // // amgcl::mpi::schur_pressure_correction<
                // // amgcl::mpi::make_solver<
                    // // amgcl::mpi::block_preconditioner<
                        // // amgcl::relaxation::as_preconditioner<PBackend, amgcl::runtime::relaxation::wrapper>
                        // // >,
                    // // amgcl::runtime::mpi::solver::wrapper<PBackend>
                    // // >,
                // // amgcl::mpi::make_solver<
                    // // amgcl::mpi::block_preconditioner<
                        // // amgcl::relaxation::as_preconditioner<PBackend, amgcl::runtime::relaxation::wrapper>
                        // // >,
                    // // amgcl::runtime::mpi::solver::wrapper<PBackend>
                    // // >
                // // >,
            // // amgcl::runtime::mpi::solver::wrapper<SBackend>
            // // > Solver;

	// // typedef amgcl::mpi::make_solver<
				// // amgcl::runtime::mpi::relaxation::wrapper
				// // <amgcl::runtime::mpi::relaxation::wrapper<PBackend>>,
				// // amgcl::runtime::mpi::solver::wrapper<PBackend>
				// // > P1;
				
	// // typedef amgcl::mpi::make_solver<
					// // amgcl::runtime::mpi::amg<
						// // PBackend,
						// // amgcl::runtime::mpi::coarsening::wrapper<PBackend>,
						// // amgcl::runtime::mpi::relaxation::wrapper<PBackend>
					// // >,
					// // amgcl::runtime::mpi::solver::wrapper<PBackend>
				// // > U1;
				
	// // typedef amgcl::runtime::mpi::solver::wrapper<SBackend> S1;
	
	// // typedef amgcl::runtime::mpi::schur_pressure_correction<P1,U1> SCHUR;
	
    // // typedef amgcl::::mpi::make_solver<SCHUR,S1> Solver;
				
        // // typedef
            // // amgcl::mpi::make_solver<
                // // amgcl::mpi::amg<
                    // // SBackend,
                    // // amgcl::runtime::mpi::coarsening::wrapper<SBackend>,
                    // // amgcl::runtime::mpi::relaxation::wrapper<SBackend>,
                    // // amgcl::runtime::mpi::direct::solver<double>,
                    // // amgcl::runtime::mpi::partition::wrapper<SBackend>
                    // // >,
                // // amgcl::runtime::mpi::solver::wrapper<SBackend>
                // // >
            // // Solver;
			



			
    // typedef
        // amgcl::mpi::make_solver<
            // amgcl::mpi::schur_pressure_correction<
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::block_preconditioner<
                        // amgcl::relaxation::as_preconditioner<SBackend, amgcl::runtime::relaxation::wrapper>
                        // >,
                    // amgcl::runtime::mpi::solver::wrapper<SBackend>
                    // >,
                // amgcl::mpi::make_solver<
                    // amgcl::mpi::amg<
                        // SBackend,
                        // amgcl::runtime::mpi::coarsening::wrapper<SBackend>,
                        // amgcl::runtime::mpi::relaxation::wrapper<SBackend>,
                        // amgcl::runtime::mpi::direct::solver<double>,
                        // amgcl::runtime::mpi::partition::wrapper<SBackend>
                        // >,
                        // amgcl::runtime::mpi::solver::wrapper<SBackend>
                    // >
                // >,
                // amgcl::runtime::mpi::solver::wrapper<SBackend>
            // > Solver;



	// // boost::property_tree::ptree prm;
	// // std::vector<char> pm(chunk, 0);
	// // for(ptrdiff_t i = 0; i < n_row; i += B_n) pm[i] = 1;
	// // prm.put("precond.pmask", static_cast<void*>(&pm[0]));
	// // prm.put("precond.pmask_size", n_row/B_n);

    // // std::function<double(ptrdiff_t,unsigned)> dv = amgcl::mpi::constant_deflation(1);
    // // prm.put("precond.psolver.num_def_vec", 1);
    // // prm.put("precond.psolver.def_vec", &dv);
			
	// // prm.put("solver.type", "cg");
	// // prm.put("solver.maxiter", 200);
	// // // prm.put("solver.tol", 1e-4);
	// // // prm.put("solver.M", 50);
	
	// // prm.put("precond.simplec_dia", false);
	
	// // prm.put("precond.usolver.solver.type", "preonly");
	// // prm.put("precond.usolver.precond.coarsening.type", "aggregation");
	// // prm.put("precond.usolver.precond.relax.type", "ilu0");
	// // // prm.put("precond.usolver.solver.tol", 1e-3);
	// // // prm.put("precond.usolver.solver.maxiter", 5);
	
	// // prm.put("precond.psolver.solver.type", "preonly");
	// // prm.put("precond.psolver.precond.class", "relaxation");
	// // // prm.put("precond.psolver.isolver.tol", 1e-2);
	// // // prm.put("precond.psolver.isolver.maxiter", 20);
	// // // prm.put("precond.psolver.local.coarse_enough", 500);
	

    // // // Solver solve(world, std::tie(chunk, ptr, col, val), prm, bprm);
    // // Solver solve(world, std::tie(chunk, ptr, col, val), prm);
    // // double tm_setup = prof.toc("setup");
	
	
    // // size_t iters; double error;
    // // std::tie(iters, error) = solve(B_vals, resiVar);
    // // double tm_solve = prof.toc("solve");

    // Solver::params prm;

	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // world, std::tie(chunk, ptr, col, val));
    // prof.tic("setup");
    // Solver solve(world, A, prm);
    // prof.toc("setup");
    // prof.tic("solve");
    // int iters; double error;
    // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // prof.toc("solve");

    // if (world.rank == 0) {
        // std::cout
            // << "Iterations: " << iters << std::endl
            // << "Error:      " << error << std::endl
            // << std::endl
            // << prof << std::endl;
	// }

    // // // Solver::params prm;
    // // // prm.solver.maxiter = 1000;
    // // // prm.solver.tol = 1.e-12;
    // // // prm.solver.replacement = true;
    // // // prm.solver.smoothing = true;
    // // // prm.solver.s = 5;
    // // // // prm.precond.simplec_dia = false;
    // // // // prm.precond.pmask.resize(n_row,0);
	// // // // for(int i=0; i<n_row/B_n; ++i){
		// // // // prm.precond.pmask[i*B_n+0] = 1;
	// // // // }
    // // // // fill(prm.precond.pmask.begin() + 58752, prm.precond.pmask.end(),1);


	// // auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // // world, std::tie(chunk, ptr, col, val));
		
    // // // Initialize the solver:
    // // prof.tic("setup");
    // // Solver solve(world, A, prm);
    // // prof.toc("setup");
	
    // // // Show the mini-report on the constructed solver:
    // // if (world.rank == 0)
        // // std::cout << solve << std::endl;

    // // // Solve the system with the zero initial approximation:
    // // prof.tic("solve");
    // // int iters; double error;
    // // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // // // std::tie(iters, error) = solve(B_vals, resiVar);
    // // prof.toc("solve");
	
    // // // Output the number of iterations, the relative error,
    // // // and the profiling data:
    // // if (world.rank == 0)
        // // std::cout
            // // << "Iters: " << iters << std::endl
            // // << "Error: " << error << std::endl
            // // << prof << std::endl;
	
// }











// #include "build.h"
// #include "mpi.h"

// #define AMGCL_NO_BOOST

// #include <amgcl/backend/builtin.hpp>
// #include <amgcl/adapter/crs_tuple.hpp>
// // #include <amgcl/io/binary.hpp>
// #include <amgcl/profiler.hpp>

// // boost
// #include "../../../boost/boost_1_77_0/boost/property_tree/ptree.hpp"
// #include "../../../boost/boost_1_77_0/boost/property_tree/json_parser.hpp"
// #include "../../../boost/boost_1_77_0/boost/foreach.hpp"
// #include "../../../boost/boost_1_77_0/boost/range/iterator_range.hpp"
// #include "../../../boost/boost_1_77_0/boost/scope_exit.hpp"

// // mpi
// #include <amgcl/mpi/distributed_matrix.hpp>
// #include <amgcl/mpi/make_solver.hpp>

// #include <amgcl/mpi/coarsening/runtime.hpp>
// #include <amgcl/mpi/relaxation/runtime.hpp>
// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/mpi/partition/runtime.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>
		
// // preconditioners
// #include <amgcl/amg.hpp>
// #include <amgcl/mpi/amg.hpp>
// #include <amgcl/mpi/relaxation/as_preconditioner.hpp>
// #include <amgcl/mpi/cpr.hpp>
// #include <amgcl/mpi/schur_pressure_correction.hpp>
// #include <amgcl/mpi/block_preconditioner.hpp>

// // relaxation
// #include <amgcl/relaxation/runtime.hpp>
// #include <amgcl/mpi/relaxation/spai0.hpp>
// #include <amgcl/mpi/relaxation/ilu0.hpp>
// #include <amgcl/relaxation/as_preconditioner.hpp>
// // #include <amgcl/relaxation/damped_jacobi.hpp>
// // #include <amgcl/relaxation/gauss_seidel.hpp>
// // #include <amgcl/relaxation/chebyshev.hpp>
// // #include <amgcl/relaxation/ilu0.hpp>
// // #include <amgcl/relaxation/iluk.hpp>
// // #include <amgcl/relaxation/ilup.hpp>
// // #include <amgcl/relaxation/ilut.hpp>
// // #include <amgcl/relaxation/spai0.hpp>
// // #include <amgcl/relaxation/spai1.hpp>

// // coarsening strategies
// #include <amgcl/mpi/coarsening/aggregation.hpp>
// #include <amgcl/mpi/coarsening/smoothed_aggregation.hpp>
// #include <amgcl/coarsening/runtime.hpp>

// // solver
// #include <amgcl/mpi/subdomain_deflation.hpp>
// #include <amgcl/mpi/solver/runtime.hpp>
// #include <amgcl/mpi/direct_solver/runtime.hpp>
// #include <amgcl/mpi/solver/idrs.hpp>
// #include <amgcl/mpi/solver/bicgstab.hpp>
// #include <amgcl/mpi/solver/preonly.hpp>
// // #include <amgcl/solver/cg.hpp>
// // #include <amgcl/solver/bicgstab.hpp>
// // #include <amgcl/solver/bicgstabl.hpp>
// // #include <amgcl/solver/gmres.hpp>
// // #include <amgcl/solver/lgmres.hpp>
// // #include <amgcl/solver/fgmres.hpp>
// // #include <amgcl/solver/preonly.hpp>


// // #include <amgcl/mpi/distributed_matrix.hpp>


// // #include <boost/program_options.hpp>
// // #include <boost/property_tree/ptree.hpp>
// // #include <boost/property_tree/json_parser.hpp>
// // #include <boost/range/iterator_range.hpp>
// // #include <boost/scope_exit.hpp>

// // #include <amgcl/io/binary.hpp>
// // #include <amgcl/io/mm.hpp>
// // #include <amgcl/adapter/crs_tuple.hpp>
// // #include <amgcl/amg.hpp>
// // #include <amgcl/coarsening/runtime.hpp>
// // #include <amgcl/relaxation/runtime.hpp>
// // #include <amgcl/relaxation/as_preconditioner.hpp>
// // #include <amgcl/mpi/make_solver.hpp>
// // #include <amgcl/mpi/schur_pressure_correction.hpp>
// // #include <amgcl/mpi/block_preconditioner.hpp>
// // #include <amgcl/mpi/subdomain_deflation.hpp>
// // #include <amgcl/mpi/solver/runtime.hpp>
// // #include <amgcl/mpi/direct_solver/runtime.hpp>
// // #include <amgcl/profiler.hpp>





// void SEMO_Solvers_Builder::solveAMGCL(
	// SEMO_Mesh_Builder& mesh,
	// int B_n, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B_vals,
	// vector<double>& resiVar){
		
	// // int rank = MPI::COMM_WORLD.Get_rank(); 
	// // int size = MPI::COMM_WORLD.Get_size(); 
	
	// // // COO to CSR
	// // int strRow = mesh.startCellGlobal * B_n;
	
    // // //compute number of non-zero entries per row of A 
	// // int n_row = B_vals.size();
	// // int nnz = A_rows.size();
	
	// // vector<double> val(nnz,0.0);
	// // vector<int> col(nnz,0);
	// // vector<int> ptr(n_row,0);

    // // for (int n = 0; n < nnz; n++){            
        // // ptr[A_rows[n]-strRow]++;
    // // }

    // // //cumsum the nnz per row to get ptr[]
    // // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // // int temp = ptr[i];
        // // ptr[i] = cumsum;
        // // cumsum += temp;
    // // }
    // // ptr[n_row] = nnz; 

    // // //write col,val into col,val
    // // for(int n = 0; n < nnz; n++){
        // // int row  = A_rows[n] - strRow;
        // // int dest = ptr[row];

        // // col[dest] = A_cols[n];
        // // val[dest] = A_vals[n];

        // // ptr[row]++;
    // // }
	// // // MPI_Barrier(MPI_COMM_WORLD);
	// // // cout << "AAAAAAAAAAAAAAA" << endl;

    // // for(int i = 0, last = 0; i <= n_row; i++){
        // // int temp = ptr[i];
        // // ptr[i]  = last;
        // // last   = temp;
    // // }

    // // //now ptr,col,val form a CSR representation (with possible duplicates)
	
	// // amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
    // // // The profiler:
    // // amgcl::profiler<> prof("AMGCL MPI");
	
	
	// // ptrdiff_t chunk = B_vals.size();
		

	// // /*
    // // // Compose the solver type
    // // typedef amgcl::backend::builtin<double> DBackend;
    // // typedef amgcl::backend::builtin<double> FBackend;
    // // typedef amgcl::mpi::make_solver<
        // // amgcl::mpi::amg<
            // // FBackend,
            // // amgcl::mpi::coarsening::smoothed_aggregation<FBackend>,
            // // amgcl::mpi::relaxation::spai0<FBackend>
            // // >,
        // // amgcl::mpi::solver::bicgstab<DBackend>
        // // > Solver;
	// // */
	
	
    // // // Compose the solver type
    // // typedef amgcl::backend::builtin<double> SBackend; // the outer iterative solver backend
    // // typedef amgcl::backend::builtin<float> PBackend;  // the PSolver backend
    // // //typedef amgcl::backend::builtin<amgcl::static_matrix<float, 3, 3>> UBackend;    // the USolver backend
    // // typedef amgcl::backend::builtin<float> UBackend;  // the USolver backend

    // // // typedef amgcl::mpi::make_solver
        // // // <
            // // // amgcl::mpi::schur_pressure_correction
            // // // <
                // // // //amgcl::make_block_solver
                // // // amgcl::mpi::make_solver
                // // // <
                    // // // amgcl::mpi::amg
                    // // // <UBackend,amgcl::mpi::coarsening::aggregation<UBackend>,amgcl::mpi::relaxation::ilu0<UBackend>>,
                    // // // // <UBackend,amgcl::mpi::coarsening::smoothed_aggregation<UBackend>,amgcl::mpi::relaxation::ilu0<UBackend>>,
                    // // // amgcl::mpi::solver::preonly<UBackend>
                // // // >,
                // // // amgcl::mpi::make_solver
                // // // <
                    // // // amgcl::mpi::relaxation::as_preconditioner
                    // // // <amgcl::mpi::relaxation::spai0<PBackend>>,
                    // // // amgcl::mpi::solver::preonly<PBackend>
                // // // >
            // // // >,
            // // // amgcl::mpi::solver::idrs<SBackend>
        // // // > Solver;
		


    // // // typedef
        // // // amgcl::mpi::make_solver<
            // // // amgcl::mpi::schur_pressure_correction<
                // // // amgcl::mpi::make_solver<
                    // // // amgcl::mpi::block_preconditioner<
                        // // // amgcl::relaxation::as_preconditioner<SBackend, amgcl::runtime::relaxation::wrapper>
                        // // // >,
                    // // // amgcl::runtime::mpi::solver::wrapper<SBackend>
                    // // // >,
                // // // amgcl::mpi::subdomain_deflation<
                    // // // amgcl::amg<SBackend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>,
                    // // // amgcl::runtime::mpi::solver::wrapper<SBackend>,
                    // // // amgcl::runtime::mpi::direct::solver<double>
                    // // // >
                // // // >,
            // // // amgcl::runtime::mpi::solver::wrapper<SBackend>
            // // // > Solver;
		
		
    // // // typedef amgcl::mpi::make_solver<
        // // // amgcl::mpi::amg<
            // // // PBackend,
            // // // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
            // // // amgcl::mpi::relaxation::spai0<PBackend>
            // // // >,
        // // // amgcl::mpi::solver::bicgstab<SBackend>
        // // // > Solver;
		
		
	// // // typedef amgcl::make_solver<
		// // // amgcl::amg<
			// // // SBackend,
			// // // amgcl::runtime::coarsening::wrapper,
			// // // amgcl::runtime::relaxation::wrapper
		// // // >,
		// // // amgcl::runtime::solver::wrapper<SBackend>
		// // // > Solver;
		
		
	// // // typedef amgcl::mpi::make_solver<
				// // // amgcl::mpi::relaxation::as_preconditioner
				// // // <amgcl::mpi::relaxation::ilu0<PBackend>>,
				// // // amgcl::mpi::solver::preonly<PBackend>
				// // // > P1;
				
	// // // typedef amgcl::mpi::make_solver<
					// // // amgcl::mpi::amg<
						// // // UBackend,
						// // // amgcl::mpi::coarsening::smoothed_aggregation<UBackend>,
						// // // amgcl::mpi::relaxation::ilu0<UBackend>
					// // // >,
					// // // amgcl::mpi::solver::bicgstab<UBackend>
				// // // > U1;
				
	// // // typedef amgcl::mpi::solver::idrs<SBackend> S1;
	
	// // // typedef amgcl::mpi::schur_pressure_correction<P1,U1> SCHUR;
		
	// // // typedef amgcl::mpi::make_solver<
					// // // amgcl::mpi::amg<
						// // // UBackend,
						// // // amgcl::mpi::coarsening::smoothed_aggregation<UBackend>,
						// // // amgcl::mpi::relaxation::ilu0<UBackend>
					// // // >,
					// // // amgcl::mpi::solver::bicgstab<UBackend>
				// // // > Solver;
				
    // // // typedef amgcl::mpi::make_solver<P1,S1> Solver;
		




    // // typedef amgcl::mpi::make_solver<
        // // amgcl::mpi::amg<
            // // PBackend,
            // // amgcl::mpi::coarsening::smoothed_aggregation<PBackend>,
            // // amgcl::mpi::relaxation::ilu0<PBackend>
            // // >,
        // // amgcl::mpi::solver::idrs<SBackend>
        // // > Solver;




    // // Solver::params prm;
    // // prm.solver.maxiter = 1000;
    // // prm.solver.tol = 1.e-12;
    // // prm.solver.replacement = true;
    // // prm.solver.smoothing = true;
    // // prm.solver.s = 5;
    // // // prm.precond.simplec_dia = false;
    // // // prm.precond.pmask.resize(n_row,0);
	// // // for(int i=0; i<n_row/B_n; ++i){
		// // // prm.precond.pmask[i*B_n+0] = 1;
	// // // }
    // // // fill(prm.precond.pmask.begin() + 58752, prm.precond.pmask.end(),1);


	// // auto A = std::make_shared<amgcl::mpi::distributed_matrix<SBackend>>(
		     // // world, std::tie(chunk, ptr, col, val));
		
    // // // Initialize the solver:
    // // prof.tic("setup");
    // // Solver solve(world, A, prm);
    // // prof.toc("setup");
	
    // // // Show the mini-report on the constructed solver:
    // // if (world.rank == 0)
        // // std::cout << solve << std::endl;

    // // // Solve the system with the zero initial approximation:
    // // prof.tic("solve");
    // // int iters; double error;
    // // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // // // std::tie(iters, error) = solve(B_vals, resiVar);
    // // prof.toc("solve");
	
    // // // Output the number of iterations, the relative error,
    // // // and the profiling data:
    // // if (world.rank == 0)
        // // std::cout
            // // << "Iters: " << iters << std::endl
            // // << "Error: " << error << std::endl
            // // << prof << std::endl;
	
// }


// void SEMO_Solvers_Builder::solveAMGCL(
	// SEMO_Mesh_Builder& mesh,
	// vector<double>& resiVar, 
	// vector<int>& A_rows, vector<int>& A_cols, vector<double>& A_vals, 
	// vector<double>& B_vals,
	// int ncells, int ncellTot,
	// string solver, double tolerance, double relTol, string preconditioner,
	// int maxIter
	// ){
		
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 
	
	// // int str_ncell = mesh.startProcCellGlobal[rank];
	// // std::vector<int>    ptr;
	// // std::vector<int>    col;
	// // // std::vector<double> val;
	// // std::vector<int>    val;
	// // int rows = mesh.cells.size();
	// // int nonzeros = 0;
	
	// // vector<vector<int>> unsorted_row_nonzeros(mesh.cells.size(),vector<int>());
	// // for(int i=0; i<mesh.cells.size(); ++i){
		// // unsorted_row_nonzeros[i].push_back(
			// // str_ncell + i);
	// // }
	
	
	// // int proc_num = 0;
	// // vector<int> rowNonzeros(mesh.cells.size(),1);
	// // for(int i=0; i<mesh.faces.size(); ++i){
		// // auto& face = mesh.faces[i];
		// // if(face.getType() == SEMO_Types::INTERNAL_FACE){
			// // ++rowNonzeros[face.owner];
			// // ++rowNonzeros[face.neighbour];
			// // unsorted_row_nonzeros[face.owner].push_back(
				// // str_ncell + face.neighbour);
			// // unsorted_row_nonzeros[face.neighbour].push_back(
				// // str_ncell + face.owner);
		// // }
		
		// // if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			// // ++rowNonzeros[face.owner];
			// // unsorted_row_nonzeros[face.owner].push_back(
				// // mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] +
				// // mesh.procNeighbCellNo[proc_num]);
			
			// // ++proc_num;
		// // }
	// // }
	
	// // ptr.push_back(0);
	// // for(auto& ncols : rowNonzeros){
		// // int strPtr = ptr.back();
		// // ptr.push_back(strPtr + ncols);
	// // }
	
	// // int non_zeros = ptr.back();
	// // col.reserve(non_zeros);
	// // val.reserve(non_zeros);
	// // for(int i=0; i<rows; ++i){
		// // std::sort(unsorted_row_nonzeros[i].begin(), 
				  // // unsorted_row_nonzeros[i].end());
		// // int nn = 0;
		// // for(int j=ptr[i]; j<ptr[i+1]; ++j){
			// // col[j] = unsorted_row_nonzeros[i][nn];
			// // val[j] = unsorted_row_nonzeros[i][nn];
			// // ++nn;
		// // }
	// // }
	
	
	
	
	// // COO to CSR
	// int strRow = mesh.startCellGlobal;
	// if(rank==0){
		// strRow = 0;
	// }
	// else if(rank==1){
		// strRow = 3;
	// }
	// else if(rank==2){
		// strRow = 5;
	// }
	// else if(rank==3){
		// strRow = 7;
	// }

	
    // //compute number of non-zero entries per row of A 
	// int n_row = B_vals.size();
	// int nnz = A_rows.size();
	
	// vector<double> val(nnz,0.0);
	// vector<int> col(nnz,0);
	// vector<int> ptr(n_row,0);

    // for (int n = 0; n < nnz; n++){            
        // ptr[A_rows[n]-strRow]++;
    // }

    // //cumsum the nnz per row to get ptr[]
    // for(int i = 0, cumsum = 0; i < n_row; i++){     
        // int temp = ptr[i];
        // ptr[i] = cumsum;
        // cumsum += temp;
    // }
    // ptr[n_row] = nnz; 

    // //write col,val into col,val
    // for(int n = 0; n < nnz; n++){
        // int row  = A_rows[n] - strRow;
        // int dest = ptr[row];

        // col[dest] = A_cols[n];
        // val[dest] = A_vals[n];

        // ptr[row]++;
    // }

    // for(int i = 0, last = 0; i <= n_row; i++){
        // int temp = ptr[i];
        // ptr[i]  = last;
        // last   = temp;
    // }

    // //now ptr,col,val form a CSR representation (with possible duplicates)
	
		
	
	// amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
	
    // // The profiler:
    // amgcl::profiler<> prof("AMGCL MPI");
	
	
	// ptrdiff_t chunk = B_vals.size();
	
    // // Compose the solver type
    // typedef amgcl::backend::builtin<double> DBackend;
    // typedef amgcl::backend::builtin<double> FBackend;
    // typedef amgcl::mpi::make_solver<
        // amgcl::mpi::amg<
            // FBackend,
            // amgcl::mpi::coarsening::smoothed_aggregation<FBackend>,
            // amgcl::mpi::relaxation::spai0<FBackend>
            // >,
        // amgcl::mpi::solver::bicgstab<DBackend>
        // > Solver;


	// // cout << "setup A!!!!!!!!" << endl;
	
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<DBackend>>(
		     // world, std::tie(chunk, ptr, col, val));
		
	// // cout << "finished setup A!!!!!!!!" << endl;

    // // Initialize the solver:
    // prof.tic("setup");
    // Solver solve(world, A);
    // prof.toc("setup");
	

	// // cout << "finished 2 setup A!!!!!!!!" << endl;

    // // Show the mini-report on the constructed solver:
    // if (world.rank == 0)
        // std::cout << solve << std::endl;

    // // Solve the system with the zero initial approximation:
    // int iters;
    // double error;
    // // std::vector<double> x(chunk, 0.0);

	// // cout << "SOLVE!!!!!!!!" << endl;

    // prof.tic("solve");
    // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // prof.toc("solve");
	
	// for(int loc=0; loc<size; ++loc){
		// if(rank==loc){
			// cout << "proc = " << loc << endl;
			// for(auto& i : resiVar){
				// cout << i << endl;
			// }
		// }
		// MPI_Barrier(MPI_COMM_WORLD);
	// }
    // // Output the number of iterations, the relative error,
    // // and the profiling data:
    // if (world.rank == 0)
        // std::cout
            // << "Iters: " << iters << std::endl
            // << "Error: " << error << std::endl
            // << prof << std::endl;
	
		
		
// }


void SEMO_Solvers_Builder::solveAMGCL_Flows(
	SEMO_Mesh_Builder& mesh,
	int B_n,
	vector<double>& A_vals, 
	vector<double>& B_vals,
	vector<double>& resiVar
	){
		
		
	// int rank = MPI::COMM_WORLD.Get_rank(); 
	// int size = MPI::COMM_WORLD.Get_size(); 
	
	
	// amgcl::mpi::communicator world(MPI_COMM_WORLD);
	
	
    // // The profiler:
    // amgcl::profiler<> prof("AMGCL MPI");
	
	// int nrow_loc = B_vals.size()/B_n;
	// int nrow = nrow_loc*B_n;
	// ptrdiff_t chunk = nrow;
	
    // // Compose the solver type
    // typedef amgcl::backend::builtin<double> DBackend;
    // typedef amgcl::backend::builtin<double> FBackend;
    // typedef amgcl::mpi::make_solver<
        // amgcl::mpi::amg<
            // FBackend,
            // amgcl::mpi::coarsening::smoothed_aggregation<FBackend>,
            // amgcl::mpi::relaxation::spai0<FBackend>
            // >,
        // amgcl::mpi::solver::bicgstab<DBackend>
        // > Solver;

	// cout << "AAAAAAAAAAAAAAA " << nrow_loc << " " << mesh.cells.size() << endl;
	// vector<int> ptr(nrow+1,0);
	// for(int i=0; i<nrow_loc; ++i){
		// for(int si=0; si<B_n; ++si){
			// int n_col_loc = (mesh.CRS_ptr[i+1]-mesh.CRS_ptr[i]);
			// ptr[B_n*i+si+1] = ptr[B_n*i+si] + n_col_loc*B_n;
		// }
	// }
	
	// int nnz_loc = mesh.non_zeros;
	// int nnz = nnz_loc*B_n*B_n;
	// vector<int> col(nnz,0);
	
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "BBBBBBBBBBBBB " << nnz << endl;
	
	// for(int i=0, tmpNum=0; i<nrow_loc; ++i){
		// int str_col_loc = mesh.CRS_ptr[i];
		// int end_col_loc = mesh.CRS_ptr[i+1];
		// int n_col_loc = end_col_loc - str_col_loc;
		// for(int sj=0; sj<B_n; ++sj){
			// for(int si=0; si<B_n; ++si){
				// for(int k=str_col_loc; k<end_col_loc; ++k){
					// col[tmpNum] = nrow_loc*si + mesh.CRS_col[k];
					// ++tmpNum;
				// }
			// }
		// }
	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// cout << "CCCCCCCCCCCCC " << ptr[nrow] << " " << nnz << " " << A_vals.size() << endl;

	// // cout << "setup A!!!!!!!!" << endl;
	
	// auto A = std::make_shared<amgcl::mpi::distributed_matrix<DBackend>>(
		     // world, std::tie(chunk, ptr, col, A_vals));
		
	// cout << "finished setup A!!!!!!!!" << endl;

    // // Initialize the solver:
    // prof.tic("setup");
    // Solver solve(world, A);
    // prof.toc("setup");
	

	// // // cout << "finished 2 setup A!!!!!!!!" << endl;

    // // // Show the mini-report on the constructed solver:
    // // if (world.rank == 0)
        // // std::cout << solve << std::endl;

    // // // Solve the system with the zero initial approximation:
    // // int iters;
    // // double error;
    // // // std::vector<double> x(chunk, 0.0);

	// // // cout << "SOLVE!!!!!!!!" << endl;

    // // prof.tic("solve");
    // // std::tie(iters, error) = solve(*A, B_vals, resiVar);
    // // prof.toc("solve");
	
    // // // Output the number of iterations, the relative error,
    // // // and the profiling data:
    // // if (world.rank == 0)
        // // std::cout
            // // << "Iters: " << iters << std::endl
            // // << "Error: " << error << std::endl
            // // << prof << std::endl;
	
		
		
}