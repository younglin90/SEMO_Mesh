#include "build.h"
#include "mpi.h"

#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/adapter/reorder.hpp>


void SEMO_Mesh_Builder::reorder(){
	
	
	int rank = MPI::COMM_WORLD.Get_rank(); 
	int size = MPI::COMM_WORLD.Get_size(); 
	
	SEMO_Mesh_Builder& mesh = *this;
	
	int str_ncell = mesh.startProcCellGlobal[rank];
	std::vector<int>    ptr;
	std::vector<int>    col;
	// std::vector<double> val;
	std::vector<int>    val;
	int rows = mesh.cells.size();
	int nonzeros = 0;
	
	vector<vector<int>> unsorted_row_nonzeros(mesh.cells.size(),vector<int>());
	for(int i=0; i<mesh.cells.size(); ++i){
		unsorted_row_nonzeros[i].push_back(
			str_ncell + i);
	}
	
	
	int proc_num = 0;
	vector<int> rowNonzeros(mesh.cells.size(),1);
	for(int i=0; i<mesh.faces.size(); ++i){
		auto& face = mesh.faces[i];
		if(face.getType() == SEMO_Types::INTERNAL_FACE){
			++rowNonzeros[face.owner];
			++rowNonzeros[face.neighbour];
			unsorted_row_nonzeros[face.owner].push_back(
				str_ncell + face.neighbour);
			unsorted_row_nonzeros[face.neighbour].push_back(
				str_ncell + face.owner);
		}
		
		if(face.getType() == SEMO_Types::PROCESSOR_FACE){
			++rowNonzeros[face.owner];
			unsorted_row_nonzeros[face.owner].push_back(
				mesh.startProcCellGlobal[mesh.neighbProcNo[proc_num]] +
				mesh.procNeighbCellNo[proc_num]);
			
			++proc_num;
		}
	}
	
	ptr.push_back(0);
	for(auto& ncols : rowNonzeros){
		int strPtr = ptr.back();
		ptr.push_back(strPtr + ncols);
	}
	
	int non_zeros = ptr.back();
	col.reserve(non_zeros);
	val.reserve(non_zeros);
	for(int i=0; i<rows; ++i){
		std::sort(unsorted_row_nonzeros[i].begin(), 
				  unsorted_row_nonzeros[i].end());
		int nn = 0;
		for(int j=ptr[i]; j<ptr[i+1]; ++j){
			col[j] = unsorted_row_nonzeros[i][nn];
			val[j] = unsorted_row_nonzeros[i][nn];
			++nn;
		}
	}
	
	auto A = std::tie(rows, ptr, col, val);
	
	
    // auto Ab = amgcl::adapter::block_matrix<int>(As);
	// const int B = rows;
	// typedef amgcl::static_matrix<double, B, B> value_type;
	// auto A = amgcl::adapter::block_matrix<value_type>(std::tie(rows, ptr, col, val));
	
	// Prepare the reordering:
	// amgcl::adapter::reorder<> perm(A);

	// // Create the solver using the reordered matrix:
	// Solver solve(perm(A), prm);

	// // Reorder the RHS and solve the system:
	// solve(perm(rhs), x_ord);

	// // Postprocess the solution vector to get the original ordering:
	// perm.inverse(x_ord, x);
	
	
}
