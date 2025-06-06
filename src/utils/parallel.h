#pragma once

#include <Eigen/Core>

#include "CompFunction.h"
#include "mpi_utils.h"
#include "trees/MultiResolutionAnalysis.h"
#include <map>
#include <vector>

// define a class for things that can be sent with MPI

using namespace Eigen;

using IntVector = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;

using IntMatrix = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

namespace mrcpp {

namespace omp {
extern int n_threads;
} // namespace omp

class Bank;
extern Bank dataBank;

namespace mpi {

extern std::vector<int> bankmaster;

void initialize();
void finalize();
void barrier(MPI_Comm comm);

bool grand_master();
bool share_master();

bool my_func(int j);
bool my_func(const CompFunction<3> &func);
bool my_func(CompFunction<3> *func);

// bool my_unique_orb(const Orbital &orb);
void free_foreign(CompFunctionVector &Phi);

void send_function(const CompFunction<3> &func, int dst, int tag, MPI_Comm comm = mpi::comm_wrk);
void recv_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm = mpi::comm_wrk);
void share_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm);

void reduce_function(double prec, CompFunction<3> &func, MPI_Comm comm);
void broadcast_function(CompFunction<3> &func, MPI_Comm comm);

template <typename T> void reduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm);
template <typename T> void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, std::vector<FunctionTree<3, T>> &Phi, MPI_Comm comm);
template <typename T> void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm);

template <typename T> void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, std::vector<CompFunction<3>> &Phi, MPI_Comm comm);

void allreduce_vector(IntVector &vec, MPI_Comm comm);
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);

} // namespace mpi

} // namespace mrcpp
