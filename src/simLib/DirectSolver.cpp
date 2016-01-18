#include "simLib\DirectSolver.h"

template<class Real>
void DirectSolver<Real>::solve(const SparseMatrix &A, const Vector &rhs, Vector &result) const
{
	Eigen::SparseLU< SparseMatrix > solver;
	solver.analyzePattern(A);
	solver.factorize(A);
	result = solver.solve(rhs);
}

template class DirectSolver<float>;
template class DirectSolver<double>;
