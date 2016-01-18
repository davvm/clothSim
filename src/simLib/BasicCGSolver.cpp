#include "simLib\BasicCGSolver.h"
#include <Eigen/IterativeLinearSolvers>

template<class Real>
void BasicCGSolver<Real>::solve(const SparseMatrix &A, const Vector &rhs, Vector &result) const
{
	Eigen::ConjugateGradient<SparseMatrix > cg;
	cg.compute(A);
	result = cg.solve(rhs);
}

template class BasicCGSolver<float>;
template class BasicCGSolver<double>;
