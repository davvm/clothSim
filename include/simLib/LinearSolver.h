#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "Eigen/Dense"
#include "Eigen/Sparse"

// abstract class for solving linear systems
template <class Real>
class LinearSolver
{

public:

	typedef Eigen::Matrix<Real, -1, 1, 0, -1, 1> Vector;
	typedef Eigen::SparseMatrix<Real> SparseMatrix;

	// perform a linear solve (A * result = rhs):
	virtual void solve( const SparseMatrix &A, const Vector &rhs, Vector &result ) const = 0;

};

#endif
