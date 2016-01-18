#ifndef CONSTRAINEDCCGSOLVER_H
#define CONSTRAINEDCCGSOLVER_H

#include "LinearSolver.h"

// Solves the linear system using conjugate gradients,
// projecting out the given constraints at every step:
template <class Real>
class ConstrainedCGSolver : public LinearSolver<Real>
{

public:

	typedef Eigen::Matrix<Real, -1, 1, 0, -1, 1> Vector;
	typedef Eigen::Matrix<Real, 3, 3, 0, 3, 3> Matrix3;

	ConstrainedCGSolver(
		const std::vector<int> &constraintIndices,
		const std::vector< Matrix3 > &constraintMatrices,
		const Vector &constraintVelocityDeltas,
		Real tol,
		int maxIterations
	);

	// perform a linear solve (A * result = rhs):
	virtual void solve(const SparseMatrix &A, const Vector &rhs, Vector &result) const;

private:
	
	void filter(Vector &x) const;
	void precondition(Vector &x, Vector &p, bool inverse) const;

	const std::vector<int> &m_constraintIndices;
	const std::vector< Matrix3 > &m_constraintMatrices;
	const Vector &m_constraintVelocityDeltas;
	Real m_tol;
	int m_maxIterations;

};


#endif
