#ifndef BASICCGSOLVER_H
#define BASICCGSOLVER_H

#include "LinearSolver.h"

namespace ClothSim
{

// Solves the linear system using a direct method
template <class Real>
class BasicCGSolver : public LinearSolver<Real>
{

public:

	// perform a linear solve (A * result = rhs):
	virtual void solve(const SparseMatrix &A, const Vector &rhs, Vector &result) const;

};

} //namespace ClothSim


#endif
