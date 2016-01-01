#ifndef BENDCONDITION_H
#define BENDCONDITION_H

#include "EnergyCondition.h"

template <class Real>
class BendCondition : public EnergyCondition<Real>
{

public:

	BendCondition(int i0, int i1, int i2, int i3);
	virtual ~BendCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Vector& forces) const;
	virtual void computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

	virtual void computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const;
	virtual void computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

private:

	// bend condition is attached to two triangles, ie four vertices like so:
	//
	//     0
	//    / \
	//   /   \
	//  1-----2
	//   \   /
	//    \ /
	//     3
	//
	int m_inds[4];

};


#endif // BENDCONDITION_H
