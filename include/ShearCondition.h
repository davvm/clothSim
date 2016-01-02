#ifndef SHEARCONDITION_H
#define SHEARCONDITION_H

#include "EnergyCondition.h"

template <class Real>
class ShearCondition : public EnergyCondition<Real>
{

public:

	ShearCondition(int i0, int i1, int i2, Real area);
	virtual ~ShearCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const;

private:

	// indices of the vertices this shear condition is attached to:
	int m_inds[3];

	// rest area of the triangle this condition is attached to:
	Real m_area;
};


#endif // SHEARCONDITION_H
