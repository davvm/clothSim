#ifndef SHEARCONDITION_H
#define SHEARCONDITION_H

#include "EnergyCondition.h"

template <class Real>
class ShearCondition : public EnergyCondition<Real>
{

public:

	ShearCondition(Real k, int i0, int i1, int i2, Real area);
	virtual ~ShearCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Vector& forces) const;
	virtual void computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

	virtual void computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const;
	virtual void computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

private:

	// stiffness
	Real m_k;

	// indices of the vertices this shear condition is attached to:
	int m_inds[3];

	// rest area of the triangle this condition is attached to:
	Real m_area;
};


#endif // SHEARCONDITION_H
