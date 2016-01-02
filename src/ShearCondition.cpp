#include "ShearCondition.h"

template <class Real>
ShearCondition<Real>::ShearCondition(int i0, int i1, int i2, Real a) : m_area(a)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
}

template <class Real>
typename EnergyCondition<Real>::Vector ShearCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	return Vector();
}

template <class Real>
void ShearCondition<Real>::computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const
{

}

template class ShearCondition<float>;
template class ShearCondition<double>;
