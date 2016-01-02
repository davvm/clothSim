#include "StretchCondition.h"

template <class Real>
StretchCondition<Real>::StretchCondition(int i0, int i1, int i2, Real a, Real restU, Real restV) : m_area(a), m_restU(restU), m_restV(restV)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
}

template <class Real>
typename EnergyCondition<Real>::Vector StretchCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	return Vector();
}

template <class Real>
void StretchCondition<Real>::computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const
{

}

template class StretchCondition<float>;
template class StretchCondition<double>;
