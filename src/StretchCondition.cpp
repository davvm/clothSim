#include "StretchCondition.h"

template <class Real>
StretchCondition<Real>::StretchCondition(Real k, int i0, int i1, int i2, Real a, Real restU, Real restV) : m_k(k)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
	m_area = a;
}

template <class Real>
typename EnergyCondition<Real>::Vector StretchCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	return Vector();
}

template <class Real>
void StretchCondition<Real>::computeForces(const Vector& x, const Vector& uv, Vector& forces) const
{

}

template <class Real>
void StretchCondition<Real>::computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{

}

template <class Real>
void StretchCondition<Real>::computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const
{

}

template <class Real>
void StretchCondition<Real>::computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{

}

template class StretchCondition<float>;
template class StretchCondition<double>;
