#include "ShearCondition.h"

template <class Real>
ShearCondition<Real>::ShearCondition(Real k, int i0, int i1, int i2, Real a) : m_k(k)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
	m_area = a;
}

template <class Real>
typename EnergyCondition<Real>::Vector ShearCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	return Vector();
}

template <class Real>
void ShearCondition<Real>::computeForces(const Vector& x, const Vector& uv, Vector& forces) const
{

}

template <class Real>
void ShearCondition<Real>::computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{

}

template <class Real>
void ShearCondition<Real>::computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const
{

}

template <class Real>
void ShearCondition<Real>::computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{

}

template class ShearCondition<float>;
template class ShearCondition<double>;
