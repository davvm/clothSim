#include "simLib\GravityField.h"

using namespace ClothSim;

template <class Real>
GravityField<Real>::GravityField(const Vector &m, Vector3 g) : m_m(m), m_g( g )
{
}

template <class Real>
void GravityField<Real>::forcesAndDerivatives(Vector& forces, SparseMatrix &dfdx) const
{
	for (int i = 0; i < forces.size(); i += 3)
	{
		forces.segment<3>(i) += m_m[i / 3] * m_g;
	}
}

template class GravityField<float>;
template class GravityField<double>;
