#include "ClothMesh.h"

template<class Real>
ClothMesh<Real>::ClothMesh(
	const Vector &x,
	const Vector &v,
	const Vector &m_uv,
	const std::vector<int> &triangleIndices,
	Real kBend, Real kStretch, Real kShear,
	Real dBend, Real dStretch, Real dShear
	) :
	m_kBend(kBend), m_kStretch(kStretch), m_kShear(kShear),
	m_dBend(dBend), m_dStretch(dStretch), m_dShear(dShear),
	m_dfdx(x.size()),
	m_implicitUpdateMatrix(x.size(), x.size()),
	m_implicitUpdateRHS(x.size())
{
	for (size_t i = 0; i < triangleIndices.size(); i += 3 )
	{
		m_shearConditions.push_back(ShearCondition<Real>((int)i, (int)i + 1, (int)i + 2));
		m_stretchConditions.push_back(StretchCondition<Real>((int)i, (int)i + 1, (int)i + 2,1,1));
	}
}

template<class Real>
const typename ClothMesh<Real>::Vector &ClothMesh<Real>::x() const
{
	return m_x;
}

// accessor for velocities:
template<class Real>
const typename ClothMesh<Real>::Vector &ClothMesh<Real>::v() const
{
	return m_v;
}

// advance the simulation:
template<class Real>
void ClothMesh<Real>::advance(Real dt, const LinearSolver<Real> &solver)
{

}

template class ClothMesh<float>;
template class ClothMesh<double>;
