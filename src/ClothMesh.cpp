#include "ClothMesh.h"

static std::vector<int> sharedVertices(const std::vector<int> &triangleIndices, size_t i, size_t j)
{
	std::vector<int> ret;
	if (triangleIndices[i + 0] == triangleIndices[j + 0]) ret.push_back(triangleIndices[i + 0]);
	if (triangleIndices[i + 0] == triangleIndices[j + 1]) ret.push_back(triangleIndices[i + 0]);
	if (triangleIndices[i + 0] == triangleIndices[j + 2]) ret.push_back(triangleIndices[i + 0]);
	if (triangleIndices[i + 1] == triangleIndices[j + 0]) ret.push_back(triangleIndices[i + 1]);
	if (triangleIndices[i + 1] == triangleIndices[j + 1]) ret.push_back(triangleIndices[i + 1]);
	if (triangleIndices[i + 1] == triangleIndices[j + 2]) ret.push_back(triangleIndices[i + 1]);
	if (triangleIndices[i + 2] == triangleIndices[j + 0]) ret.push_back(triangleIndices[i + 2]);
	if (triangleIndices[i + 2] == triangleIndices[j + 1]) ret.push_back(triangleIndices[i + 2]);
	if (triangleIndices[i + 2] == triangleIndices[j + 2]) ret.push_back(triangleIndices[i + 2]);
	return ret;
}

template<class Real>
ClothMesh<Real>::ClothMesh(
	const Vector &x,
	const Vector &v,
	const Vector &uv,
	const std::vector<int> &triangleIndices,
	Real kBend, Real kStretch, Real kShear,
	Real dBend, Real dStretch, Real dShear,
	Real density
	) :
	m_kBend(kBend), m_kStretch(kStretch), m_kShear(kShear),
	m_dBend(dBend), m_dStretch(dStretch), m_dShear(dShear),
	m_x(x),
	m_v(v),
	m_uv(uv),
	m_m((int)m_uv.size() / 2),
	m_dfdx((int)x.size(), (int)x.size()),
	m_dfdv((int)x.size(), (int)x.size()),
	m_implicitUpdateMatrix((int)x.size(), (int)x.size()),
	m_implicitUpdateRHS((int)x.size()),
	m_forces((int)x.size())
{
	for(size_t i = 0; i < triangleIndices.size(); i += 3 )
	{
		// add stretch and shear terms to the equations of motion:
		m_shearConditions.push_back(ShearCondition<Real>((int)i, (int)i + 1, (int)i + 2));
		m_stretchConditions.push_back(StretchCondition<Real>((int)i, (int)i + 1, (int)i + 2,1,1));

		// find triangle area:
		const Vector3 &p0 = m_x.segment<3>(3 * triangleIndices[i+0]);
		const Vector3 &p1 = m_x.segment<3>(3 * triangleIndices[i + 1]);
		const Vector3 &p2 = m_x.segment<3>(3 * triangleIndices[i + 2]);
		Real area = (p1 - p0).cross(p2 - p0).norm() / 2;

		// calculate triangle mass and distribute it equally over the 3 vertices:
		Real triangleMass = area * density;
		m_m[triangleIndices[i + 0]] += triangleMass / 3;
		m_m[triangleIndices[i + 1]] += triangleMass / 3;
		m_m[triangleIndices[i + 2]] += triangleMass / 3;
	}

	// find shared edges and instantiate BendConditions for them.
	// Doing it in a dumb slow way right now just so it works:
	for(size_t i = 0; i < triangleIndices.size(); i += 3)
	{
		for(size_t j = i + 1; j < triangleIndices.size(); j += 3)
		{
			std::vector<int> shared = sharedVertices(triangleIndices, i, j);
			if(shared.size() == 2)
			{
				// sharing an edge!
				int i1 = shared[0];
				int i2 = shared[1];
				int i0;
				for(size_t n = 0; n < 3; ++n)
				{
					if (triangleIndices[i + n] != i1 && triangleIndices[i + n] != i2)
					{
						i0 = triangleIndices[i + n];
						break;
					}
				}
				
				int i3;
				for (size_t n = 0; n < 3; ++n)
				{
					if (triangleIndices[i + n] != i1 && triangleIndices[i + n] != i2)
					{
						i3 = triangleIndices[i + n];
						break;
					}
				}
				m_bendConditions.push_back( BendCondition<Real>( i0, i1, i2, i3 ) );
			}
		}
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
	// reset existing sim data to zero:
	m_implicitUpdateRHS.setConstant(0);
	m_forces.setConstant(0);
	m_dv.setConstant(0);
	for (int i = 0; i < m_dfdx.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(m_dfdx, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}
	for (int i = 0; i < m_dfdv.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(m_dfdv, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}
	for (int i = 0; i < m_implicitUpdateMatrix.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(m_implicitUpdateMatrix, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}

	// build dfdx, dfdv and the forces:

	// add the bend terms:
	for (size_t i = 0; i < m_bendConditions.size(); ++i)
	{
		m_bendConditions[i].computeForces(m_x, m_uv, m_kShear, m_forces, m_dfdx, m_v, m_dShear, m_forces, m_dfdv);
	}

	// add the stretch terms:
	for (size_t i = 0; i < m_stretchConditions.size(); ++i)
	{
		m_stretchConditions[i].computeForces(m_x, m_uv, m_kShear, m_forces, m_dfdx, m_v, m_dShear, m_forces, m_dfdv);
	}

	// add the shear terms:
	for (size_t i = 0; i < m_shearConditions.size(); ++i)
	{
		m_shearConditions[i].computeForces(m_x, m_uv, m_kShear, m_forces, m_dfdx, m_v, m_dShear, m_forces, m_dfdv);
	}

	// The implicit update system we want to solve is this:
	// ( M - dt * dfdv - dt * dt * dfdx ) * dv = dt * ( f + dt * dfdx * v )

	// build the implicit update matrix:
	m_implicitUpdateMatrix += - dt * m_dfdv - dt * dt * m_dfdv;
	
	// add masses onto the diagonal:
	for (int i = 0; i < m_x.size(); ++i)
	{
		// unit masses for now...
		m_implicitUpdateMatrix.coeffRef(i, i) += m_m[i/3];
	}

	// build the right hand side:
	m_implicitUpdateRHS += dt * (m_forces + dt * m_dfdx * m_v);

	// solve the linear system:
	solver.solve(m_implicitUpdateMatrix, m_implicitUpdateRHS, m_dv);

	// update:
	m_v += m_dv;
	m_x += dt * m_v;
}

template class ClothMesh<float>;
template class ClothMesh<double>;
