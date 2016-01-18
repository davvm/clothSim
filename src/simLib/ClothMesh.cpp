#include "simLib\ClothMesh.h"

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
	m_x((int)x.size()),
	m_v((int)v.size()),
	m_uv((int)uv.size()),
	m_implicitUpdateMatrix((int)x.size(), (int)x.size()),
	m_dfdx((int)x.size(), (int)x.size()),
	m_dfdv((int)x.size(), (int)x.size()),
	m_triangleIndices( triangleIndices ),
	m_m((int)uv.size() / 2)
{
	for (int i = 0; i < x.size(); ++i)
	{
		m_x[i] = x[i];
	}
	for (int i = 0; i < v.size(); ++i)
	{
		m_v[i] = v[i];
	}
	for (int i = 0; i < uv.size(); ++i)
	{
		m_uv[i] = uv[i];
	}

	m_m.setConstant(0);
	for(size_t i = 0; i < triangleIndices.size(); i += 3 )
	{
		// add stretch and shear terms to the equations of motion:
		m_shearConditions.push_back(ShearCondition<Real>(triangleIndices[i + 0], triangleIndices[i + 1], triangleIndices[i + 2]));
		m_stretchConditions.push_back(StretchCondition<Real>(triangleIndices[i + 0], triangleIndices[i + 1], triangleIndices[i + 2], 1, 1));

		// find triangle rest area:
		const Vector2 &uv0 = m_uv.segment<2>(2 * triangleIndices[i+0]);
		const Vector2 &uv1 = m_uv.segment<2>(2 * triangleIndices[i + 1]);
		const Vector2 &uv2 = m_uv.segment<2>(2 * triangleIndices[i + 2]);
		
		Vector2 duv1 = uv1 - uv0;
		Vector2 duv2 = uv2 - uv0;

		Real du1 = duv1[0];
		Real dv1 = duv1[1];
		Real du2 = duv2[0];
		Real dv2 = duv2[1];

		// triangle area in reference pose:
		Real area = Real(0.5) * abs(du1 * dv2 - du2 * dv1);

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
		for(size_t j = i + 3; j < triangleIndices.size(); j += 3)
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
					if (triangleIndices[j + n] != i1 && triangleIndices[j + n] != i2)
					{
						i3 = triangleIndices[j + n];
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

// accessor for masses:
template<class Real>
const typename ClothMesh<Real>::Vector &ClothMesh<Real>::m() const
{
	return m_m;
}
template<class Real>
const std::vector< BendCondition<Real> > &ClothMesh<Real>::bendConditions() const
{
	return m_bendConditions;
}

template<class Real>
const std::vector< ShearCondition<Real> > &ClothMesh<Real>::shearConditions() const
{
	return m_shearConditions;
}

template<class Real>
const std::vector< StretchCondition<Real> > &ClothMesh<Real>::stretchConditions() const
{
	return m_stretchConditions;
}

template<class Real>
const std::vector<int> &ClothMesh<Real>::triangleIndices()
{
	return m_triangleIndices;
}

// advance the simulation:
template<class Real>
void ClothMesh<Real>::advance(std::vector< ForceField<Real>* > &forceFields, Real dt, const LinearSolver<Real> &solver)
{
	Vector forces((int)m_x.size());
	forcesAndDerivatives(m_x, m_uv, m_v, forces, forces, m_dfdx, m_dfdx, m_dfdv);
	for (size_t i = 0; i < forceFields.size(); ++i)
	{
		forceFields[i]->forcesAndDerivatives(forces, m_dfdx);
	}

	Vector implicitUpdateRHS((int)m_x.size());
	assembleImplicitUpdateEquations(dt, forces, m_v, m_dfdx, m_dfdv, m_implicitUpdateMatrix, implicitUpdateRHS);

	// solve the linear system:
	Vector dv((int)m_x.size());
	solver.solve(m_implicitUpdateMatrix, implicitUpdateRHS, dv);

	// update:
	m_v += dv;
	m_x += dt * m_v;
}

template<class Real>
void ClothMesh<Real>::forcesAndDerivatives(const Vector &x, const Vector &uv, const Vector &v, Vector &f, Vector &d, SparseMatrix &dfdx, SparseMatrix &dddx, SparseMatrix &dddv) const
{
	d.setConstant(0);
	f.setConstant(0);
	for (int i = 0; i < dfdx.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(dfdx, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}
	for (int i = 0; i < dddx.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(dddx, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}
	for (int i = 0; i < dddv.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(dddv, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}

	// add the bend terms:
	for (size_t i = 0; i < m_bendConditions.size(); ++i)
	{
		m_bendConditions[i].computeForces(x, uv, m_kBend, f, dfdx, v, m_dBend, d, dddx, dddv);
	}

	// add the stretch terms:
	for (size_t i = 0; i < m_stretchConditions.size(); ++i)
	{
		m_stretchConditions[i].computeForces(x, uv, m_kStretch, f, dfdx, v, m_dStretch, d, dddx, dddv);
	}

	// add the shear terms:
	for (size_t i = 0; i < m_shearConditions.size(); ++i)
	{
		m_shearConditions[i].computeForces(x, uv, m_kShear, f, dfdx, v, m_dShear, d, dddx, dddv);
	}
}


template<class Real>
void ClothMesh<Real>::assembleImplicitUpdateEquations(Real dt, const Vector &f, const Vector &v, const SparseMatrix &dfdx, const SparseMatrix &dfdv, SparseMatrix &implicitUpdateMatrix, Vector &implicitUpdateRHS) const
{

	// reset existing sim data to zero:
	for (int i = 0; i < implicitUpdateMatrix.outerSize(); ++i)
	{
		for (SparseMatrix::InnerIterator it(implicitUpdateMatrix, i); it; ++it)
		{
			it.valueRef() = 0;
		}
	}

	// right hand side:
	implicitUpdateRHS = f * dt + dfdx * v * dt * dt;

	// matrix:
	implicitUpdateMatrix = -dfdx * dt * dt - dfdv * dt;

	// add masses onto the diagonal:
	for (int i = 0; i < f.size(); ++i)
	{
		// unit masses for now...
		implicitUpdateMatrix.coeffRef(i, i) += m_m[i / 3];
	}
	implicitUpdateMatrix.makeCompressed();
}

template<class Real>
Real ClothMesh<Real>::energy(const Vector& x, const Vector& uv) const
{
	Real e = 0;
	for (size_t i = 0; i < m_bendConditions.size(); ++i)
	{
		Vector c = m_bendConditions[i].C(x, uv);
		e += m_kBend * Real(0.5) * c.dot(c);
	}
	for (size_t i = 0; i < m_shearConditions.size(); ++i)
	{
		Vector c = m_shearConditions[i].C(x, uv);
		e += m_kShear * Real(0.5) * c.dot(c);
	}
	for (size_t i = 0; i < m_stretchConditions.size(); ++i)
	{
		Vector c = m_stretchConditions[i].C(x, uv);
		e += m_kStretch * Real(0.5) * c.dot(c);
	}
	return e;
}

template<class Real>
void ClothMesh<Real>::C(const Vector& x, const Vector& uv, Vector& c) const
{
	c.resize(m_bendConditions.size() + m_shearConditions.size() + 2 * m_stretchConditions.size());
	c.setConstant(0);

	size_t n = 0;
	for (size_t i = 0; i < m_bendConditions.size(); ++i)
	{
		c[n++] = m_bendConditions[i].C(x, uv)[0];
	}
	for (size_t i = 0; i < m_shearConditions.size(); ++i)
	{
		c[n++] = m_shearConditions[i].C(x, uv)[0];
	}
	for (size_t i = 0; i < m_stretchConditions.size(); ++i)
	{
		Vector cc = m_stretchConditions[i].C(x, uv);
		c[n++] = cc[0];
		c[n++] = cc[1];
	}
}

template class ClothMesh<float>;
template class ClothMesh<double>;
