#include "simLib\ConstrainedCGSolver.h"

using namespace ClothSim;

template<class Real>
ConstrainedCGSolver<Real>::ConstrainedCGSolver(
	const std::vector<int> &constraintIndices,
	const std::vector< Matrix3 > &constraintMatrices,
	const Vector &constraintVelocityDeltas,
	Real tol,
	int maxIterations
) :
	m_constraintIndices(constraintIndices),
	m_constraintMatrices(constraintMatrices),
	m_constraintVelocityDeltas(constraintVelocityDeltas),
	m_tol(tol),
	m_maxIterations(maxIterations)
{
}

template<class Real>
void ConstrainedCGSolver<Real>::solve(const SparseMatrix &A, const Vector &b, Vector &deltaV) const
{
	//\todo: add preconditioner
	Vector P(b.size());
	for (int i = 0; i < b.size(); ++i)
	{
		P[i] = A.coeff(i,i);
	}

	deltaV = m_constraintVelocityDeltas;

	Real delta0;
	{
		Vector filterB = b * Real(1.0);
		filter(filterB);
		Vector filterBPrecond = filterB * Real(1.0);
		precondition(filterBPrecond, P, false);
		delta0 = filterB.dot(filterBPrecond);
	}

	Vector r = b - A * deltaV;
	filter(r);
	Vector c = r * Real(1.0);
	precondition(c, P, true);
	filter(c);

	Real deltaNew = r.dot(c);
	
	Real deltaOld;
	Vector q, s;

	for (int i = 0; i < m_maxIterations && deltaNew > m_tol * m_tol * delta0; ++i)
	{
		q = A * c;
		filter(q);

		Real alpha = deltaNew / c.dot(q);
		deltaV += alpha * c;
		r -= alpha * q;
		s = r * Real(1.0);
		precondition(s, P, true);

		deltaOld = deltaNew;
		deltaNew = r.dot(s);

		c = s + (deltaNew / deltaOld) * c;
		filter(c);
	}
}

template<class Real>
void ConstrainedCGSolver<Real>::filter(Vector &x) const
{
	std::vector<int>::const_iterator it = m_constraintIndices.begin();
	std::vector<Matrix3>::const_iterator matrixIt = m_constraintMatrices.begin();
	for (; it != m_constraintIndices.end(); ++it, ++matrixIt)
	{
		int idx = *it;
		x.segment<3>(3 * idx) = (*matrixIt) * x.segment<3>(3 * idx);
	}
}

template<class Real>
void ConstrainedCGSolver<Real>::precondition(Vector &x, Vector &p, bool inverse) const
{
	if (inverse)
	{
		x.array() /= p.array();
	}
	else
	{
		x.array() *= p.array();
	}
}


template class ConstrainedCGSolver<float>;
template class ConstrainedCGSolver<double>;
