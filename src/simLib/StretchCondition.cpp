#include "simLib\StretchCondition.h"

template <class Real>
StretchCondition<Real>::StretchCondition(int i0, int i1, int i2, Real restU, Real restV) : m_restU(restU), m_restV(restV)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
}

template<class Real>
StretchCondition<Real>::TriangleQuantities::TriangleQuantities() : TangentTriangleQuantities()
{
}

template<class Real>
StretchCondition<Real>::TriangleQuantities::TriangleQuantities(
	const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
	const Vector2 &uv0, const Vector2 &uv1, const Vector2 &uv2,
	const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
	Real bu, Real bv
) : TangentTriangleQuantities<Real>( p0, p1, p2, uv0, uv1, uv2 )
{

	// first derivatives of condition quantities:
	Real wuNorm = wu.norm();
	dC0dP0 = a * dwudP0 * wu / wuNorm;
	dC0dP1 = a * dwudP1 * wu / wuNorm;
	dC0dP2 = a * dwudP2 * wu / wuNorm;

	Real wvNorm = wv.norm();
	dC1dP0 = a * dwvdP0 * wv / wvNorm;
	dC1dP1 = a * dwvdP1 * wv / wvNorm;
	dC1dP2 = a * dwvdP2 * wv / wvNorm;

	// condition quantities:
	C0 = a * (wuNorm - bu);
	C1 = a * (wvNorm - bv);

	// time derivative of condition quantities:
	dC0dt = dC0dP0.dot(v0) + dC0dP1.dot(v1) + dC0dP2.dot(v2);
	dC1dt = dC1dP0.dot(v0) + dC1dP1.dot(v1) + dC1dP2.dot(v2);

	// second derivatives of condition quantities:
	Matrix3 wuMatrix = (Matrix3::Identity() - wu * wu.transpose() / (wuNorm * wuNorm));
	d2C0dP0dP0 = (a / wuNorm) * dwudP0 * dwudP0 * wuMatrix;
	d2C0dP0dP1 = (a / wuNorm) * dwudP0 * dwudP1 * wuMatrix;
	d2C0dP0dP2 = (a / wuNorm) * dwudP0 * dwudP2 * wuMatrix;

	d2C0dP1dP0 = (a / wuNorm) * dwudP1 * dwudP0 * wuMatrix;
	d2C0dP1dP1 = (a / wuNorm) * dwudP1 * dwudP1 * wuMatrix;
	d2C0dP1dP2 = (a / wuNorm) * dwudP1 * dwudP2 * wuMatrix;

	d2C0dP2dP0 = (a / wuNorm) * dwudP2 * dwudP0 * wuMatrix;
	d2C0dP2dP1 = (a / wuNorm) * dwudP2 * dwudP1 * wuMatrix;
	d2C0dP2dP2 = (a / wuNorm) * dwudP2 * dwudP2 * wuMatrix;

	Matrix3 wvMatrix = (Matrix3::Identity() - wv * wv.transpose() / (wvNorm * wvNorm));
	d2C1dP0dP0 = (a / wvNorm) * dwvdP0 * dwvdP0 * wvMatrix;
	d2C1dP0dP1 = (a / wvNorm) * dwvdP0 * dwvdP1 * wvMatrix;
	d2C1dP0dP2 = (a / wvNorm) * dwvdP0 * dwvdP2 * wvMatrix;

	d2C1dP1dP0 = (a / wvNorm) * dwvdP1 * dwvdP0 * wvMatrix;
	d2C1dP1dP1 = (a / wvNorm) * dwvdP1 * dwvdP1 * wvMatrix;
	d2C1dP1dP2 = (a / wvNorm) * dwvdP1 * dwvdP2 * wvMatrix;

	d2C1dP2dP0 = (a / wvNorm) * dwvdP2 * dwvdP0 * wvMatrix;
	d2C1dP2dP1 = (a / wvNorm) * dwvdP2 * dwvdP1 * wvMatrix;
	d2C1dP2dP2 = (a / wvNorm) * dwvdP2 * dwvdP2 * wvMatrix;

}

template <class Real>
typename EnergyCondition<Real>::Vector StretchCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]),
		uv.segment<2>(2 * m_inds[0]), uv.segment<2>(2 * m_inds[1]), uv.segment<2>(2 * m_inds[2]),
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]),
		m_restU, m_restV
	);

	Vector ret(2);
	ret[0] = q.C0;
	ret[1] = q.C1;
	return ret;
}

template <class Real>
void StretchCondition<Real>::computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoXDerivatives, SparseMatrix &dddv) const
{
	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]),
		uv.segment<2>(2 * m_inds[0]), uv.segment<2>(2 * m_inds[1]), uv.segment<2>(2 * m_inds[2]),
		v.segment<3>(3 * m_inds[0]), v.segment<3>(3 * m_inds[1]), v.segment<3>(3 * m_inds[2]),
		m_restU, m_restV
	);

	// E = 0.5 * ( C0*C0 + C1*C1 )
	// f = -dE/dx = C0 dC0/dx + C1 dC1/dx
	forces.segment<3>(3 * m_inds[0]) -= k * (q.C0 * q.dC0dP0 + q.C1 * q.dC1dP0);
	forces.segment<3>(3 * m_inds[1]) -= k * (q.C0 * q.dC0dP1 + q.C1 * q.dC1dP1);
	forces.segment<3>(3 * m_inds[2]) -= k * (q.C0 * q.dC0dP2 + q.C1 * q.dC1dP2);

	// compute force derivatives and insert them into the sparse matrix:
	Matrix3 df0dP0 = -k * (q.dC0dP0 * q.dC0dP0.transpose() + q.C0 * q.d2C0dP0dP0 + q.dC1dP0 * q.dC1dP0.transpose() + q.C1 * q.d2C1dP0dP0);
	Matrix3 df0dP1 = -k * (q.dC0dP0 * q.dC0dP1.transpose() + q.C0 * q.d2C0dP0dP1 + q.dC1dP0 * q.dC1dP1.transpose() + q.C1 * q.d2C1dP0dP1);
	Matrix3 df0dP2 = -k * (q.dC0dP0 * q.dC0dP2.transpose() + q.C0 * q.d2C0dP0dP2 + q.dC1dP0 * q.dC1dP2.transpose() + q.C1 * q.d2C1dP0dP2);

	Matrix3 df1dP0 = -k * (q.dC0dP1 * q.dC0dP0.transpose() + q.C0 * q.d2C0dP1dP0 + q.dC1dP1 * q.dC1dP0.transpose() + q.C1 * q.d2C1dP1dP0);
	Matrix3 df1dP1 = -k * (q.dC0dP1 * q.dC0dP1.transpose() + q.C0 * q.d2C0dP1dP1 + q.dC1dP1 * q.dC1dP1.transpose() + q.C1 * q.d2C1dP1dP1);
	Matrix3 df1dP2 = -k * (q.dC0dP1 * q.dC0dP2.transpose() + q.C0 * q.d2C0dP1dP2 + q.dC1dP1 * q.dC1dP2.transpose() + q.C1 * q.d2C1dP1dP2);

	Matrix3 df2dP0 = -k * (q.dC0dP2 * q.dC0dP0.transpose() + q.C0 * q.d2C0dP2dP0 + q.dC1dP2 * q.dC1dP0.transpose() + q.C1 * q.d2C1dP2dP0);
	Matrix3 df2dP1 = -k * (q.dC0dP2 * q.dC0dP1.transpose() + q.C0 * q.d2C0dP2dP1 + q.dC1dP2 * q.dC1dP1.transpose() + q.C1 * q.d2C1dP2dP1);
	Matrix3 df2dP2 = -k * (q.dC0dP2 * q.dC0dP2.transpose() + q.C0 * q.d2C0dP2dP2 + q.dC1dP2 * q.dC1dP2.transpose() + q.C1 * q.d2C1dP2dP2);

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[0] + j) += df0dP0(i, j);
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[1] + j) += df0dP1(i, j);
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[2] + j) += df0dP2(i, j);

			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[0] + j) += df1dP0(i, j);
			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[1] + j) += df1dP1(i, j);
			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[2] + j) += df1dP2(i, j);

			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[0] + j) += df2dP0(i, j);
			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[1] + j) += df2dP1(i, j);
			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[2] + j) += df2dP2(i, j);
		}
	}

	// compute damping forces and v derivatives:

	// fd = -d * ( dC0/dt * dC0/dx + dC1/dt * dC1/dx ):

	dampingForces.segment<3>(3 * m_inds[0]) -= d * ( q.dC0dt * q.dC0dP0 + q.dC1dt * q.dC1dP0 );
	dampingForces.segment<3>(3 * m_inds[1]) -= d * ( q.dC0dt * q.dC0dP1 + q.dC1dt * q.dC1dP1 );
	dampingForces.segment<3>(3 * m_inds[2]) -= d * ( q.dC0dt * q.dC0dP2 + q.dC1dt * q.dC1dP2 );

	Matrix3 dfd0dV0 = -d * (q.dC0dP0 * q.dC0dP0.transpose() + q.dC1dP0 * q.dC1dP0.transpose());
	Matrix3 dfd0dV1 = -d * (q.dC0dP0 * q.dC0dP1.transpose() + q.dC1dP0 * q.dC1dP1.transpose());
	Matrix3 dfd0dV2 = -d * (q.dC0dP0 * q.dC0dP2.transpose() + q.dC1dP0 * q.dC1dP2.transpose());

	Matrix3 dfd1dV0 = -d * (q.dC0dP1 * q.dC0dP0.transpose() + q.dC1dP1 * q.dC1dP0.transpose());
	Matrix3 dfd1dV1 = -d * (q.dC0dP1 * q.dC0dP1.transpose() + q.dC1dP1 * q.dC1dP1.transpose());
	Matrix3 dfd1dV2 = -d * (q.dC0dP1 * q.dC0dP2.transpose() + q.dC1dP1 * q.dC1dP2.transpose());

	Matrix3 dfd2dV0 = -d * (q.dC0dP2 * q.dC0dP0.transpose() + q.dC1dP2 * q.dC1dP0.transpose());
	Matrix3 dfd2dV1 = -d * (q.dC0dP2 * q.dC0dP1.transpose() + q.dC1dP2 * q.dC1dP1.transpose());
	Matrix3 dfd2dV2 = -d * (q.dC0dP2 * q.dC0dP2.transpose() + q.dC1dP2 * q.dC1dP2.transpose());

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			dddv.coeffRef(3 * m_inds[0] + i, 3 * m_inds[0] + j) += dfd0dV0(i, j);
			dddv.coeffRef(3 * m_inds[0] + i, 3 * m_inds[1] + j) += dfd0dV1(i, j);
			dddv.coeffRef(3 * m_inds[0] + i, 3 * m_inds[2] + j) += dfd0dV2(i, j);

			dddv.coeffRef(3 * m_inds[1] + i, 3 * m_inds[0] + j) += dfd1dV0(i, j);
			dddv.coeffRef(3 * m_inds[1] + i, 3 * m_inds[1] + j) += dfd1dV1(i, j);
			dddv.coeffRef(3 * m_inds[1] + i, 3 * m_inds[2] + j) += dfd1dV2(i, j);

			dddv.coeffRef(3 * m_inds[2] + i, 3 * m_inds[0] + j) += dfd2dV0(i, j);
			dddv.coeffRef(3 * m_inds[2] + i, 3 * m_inds[1] + j) += dfd2dV1(i, j);
			dddv.coeffRef(3 * m_inds[2] + i, 3 * m_inds[2] + j) += dfd2dV2(i, j);
		}
	}

	// Damping force x derivatives are given by:

	// -kd ( dc/dpi pc/dpi^T + d2c/dpidpj dc/dt )

	// Unfortunately, the first term isn't symmetric, so it messes up the
	// linear solver. However, apparently it's quite small in practice so
	// you can leave it out, hence the name dampingForcePseudoXDerivatives.
	// see page 6 of http://www.cs.cmu.edu/~baraff/papers/sig98.pdf for
	// details.

	// Therefore, lets compute -kd * d2c/dpidpj dc/dt:
	Matrix3 dD0dP0Pseudo = -d * (q.d2C0dP0dP0 * q.dC0dt + q.d2C1dP0dP0 * q.dC1dt);
	Matrix3 dD1dP0Pseudo = -d * (q.d2C0dP1dP0 * q.dC0dt + q.d2C1dP1dP0 * q.dC1dt);
	Matrix3 dD2dP0Pseudo = -d * (q.d2C0dP2dP0 * q.dC0dt + q.d2C1dP2dP0 * q.dC1dt);

	Matrix3 dD0dP1Pseudo = -d * (q.d2C0dP0dP1 * q.dC0dt + q.d2C1dP0dP1 * q.dC1dt);
	Matrix3 dD1dP1Pseudo = -d * (q.d2C0dP1dP1 * q.dC0dt + q.d2C1dP1dP1 * q.dC1dt);
	Matrix3 dD2dP1Pseudo = -d * (q.d2C0dP2dP1 * q.dC0dt + q.d2C1dP2dP1 * q.dC1dt);

	Matrix3 dD0dP2Pseudo = -d * (q.d2C0dP0dP2 * q.dC0dt + q.d2C1dP0dP2 * q.dC1dt);
	Matrix3 dD1dP2Pseudo = -d * (q.d2C0dP1dP2 * q.dC0dt + q.d2C1dP1dP2 * q.dC1dt);
	Matrix3 dD2dP2Pseudo = -d * (q.d2C0dP2dP2 * q.dC0dt + q.d2C1dP2dP2 * q.dC1dt);

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[0] + j) += dD0dP0Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[1] + j) += dD1dP0Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[2] + j) += dD2dP0Pseudo(i, j);

			dampingPseudoXDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[0] + j) += dD0dP1Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[1] + j) += dD1dP1Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[2] + j) += dD2dP1Pseudo(i, j);

			dampingPseudoXDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[0] + j) += dD0dP2Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[1] + j) += dD1dP2Pseudo(i, j);
			dampingPseudoXDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[2] + j) += dD2dP2Pseudo(i, j);
		}
	}


}

template class StretchCondition<float>;
template class StretchCondition<double>;
