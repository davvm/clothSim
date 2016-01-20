#include "simLib\ShearCondition.h"

using namespace ClothSim;

template <class Real>
ShearCondition<Real>::ShearCondition(int i0, int i1, int i2)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
}

template<class Real>
ShearCondition<Real>::TriangleQuantities::TriangleQuantities() : TangentTriangleQuantities()
{
}

template<class Real>
ShearCondition<Real>::TriangleQuantities::TriangleQuantities(
	const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
	const Vector2 &uv0, const Vector2 &uv1, const Vector2 &uv2,
	const Vector3 &v0, const Vector3 &v1, const Vector3 &v2
	) : TangentTriangleQuantities<Real>(p0, p1, p2, uv0, uv1, uv2)
{

	Real wuNorm = wu.norm();
	Real wvNorm = wv.norm();

	wuHat = wu / wuNorm;
	wvHat = wv / wvNorm;

	// energy condition:
	C = a * wu.dot(wv);

	// first derivatives of condition quantity:
	dCdP0 = a * (dwudP0Scalar * wv + dwvdP0Scalar * wu);
	dCdP1 = a * (dwudP1Scalar * wv + dwvdP1Scalar * wu);
	dCdP2 = a * (dwudP2Scalar * wv + dwvdP2Scalar * wu);

	// d/dP1  a * (dwudP0Scalar * wv + dwvdP0Scalar * wu)
	// == a * (dwudP0Scalar * dwvdP1 + dwvdP0Scalar * dwudP1)

	// time derivative of energy condition:
	dCdt = dCdP0.dot(v0) + dCdP1.dot(v1) + dCdP2.dot(v2);
	
	d2CdP0dP0 = 2 * a * dwudP0Scalar * dwvdP0Scalar * Matrix3::Identity();
	d2CdP0dP1 = a * (dwudP0Scalar * dwvdP1Scalar + dwvdP0Scalar * dwudP1Scalar) * Matrix3::Identity();
	d2CdP0dP2 = a * (dwudP0Scalar * dwvdP2Scalar + dwvdP0Scalar * dwudP2Scalar) * Matrix3::Identity();

	d2CdP1dP0 = a * (dwudP1Scalar * dwvdP0Scalar + dwvdP1Scalar * dwudP0Scalar) * Matrix3::Identity();
	d2CdP1dP1 = 2 * a * dwvdP1Scalar * dwudP1Scalar * Matrix3::Identity();
	d2CdP1dP2 = a * (dwudP1Scalar * dwvdP2Scalar + dwvdP1Scalar * dwudP2Scalar) * Matrix3::Identity();

	d2CdP2dP0 = a * (dwudP2Scalar * dwvdP0Scalar + dwvdP2Scalar * dwudP0Scalar) * Matrix3::Identity();
	d2CdP2dP1 = a * (dwudP2Scalar * dwvdP1Scalar + dwvdP2Scalar * dwudP1Scalar) * Matrix3::Identity();
	d2CdP2dP2 = 2 * a * dwvdP2Scalar * dwudP2Scalar * Matrix3::Identity();
	
}

template <class Real>
typename EnergyCondition<Real>::Vector ShearCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]),
		uv.segment<2>(2 * m_inds[0]), uv.segment<2>(2 * m_inds[1]), uv.segment<2>(2 * m_inds[2]),
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2])
		);

	Vector ret(1);
	ret[0] = q.C;
	return ret;
}

template <class Real>
void ShearCondition<Real>::computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoXDerivatives, SparseMatrix &dddv) const
{
	// I can't be bothered doing C = a * wuHat . wvHat, so I'm doing C = a * wu * wv instead.
	// I guess this term will have the effect of trying to shrink the cloth. Hopefully that'll
	// get counteracted by the stretch condition.

	// Actually, it won't try and shrink it at all if the tangents are orthogonal will it...
	// probably won't be a promlem then.

	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]),
		uv.segment<2>(2 * m_inds[0]), uv.segment<2>(2 * m_inds[1]), uv.segment<2>(2 * m_inds[2]),
		v.segment<3>(3 * m_inds[0]), v.segment<3>(3 * m_inds[1]), v.segment<3>(3 * m_inds[2])
		);

	// E = 0.5 * C * C
	// f = -dE/dx = C dC/dx
	forces.segment<3>(3 * m_inds[0]) -= k * q.C * q.dCdP0;
	forces.segment<3>(3 * m_inds[1]) -= k * q.C * q.dCdP1;
	forces.segment<3>(3 * m_inds[2]) -= k * q.C * q.dCdP2;


	// compute force derivatives and insert them into the sparse matrix:
	Matrix3 df0dP0 = -k * (q.dCdP0 * q.dCdP0.transpose() + q.C * q.d2CdP0dP0);
	Matrix3 df0dP1 = -k * (q.dCdP0 * q.dCdP1.transpose() + q.C * q.d2CdP0dP1);
	Matrix3 df0dP2 = -k * (q.dCdP0 * q.dCdP2.transpose() + q.C * q.d2CdP0dP2);

	Matrix3 df1dP0 = -k * (q.dCdP1 * q.dCdP0.transpose() + q.C * q.d2CdP1dP0);
	Matrix3 df1dP1 = -k * (q.dCdP1 * q.dCdP1.transpose() + q.C * q.d2CdP1dP1);
	Matrix3 df1dP2 = -k * (q.dCdP1 * q.dCdP2.transpose() + q.C * q.d2CdP1dP2);

	Matrix3 df2dP0 = -k * (q.dCdP2 * q.dCdP0.transpose() + q.C * q.d2CdP2dP0);
	Matrix3 df2dP1 = -k * (q.dCdP2 * q.dCdP1.transpose() + q.C * q.d2CdP2dP1);
	Matrix3 df2dP2 = -k * (q.dCdP2 * q.dCdP2.transpose() + q.C * q.d2CdP2dP2);


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

	// fd = -d * dC/dt * dC/dx:

	dampingForces.segment<3>(3 * m_inds[0]) -= d * q.dCdt * q.dCdP0;
	dampingForces.segment<3>(3 * m_inds[1]) -= d * q.dCdt * q.dCdP1;
	dampingForces.segment<3>(3 * m_inds[2]) -= d * q.dCdt * q.dCdP2;

	Matrix3 dfd0dV0 = -d * (q.dCdP0 * q.dCdP0.transpose());
	Matrix3 dfd0dV1 = -d * (q.dCdP0 * q.dCdP1.transpose());
	Matrix3 dfd0dV2 = -d * (q.dCdP0 * q.dCdP2.transpose());

	Matrix3 dfd1dV0 = -d * (q.dCdP1 * q.dCdP0.transpose());
	Matrix3 dfd1dV1 = -d * (q.dCdP1 * q.dCdP1.transpose());
	Matrix3 dfd1dV2 = -d * (q.dCdP1 * q.dCdP2.transpose());

	Matrix3 dfd2dV0 = -d * (q.dCdP2 * q.dCdP0.transpose());
	Matrix3 dfd2dV1 = -d * (q.dCdP2 * q.dCdP1.transpose());
	Matrix3 dfd2dV2 = -d * (q.dCdP2 * q.dCdP2.transpose());

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
	Matrix3 dD0dP0Pseudo = -d * (q.d2CdP0dP0 * q.dCdt);
	Matrix3 dD1dP0Pseudo = -d * (q.d2CdP1dP0 * q.dCdt);
	Matrix3 dD2dP0Pseudo = -d * (q.d2CdP2dP0 * q.dCdt);

	Matrix3 dD0dP1Pseudo = -d * (q.d2CdP0dP1 * q.dCdt);
	Matrix3 dD1dP1Pseudo = -d * (q.d2CdP1dP1 * q.dCdt);
	Matrix3 dD2dP1Pseudo = -d * (q.d2CdP2dP1 * q.dCdt);

	Matrix3 dD0dP2Pseudo = -d * (q.d2CdP0dP2 * q.dCdt);
	Matrix3 dD1dP2Pseudo = -d * (q.d2CdP1dP2 * q.dCdt);
	Matrix3 dD2dP2Pseudo = -d * (q.d2CdP2dP2 * q.dCdt);

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

template class ShearCondition<float>;
template class ShearCondition<double>;
