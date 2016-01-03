#ifndef TANGENTTRIANGLEQUANTITIES_H
#define TANGENTTRIANGLEQUANTITIES_H

#include "EnergyCondition.h"

template<class Real>
struct TangentTriangleQuantities
{
	typedef typename EnergyCondition<Real>::Vector2 Vector2;
	typedef typename EnergyCondition<Real>::Vector3 Vector3;
	typedef typename EnergyCondition<Real>::Matrix3 Matrix3;

	// area:
	Real a;

	// tangent vectors:
	Vector3 wu, wv;

	// derivatives of the tangent vectors:
	Matrix3 dwudP0, dwudP1, dwudP2;
	Matrix3 dwvdP0, dwvdP1, dwvdP2;

	Real dwudP0Scalar, dwudP1Scalar, dwudP2Scalar;
	Real dwvdP0Scalar, dwvdP1Scalar, dwvdP2Scalar;

	TangentTriangleQuantities() {}
	TangentTriangleQuantities(
		const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
		const Vector2 &uv0, const Vector2 &uv1, const Vector2 &uv2
		)
	{
		Vector2 duv1 = uv1 - uv0;
		Vector2 duv2 = uv2 - uv0;

		Real du1 = duv1[0];
		Real dv1 = duv1[1];
		Real du2 = duv2[0];
		Real dv2 = duv2[1];

		// triangle area in reference pose:
		a = Real(0.5) * (du1 * dv2 - du2 * dv1);

		// trangle tangents in reference directions:
		wu = ((p1 - p0) * dv2 - (p2 - p0) * dv1) / (2 * a);
		wv = (-(p1 - p0) * du2 + (p2 - p0) * du1) / (2 * a);

		// first derivatives of uv tangents:
		dwudP0Scalar = (dv1 - dv2) / (2 * a);
		dwudP1Scalar = dv2 / (2 * a);
		dwudP2Scalar = -dv1 / (2 * a);

		dwvdP0Scalar = (du2 - du1) / (2 * a);
		dwvdP1Scalar = -du2 / (2 * a);
		dwvdP2Scalar = du1 / (2 * a);


		dwudP0 = Matrix3::Identity() * dwudP0Scalar;
		dwudP1 = Matrix3::Identity() * dwudP1Scalar;
		dwudP2 = Matrix3::Identity() * dwudP2Scalar;

		dwvdP0 = Matrix3::Identity() * dwvdP0Scalar;
		dwvdP1 = Matrix3::Identity() * dwvdP1Scalar;
		dwvdP2 = Matrix3::Identity() * dwvdP2Scalar;

	}
};

#endif
