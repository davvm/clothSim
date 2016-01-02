#ifndef STRETCHCONDITION_H
#define STRETCHCONDITION_H

#include "EnergyCondition.h"
#include "TangentTriangleQuantities.h"

template <class Real>
class StretchCondition : public EnergyCondition<Real>
{

public:

	StretchCondition(int i0, int i1, int i2, Real restU, Real restV);
	virtual ~StretchCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const;

	struct TriangleQuantities : public TangentTriangleQuantities<Real>
	{
		TriangleQuantities();
		TriangleQuantities(
			const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
			const Vector2 &uv0, const Vector2 &uv1, const Vector2 &uv2,
			const Vector3 &v0, const Vector3 &v1, const Vector3 &v2,
			Real bu, Real bv
		);

		// energy condition:
		Real C0, C1;

		// time derivatives of energy condition:
		Real dC0dt, dC1dt;

		// derivatives of the energy conditions:
		Vector3 dC0dP0, dC0dP1, dC0dP2;
		Vector3 dC1dP0, dC1dP1, dC1dP2;

		// second derivatives of the energy conditions:
		Matrix3 d2C0dP0dP0, d2C0dP0dP1, d2C0dP0dP2;
		Matrix3 d2C0dP1dP0, d2C0dP1dP1, d2C0dP1dP2;
		Matrix3 d2C0dP2dP0, d2C0dP2dP1, d2C0dP2dP2;

		Matrix3 d2C1dP0dP0, d2C1dP0dP1, d2C1dP0dP2;
		Matrix3 d2C1dP1dP0, d2C1dP1dP1, d2C1dP1dP2;
		Matrix3 d2C1dP2dP0, d2C1dP2dP1, d2C1dP2dP2;

	};

private:

	// indices of the vertices this stretch condition is attached to:
	int m_inds[3];

	// rest stretches:
	Real m_restU;
	Real m_restV;
};


#endif // STRETCHCONDITION_H
