#ifndef SHEARCONDITION_H
#define SHEARCONDITION_H

#include "EnergyCondition.h"
#include "TangentTriangleQuantities.h"

template <class Real>
class ShearCondition : public EnergyCondition<Real>
{

public:

	ShearCondition(int i0, int i1, int i2);
	virtual ~ShearCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const;

	struct TriangleQuantities : public TangentTriangleQuantities<Real>
	{
		TriangleQuantities();
		TriangleQuantities(
			const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
			const Vector2 &uv0, const Vector2 &uv1, const Vector2 &uv2,
			const Vector3 &v0, const Vector3 &v1, const Vector3 &v2
		);

		// normalized tangents:
		Vector3 wuHat, wvHat;

		// energy condition:
		Real C;

		// time derivative of energy condition:
		Real dCdt;

		// derivatives of the energy condition:
		Vector3 dCdP0, dCdP1, dCdP2;

		// second derivatives of the energy condition:
		Matrix3 d2CdP0dP0, d2CdP0dP1, d2CdP0dP2;
		Matrix3 d2CdP1dP0, d2CdP1dP1, d2CdP1dP2;
		Matrix3 d2CdP2dP0, d2CdP2dP1, d2CdP2dP2;

	};

private:

	// indices of the vertices this shear condition is attached to:
	int m_inds[3];

};


#endif // SHEARCONDITION_H
