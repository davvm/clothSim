#ifndef BENDCONDITION_H
#define BENDCONDITION_H

#include "EnergyCondition.h"

template <class Real>
class BendCondition : public EnergyCondition<Real>
{

public:

	BendCondition(int i0, int i1, int i2, int i3);
	virtual ~BendCondition() {}

	virtual Vector C(const Vector& x, const Vector& uv) const;

	virtual void computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const;

	struct TriangleQuantities
	{
		TriangleQuantities();
		TriangleQuantities(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, const Vector3 &v3);

		// edge vectors:
		Vector3 v01, v02, v32, v31, v21;

		// normalized edge vectors:
		Vector3 e01, e02, e32, e31, e21;

		// cosines of the angles at each of the points:
		Real c00, c01, c02, c13, c11, c12;

		// normalized triangle normals:
		Vector3 n0, n1;

		// normalized binormals for each triangle, pointing from a vertex to its opposite edge:
		Vector3 b00, b01, b02, b13, b12, b11;

		// vertex distances to opposite edges:
		Real d00, d01, d02, d11, d12, d13;

		// compute angle between triangles:
		Real theta;

		// derivatives of theta with respect to the different vertex positions:
		Vector3 dThetadP0, dThetadP1, dThetadP2, dThetadP3;

		// time derivative of angle:
		Real dThetadt;

		// derivatives of the normals:
		Matrix3 dn0dP0, dn0dP1, dn0dP2, dn0dP3, dn1dP0, dn1dP1, dn1dP2, dn1dP3;

		// derivatives of the cosines:
		Vector3 dc01dP0, dc01dP1, dc01dP2, dc01dP3, dc02dP0, dc02dP1, dc02dP2, dc02dP3, dc11dP0, dc11dP1, dc11dP2, dc11dP3, dc12dP0, dc12dP1, dc12dP2, dc12dP3;

		// derivatives of the perpendicular distances:
		Vector3 dd00dP0, dd00dP1, dd00dP2, dd00dP3, dd01dP0, dd01dP1, dd01dP2, dd01dP3, dd02dP0, dd02dP1, dd02dP2, dd02dP3, dd11dP0, dd11dP1, dd11dP2, dd11dP3, dd12dP0, dd12dP1, dd12dP2, dd12dP3, dd13dP0, dd13dP1, dd13dP2, dd13dP3;

		// second derivatives of theta with respect to the different vertex positions:
		Matrix3 d2ThetadP0dP0, d2ThetadP0dP1, d2ThetadP0dP2, d2ThetadP0dP3;
		Matrix3 d2ThetadP1dP0, d2ThetadP1dP1, d2ThetadP1dP2, d2ThetadP1dP3;
		Matrix3 d2ThetadP2dP0, d2ThetadP2dP1, d2ThetadP2dP2, d2ThetadP2dP3;
		Matrix3 d2ThetadP3dP0, d2ThetadP3dP1, d2ThetadP3dP2, d2ThetadP3dP3;
	};

private:

	// bend condition is attached to two triangles, ie four vertices like so:
	//
	//     0
	//    / \
	//   /   \
	//  1-----2
	//   \   /
	//    \ /
	//     3
	//
	int m_inds[4];

};


#endif // BENDCONDITION_H
