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

	virtual void computeForces(const Vector& x, const Vector& uv, Vector& forces) const;
	virtual void computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

	virtual void computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const;
	virtual void computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const;

	typedef Eigen::Matrix<Real, 3, 3, 0, 3, 3> Matrix3;
	static Matrix3 dndP(const Vector3& p0, const Vector3& p1, const Vector3& p2);
	static Vector3 dNormalCrossProductdP(const Vector3& p, const Vector3& pOpposite, const Vector3& p0, const Vector3& p1);

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
