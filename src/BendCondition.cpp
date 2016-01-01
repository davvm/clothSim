#include "BendCondition.h"
#include <iostream>
template<class Real>
BendCondition<Real>::BendCondition(int i0, int i1, int i2, int i3)
{
	m_inds[0] = i0;
	m_inds[1] = i1;
	m_inds[2] = i2;
	m_inds[3] = i3;
}

template<class Real>
typename EnergyCondition<Real>::Vector BendCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	// triangles are laid out like this:
	//     0
	//    / \
	//   /   \
	//  1-----2
	//   \   /
	//    \ /
	//     3

	const Vector3 &p0 = x.segment<3>(3 * m_inds[0]);
	const Vector3 &p1 = x.segment<3>(3 * m_inds[1]);
	const Vector3 &p2 = x.segment<3>(3 * m_inds[2]);
	const Vector3 &p3 = x.segment<3>(3 * m_inds[3]);

	Vector3 v01 = p1 - p0;
	Vector3 v32 = p2 - p3;
	Vector3 v21 = p1 - p2;

	// normalized vector along shared edge:
	Vector3 e = -v21.normalized();

	// compute normalized triangle normals:
	Vector3 n0 = (v21.cross(v01)).normalized();
	Vector3 n1 = (v32.cross(v21)).normalized();
	Vector3 b0 = (v01 - e * (e.dot(v01))).normalized();

	Real cosTheta = n0.dot(n1);
	Real sinTheta = n1.dot(b0);
	Real theta = atan2(sinTheta, cosTheta);

	Vector ret(1);
	ret[0] = theta;

	return ret;
}

template<class Real>
void BendCondition<Real>::computeForces(const Vector& x, const Vector& uv, Vector& forces) const
{
	// triangles are laid out like this:
	//     0
	//    / \
	//   / 0 \
	//  1-----2
	//   \ 1 /
	//    \ /
	//     3

	const Vector3 &p0 = x.segment<3>(3 * m_inds[0]);
	const Vector3 &p1 = x.segment<3>(3 * m_inds[1]);
	const Vector3 &p2 = x.segment<3>(3 * m_inds[2]);
	const Vector3 &p3 = x.segment<3>(3 * m_inds[3]);

	// edge vectors:
	Vector3 v01 = p1 - p0;
	Vector3 v02 = p2 - p0;
	Vector3 v32 = p2 - p3;
	Vector3 v31 = p1 - p3;
	Vector3 v21 = p1 - p2;

	// normalized edge vectors:
	Vector3 e01 = v01.normalized();
	Vector3 e02 = v02.normalized();
	Vector3 e32 = v32.normalized();
	Vector3 e31 = v31.normalized();
	Vector3 e21 = v21.normalized();

	// normalized triangle normals:
	Vector3 n0 = (e21.cross(e01)).normalized();
	Vector3 n1 = (e32.cross(e21)).normalized();

	// normalized binormals for each triangle, pointing from a vertex to its opposite edge:
	Vector3 b00 = (e01 - e21 * (e21.dot(e01))).normalized();
	Vector3 b01 = (-e01 - e02 * (e02.dot(-e01))).normalized();
	Vector3 b02 = (-e02 - e01 * (e01.dot(-e02))).normalized();

	Vector3 b13 = (e32 - e21 * (e21.dot(e32))).normalized();
	Vector3 b12 = (-e32 - e31 * (e31.dot(-e32))).normalized();
	Vector3 b11 = (-e31 - e32 * (e32.dot(-e31))).normalized();

	// compute angle between triangles, using a sin as that can give us negative answers as well as positive ones:
	Real cosTheta = n0.dot(n1);
	Real sinTheta = n1.dot(b00);
	Real theta = atan2(sinTheta,cosTheta);

	// E = 1/2 * theta^2
	// dE/dx = dtheta/dx * theta

	// now we can fill in the derivatives for dE/dp0 and dE/dp3, which are just their
	// respective normals divided by the perpendicular distances times the angle:
	forces.segment<3>(3 * m_inds[0]) += (theta / b00.dot(v01) ) * n0;
	forces.segment<3>(3 * m_inds[3]) += (theta / b13.dot(v31) ) * n1;

	// d/dP1 cos theta = d/dP1(n0.n1) = -( sin theta ) dtheta/dP1
	Real invSinc;
	if (abs(theta) < Real(1.e-4))
	{
		invSinc = Real(1.0) / (Real(1.0) - theta * theta / Real(6.0));
	}
	else
	{
		invSinc = theta / sinTheta;
	}

	forces.segment<3>(3 * m_inds[1]) += invSinc * (n1.dot( b01 ) * n0 / (p2 - p1).dot(b01) + n0.dot( b11 ) * n1 / (p3 - p1).dot(b11));
	forces.segment<3>(3 * m_inds[2]) += invSinc * (n1.dot( b02 ) * n0 / (p0 - p2).dot(b02) + n0.dot( b12 ) * n1 / (p1 - p2).dot(b12));
}

template<class Real>
void BendCondition<Real>::computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{
}

template<class Real>
void BendCondition<Real>::computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const
{
}

template<class Real>
void BendCondition<Real>::computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const
{
}

template class BendCondition<float>;
template class BendCondition<double>;
