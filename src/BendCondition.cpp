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
BendCondition<Real>::TriangleQuantities::TriangleQuantities()
{
}

template<class Real>
BendCondition<Real>::TriangleQuantities::TriangleQuantities(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3)
{
	// triangles are laid out like this:
	//     0
	//    / \
	//   / 0 \
	//  1-----2
	//   \ 1 /
	//    \ /
	//     3

	// edge vectors:
	v01 = p1 - p0;
	v02 = p2 - p0;
	v32 = p2 - p3;
	v31 = p1 - p3;
	v21 = p1 - p2;

	// normalized edge vectors:
	e01 = v01.normalized();
	e02 = v02.normalized();
	e32 = v32.normalized();
	e31 = v31.normalized();
	e21 = v21.normalized();

	// cosines of the angles at each of the points:
	c00 = e01.dot(e02);
	c01 = e01.dot(e21);
	c02 = -e02.dot(e21);

	c13 = e31.dot(e32);
	c11 = e21.dot(e31);
	c12 = -e32.dot(e21);

	// normalized triangle normals:
	n0 = (e21.cross(e01)).normalized();
	n1 = (e32.cross(e21)).normalized();

	// normalized binormals for each triangle, pointing from a vertex to its opposite edge:
	b00 = (e01 - e21 * (e21.dot(e01))).normalized();
	b01 = (-e01 - e02 * (e02.dot(-e01))).normalized();
	b02 = (-e02 - e01 * (e01.dot(-e02))).normalized();

	b13 = (e32 - e21 * (e21.dot(e32))).normalized();
	b12 = (-e32 - e31 * (e31.dot(-e32))).normalized();
	b11 = (-e31 - e32 * (e32.dot(-e31))).normalized();

	// vertex distances to opposite edges:
	d00 = (v01 - v01.dot(v21) / v21.dot(v21) * v21).norm();//b00.dot(v01);
	d01 = b01.dot(-v21);
	d02 = b02.dot(-v02);

	d11 = b11.dot(-v31);
	d12 = b12.dot(v21);
	d13 = b13.dot(v31);

	// compute angle between triangles, using a sin as that can give us negative answers as well as positive ones:
	Real cosTheta = n0.dot(n1);
	Real sinTheta = n1.dot(b00);
	theta = atan2(sinTheta, cosTheta);

	// derivatives of theta with respect to the different vertex positions:
	dThetadP0 = -n0 / d00;
	dThetadP1 = c02 * n0 / d01 + c12 * n1 / d11;
	dThetadP2 = c01 * n0 / d02 + c11 * n1 / d12;
	dThetadP3 = -n1 / d13;

	// now compute the derivatives of some terms in those formulae, so we can get second derivatives:

	// derivatives of the normals:
	dn0dP0 = b00 * n0.transpose() / d00;
	dn0dP1 = b01 * n0.transpose() / d01;
	dn0dP2 = b02 * n0.transpose() / d02;
	dn0dP3 = Matrix3::Zero();

	dn1dP0 = Matrix3::Zero();
	dn1dP1 = b11 * n1.transpose() / d11;
	dn1dP2 = b12 * n1.transpose() / d12;
	dn1dP3 = b13 * n1.transpose() / d13;

	// derivatives of the cosines:
	dc01dP0 = -b02 * b00.dot(v01) / v01.dot(v01);
	dc01dP2 = -b00 * b02.dot(v21) / v21.dot(v21);
	dc01dP1 = -dc01dP0 - dc01dP2;
	dc01dP3 = Vector3::Zero();

	dc02dP0 = -b01 * b00.dot(v02) / v02.dot(v02);
	dc02dP1 = b00 * b01.dot(v21) / v21.dot(v21);
	dc02dP2 = -dc02dP0 - dc02dP1;
	dc02dP3 = Vector3::Zero();

	dc11dP0 = Vector3::Zero();
	dc11dP2 = -b13 * b12.dot(v21) / v21.dot(v21);
	dc11dP3 = -b12 * b13.dot(v31) / v31.dot(v31);
	dc11dP1 = -dc11dP2 - dc11dP3;

	dc12dP0 = Vector3::Zero();
	dc12dP1 = b13 * b11.dot(v21) / v21.dot(v21);
	dc12dP3 = -b11 * b13.dot(v32) / v32.dot(v32);
	dc12dP2 = -dc12dP1 - dc12dP3;

	// derivatives of the perpendicular distances:
	dd00dP0 = -b00;
	dd00dP1 = b00 * -v21.dot(v02) / v21.dot(v21);
	dd00dP2 = b00 * v21.dot(v01) / v21.dot(v21);
	dd00dP3 = Vector3::Zero();

	dd01dP0 = b01 * v02.dot(-v21) / v02.dot(v02);
	dd01dP1 = -b01;
	dd01dP2 = b01 * v02.dot(v01) / v02.dot(v02);
	dd01dP3 = Vector3::Zero();

	dd02dP0 = b02 * v01.dot(v21) / v01.dot(v01);
	dd02dP1 = b02 * v01.dot(v02) / v01.dot(v01);
	dd02dP2 = -b02;
	dd02dP3 = Vector3::Zero();

	dd11dP0 = Vector3::Zero();
	dd11dP1 = -b11;
	dd11dP2 = b11 * v32.dot(v31) / v32.dot(v32);
	dd11dP3 = b11 * v32.dot(-v21) / v32.dot(v32);

	dd12dP0 = Vector3::Zero();
	dd12dP1 = b12 * v31.dot(v32) / v31.dot(v31);
	dd12dP2 = -b12;
	dd12dP3 = b12 * v31.dot(v21) / v31.dot(v31);

	dd13dP0 = Vector3::Zero();
	dd13dP1 = b13 * v21.dot(-v32) / v21.dot(v21);
	dd13dP2 = b13 * v21.dot(v31) / v21.dot(v21);
	dd13dP3 = -b13;
}

template<class Real>
typename EnergyCondition<Real>::Vector BendCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	TriangleQuantities q(x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]), x.segment<3>(3 * m_inds[3]));

	Vector ret(1);
	ret[0] = q.theta;

	return ret;
}

template<class Real>
void BendCondition<Real>::computeForces(const Vector& x, const Vector& uv, Vector& forces) const
{
	
	TriangleQuantities q(x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]), x.segment<3>(3 * m_inds[3]));

	// Compute forces:

	// E = 1/2 C^t C
	// dE/dx = theta * dThetadX
	// f = -dE/dx
	// f = - theta * dThetadX

	forces.segment<3>(3 * m_inds[0]) -= q.theta * q.dThetadP0;
	forces.segment<3>(3 * m_inds[1]) -= q.theta * q.dThetadP1;
	forces.segment<3>(3 * m_inds[2]) -= q.theta * q.dThetadP2;
	forces.segment<3>(3 * m_inds[3]) -= q.theta * q.dThetadP3;
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
