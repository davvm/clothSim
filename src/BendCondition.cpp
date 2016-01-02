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
BendCondition<Real>::TriangleQuantities::TriangleQuantities(
	const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3,
	const Vector3 &v0, const Vector3 &v1, const Vector3 &v2, const Vector3 &v3
)
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

	// cosines of the angles at each of the points for each of the triangles:
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

	// time derivative of theta:
	dThetadt = dThetadP0.dot(v0) + dThetadP1.dot(v1) + dThetadP2.dot(v2) + dThetadP3.dot(v3);

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


	// second derivatives of theta:
	d2ThetadP0dP0 = -dn0dP0 / d00 + n0 * dd00dP0.transpose() / (d00 * d00);
	d2ThetadP0dP1 = -dn0dP1 / d00 + n0 * dd00dP1.transpose() / (d00 * d00);
	d2ThetadP0dP2 = -dn0dP2 / d00 + n0 * dd00dP2.transpose() / (d00 * d00);
	d2ThetadP0dP3 = -dn0dP3 / d00 + n0 * dd00dP3.transpose() / (d00 * d00);

	d2ThetadP1dP0 = ((c02 / d01) * dn0dP0 + n0 * (d01 * dc02dP0 - c02 * dd01dP0).transpose() / (d01 * d01)) + ((c12 / d11) * dn1dP0 + n1 * (d11 * dc12dP0 - c12 * dd11dP0).transpose() / (d11 * d11));
	d2ThetadP1dP1 = ((c02 / d01) * dn0dP1 + n0 * (d01 * dc02dP1 - c02 * dd01dP1).transpose() / (d01 * d01)) + ((c12 / d11) * dn1dP1 + n1 * (d11 * dc12dP1 - c12 * dd11dP1).transpose() / (d11 * d11));
	d2ThetadP1dP2 = ((c02 / d01) * dn0dP2 + n0 * (d01 * dc02dP2 - c02 * dd01dP2).transpose() / (d01 * d01)) + ((c12 / d11) * dn1dP2 + n1 * (d11 * dc12dP2 - c12 * dd11dP2).transpose() / (d11 * d11));
	d2ThetadP1dP3 = ((c02 / d01) * dn0dP3 + n0 * (d01 * dc02dP3 - c02 * dd01dP3).transpose() / (d01 * d01)) + ((c12 / d11) * dn1dP3 + n1 * (d11 * dc12dP3 - c12 * dd11dP3).transpose() / (d11 * d11));

	d2ThetadP2dP0 = ((c01 / d02) * dn0dP0 + n0 * (d02 * dc01dP0 - c01 * dd02dP0).transpose() / (d02 * d02)) + ((c11 / d12) * dn1dP0 + n1 * (d12 * dc11dP0 - c11 * dd12dP0).transpose() / (d12 * d12));
	d2ThetadP2dP1 = ((c01 / d02) * dn0dP1 + n0 * (d02 * dc01dP1 - c01 * dd02dP1).transpose() / (d02 * d02)) + ((c11 / d12) * dn1dP1 + n1 * (d12 * dc11dP1 - c11 * dd12dP1).transpose() / (d12 * d12));
	d2ThetadP2dP2 = ((c01 / d02) * dn0dP2 + n0 * (d02 * dc01dP2 - c01 * dd02dP2).transpose() / (d02 * d02)) + ((c11 / d12) * dn1dP2 + n1 * (d12 * dc11dP2 - c11 * dd12dP2).transpose() / (d12 * d12));
	d2ThetadP2dP3 = ((c01 / d02) * dn0dP3 + n0 * (d02 * dc01dP3 - c01 * dd02dP3).transpose() / (d02 * d02)) + ((c11 / d12) * dn1dP3 + n1 * (d12 * dc11dP3 - c11 * dd12dP3).transpose() / (d12 * d12));

	d2ThetadP3dP0 = -dn1dP0 / d13 + n1 * dd13dP0.transpose() / (d13 * d13);
	d2ThetadP3dP1 = -dn1dP1 / d13 + n1 * dd13dP1.transpose() / (d13 * d13);
	d2ThetadP3dP2 = -dn1dP2 / d13 + n1 * dd13dP2.transpose() / (d13 * d13);
	d2ThetadP3dP3 = -dn1dP3 / d13 + n1 * dd13dP3.transpose() / (d13 * d13);


}

template<class Real>
typename EnergyCondition<Real>::Vector BendCondition<Real>::C(const Vector& x, const Vector& uv) const
{
	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]), x.segment<3>(3 * m_inds[3]),
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]), x.segment<3>(3 * m_inds[3])
	);

	Vector ret(1);
	ret[0] = q.theta;

	return ret;
}

template<class Real>
void BendCondition<Real>::computeForces(
	const Vector& x,
	const Vector& uv,
	Real k,
	Vector& forces,
	SparseMatrix &dfdx,
	const Vector& v,
	Real d,
	Vector &dampingForces,
	SparseMatrix &dampingPseudoDerivatives
) const
{
	
	TriangleQuantities q(
		x.segment<3>(3 * m_inds[0]), x.segment<3>(3 * m_inds[1]), x.segment<3>(3 * m_inds[2]), x.segment<3>(3 * m_inds[3]),
		v.segment<3>(3 * m_inds[0]), v.segment<3>(3 * m_inds[1]), v.segment<3>(3 * m_inds[2]), v.segment<3>(3 * m_inds[3])
	);

	// Compute forces:

	// E = 1/2 C^t C
	// dE/dx = theta * dThetadX
	// f = -dE/dx
	// f = - theta * dThetadX

	forces.segment<3>(3 * m_inds[0]) -= k * q.theta * q.dThetadP0;
	forces.segment<3>(3 * m_inds[1]) -= k * q.theta * q.dThetadP1;
	forces.segment<3>(3 * m_inds[2]) -= k * q.theta * q.dThetadP2;
	forces.segment<3>(3 * m_inds[3]) -= k * q.theta * q.dThetadP3;

	// d/dP0( - q.theta * q.dThetadP0 )
	// = - q.dThetadP0 * q.dThetadP0.transpose() - q.theta * q.d2ThetadP0dP0

	// compute force derivatives and insert them into the sparse matrix:
	Matrix3 df0dP0 = -k * (q.dThetadP0 * q.dThetadP0.transpose() + q.theta * q.d2ThetadP0dP0);
	Matrix3 df0dP1 = -k * (q.dThetadP0 * q.dThetadP1.transpose() + q.theta * q.d2ThetadP0dP1);
	Matrix3 df0dP2 = -k * (q.dThetadP0 * q.dThetadP2.transpose() + q.theta * q.d2ThetadP0dP2);
	Matrix3 df0dP3 = -k * (q.dThetadP0 * q.dThetadP3.transpose() + q.theta * q.d2ThetadP0dP3);

	Matrix3 df1dP0 = -k * (q.dThetadP1 * q.dThetadP0.transpose() + q.theta * q.d2ThetadP1dP0);
	Matrix3 df1dP1 = -k * (q.dThetadP1 * q.dThetadP1.transpose() + q.theta * q.d2ThetadP1dP1);
	Matrix3 df1dP2 = -k * (q.dThetadP1 * q.dThetadP2.transpose() + q.theta * q.d2ThetadP1dP2);
	Matrix3 df1dP3 = -k * (q.dThetadP1 * q.dThetadP3.transpose() + q.theta * q.d2ThetadP1dP3);

	Matrix3 df2dP0 = -k * (q.dThetadP2 * q.dThetadP0.transpose() + q.theta * q.d2ThetadP2dP0);
	Matrix3 df2dP1 = -k * (q.dThetadP2 * q.dThetadP1.transpose() + q.theta * q.d2ThetadP2dP1);
	Matrix3 df2dP2 = -k * (q.dThetadP2 * q.dThetadP2.transpose() + q.theta * q.d2ThetadP2dP2);
	Matrix3 df2dP3 = -k * (q.dThetadP2 * q.dThetadP3.transpose() + q.theta * q.d2ThetadP2dP3);

	Matrix3 df3dP0 = -k * (q.dThetadP3 * q.dThetadP0.transpose() + q.theta * q.d2ThetadP3dP0);
	Matrix3 df3dP1 = -k * (q.dThetadP3 * q.dThetadP1.transpose() + q.theta * q.d2ThetadP3dP1);
	Matrix3 df3dP2 = -k * (q.dThetadP3 * q.dThetadP2.transpose() + q.theta * q.d2ThetadP3dP2);
	Matrix3 df3dP3 = -k * (q.dThetadP3 * q.dThetadP3.transpose() + q.theta * q.d2ThetadP3dP3);

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[0] + j) += df0dP0(i, j);
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[1] + j) += df0dP1(i, j);
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[2] + j) += df0dP2(i, j);
			dfdx.coeffRef(3 * m_inds[0] + i, 3 * m_inds[3] + j) += df0dP3(i, j);

			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[0] + j) += df1dP0(i, j);
			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[1] + j) += df1dP1(i, j);
			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[2] + j) += df1dP2(i, j);
			dfdx.coeffRef(3 * m_inds[1] + i, 3 * m_inds[3] + j) += df1dP3(i, j);

			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[0] + j) += df2dP0(i, j);
			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[1] + j) += df2dP1(i, j);
			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[2] + j) += df2dP2(i, j);
			dfdx.coeffRef(3 * m_inds[2] + i, 3 * m_inds[3] + j) += df2dP3(i, j);

			dfdx.coeffRef(3 * m_inds[3] + i, 3 * m_inds[0] + j) += df3dP0(i, j);
			dfdx.coeffRef(3 * m_inds[3] + i, 3 * m_inds[1] + j) += df3dP1(i, j);
			dfdx.coeffRef(3 * m_inds[3] + i, 3 * m_inds[2] + j) += df3dP2(i, j);
			dfdx.coeffRef(3 * m_inds[3] + i, 3 * m_inds[3] + j) += df3dP3(i, j);
		}
	}

	// compute damping forces and pseudo derivatives:

	// fd = -d * dTheta/dt * dTheta/dx:

	dampingForces.segment<3>(3 * m_inds[0]) -= d * q.dThetadt * q.dThetadP0;
	dampingForces.segment<3>(3 * m_inds[1]) -= d * q.dThetadt * q.dThetadP1;
	dampingForces.segment<3>(3 * m_inds[2]) -= d * q.dThetadt * q.dThetadP2;
	dampingForces.segment<3>(3 * m_inds[3]) -= d * q.dThetadt * q.dThetadP3;

	Matrix3 dfd0dP0 = -d * (q.dThetadP0 * q.dThetadP0.transpose());
	Matrix3 dfd0dP1 = -d * (q.dThetadP0 * q.dThetadP1.transpose());
	Matrix3 dfd0dP2 = -d * (q.dThetadP0 * q.dThetadP2.transpose());
	Matrix3 dfd0dP3 = -d * (q.dThetadP0 * q.dThetadP3.transpose());

	Matrix3 dfd1dP0 = -d * (q.dThetadP1 * q.dThetadP0.transpose());
	Matrix3 dfd1dP1 = -d * (q.dThetadP1 * q.dThetadP1.transpose());
	Matrix3 dfd1dP2 = -d * (q.dThetadP1 * q.dThetadP2.transpose());
	Matrix3 dfd1dP3 = -d * (q.dThetadP1 * q.dThetadP3.transpose());

	Matrix3 dfd2dP0 = -d * (q.dThetadP2 * q.dThetadP0.transpose());
	Matrix3 dfd2dP1 = -d * (q.dThetadP2 * q.dThetadP1.transpose());
	Matrix3 dfd2dP2 = -d * (q.dThetadP2 * q.dThetadP2.transpose());
	Matrix3 dfd2dP3 = -d * (q.dThetadP2 * q.dThetadP3.transpose());

	Matrix3 dfd3dP0 = -d * (q.dThetadP3 * q.dThetadP0.transpose());
	Matrix3 dfd3dP1 = -d * (q.dThetadP3 * q.dThetadP1.transpose());
	Matrix3 dfd3dP2 = -d * (q.dThetadP3 * q.dThetadP2.transpose());
	Matrix3 dfd3dP3 = -d * (q.dThetadP3 * q.dThetadP3.transpose());

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			dampingPseudoDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[0] + j) += dfd0dP0(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[1] + j) += dfd0dP1(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[2] + j) += dfd0dP2(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[0] + i, 3 * m_inds[3] + j) += dfd0dP3(i, j);

			dampingPseudoDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[0] + j) += dfd1dP0(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[1] + j) += dfd1dP1(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[2] + j) += dfd1dP2(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[1] + i, 3 * m_inds[3] + j) += dfd1dP3(i, j);

			dampingPseudoDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[0] + j) += dfd2dP0(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[1] + j) += dfd2dP1(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[2] + j) += dfd2dP2(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[2] + i, 3 * m_inds[3] + j) += dfd2dP3(i, j);

			dampingPseudoDerivatives.coeffRef(3 * m_inds[3] + i, 3 * m_inds[0] + j) += dfd3dP0(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[3] + i, 3 * m_inds[1] + j) += dfd3dP1(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[3] + i, 3 * m_inds[2] + j) += dfd3dP2(i, j);
			dampingPseudoDerivatives.coeffRef(3 * m_inds[3] + i, 3 * m_inds[3] + j) += dfd3dP3(i, j);
		}
	}
}

template class BendCondition<float>;
template class BendCondition<double>;
