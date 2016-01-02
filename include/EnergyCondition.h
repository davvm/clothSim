#ifndef ENERGYCONDITION_H
#define ENERGYCONDITION_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

// abstract class for defining energy conditions and adding to a vector of
// forces/matrix of force derivatives:
template <class Real>
class EnergyCondition
{
public:

	virtual ~EnergyCondition() {}

	typedef Eigen::Matrix<Real, -1, 1, 0, -1, 1> Vector;
	typedef Eigen::Matrix<Real, 3, 1, 0, 3, 1> Vector3;
	typedef Eigen::Matrix<Real, 3, 3, 0, 3, 3> Matrix3;
	typedef Eigen::SparseMatrix<Real> SparseMatrix;

	// compute the vector condition function C for testing purposes:
	virtual Vector C(const Vector& x, const Vector& uv) const = 0;

	// computes forces and their derivatives:
	virtual void computeForces(const Vector& x, const Vector& uv, Real k, Vector& forces, SparseMatrix &dfdx, const Vector& v, Real d, Vector &dampingForces, SparseMatrix &dampingPseudoDerivatives) const = 0;

};

#endif // ENERGYCONDITION_H
