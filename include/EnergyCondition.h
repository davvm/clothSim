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

	// compute the vector condition function C for testing purposes:
	virtual Vector C(const Vector& x, const Vector& uv) const = 0;

	// computes -dE/dX
	virtual void computeForces(const Vector& x, const Vector& uv, Vector& forces) const = 0;

	// computes -d2E/dXidXj
	virtual void computeForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const = 0;

	virtual void computeDampingForces(const Vector& x, const Vector& v, const Vector& uv, Vector& forces) const = 0;
	virtual void computeDampingForceDerivatives(const Vector& x, const Vector& uv, Eigen::SparseMatrix<Real>& dfdx) const = 0;

};

#endif // ENERGYCONDITION_H
