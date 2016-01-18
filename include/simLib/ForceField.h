#ifndef FORCEFIELD_H
#define FORCEFIELD_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

// abstract class for defining force fields and adding to a vector of
// forces/matrix of force derivatives:
template <class Real>
class ForceField
{
public:

	virtual ~ForceField() {}

	typedef Eigen::Matrix<Real, -1, 1, 0, -1, 1> Vector;
	typedef Eigen::Matrix<Real, 2, 1, 0, 2, 1> Vector2;
	typedef Eigen::Matrix<Real, 3, 1, 0, 3, 1> Vector3;
	typedef Eigen::Matrix<Real, 3, 3, 0, 3, 3> Matrix3;
	typedef Eigen::SparseMatrix<Real> SparseMatrix;

	// computes forces and their derivatives:
	virtual void forcesAndDerivatives(Vector& forces, SparseMatrix &dfdx) const = 0;

};

#endif // FORCEFIELD_H
