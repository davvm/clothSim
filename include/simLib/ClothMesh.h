#ifndef CLOTHMESH_H
#define CLOTHMESH_H


#include "LinearSolver.h"

#include "BendCondition.h"
#include "ShearCondition.h"
#include "StretchCondition.h"

#include "ForceField.h"

template <class Real>
class ClothMesh
{

public:

	typedef Eigen::Matrix<Real, -1, 1, 0, -1, 1> Vector;
	typedef Eigen::Matrix<Real, 2, 1, 0, 2, 1> Vector2;
	typedef Eigen::Matrix<Real, 3, 1, 0, 3, 1> Vector3;
	typedef Eigen::Matrix<Real, 3, 3, 0, 3, 3> Matrix3;
	typedef Eigen::SparseMatrix<Real> SparseMatrix;

	ClothMesh(
		const Vector &x,
		const Vector &v,
		const Vector &uv,
		const std::vector<int> &triangleIndices,
		Real kBend, Real kStretch, Real kShear,
		Real dBend, Real dStretch, Real dShear,
		Real density
	);

	// accessor for positions:
	const Vector &x() const;

	// accessor for velocities:
	const Vector &v() const;

	// accessor for masses:
	const Vector &m() const;

	// accessor for indices:
	const std::vector<int> &triangleIndices();

	// accessors for the energy terms:
	const std::vector< BendCondition<Real> > &bendConditions() const;
	const std::vector< ShearCondition<Real> > &shearConditions() const;
	const std::vector< StretchCondition<Real> > &stretchConditions() const;

	// advance the simulation:
	void advance( std::vector< ForceField<Real>* > &forceFields, Real dt, const LinearSolver<Real> &solver);

	// calculate forces and derivatives:
	void forcesAndDerivatives(const Vector &x, const Vector &uv, const Vector &v, Vector &f, Vector &d, SparseMatrix &dfdx, SparseMatrix &dddx, SparseMatrix &dddv) const;

	// get the implicit update linear system ready:
	void assembleImplicitUpdateEquations(Real dt, const Vector &f, const Vector &v, const SparseMatrix &dfdx, const SparseMatrix &dfdv, SparseMatrix &implicitUpdateMatrix, Vector &implicitUpdateRHS) const;

	// energy function for testing:
	Real energy(const Vector& x, const Vector& uv) const;

	// vector containing all the constraint functions for testing:
	void C(const Vector& x, const Vector& uv, Vector& c) const;

private:

	// positions:
	Vector m_x;

	// velocities:
	Vector m_v;

	// masses:
	Vector m_m;

	// reference pose:
	Vector m_uv;

	SparseMatrix m_implicitUpdateMatrix;
	SparseMatrix m_dfdx;
	SparseMatrix m_dfdv;

	const std::vector<int> m_triangleIndices;

	// material properties:
	Real m_kBend;
	Real m_kStretch;
	Real m_kShear;

	Real m_dBend;
	Real m_dStretch;
	Real m_dShear;

	// energy terms:
	std::vector< BendCondition<Real> > m_bendConditions;
	std::vector< ShearCondition<Real> > m_shearConditions;
	std::vector< StretchCondition<Real> > m_stretchConditions;

};

#endif
