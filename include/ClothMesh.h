#ifndef CLOTHMESH_H
#define CLOTHMESH_H

#include "LinearSolver.h"

#include "BendCondition.h"
#include "ShearCondition.h"
#include "StretchCondition.h"

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
		const Vector &m_uv,
		const std::vector<int> &triangleIndices,
		Real kBend, Real kStretch, Real kShear,
		Real dBend, Real dStretch, Real dShear
	);

	// accessor for positions:
	const Vector &x() const;

	// accessor for velocities:
	const Vector &v() const;

	// advance the simulation:
	void advance(Real dt, const LinearSolver<Real> &solver);

private:

	// positions:
	Vector m_x;

	// velocities:
	Vector m_v;

	// reference pose:
	const Vector m_uv;

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

	// solver data:
	SparseMatrix m_dfdx;
	SparseMatrix m_implicitUpdateMatrix;
	Vector m_implicitUpdateRHS;

};

#endif
