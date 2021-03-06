#ifndef GRAVITYFIELD_H
#define GRAVITYFIELD_H

#include "ForceField.h"

namespace ClothSim
{

template <class Real>
class GravityField : public ForceField<Real>
{
public:

	GravityField(const Vector &m, Vector3 g);

	// computes forces and their derivatives:
	virtual void forcesAndDerivatives(Vector& forces, SparseMatrix &dfdx) const;

private:

	const Vector &m_m;
	Vector3 m_g;

};

} //namespace ClothSim


#endif // GRAVITYFIELD_H
