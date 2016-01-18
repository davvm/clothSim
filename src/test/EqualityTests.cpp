#include "test\EqualityTests.h"

#include "CppUnitTest.h"

void checkVectorEquality(Eigen::Vector3d v0, Eigen::Vector3d v1, double tol)
{
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[0], v1[0], tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[1], v1[1], tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[2], v1[2], tol);
}

void checkVectorEquality(Eigen::VectorXd v0, Eigen::VectorXd v1, double tol, bool relative)
{
	if (relative)
	{
		for (int i = 0; i < v0.size(); ++i)
		{
			if (fabs(v0[i]) > 0.01)
			{
				Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[i], v1[i], tol*v0[i]);
			}
			else
			{
				Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[i], v1[i], tol);
			}
		}
	}
	else
	{
		for (int i = 0; i < v0.size(); ++i)
		{
			Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(v0[i], v1[i], tol);
		}
	}
}

void checkMatrixEquality(Eigen::Matrix3d m0, Eigen::Matrix3d m1, double tol)
{
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(0, 0), m1(0, 0), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(0, 1), m1(0, 1), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(0, 2), m1(0, 2), tol);

	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(1, 0), m1(1, 0), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(1, 1), m1(1, 1), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(1, 2), m1(1, 2), tol);

	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(2, 0), m1(2, 0), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(2, 1), m1(2, 1), tol);
	Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(m0(2, 2), m1(2, 2), tol);
}

double numericalForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, const Eigen::VectorXd &x, double k, int i, double dx)
{
	Eigen::VectorXd xtest = x;
	xtest[i] = x[i] + dx;
	double ePlus = 0.5 * c.C(xtest, uv).squaredNorm();
	xtest[i] = x[i] - dx;
	double eMinus = 0.5 * c.C(xtest, uv).squaredNorm();

	// f = -k * dE/dx
	return -k * (ePlus - eMinus) / (2 * dx);
}

static Eigen::VectorXd numericalCTimeDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &x, const Eigen::VectorXd &v, const Eigen::VectorXd &uv, double dx)
{
	Eigen::VectorXd xtest;
	xtest = x + dx * v;
	Eigen::VectorXd cPlus = c.C(xtest, uv);
	xtest = x - dx * v;
	Eigen::VectorXd cMinus = c.C(xtest, uv);

	return (cPlus - cMinus) / (2 * dx);
}

static Eigen::VectorXd numericalFirstCDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &x, const Eigen::VectorXd &uv, int i, double dx)
{
	Eigen::VectorXd xtest = x;

	xtest[i] = x[i] + dx;
	Eigen::VectorXd cPlus = c.C(xtest, uv);
	xtest[i] = x[i] - dx;
	Eigen::VectorXd cMinus = c.C(xtest, uv);

	return (cPlus - cMinus) / (2 * dx);
}

static Eigen::VectorXd numericalSecondCDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &x, const Eigen::VectorXd &uv, int i, int j, double dx)
{
	Eigen::VectorXd xtest = x;

	xtest[j] = x[j] + dx;
	Eigen::VectorXd dCdxPlus = numericalFirstCDerivative(c, xtest, uv, i, dx);
	xtest[j] = x[j] - dx;
	Eigen::VectorXd dCdxcMinus = numericalFirstCDerivative(c, xtest, uv, i, dx);

	return (dCdxPlus - dCdxcMinus) / (2 * dx);
}

double numericalDampingForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, const Eigen::VectorXd &x, const Eigen::VectorXd &v, double d, int i, double dx)
{
	// find dC/dt
	Eigen::VectorXd dCdt = numericalCTimeDerivative(c, x, v, uv, dx);

	// find dC/dx
	Eigen::VectorXd dCdx = numericalFirstCDerivative(c, x, uv, i, dx);
	
	// fd = -d * sum( i, dC_i/dx * dC_i/dt )
	return -d * (dCdx.array() * dCdt.array()).sum();
}

Eigen::VectorXd numericalForceDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, double d, int i, double dx)
{
	Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
	Eigen::SparseMatrix<double> dfdv((int)x.size(), (int)x.size());
	Eigen::VectorXd dampingForces((int)x.size());
	Eigen::SparseMatrix<double> dampingPseudoDerivatives((int)x.size(), (int)x.size());

	double xOrig = x[i];

	x[i] = xOrig + dx;

	Eigen::VectorXd fPlus(x.size());
	fPlus.setConstant(0);
	c.computeForces(x, uv, k, fPlus, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

	x[i] = xOrig - dx;

	Eigen::VectorXd fMinus(x.size());
	fMinus.setConstant(0);
	c.computeForces(x, uv, k, fMinus, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

	x[i] = xOrig;

	// df/dx
	return (fPlus - fMinus) / (2 * dx);
}

Eigen::VectorXd numericalDampingForceDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, double d, int i, double dx)
{
	Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
	Eigen::SparseMatrix<double> dfdv((int)x.size(), (int)x.size());
	Eigen::VectorXd f((int)x.size());
	Eigen::SparseMatrix<double> dampingPseudoDerivatives((int)x.size(), (int)x.size());

	double vOrig = v[i];

	v[i] = vOrig + dx;

	Eigen::VectorXd dampingForcePlus(x.size());
	dampingForcePlus.setConstant(0);
	c.computeForces(x, uv, k, f, dfdx, v, d, dampingForcePlus, dampingPseudoDerivatives, dfdv);

	v[i] = vOrig - dx;

	Eigen::VectorXd dampingForceMinus(x.size());
	dampingForceMinus.setConstant(0);
	c.computeForces(x, uv, k, f, dfdx, v, d, dampingForceMinus, dampingPseudoDerivatives, dfdv);

	v[i] = vOrig;

	// df/dx
	return (dampingForcePlus - dampingForceMinus) / (2 * dx);
}

void checkPseudoDerivatives(const Eigen::SparseMatrix<double> &dampingPseudoDerivatives, const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, double d, double dx, double tol)
{
	Eigen::VectorXd dCdt = numericalCTimeDerivative(c, x, v, uv, dx);

	for (int i = 0; i < x.size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			Eigen::VectorXd d2Cdidj = numericalSecondCDerivative(c, x, uv, i, j, dx);
			double expectedCoeff = -d * (d2Cdidj.array() * dCdt.array()).sum();
			double actualCoeff = dampingPseudoDerivatives.coeff(i, j);
			Microsoft::VisualStudio::CppUnitTestFramework::Assert::AreEqual(actualCoeff, expectedCoeff, tol);
		}
	}
}
