#include "EqualityTests.h"

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

double numericalForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, double k, int i, double dx)
{
	double xOrig = x[i];
	x[i] = xOrig + dx;
	double ePlus = 0.5 * c.C(x, uv).squaredNorm();
	x[i] = xOrig - dx;
	double eMinus = 0.5 * c.C(x, uv).squaredNorm();
	x[i] = xOrig;

	// f = -dE/dx
	return -k * (ePlus - eMinus) / (2 * dx);
}

Eigen::VectorXd numericalForceDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, int i, double dx)
{
	Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
	Eigen::SparseMatrix<double> dfdv((int)x.size(), (int)x.size());
	double d = 1.0;
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
