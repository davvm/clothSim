#include "CppUnitTest.h"

#include "ShearCondition.h"
#include "StretchCondition.h"
#include "BendCondition.h"

#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{		
	TEST_CLASS(UnitTest1)
	{
	public:

		double numericalForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, int i, double dx)
		{
			double xOrig = x[i];
			x[i] = xOrig + dx;
			double ePlus = 0.5 * c.C(x, uv).squaredNorm();
			x[i] = xOrig - dx;
			double eMinus = 0.5 * c.C(x, uv).squaredNorm();
			x[i] = xOrig;

			// f = -dE/dx
			return -(ePlus - eMinus) / (2 * dx);
		}

		TEST_METHOD(TestBendCondition)
		{
			Eigen::VectorXd x(5 * 3);
			
			x[3 * 0] = -2.0; 
			x[3 * 0 + 1] = 0;
			x[3 * 0 + 2] = 0.4;

			x[3 * 1] = 0;
			x[3 * 1 + 1] = 1.0;
			x[3 * 1 + 2] = 0;

			x[3 * 2] = 0;
			x[3 * 2 + 1] = -1.0;
			x[3 * 2 + 2] = 0;

			x[3 * 3] = 1.0;
			x[3 * 3 + 1] = 0;
			x[3 * 3 + 2] = 0.2;

			x[3 * 4] = 0.0;
			x[3 * 4 + 1] = 0;
			x[3 * 4 + 2] = 0;

			Eigen::VectorXd uv(5 * 2);
			Eigen::VectorXd f(5 * 3);
			f.setConstant(0);
			BendCondition<double> bc(0, 1, 2, 3);

			// check we've got the energy condition right, and we're measuring the angle between the normals:
			Eigen::VectorXd c = bc.C(x, uv);
			double expectedAngle = 2 * atan( 0.2 );
			Assert::AreEqual(c[0], expectedAngle, 1.e-4);

			// compute forces analytically:
			bc.computeForces(x, uv, f);

			// compare to numerically computed forces:
			double dx = 0.0001;
			double tol = 1.e-7;
			Assert::AreEqual(numericalForce(bc, uv, x, 0, dx), f[0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 1, dx), f[1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 2, dx), f[2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 0, dx), f[3 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 1, dx), f[3 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 2, dx), f[3 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 0, dx), f[6 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 1, dx), f[6 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 2, dx), f[6 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 0, dx), f[9 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 1, dx), f[9 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 2, dx), f[9 + 2], tol);

			// turn the angle the other way:
			x[3 * 0 + 2] = -0.4;
			x[3 * 3 + 2] = -0.2;
			f.setConstant(0);
			bc.computeForces(x, uv, f);
			
			Assert::AreEqual(numericalForce(bc, uv, x, 0, dx), f[0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 1, dx), f[1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 2, dx), f[2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 0, dx), f[3 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 1, dx), f[3 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 2, dx), f[3 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 0, dx), f[6 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 1, dx), f[6 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 2, dx), f[6 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 0, dx), f[9 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 1, dx), f[9 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 2, dx), f[9 + 2], tol);

			// randomize positions:
			x = Eigen::VectorXd::Random(5 * 3);
			f.setConstant(0);
			bc.computeForces(x, uv, f);

			Assert::AreEqual(numericalForce(bc, uv, x, 0, dx), f[0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 1, dx), f[1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 2, dx), f[2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 0, dx), f[3 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 1, dx), f[3 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 2, dx), f[3 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 0, dx), f[6 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 1, dx), f[6 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 2, dx), f[6 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 0, dx), f[9 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 1, dx), f[9 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 2, dx), f[9 + 2], tol);

			// randomize positions again:
			x = Eigen::VectorXd::Random(5 * 3);
			f.setConstant(0);
			bc.computeForces(x, uv, f);

			Assert::AreEqual(numericalForce(bc, uv, x, 0, dx), f[0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 1, dx), f[1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 2, dx), f[2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 0, dx), f[3 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 1, dx), f[3 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 3 + 2, dx), f[3 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 0, dx), f[6 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 1, dx), f[6 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 6 + 2, dx), f[6 + 2], tol);

			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 0, dx), f[9 + 0], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 1, dx), f[9 + 1], tol);
			Assert::AreEqual(numericalForce(bc, uv, x, 9 + 2, dx), f[9 + 2], tol);
		}

		TEST_METHOD(TestShearCondition)
		{
			Eigen::VectorXf x(30);
			ShearCondition<double> s(1.0, 0, 1, 2, 1.0);
		}

		TEST_METHOD(TestStretchCondition)
		{
			Eigen::VectorXf x(30);
			StretchCondition<double> s(1.0, 0, 1, 2, 1.0, 1.0, 1.0);
		}

	};
}
