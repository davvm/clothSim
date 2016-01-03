#include "CppUnitTest.h"

#include "EqualityTests.h"

#include "ShearCondition.h"

#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{
	TEST_CLASS(ShearConditionTest)
	{
	public:
		
		void shearTriangleQuantitiesForDerivative(Eigen::VectorXd x, Eigen::VectorXd uv, Eigen::VectorXd v, int idx, double dx, ShearCondition<double>::TriangleQuantities &qPlus, ShearCondition<double>::TriangleQuantities &qMinus) const
		{
			double component = x[idx];
			x[idx] = component + dx;
			qPlus = ShearCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);
			x[idx] = component - dx;
			qMinus = ShearCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);
			x[idx] = component;
		}

		void checkShearTriangleQuantities(const Eigen::VectorXd &x, const Eigen::VectorXd &uv, const Eigen::VectorXd &v, double dx, double tol)
		{
			ShearCondition<double>::TriangleQuantities q(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);

			Eigen::VectorXd xtest;
			xtest = x + v * dx;
			ShearCondition<double>::TriangleQuantities qtPlus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);
			xtest = x - v * dx;
			ShearCondition<double>::TriangleQuantities qtMinus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);

			double dCdt = (qtPlus.C - qtMinus.C) / (2 * dx);
			Assert::AreEqual(dCdt, q.dCdt, tol);


			ShearCondition<double>::TriangleQuantities q0xPlus, q0xMinus, q0yPlus, q0yMinus, q0zPlus, q0zMinus;
			ShearCondition<double>::TriangleQuantities q1xPlus, q1xMinus, q1yPlus, q1yMinus, q1zPlus, q1zMinus;
			ShearCondition<double>::TriangleQuantities q2xPlus, q2xMinus, q2yPlus, q2yMinus, q2zPlus, q2zMinus;

			shearTriangleQuantitiesForDerivative(x, uv, v, 0, dx, q0xPlus, q0xMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 1, dx, q0yPlus, q0yMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 2, dx, q0zPlus, q0zMinus);

			shearTriangleQuantitiesForDerivative(x, uv, v, 3, dx, q1xPlus, q1xMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 4, dx, q1yPlus, q1yMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 5, dx, q1zPlus, q1zMinus);

			shearTriangleQuantitiesForDerivative(x, uv, v, 6, dx, q2xPlus, q2xMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 7, dx, q2yPlus, q2yMinus);
			shearTriangleQuantitiesForDerivative(x, uv, v, 8, dx, q2zPlus, q2zMinus);


			// test first normalized wu and wv derivatives:
			Eigen::Matrix3d dwudP0;
			dwudP0.col(0) = (q0xPlus.wu - q0xMinus.wu) / (2 * dx);
			dwudP0.col(1) = (q0yPlus.wu - q0yMinus.wu) / (2 * dx);
			dwudP0.col(2) = (q0zPlus.wu - q0zMinus.wu) / (2 * dx);
			checkMatrixEquality(dwudP0, q.dwudP0, tol);

			Eigen::Matrix3d dwudP1;
			dwudP1.col(0) = (q1xPlus.wu - q1xMinus.wu) / (2 * dx);
			dwudP1.col(1) = (q1yPlus.wu - q1yMinus.wu) / (2 * dx);
			dwudP1.col(2) = (q1zPlus.wu - q1zMinus.wu) / (2 * dx);
			checkMatrixEquality(dwudP1, q.dwudP1, tol);

			Eigen::Matrix3d dwudP2;
			dwudP2.col(0) = (q2xPlus.wu - q2xMinus.wu) / (2 * dx);
			dwudP2.col(1) = (q2yPlus.wu - q2yMinus.wu) / (2 * dx);
			dwudP2.col(2) = (q2zPlus.wu - q2zMinus.wu) / (2 * dx);
			checkMatrixEquality(dwudP2, q.dwudP2, tol);

			Eigen::Matrix3d dwvdP0;
			dwvdP0.col(0) = (q0xPlus.wv - q0xMinus.wv) / (2 * dx);
			dwvdP0.col(1) = (q0yPlus.wv - q0yMinus.wv) / (2 * dx);
			dwvdP0.col(2) = (q0zPlus.wv - q0zMinus.wv) / (2 * dx);
			checkMatrixEquality(dwvdP0, q.dwvdP0, tol);

			Eigen::Matrix3d dwvdP1;
			dwvdP1.col(0) = (q1xPlus.wv - q1xMinus.wv) / (2 * dx);
			dwvdP1.col(1) = (q1yPlus.wv - q1yMinus.wv) / (2 * dx);
			dwvdP1.col(2) = (q1zPlus.wv - q1zMinus.wv) / (2 * dx);
			checkMatrixEquality(dwvdP1, q.dwvdP1, tol);

			Eigen::Matrix3d dwvdP2;
			dwvdP2.col(0) = (q2xPlus.wv - q2xMinus.wv) / (2 * dx);
			dwvdP2.col(1) = (q2yPlus.wv - q2yMinus.wv) / (2 * dx);
			dwvdP2.col(2) = (q2zPlus.wv - q2zMinus.wv) / (2 * dx);
			checkMatrixEquality(dwvdP2, q.dwvdP2, tol);


			// test first C derivatives:
			Eigen::Vector3d dCdP0;
			dCdP0[0] = (q0xPlus.C - q0xMinus.C) / (2 * dx);
			dCdP0[1] = (q0yPlus.C - q0yMinus.C) / (2 * dx);
			dCdP0[2] = (q0zPlus.C - q0zMinus.C) / (2 * dx);
			checkVectorEquality(dCdP0, q.dCdP0, tol);

			Eigen::Vector3d dCdP1;
			dCdP1[0] = (q1xPlus.C - q1xMinus.C) / (2 * dx);
			dCdP1[1] = (q1yPlus.C - q1yMinus.C) / (2 * dx);
			dCdP1[2] = (q1zPlus.C - q1zMinus.C) / (2 * dx);
			checkVectorEquality(dCdP1, q.dCdP1, tol);

			Eigen::Vector3d dCdP2;
			dCdP2[0] = (q2xPlus.C - q2xMinus.C) / (2 * dx);
			dCdP2[1] = (q2yPlus.C - q2yMinus.C) / (2 * dx);
			dCdP2[2] = (q2zPlus.C - q2zMinus.C) / (2 * dx);
			checkVectorEquality(dCdP2, q.dCdP2, tol);


			// test second C derivatives:
			Eigen::Matrix3d d2CdP0dP0;
			d2CdP0dP0.col(0) = (q0xPlus.dCdP0 - q0xMinus.dCdP0) / (2 * dx);
			d2CdP0dP0.col(1) = (q0yPlus.dCdP0 - q0yMinus.dCdP0) / (2 * dx);
			d2CdP0dP0.col(2) = (q0zPlus.dCdP0 - q0zMinus.dCdP0) / (2 * dx);
			checkMatrixEquality(d2CdP0dP0, q.d2CdP0dP0, tol);

			Eigen::Matrix3d d2CdP0dP1;
			d2CdP0dP1.col(0) = (q1xPlus.dCdP0 - q1xMinus.dCdP0) / (2 * dx);
			d2CdP0dP1.col(1) = (q1yPlus.dCdP0 - q1yMinus.dCdP0) / (2 * dx);
			d2CdP0dP1.col(2) = (q1zPlus.dCdP0 - q1zMinus.dCdP0) / (2 * dx);
			checkMatrixEquality(d2CdP0dP1, q.d2CdP0dP1, tol);

			Eigen::Matrix3d d2CdP0dP2;
			d2CdP0dP2.col(0) = (q2xPlus.dCdP0 - q2xMinus.dCdP0) / (2 * dx);
			d2CdP0dP2.col(1) = (q2yPlus.dCdP0 - q2yMinus.dCdP0) / (2 * dx);
			d2CdP0dP2.col(2) = (q2zPlus.dCdP0 - q2zMinus.dCdP0) / (2 * dx);
			checkMatrixEquality(d2CdP0dP2, q.d2CdP0dP2, tol);
			


			Eigen::Matrix3d d2CdP1dP0;
			d2CdP1dP0.col(0) = (q0xPlus.dCdP1 - q0xMinus.dCdP1) / (2 * dx);
			d2CdP1dP0.col(1) = (q0yPlus.dCdP1 - q0yMinus.dCdP1) / (2 * dx);
			d2CdP1dP0.col(2) = (q0zPlus.dCdP1 - q0zMinus.dCdP1) / (2 * dx);
			checkMatrixEquality(d2CdP1dP0, q.d2CdP1dP0, tol);
			
			Eigen::Matrix3d d2CdP1dP1;
			d2CdP1dP1.col(0) = (q1xPlus.dCdP1 - q1xMinus.dCdP1) / (2 * dx);
			d2CdP1dP1.col(1) = (q1yPlus.dCdP1 - q1yMinus.dCdP1) / (2 * dx);
			d2CdP1dP1.col(2) = (q1zPlus.dCdP1 - q1zMinus.dCdP1) / (2 * dx);
			checkMatrixEquality(d2CdP1dP1, q.d2CdP1dP1, tol);
			
			Eigen::Matrix3d d2CdP1dP2;
			d2CdP1dP2.col(0) = (q2xPlus.dCdP1 - q2xMinus.dCdP1) / (2 * dx);
			d2CdP1dP2.col(1) = (q2yPlus.dCdP1 - q2yMinus.dCdP1) / (2 * dx);
			d2CdP1dP2.col(2) = (q2zPlus.dCdP1 - q2zMinus.dCdP1) / (2 * dx);
			checkMatrixEquality(d2CdP1dP2, q.d2CdP1dP2, tol);


			Eigen::Matrix3d d2CdP2dP0;
			d2CdP2dP0.col(0) = (q0xPlus.dCdP2 - q0xMinus.dCdP2) / (2 * dx);
			d2CdP2dP0.col(1) = (q0yPlus.dCdP2 - q0yMinus.dCdP2) / (2 * dx);
			d2CdP2dP0.col(2) = (q0zPlus.dCdP2 - q0zMinus.dCdP2) / (2 * dx);
			checkMatrixEquality(d2CdP2dP0, q.d2CdP2dP0, tol);

			Eigen::Matrix3d d2CdP2dP1;
			d2CdP2dP1.col(0) = (q1xPlus.dCdP2 - q1xMinus.dCdP2) / (2 * dx);
			d2CdP2dP1.col(1) = (q1yPlus.dCdP2 - q1yMinus.dCdP2) / (2 * dx);
			d2CdP2dP1.col(2) = (q1zPlus.dCdP2 - q1zMinus.dCdP2) / (2 * dx);
			checkMatrixEquality(d2CdP2dP1, q.d2CdP2dP1, tol);
			
			Eigen::Matrix3d d2CdP2dP2;
			d2CdP2dP2.col(0) = (q2xPlus.dCdP2 - q2xMinus.dCdP2) / (2 * dx);
			d2CdP2dP2.col(1) = (q2yPlus.dCdP2 - q2yMinus.dCdP2) / (2 * dx);
			d2CdP2dP2.col(2) = (q2zPlus.dCdP2 - q2zMinus.dCdP2) / (2 * dx);
			checkMatrixEquality(d2CdP2dP2, q.d2CdP2dP2, tol);

		}


		TEST_METHOD(TestShearConditionTriangleQuantities)
		{
			Eigen::VectorXd uv(3 * 2);

			uv[2 * 0] = 1.0;
			uv[2 * 0 + 1] = 2.0;

			uv[2 * 1] = 6.0;
			uv[2 * 1 + 1] = 3.0;

			uv[2 * 2] = 2.0;
			uv[2 * 2 + 1] = 5.0;

			Eigen::Vector3d p0(1, 2, 3);
			Eigen::Vector3d wu(1.3, 0.2, -0.1);
			Eigen::Vector3d wv(-0.5, 1.1, 0.3);

			Eigen::VectorXd x(3 * 3);
			x.segment<3>(3 * 0) = p0 + wu * uv[2 * 0 + 0] + wv * uv[2 * 0 + 1];
			x.segment<3>(3 * 1) = p0 + wu * uv[2 * 1 + 0] + wv * uv[2 * 1 + 1];
			x.segment<3>(3 * 2) = p0 + wu * uv[2 * 2 + 0] + wv * uv[2 * 2 + 1];

			Eigen::VectorXd v = Eigen::VectorXd::Random(3 * 3);

			ShearCondition<double>::TriangleQuantities q(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2)
				);

			double dx = 0.0001;
			double tol = 1.e-5;

			// make sure we're calculating wu, wv and a as expected:
			checkVectorEquality(q.wu, wu, tol);
			checkVectorEquality(q.wv, wv, tol);
			Eigen::Vector3d uv0(uv[2 * 0], uv[2 * 0 + 1], 0);
			Eigen::Vector3d uv1(uv[2 * 1], uv[2 * 1 + 1], 0);
			Eigen::Vector3d uv2(uv[2 * 2], uv[2 * 2 + 1], 0);
			Assert::AreEqual(q.a, 0.5 * ((uv1 - uv0).cross(uv2 - uv0)).norm(), tol);

			checkShearTriangleQuantities(x, uv, v, dx, tol);

			// test for 10 random configurations:
			for (int i = 0; i < 10; ++i)
			{
				// randomize the positions:
				x = Eigen::VectorXd::Random(3 * 3);
				uv = Eigen::VectorXd::Random(3 * 2);
				v = Eigen::VectorXd::Random(3 * 3);
				checkShearTriangleQuantities(x, uv, v, dx, tol);
			}

		}

		TEST_METHOD(TestShearCondition)
		{
			Eigen::VectorXd uv(3 * 2);

			uv[2 * 0] = 1.0;
			uv[2 * 0 + 1] = 2.0;

			uv[2 * 1] = 6.0;
			uv[2 * 1 + 1] = 3.0;

			uv[2 * 2] = 2.0;
			uv[2 * 2 + 1] = 5.0;

			Eigen::Vector3d p0(1, 2, 3);
			Eigen::Vector3d wu(1.3, 0.2, -0.1);
			Eigen::Vector3d wv(-0.5, 1.1, 0.3);

			Eigen::VectorXd x(3 * 3);
			x.segment<3>(3 * 0) = p0 + wu * uv[2 * 0 + 0] + wv * uv[2 * 0 + 1];
			x.segment<3>(3 * 1) = p0 + wu * uv[2 * 1 + 0] + wv * uv[2 * 1 + 1];
			x.segment<3>(3 * 2) = p0 + wu * uv[2 * 2 + 0] + wv * uv[2 * 2 + 1];

			Eigen::VectorXd v(3 * 3);

			v[3 * 0] = 0.1;
			v[3 * 0 + 1] = 0.2;
			v[3 * 0 + 2] = 0.3;

			v[3 * 1] = -.1;
			v[3 * 1 + 1] = -1;
			v[3 * 1 + 2] = 0.4;

			v[3 * 2] = 0;
			v[3 * 2 + 1] = -0.2;
			v[3 * 2 + 2] = 0.1;

			double bu = 0.9;
			double bv = 1.1;

			Eigen::VectorXd f(3 * 3);
			f.setConstant(0);
			Eigen::VectorXd dampingForces(3 * 3);
			dampingForces.setConstant(0);
			ShearCondition<double> sc(0, 1, 2);

			double k = 1.5;
			Eigen::SparseMatrix<double> dfdx(3 * 3, 3 * 3);
			double d = 1.0;
			Eigen::SparseMatrix<double> dampingPseudoDerivatives(3 * 3, 3 * 3);
			Eigen::SparseMatrix<double> dfdv(3 * 3, 3 * 3);

			// compute forces analytically:
			sc.computeForces(x, uv, k, f, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

			// compare to numerically computed forces:
			double dx = 0.0001;
			double tol = 1.e-4;
			Assert::AreEqual(numericalForce(sc, uv, x, k, 0, dx), f[0], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 1, dx), f[1], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 2, dx), f[2], tol);

			Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 0, dx), f[3 + 0], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 1, dx), f[3 + 1], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 2, dx), f[3 + 2], tol);

			Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 0, dx), f[6 + 0], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 1, dx), f[6 + 1], tol);
			Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 2, dx), f[6 + 2], tol);

			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 0, dx), dampingForces[0], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 1, dx), dampingForces[1], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 2, dx), dampingForces[2], tol);

			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 0, dx), dampingForces[3 + 0], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 1, dx), dampingForces[3 + 1], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 2, dx), dampingForces[3 + 2], tol);

			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 0, dx), dampingForces[6 + 0], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 1, dx), dampingForces[6 + 1], tol);
			Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 2, dx), dampingForces[6 + 2], tol);

			for (int i = 0; i < x.size(); ++i)
			{
				Eigen::VectorXd fdN = dfdx.row(i);
				checkVectorEquality(fdN, numericalForceDerivative(sc, uv, x, v, k, d, i, dx), tol, true);
				fdN = dfdv.row(i);
				checkVectorEquality(fdN, numericalDampingForceDerivative(sc, uv, x, v, k, d, i, dx), tol, true);
			}

			// test on 10 randomized configurations:
			for (size_t n = 0; n < 10; ++n)
			{
				// randomize positions:
				x = Eigen::VectorXd::Random(3 * 3);
				f.setConstant(0);
				for (int i = 0; i < dfdx.outerSize(); ++i)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(dfdx, i); it; ++it)
					{
						it.valueRef() = 0;
					}
				}
				dampingForces.setConstant(0);
				for (int i = 0; i < dampingPseudoDerivatives.outerSize(); ++i)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(dampingPseudoDerivatives, i); it; ++it)
					{
						it.valueRef() = 0;
					}
				}
				for (int i = 0; i < dfdv.outerSize(); ++i)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(dfdv, i); it; ++it)
					{
						it.valueRef() = 0;
					}
				}
				sc.computeForces(x, uv, k, f, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

				Assert::AreEqual(numericalForce(sc, uv, x, k, 0, dx), f[0], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 1, dx), f[1], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 2, dx), f[2], tol);

				Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 0, dx), f[3 + 0], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 1, dx), f[3 + 1], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 3 + 2, dx), f[3 + 2], tol);

				Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 0, dx), f[6 + 0], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 1, dx), f[6 + 1], tol);
				Assert::AreEqual(numericalForce(sc, uv, x, k, 6 + 2, dx), f[6 + 2], tol);

				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 0, dx), dampingForces[0], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 1, dx), dampingForces[1], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 2, dx), dampingForces[2], tol);

				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 0, dx), dampingForces[3 + 0], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 1, dx), dampingForces[3 + 1], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 3 + 2, dx), dampingForces[3 + 2], tol);

				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 0, dx), dampingForces[6 + 0], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 1, dx), dampingForces[6 + 1], tol);
				Assert::AreEqual(numericalDampingForce(sc, uv, x, v, d, 6 + 2, dx), dampingForces[6 + 2], tol);

				for (int i = 0; i < x.size(); ++i)
				{
					Eigen::VectorXd fdN = dfdx.row(i);
					checkVectorEquality(fdN, numericalForceDerivative(sc, uv, x, v, k, d, i, dx), tol, true);
					fdN = dfdv.row(i);
					checkVectorEquality(fdN, numericalDampingForceDerivative(sc, uv, x, v, k, d, i, dx), tol, true);
				}
			}
		}

	};
}
