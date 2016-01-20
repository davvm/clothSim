#include "CppUnitTest.h"

#include "test\EqualityTests.h"

#include "simLib\BendCondition.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace ClothSim;

namespace UnitTest1
{
	TEST_CLASS(BendConditionTest)
	{
	public:

		void bendTriangleQuantitiesForDerivative(Eigen::VectorXd x, Eigen::VectorXd v, int idx, double dx, BendCondition<double>::TriangleQuantities &qPlus, BendCondition<double>::TriangleQuantities &qMinus) const
		{
			double component = x[idx];
			x[idx] = component + dx;
			qPlus = BendCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2), x.segment<3>(3 * 3),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2), v.segment<3>(3 * 3)
				);
			x[idx] = component - dx;
			qMinus = BendCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2), x.segment<3>(3 * 3),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2), v.segment<3>(3 * 3)
				);
			x[idx] = component;
		}

		void checkBendTriangleQuantities(const Eigen::VectorXd &x, const Eigen::VectorXd &v, double dx, double tol)
		{
			BendCondition<double>::TriangleQuantities q(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2), x.segment<3>(3 * 3),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2), v.segment<3>(3 * 3)
				);

			Eigen::VectorXd xtest;
			xtest = x + v * dx;
			BendCondition<double>::TriangleQuantities qtPlus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2), xtest.segment<3>(3 * 3),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2), v.segment<3>(3 * 3)
				);
			xtest = x - v * dx;
			BendCondition<double>::TriangleQuantities qtMinus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2), xtest.segment<3>(3 * 3),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2), v.segment<3>(3 * 3)
				);

			double dThetadt = (qtPlus.theta - qtMinus.theta) / (2 * dx);
			Assert::AreEqual(dThetadt, q.dThetadt, tol);

			BendCondition<double>::TriangleQuantities q0xPlus, q0xMinus, q0yPlus, q0yMinus, q0zPlus, q0zMinus;
			BendCondition<double>::TriangleQuantities q1xPlus, q1xMinus, q1yPlus, q1yMinus, q1zPlus, q1zMinus;
			BendCondition<double>::TriangleQuantities q2xPlus, q2xMinus, q2yPlus, q2yMinus, q2zPlus, q2zMinus;
			BendCondition<double>::TriangleQuantities q3xPlus, q3xMinus, q3yPlus, q3yMinus, q3zPlus, q3zMinus;

			bendTriangleQuantitiesForDerivative(x, v, 0, dx, q0xPlus, q0xMinus);
			bendTriangleQuantitiesForDerivative(x, v, 1, dx, q0yPlus, q0yMinus);
			bendTriangleQuantitiesForDerivative(x, v, 2, dx, q0zPlus, q0zMinus);

			bendTriangleQuantitiesForDerivative(x, v, 3, dx, q1xPlus, q1xMinus);
			bendTriangleQuantitiesForDerivative(x, v, 4, dx, q1yPlus, q1yMinus);
			bendTriangleQuantitiesForDerivative(x, v, 5, dx, q1zPlus, q1zMinus);

			bendTriangleQuantitiesForDerivative(x, v, 6, dx, q2xPlus, q2xMinus);
			bendTriangleQuantitiesForDerivative(x, v, 7, dx, q2yPlus, q2yMinus);
			bendTriangleQuantitiesForDerivative(x, v, 8, dx, q2zPlus, q2zMinus);

			bendTriangleQuantitiesForDerivative(x, v, 9, dx, q3xPlus, q3xMinus);
			bendTriangleQuantitiesForDerivative(x, v, 10, dx, q3yPlus, q3yMinus);
			bendTriangleQuantitiesForDerivative(x, v, 11, dx, q3zPlus, q3zMinus);


			// test theta derivatives:
			BendCondition<double>::Vector3 dThetadP0;
			dThetadP0[0] = (q0xPlus.theta - q0xMinus.theta) / (2 * dx);
			dThetadP0[1] = (q0yPlus.theta - q0yMinus.theta) / (2 * dx);
			dThetadP0[2] = (q0zPlus.theta - q0zMinus.theta) / (2 * dx);
			checkVectorEquality(dThetadP0, q.dThetadP0, tol);

			BendCondition<double>::Vector3 dThetadP1;
			dThetadP1[0] = (q1xPlus.theta - q1xMinus.theta) / (2 * dx);
			dThetadP1[1] = (q1yPlus.theta - q1yMinus.theta) / (2 * dx);
			dThetadP1[2] = (q1zPlus.theta - q1zMinus.theta) / (2 * dx);
			checkVectorEquality(dThetadP1, q.dThetadP1, tol);

			BendCondition<double>::Vector3 dThetadP2;
			dThetadP2[0] = (q2xPlus.theta - q2xMinus.theta) / (2 * dx);
			dThetadP2[1] = (q2yPlus.theta - q2yMinus.theta) / (2 * dx);
			dThetadP2[2] = (q2zPlus.theta - q2zMinus.theta) / (2 * dx);
			checkVectorEquality(dThetadP2, q.dThetadP2, tol);

			BendCondition<double>::Vector3 dThetadP3;
			dThetadP3[0] = (q3xPlus.theta - q3xMinus.theta) / (2 * dx);
			dThetadP3[1] = (q3yPlus.theta - q3yMinus.theta) / (2 * dx);
			dThetadP3[2] = (q3zPlus.theta - q3zMinus.theta) / (2 * dx);
			checkVectorEquality(dThetadP3, q.dThetadP3, tol);


			// test normal derivatives:
			BendCondition<double>::Matrix3 dn0dP0N;
			dn0dP0N.col(0) = (q0xPlus.n0 - q0xMinus.n0) / (2 * dx);
			dn0dP0N.col(1) = (q0yPlus.n0 - q0yMinus.n0) / (2 * dx);
			dn0dP0N.col(2) = (q0zPlus.n0 - q0zMinus.n0) / (2 * dx);
			checkMatrixEquality(dn0dP0N, q.dn0dP0, tol);

			BendCondition<double>::Matrix3 dn0dP1N;
			dn0dP1N.col(0) = (q1xPlus.n0 - q1xMinus.n0) / (2 * dx);
			dn0dP1N.col(1) = (q1yPlus.n0 - q1yMinus.n0) / (2 * dx);
			dn0dP1N.col(2) = (q1zPlus.n0 - q1zMinus.n0) / (2 * dx);
			checkMatrixEquality(dn0dP1N, q.dn0dP1, tol);

			BendCondition<double>::Matrix3 dn0dP2N;
			dn0dP2N.col(0) = (q2xPlus.n0 - q2xMinus.n0) / (2 * dx);
			dn0dP2N.col(1) = (q2yPlus.n0 - q2yMinus.n0) / (2 * dx);
			dn0dP2N.col(2) = (q2zPlus.n0 - q2zMinus.n0) / (2 * dx);
			checkMatrixEquality(dn0dP2N, q.dn0dP2, tol);

			BendCondition<double>::Matrix3 dn0dP3N;
			dn0dP3N.col(0) = (q3xPlus.n0 - q3xMinus.n0) / (2 * dx);
			dn0dP3N.col(1) = (q3yPlus.n0 - q3yMinus.n0) / (2 * dx);
			dn0dP3N.col(2) = (q3zPlus.n0 - q3zMinus.n0) / (2 * dx);
			checkMatrixEquality(dn0dP3N, q.dn0dP3, tol);


			// test perpendicular distance derivatives:
			BendCondition<double>::Vector3 dd00dP0;
			dd00dP0[0] = (q0xPlus.d00 - q0xMinus.d00) / (2 * dx);
			dd00dP0[1] = (q0yPlus.d00 - q0yMinus.d00) / (2 * dx);
			dd00dP0[2] = (q0zPlus.d00 - q0zMinus.d00) / (2 * dx);
			checkVectorEquality(dd00dP0, q.dd00dP0, tol);

			BendCondition<double>::Vector3 dd00dP1;
			dd00dP1[0] = (q1xPlus.d00 - q1xMinus.d00) / (2 * dx);
			dd00dP1[1] = (q1yPlus.d00 - q1yMinus.d00) / (2 * dx);
			dd00dP1[2] = (q1zPlus.d00 - q1zMinus.d00) / (2 * dx);
			checkVectorEquality(dd00dP1, q.dd00dP1, tol);

			BendCondition<double>::Vector3 dd00dP2;
			dd00dP2[0] = (q2xPlus.d00 - q2xMinus.d00) / (2 * dx);
			dd00dP2[1] = (q2yPlus.d00 - q2yMinus.d00) / (2 * dx);
			dd00dP2[2] = (q2zPlus.d00 - q2zMinus.d00) / (2 * dx);
			checkVectorEquality(dd00dP2, q.dd00dP2, tol);

			BendCondition<double>::Vector3 dd00dP3;
			dd00dP3[0] = (q3xPlus.d00 - q3xMinus.d00) / (2 * dx);
			dd00dP3[1] = (q3yPlus.d00 - q3yMinus.d00) / (2 * dx);
			dd00dP3[2] = (q3zPlus.d00 - q3zMinus.d00) / (2 * dx);
			checkVectorEquality(dd00dP3, q.dd00dP3, tol);

			BendCondition<double>::Vector3 dd01dP0;
			dd01dP0[0] = (q0xPlus.d01 - q0xMinus.d01) / (2 * dx);
			dd01dP0[1] = (q0yPlus.d01 - q0yMinus.d01) / (2 * dx);
			dd01dP0[2] = (q0zPlus.d01 - q0zMinus.d01) / (2 * dx);
			checkVectorEquality(dd01dP0, q.dd01dP0, tol);

			BendCondition<double>::Vector3 dd01dP1;
			dd01dP1[0] = (q1xPlus.d01 - q1xMinus.d01) / (2 * dx);
			dd01dP1[1] = (q1yPlus.d01 - q1yMinus.d01) / (2 * dx);
			dd01dP1[2] = (q1zPlus.d01 - q1zMinus.d01) / (2 * dx);
			checkVectorEquality(dd01dP1, q.dd01dP1, tol);

			BendCondition<double>::Vector3 dd01dP2;
			dd01dP2[0] = (q2xPlus.d01 - q2xMinus.d01) / (2 * dx);
			dd01dP2[1] = (q2yPlus.d01 - q2yMinus.d01) / (2 * dx);
			dd01dP2[2] = (q2zPlus.d01 - q2zMinus.d01) / (2 * dx);
			checkVectorEquality(dd01dP2, q.dd01dP2, tol);

			BendCondition<double>::Vector3 dd01dP3;
			dd01dP3[0] = (q3xPlus.d01 - q3xMinus.d01) / (2 * dx);
			dd01dP3[1] = (q3yPlus.d01 - q3yMinus.d01) / (2 * dx);
			dd01dP3[2] = (q3zPlus.d01 - q3zMinus.d01) / (2 * dx);
			checkVectorEquality(dd01dP3, q.dd01dP3, tol);

			BendCondition<double>::Vector3 dd02dP0;
			dd02dP0[0] = (q0xPlus.d02 - q0xMinus.d02) / (2 * dx);
			dd02dP0[1] = (q0yPlus.d02 - q0yMinus.d02) / (2 * dx);
			dd02dP0[2] = (q0zPlus.d02 - q0zMinus.d02) / (2 * dx);
			checkVectorEquality(dd02dP0, q.dd02dP0, tol);

			BendCondition<double>::Vector3 dd02dP1;
			dd02dP1[0] = (q1xPlus.d02 - q1xMinus.d02) / (2 * dx);
			dd02dP1[1] = (q1yPlus.d02 - q1yMinus.d02) / (2 * dx);
			dd02dP1[2] = (q1zPlus.d02 - q1zMinus.d02) / (2 * dx);
			checkVectorEquality(dd02dP1, q.dd02dP1, tol);

			BendCondition<double>::Vector3 dd02dP2;
			dd02dP2[0] = (q2xPlus.d02 - q2xMinus.d02) / (2 * dx);
			dd02dP2[1] = (q2yPlus.d02 - q2yMinus.d02) / (2 * dx);
			dd02dP2[2] = (q2zPlus.d02 - q2zMinus.d02) / (2 * dx);
			checkVectorEquality(dd02dP2, q.dd02dP2, tol);

			BendCondition<double>::Vector3 dd02dP3;
			dd02dP3[0] = (q3xPlus.d02 - q3xMinus.d02) / (2 * dx);
			dd02dP3[1] = (q3yPlus.d02 - q3yMinus.d02) / (2 * dx);
			dd02dP3[2] = (q3zPlus.d02 - q3zMinus.d02) / (2 * dx);
			checkVectorEquality(dd02dP3, q.dd02dP3, tol);

			BendCondition<double>::Vector3 dd11dP0;
			dd11dP0[0] = (q0xPlus.d11 - q0xMinus.d11) / (2 * dx);
			dd11dP0[1] = (q0yPlus.d11 - q0yMinus.d11) / (2 * dx);
			dd11dP0[2] = (q0zPlus.d11 - q0zMinus.d11) / (2 * dx);
			checkVectorEquality(dd11dP0, q.dd11dP0, tol);

			BendCondition<double>::Vector3 dd11dP1;
			dd11dP1[0] = (q1xPlus.d11 - q1xMinus.d11) / (2 * dx);
			dd11dP1[1] = (q1yPlus.d11 - q1yMinus.d11) / (2 * dx);
			dd11dP1[2] = (q1zPlus.d11 - q1zMinus.d11) / (2 * dx);
			checkVectorEquality(dd11dP1, q.dd11dP1, tol);

			BendCondition<double>::Vector3 dd11dP2;
			dd11dP2[0] = (q2xPlus.d11 - q2xMinus.d11) / (2 * dx);
			dd11dP2[1] = (q2yPlus.d11 - q2yMinus.d11) / (2 * dx);
			dd11dP2[2] = (q2zPlus.d11 - q2zMinus.d11) / (2 * dx);
			checkVectorEquality(dd11dP2, q.dd11dP2, tol);

			BendCondition<double>::Vector3 dd11dP3;
			dd11dP3[0] = (q3xPlus.d11 - q3xMinus.d11) / (2 * dx);
			dd11dP3[1] = (q3yPlus.d11 - q3yMinus.d11) / (2 * dx);
			dd11dP3[2] = (q3zPlus.d11 - q3zMinus.d11) / (2 * dx);
			checkVectorEquality(dd11dP3, q.dd11dP3, tol);

			BendCondition<double>::Vector3 dd12dP0;
			dd12dP0[0] = (q0xPlus.d12 - q0xMinus.d12) / (2 * dx);
			dd12dP0[1] = (q0yPlus.d12 - q0yMinus.d12) / (2 * dx);
			dd12dP0[2] = (q0zPlus.d12 - q0zMinus.d12) / (2 * dx);
			checkVectorEquality(dd12dP0, q.dd12dP0, tol);

			BendCondition<double>::Vector3 dd12dP1;
			dd12dP1[0] = (q1xPlus.d12 - q1xMinus.d12) / (2 * dx);
			dd12dP1[1] = (q1yPlus.d12 - q1yMinus.d12) / (2 * dx);
			dd12dP1[2] = (q1zPlus.d12 - q1zMinus.d12) / (2 * dx);
			checkVectorEquality(dd12dP1, q.dd12dP1, tol);

			BendCondition<double>::Vector3 dd12dP2;
			dd12dP2[0] = (q2xPlus.d12 - q2xMinus.d12) / (2 * dx);
			dd12dP2[1] = (q2yPlus.d12 - q2yMinus.d12) / (2 * dx);
			dd12dP2[2] = (q2zPlus.d12 - q2zMinus.d12) / (2 * dx);
			checkVectorEquality(dd12dP2, q.dd12dP2, tol);

			BendCondition<double>::Vector3 dd12dP3;
			dd12dP3[0] = (q3xPlus.d12 - q3xMinus.d12) / (2 * dx);
			dd12dP3[1] = (q3yPlus.d12 - q3yMinus.d12) / (2 * dx);
			dd12dP3[2] = (q3zPlus.d12 - q3zMinus.d12) / (2 * dx);
			checkVectorEquality(dd12dP3, q.dd12dP3, tol);

			BendCondition<double>::Vector3 dd13dP0;
			dd13dP0[0] = (q0xPlus.d13 - q0xMinus.d13) / (2 * dx);
			dd13dP0[1] = (q0yPlus.d13 - q0yMinus.d13) / (2 * dx);
			dd13dP0[2] = (q0zPlus.d13 - q0zMinus.d13) / (2 * dx);
			checkVectorEquality(dd13dP0, q.dd13dP0, tol);

			BendCondition<double>::Vector3 dd13dP1;
			dd13dP1[0] = (q1xPlus.d13 - q1xMinus.d13) / (2 * dx);
			dd13dP1[1] = (q1yPlus.d13 - q1yMinus.d13) / (2 * dx);
			dd13dP1[2] = (q1zPlus.d13 - q1zMinus.d13) / (2 * dx);
			checkVectorEquality(dd13dP1, q.dd13dP1, tol);

			BendCondition<double>::Vector3 dd13dP2;
			dd13dP2[0] = (q2xPlus.d13 - q2xMinus.d13) / (2 * dx);
			dd13dP2[1] = (q2yPlus.d13 - q2yMinus.d13) / (2 * dx);
			dd13dP2[2] = (q2zPlus.d13 - q2zMinus.d13) / (2 * dx);
			checkVectorEquality(dd13dP2, q.dd13dP2, tol);

			BendCondition<double>::Vector3 dd13dP3;
			dd13dP3[0] = (q3xPlus.d13 - q3xMinus.d13) / (2 * dx);
			dd13dP3[1] = (q3yPlus.d13 - q3yMinus.d13) / (2 * dx);
			dd13dP3[2] = (q3zPlus.d13 - q3zMinus.d13) / (2 * dx);
			checkVectorEquality(dd13dP3, q.dd13dP3, tol);


			// test cosine derivatives:
			BendCondition<double>::Vector3 dc01dP0;
			dc01dP0[0] = (q0xPlus.c01 - q0xMinus.c01) / (2 * dx);
			dc01dP0[1] = (q0yPlus.c01 - q0yMinus.c01) / (2 * dx);
			dc01dP0[2] = (q0zPlus.c01 - q0zMinus.c01) / (2 * dx);
			checkVectorEquality(dc01dP0, q.dc01dP0, tol);

			BendCondition<double>::Vector3 dc01dP1;
			dc01dP1[0] = (q1xPlus.c01 - q1xMinus.c01) / (2 * dx);
			dc01dP1[1] = (q1yPlus.c01 - q1yMinus.c01) / (2 * dx);
			dc01dP1[2] = (q1zPlus.c01 - q1zMinus.c01) / (2 * dx);
			checkVectorEquality(dc01dP1, q.dc01dP1, tol);

			BendCondition<double>::Vector3 dc01dP2;
			dc01dP2[0] = (q2xPlus.c01 - q2xMinus.c01) / (2 * dx);
			dc01dP2[1] = (q2yPlus.c01 - q2yMinus.c01) / (2 * dx);
			dc01dP2[2] = (q2zPlus.c01 - q2zMinus.c01) / (2 * dx);
			checkVectorEquality(dc01dP2, q.dc01dP2, tol);

			BendCondition<double>::Vector3 dc01dP3;
			dc01dP3[0] = (q3xPlus.c01 - q3xMinus.c01) / (2 * dx);
			dc01dP3[1] = (q3yPlus.c01 - q3yMinus.c01) / (2 * dx);
			dc01dP3[2] = (q3zPlus.c01 - q3zMinus.c01) / (2 * dx);
			checkVectorEquality(dc01dP3, q.dc01dP3, tol);

			BendCondition<double>::Vector3 dc02dP0;
			dc02dP0[0] = (q0xPlus.c02 - q0xMinus.c02) / (2 * dx);
			dc02dP0[1] = (q0yPlus.c02 - q0yMinus.c02) / (2 * dx);
			dc02dP0[2] = (q0zPlus.c02 - q0zMinus.c02) / (2 * dx);
			checkVectorEquality(dc02dP0, q.dc02dP0, tol);

			BendCondition<double>::Vector3 dc02dP1;
			dc02dP1[0] = (q1xPlus.c02 - q1xMinus.c02) / (2 * dx);
			dc02dP1[1] = (q1yPlus.c02 - q1yMinus.c02) / (2 * dx);
			dc02dP1[2] = (q1zPlus.c02 - q1zMinus.c02) / (2 * dx);
			checkVectorEquality(dc02dP1, q.dc02dP1, tol);

			BendCondition<double>::Vector3 dc02dP2;
			dc02dP2[0] = (q2xPlus.c02 - q2xMinus.c02) / (2 * dx);
			dc02dP2[1] = (q2yPlus.c02 - q2yMinus.c02) / (2 * dx);
			dc02dP2[2] = (q2zPlus.c02 - q2zMinus.c02) / (2 * dx);
			checkVectorEquality(dc02dP2, q.dc02dP2, tol);

			BendCondition<double>::Vector3 dc02dP3;
			dc02dP3[0] = (q3xPlus.c02 - q3xMinus.c02) / (2 * dx);
			dc02dP3[1] = (q3yPlus.c02 - q3yMinus.c02) / (2 * dx);
			dc02dP3[2] = (q3zPlus.c02 - q3zMinus.c02) / (2 * dx);
			checkVectorEquality(dc02dP3, q.dc02dP3, tol);

			BendCondition<double>::Vector3 dc11dP0;
			dc11dP0[0] = (q0xPlus.c11 - q0xMinus.c11) / (2 * dx);
			dc11dP0[1] = (q0yPlus.c11 - q0yMinus.c11) / (2 * dx);
			dc11dP0[2] = (q0zPlus.c11 - q0zMinus.c11) / (2 * dx);
			checkVectorEquality(dc11dP0, q.dc11dP0, tol);

			BendCondition<double>::Vector3 dc11dP1;
			dc11dP1[0] = (q1xPlus.c11 - q1xMinus.c11) / (2 * dx);
			dc11dP1[1] = (q1yPlus.c11 - q1yMinus.c11) / (2 * dx);
			dc11dP1[2] = (q1zPlus.c11 - q1zMinus.c11) / (2 * dx);
			checkVectorEquality(dc11dP1, q.dc11dP1, tol);

			BendCondition<double>::Vector3 dc11dP2;
			dc11dP2[0] = (q2xPlus.c11 - q2xMinus.c11) / (2 * dx);
			dc11dP2[1] = (q2yPlus.c11 - q2yMinus.c11) / (2 * dx);
			dc11dP2[2] = (q2zPlus.c11 - q2zMinus.c11) / (2 * dx);
			checkVectorEquality(dc11dP2, q.dc11dP2, tol);

			BendCondition<double>::Vector3 dc11dP3;
			dc11dP3[0] = (q3xPlus.c11 - q3xMinus.c11) / (2 * dx);
			dc11dP3[1] = (q3yPlus.c11 - q3yMinus.c11) / (2 * dx);
			dc11dP3[2] = (q3zPlus.c11 - q3zMinus.c11) / (2 * dx);
			checkVectorEquality(dc11dP3, q.dc11dP3, tol);

			BendCondition<double>::Vector3 dc12dP0;
			dc12dP0[0] = (q0xPlus.c12 - q0xMinus.c12) / (2 * dx);
			dc12dP0[1] = (q0yPlus.c12 - q0yMinus.c12) / (2 * dx);
			dc12dP0[2] = (q0zPlus.c12 - q0zMinus.c12) / (2 * dx);
			checkVectorEquality(dc12dP0, q.dc12dP0, tol);

			BendCondition<double>::Vector3 dc12dP1;
			dc12dP1[0] = (q1xPlus.c12 - q1xMinus.c12) / (2 * dx);
			dc12dP1[1] = (q1yPlus.c12 - q1yMinus.c12) / (2 * dx);
			dc12dP1[2] = (q1zPlus.c12 - q1zMinus.c12) / (2 * dx);
			checkVectorEquality(dc12dP1, q.dc12dP1, tol);

			BendCondition<double>::Vector3 dc12dP2;
			dc12dP2[0] = (q2xPlus.c12 - q2xMinus.c12) / (2 * dx);
			dc12dP2[1] = (q2yPlus.c12 - q2yMinus.c12) / (2 * dx);
			dc12dP2[2] = (q2zPlus.c12 - q2zMinus.c12) / (2 * dx);
			checkVectorEquality(dc12dP2, q.dc12dP2, tol);

			BendCondition<double>::Vector3 dc12dP3;
			dc12dP3[0] = (q3xPlus.c12 - q3xMinus.c12) / (2 * dx);
			dc12dP3[1] = (q3yPlus.c12 - q3yMinus.c12) / (2 * dx);
			dc12dP3[2] = (q3zPlus.c12 - q3zMinus.c12) / (2 * dx);
			checkVectorEquality(dc12dP3, q.dc12dP3, tol);


			// second derivatives of theta:
			BendCondition<double>::Matrix3 d2ThetadP0dP0;
			d2ThetadP0dP0.col(0) = (q0xPlus.dThetadP0 - q0xMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP0.col(1) = (q0yPlus.dThetadP0 - q0yMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP0.col(2) = (q0zPlus.dThetadP0 - q0zMinus.dThetadP0) / (2 * dx);
			checkMatrixEquality(d2ThetadP0dP0, q.d2ThetadP0dP0, tol);

			BendCondition<double>::Matrix3 d2ThetadP0dP1;
			d2ThetadP0dP1.col(0) = (q1xPlus.dThetadP0 - q1xMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP1.col(1) = (q1yPlus.dThetadP0 - q1yMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP1.col(2) = (q1zPlus.dThetadP0 - q1zMinus.dThetadP0) / (2 * dx);
			checkMatrixEquality(d2ThetadP0dP1, q.d2ThetadP0dP1, tol);

			BendCondition<double>::Matrix3 d2ThetadP0dP2;
			d2ThetadP0dP2.col(0) = (q2xPlus.dThetadP0 - q2xMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP2.col(1) = (q2yPlus.dThetadP0 - q2yMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP2.col(2) = (q2zPlus.dThetadP0 - q2zMinus.dThetadP0) / (2 * dx);
			checkMatrixEquality(d2ThetadP0dP2, q.d2ThetadP0dP2, tol);

			BendCondition<double>::Matrix3 d2ThetadP0dP3;
			d2ThetadP0dP3.col(0) = (q3xPlus.dThetadP0 - q3xMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP3.col(1) = (q3yPlus.dThetadP0 - q3yMinus.dThetadP0) / (2 * dx);
			d2ThetadP0dP3.col(2) = (q3zPlus.dThetadP0 - q3zMinus.dThetadP0) / (2 * dx);
			checkMatrixEquality(d2ThetadP0dP3, q.d2ThetadP0dP3, tol);

			BendCondition<double>::Matrix3 d2ThetadP1dP0;
			d2ThetadP1dP0.col(0) = (q0xPlus.dThetadP1 - q0xMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP0.col(1) = (q0yPlus.dThetadP1 - q0yMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP0.col(2) = (q0zPlus.dThetadP1 - q0zMinus.dThetadP1) / (2 * dx);
			checkMatrixEquality(d2ThetadP1dP0, q.d2ThetadP1dP0, tol);

			BendCondition<double>::Matrix3 d2ThetadP1dP1;
			d2ThetadP1dP1.col(0) = (q1xPlus.dThetadP1 - q1xMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP1.col(1) = (q1yPlus.dThetadP1 - q1yMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP1.col(2) = (q1zPlus.dThetadP1 - q1zMinus.dThetadP1) / (2 * dx);
			checkMatrixEquality(d2ThetadP1dP1, q.d2ThetadP1dP1, tol);

			BendCondition<double>::Matrix3 d2ThetadP1dP2;
			d2ThetadP1dP2.col(0) = (q2xPlus.dThetadP1 - q2xMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP2.col(1) = (q2yPlus.dThetadP1 - q2yMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP2.col(2) = (q2zPlus.dThetadP1 - q2zMinus.dThetadP1) / (2 * dx);
			checkMatrixEquality(d2ThetadP1dP2, q.d2ThetadP1dP2, tol);

			BendCondition<double>::Matrix3 d2ThetadP1dP3;
			d2ThetadP1dP3.col(0) = (q3xPlus.dThetadP1 - q3xMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP3.col(1) = (q3yPlus.dThetadP1 - q3yMinus.dThetadP1) / (2 * dx);
			d2ThetadP1dP3.col(2) = (q3zPlus.dThetadP1 - q3zMinus.dThetadP1) / (2 * dx);
			checkMatrixEquality(d2ThetadP1dP3, q.d2ThetadP1dP3, tol);

			BendCondition<double>::Matrix3 d2ThetadP2dP0;
			d2ThetadP2dP0.col(0) = (q0xPlus.dThetadP2 - q0xMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP0.col(1) = (q0yPlus.dThetadP2 - q0yMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP0.col(2) = (q0zPlus.dThetadP2 - q0zMinus.dThetadP2) / (2 * dx);
			checkMatrixEquality(d2ThetadP2dP0, q.d2ThetadP2dP0, tol);

			BendCondition<double>::Matrix3 d2ThetadP2dP1;
			d2ThetadP2dP1.col(0) = (q1xPlus.dThetadP2 - q1xMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP1.col(1) = (q1yPlus.dThetadP2 - q1yMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP1.col(2) = (q1zPlus.dThetadP2 - q1zMinus.dThetadP2) / (2 * dx);
			checkMatrixEquality(d2ThetadP2dP1, q.d2ThetadP2dP1, tol);

			BendCondition<double>::Matrix3 d2ThetadP2dP2;
			d2ThetadP2dP2.col(0) = (q2xPlus.dThetadP2 - q2xMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP2.col(1) = (q2yPlus.dThetadP2 - q2yMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP2.col(2) = (q2zPlus.dThetadP2 - q2zMinus.dThetadP2) / (2 * dx);
			checkMatrixEquality(d2ThetadP2dP2, q.d2ThetadP2dP2, tol);

			BendCondition<double>::Matrix3 d2ThetadP2dP3;
			d2ThetadP2dP3.col(0) = (q3xPlus.dThetadP2 - q3xMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP3.col(1) = (q3yPlus.dThetadP2 - q3yMinus.dThetadP2) / (2 * dx);
			d2ThetadP2dP3.col(2) = (q3zPlus.dThetadP2 - q3zMinus.dThetadP2) / (2 * dx);
			checkMatrixEquality(d2ThetadP2dP3, q.d2ThetadP2dP3, tol);

			BendCondition<double>::Matrix3 d2ThetadP3dP0;
			d2ThetadP3dP0.col(0) = (q0xPlus.dThetadP3 - q0xMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP0.col(1) = (q0yPlus.dThetadP3 - q0yMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP0.col(2) = (q0zPlus.dThetadP3 - q0zMinus.dThetadP3) / (2 * dx);
			checkMatrixEquality(d2ThetadP3dP0, q.d2ThetadP3dP0, tol);

			BendCondition<double>::Matrix3 d2ThetadP3dP1;
			d2ThetadP3dP1.col(0) = (q1xPlus.dThetadP3 - q1xMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP1.col(1) = (q1yPlus.dThetadP3 - q1yMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP1.col(2) = (q1zPlus.dThetadP3 - q1zMinus.dThetadP3) / (2 * dx);
			checkMatrixEquality(d2ThetadP3dP1, q.d2ThetadP3dP1, tol);

			BendCondition<double>::Matrix3 d2ThetadP3dP2;
			d2ThetadP3dP2.col(0) = (q2xPlus.dThetadP3 - q2xMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP2.col(1) = (q2yPlus.dThetadP3 - q2yMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP2.col(2) = (q2zPlus.dThetadP3 - q2zMinus.dThetadP3) / (2 * dx);
			checkMatrixEquality(d2ThetadP3dP2, q.d2ThetadP3dP2, tol);

			BendCondition<double>::Matrix3 d2ThetadP3dP3;
			d2ThetadP3dP3.col(0) = (q3xPlus.dThetadP3 - q3xMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP3.col(1) = (q3yPlus.dThetadP3 - q3yMinus.dThetadP3) / (2 * dx);
			d2ThetadP3dP3.col(2) = (q3zPlus.dThetadP3 - q3zMinus.dThetadP3) / (2 * dx);
			checkMatrixEquality(d2ThetadP3dP3, q.d2ThetadP3dP3, tol);
		}

		TEST_METHOD(TestBendConditionTriangleQuantities)
		{
			Eigen::VectorXd v(5 * 3);

			v[3 * 0] = 0;
			v[3 * 0 + 1] = 0;
			v[3 * 0 + 2] = -0.4;

			v[3 * 1] = 0;
			v[3 * 1 + 1] = 0;
			v[3 * 1 + 2] = 0.1;

			v[3 * 2] = 0;
			v[3 * 2 + 1] = 0;
			v[3 * 2 + 2] = 0.1;

			v[3 * 3] = 1.0;
			v[3 * 3 + 1] = 0;
			v[3 * 3 + 2] = -0.2;

			v[3 * 4] = 0.0;
			v[3 * 4 + 1] = 0;
			v[3 * 4 + 2] = 0;

			Eigen::VectorXd x(5 * 3);

			x[3 * 0] = -2.0;
			x[3 * 0 + 1] = 0;
			x[3 * 0 + 2] = 0.4;

			x[3 * 1] = 0;
			x[3 * 1 + 1] = 1.0;
			x[3 * 1 + 2] = 0;

			x[3 * 2] = 0;
			x[3 * 2 + 1] = -0.5;
			x[3 * 2 + 2] = 0;

			x[3 * 3] = 1.0;
			x[3 * 3 + 1] = 0;
			x[3 * 3 + 2] = 0.2;

			x[3 * 4] = 0.0;
			x[3 * 4 + 1] = 0;
			x[3 * 4 + 2] = 0;

			double dx = 0.0001;
			double tol = 1.e-5;

			// test with positive angle:
			checkBendTriangleQuantities(x, v, dx, tol);

			// test with negative angle:
			x[3 * 0 + 2] = -0.4;
			x[3 * 3 + 2] = -0.2;
			checkBendTriangleQuantities(x, v, dx, tol);

			// test for 10 random configurations:
			for (int i = 0; i < 10; ++i)
			{
				// randomize the positions:
				x = Eigen::VectorXd::Random(5 * 3);
				v = Eigen::VectorXd::Random(5 * 3);
				checkBendTriangleQuantities(x, v, dx, tol);
			}
		}

		TEST_METHOD(TestBendCondition)
		{

			Eigen::VectorXd v(5 * 3);

			v[3 * 0] = 0;
			v[3 * 0 + 1] = 0;
			v[3 * 0 + 2] = -0.4;

			v[3 * 1] = 0;
			v[3 * 1 + 1] = 0;
			v[3 * 1 + 2] = 0.1;

			v[3 * 2] = 0;
			v[3 * 2 + 1] = 0;
			v[3 * 2 + 2] = 0.1;

			v[3 * 3] = 1.0;
			v[3 * 3 + 1] = 0;
			v[3 * 3 + 2] = -0.2;

			v[3 * 4] = 0.0;
			v[3 * 4 + 1] = 0;
			v[3 * 4 + 2] = 0;

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

			uv[2 * 0] = -2.0;
			uv[2 * 0 + 1] = 0;

			uv[2 * 1] = 0;
			uv[2 * 1 + 1] = 1.0;

			uv[2 * 2] = 0;
			uv[2 * 2 + 1] = -1.0;

			uv[2 * 3] = 1.0;
			uv[2 * 3 + 1] = 0;

			uv[2 * 4] = 0.0;
			uv[2 * 4 + 1] = 0;

			Eigen::VectorXd f(5 * 3);
			f.setConstant(0);
			Eigen::VectorXd dampingForces(5 * 3);
			dampingForces.setConstant(0);
			BendCondition<double> bc(0, 1, 2, 3);

			double k = 1.5;
			Eigen::SparseMatrix<double> dfdx(5 * 3, 5 * 3);
			double d = 1.0;
			Eigen::SparseMatrix<double> dampingPseudoDerivatives(5 * 3, 5 * 3);
			Eigen::SparseMatrix<double> dfdv(5 * 3, 5 * 3);

			// check we've got the energy condition right, and we're measuring the angle between the normals:
			Eigen::VectorXd c;
			c = bc.C(x, uv);

			// sqrt(2) comes from the length of the bend edge:
			Assert::AreEqual(c[0], sqrt(2.0) * 2 * atan(0.2), 1.e-4);

			// compute forces analytically:
			bc.computeForces(x, uv, k, f, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

			// sanity check: forces should be opposing the bend, so the ones on x0,x3 should be pointing down,
			// the ones on x1,x2 should be pointing up:
			Assert::IsTrue(f[3 * 0 + 2] < 0);
			Assert::IsTrue(f[3 * 3 + 2] < 0);
			Assert::IsTrue(f[3 * 1 + 2] > 0);
			Assert::IsTrue(f[3 * 2 + 2] > 0);

			// you'd expect the damping forces to oppose the motion:
			Assert::IsTrue(dampingForces[3 * 0 + 2] > 0);
			Assert::IsTrue(dampingForces[3 * 3 + 2] > 0);
			Assert::IsTrue(dampingForces[3 * 1 + 2] < 0);
			Assert::IsTrue(dampingForces[3 * 2 + 2] < 0);

			// compare to numerically computed forces:
			double dx = 0.0001;
			double tol = 1.e-3;
			for (int i = 0; i < x.size(); ++i)
			{
				Assert::AreEqual(numericalForce(bc, uv, x, k, i, dx), f[i], tol);
				Assert::AreEqual(numericalDampingForce(bc, uv, x, v, d, i, dx), dampingForces[i], tol);
			}

			for (int i = 0; i < x.size(); ++i)
			{
				Eigen::VectorXd fdN = dfdx.row(i);
				checkVectorEquality(fdN, numericalForceDerivative(bc, uv, x, v, k, d, i, dx), tol, true);
			}

			checkPseudoDerivatives(dampingPseudoDerivatives, bc, uv, x, v, k, d, dx*0.1, tol);

			// turn the angle the other way to make sure we're handling bends in both directions:
			x[3 * 0 + 2] = -0.4;
			x[3 * 3 + 2] = -0.2;

			// test we're measuring a negative angle:
			c = bc.C(x, uv);
			Assert::AreEqual(c[0], -sqrt(2) * 2 * atan(0.2), 1.e-4);

			// recompute the forces and check the derivatives:
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
			bc.computeForces(x, uv, k, f, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

			// sanity check: forces should be opposing the bend, so now the ones on x0,x3 should be pointing up,
			// and the ones on x1,x2 should be pointing down:
			Assert::IsTrue(f[3 * 0 + 2] > 0);
			Assert::IsTrue(f[3 * 3 + 2] > 0);
			Assert::IsTrue(f[3 * 1 + 2] < 0);
			Assert::IsTrue(f[3 * 2 + 2] < 0);

			// you'd expect the damping forces to oppose the motion:
			Assert::IsTrue(dampingForces[3 * 0 + 2] > 0);
			Assert::IsTrue(dampingForces[3 * 3 + 2] > 0);
			Assert::IsTrue(dampingForces[3 * 1 + 2] < 0);
			Assert::IsTrue(dampingForces[3 * 2 + 2] < 0);

			for (int i = 0; i < x.size(); ++i)
			{
				Assert::AreEqual(numericalForce(bc, uv, x, k, i, dx), f[i], tol);
				Assert::AreEqual(numericalDampingForce(bc, uv, x, v, d, i, dx), dampingForces[i], tol);
			}


			for (int i = 0; i < x.size(); ++i)
			{
				Eigen::VectorXd fdN = dfdx.row(i);
				checkVectorEquality(fdN, numericalForceDerivative(bc, uv, x, v, k, d, i, dx), tol, true);
				fdN = dfdv.row(i);
				checkVectorEquality(fdN, numericalDampingForceDerivative(bc, uv, x, v, k, d, i, dx), tol, true);
			}

			// check the terms that are almost dfd/dx but not quite:
			checkPseudoDerivatives(dampingPseudoDerivatives, bc, uv, x, v, k, d, dx*0.1, tol);

			// test on 10 randomized configurations:
			for (size_t n = 0; n < 10; ++n)
			{
				// randomize positions:
				x = Eigen::VectorXd::Random(5 * 3);
				v = Eigen::VectorXd::Random(5 * 3);
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
				bc.computeForces(x, uv, k, f, dfdx, v, d, dampingForces, dampingPseudoDerivatives, dfdv);

				for (int i = 0; i < x.size(); ++i)
				{
					Assert::AreEqual(numericalForce(bc, uv, x, k, i, dx), f[i], tol);
					Assert::AreEqual(numericalDampingForce(bc, uv, x, v, d, i, dx), dampingForces[i], tol);
				}

				for (int i = 0; i < x.size(); ++i)
				{
					Eigen::VectorXd fdN = dfdx.row(i);
					checkVectorEquality(fdN, numericalForceDerivative(bc, uv, x, v, k, d, i, dx), tol, true);
					fdN = dfdv.row(i);
					checkVectorEquality(fdN, numericalDampingForceDerivative(bc, uv, x, v, k, d, i, dx), tol, true);
				}

				// check the terms that are almost dfd/dx but not quite:
				checkPseudoDerivatives(dampingPseudoDerivatives, bc, uv, x, v, k, d, dx*0.1, tol);
			}
		}

	};
}
