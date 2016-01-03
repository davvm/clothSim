#include "CppUnitTest.h"

#include "EqualityTests.h"

#include "StretchCondition.h"

#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{
	TEST_CLASS(StretchConditionTest)
	{
	public:

		void stretchTriangleQuantitiesForDerivative(Eigen::VectorXd x, Eigen::VectorXd uv, Eigen::VectorXd v, double bu, double bv, int idx, double dx, StretchCondition<double>::TriangleQuantities &qPlus, StretchCondition<double>::TriangleQuantities &qMinus) const
		{
			double component = x[idx];
			x[idx] = component + dx;
			qPlus = StretchCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
				);
			x[idx] = component - dx;
			qMinus = StretchCondition<double>::TriangleQuantities(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
				);
			x[idx] = component;
		}

		void checkStretchTriangleQuantities(const Eigen::VectorXd &x, const Eigen::VectorXd &uv, const Eigen::VectorXd &v, double bu, double bv, double dx, double tol)
		{
			StretchCondition<double>::TriangleQuantities q(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
				);

			Eigen::VectorXd xtest;
			xtest = x + v * dx;
			StretchCondition<double>::TriangleQuantities qtPlus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
				);
			xtest = x - v * dx;
			StretchCondition<double>::TriangleQuantities qtMinus(
				xtest.segment<3>(3 * 0), xtest.segment<3>(3 * 1), xtest.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
				);

			double dC0dt = (qtPlus.C0 - qtMinus.C0) / (2 * dx);
			double dC1dt = (qtPlus.C1 - qtMinus.C1) / (2 * dx);
			Assert::AreEqual(dC0dt, q.dC0dt, tol);
			Assert::AreEqual(dC1dt, q.dC1dt, tol);


			StretchCondition<double>::TriangleQuantities q0xPlus, q0xMinus, q0yPlus, q0yMinus, q0zPlus, q0zMinus;
			StretchCondition<double>::TriangleQuantities q1xPlus, q1xMinus, q1yPlus, q1yMinus, q1zPlus, q1zMinus;
			StretchCondition<double>::TriangleQuantities q2xPlus, q2xMinus, q2yPlus, q2yMinus, q2zPlus, q2zMinus;

			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 0, dx, q0xPlus, q0xMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 1, dx, q0yPlus, q0yMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 2, dx, q0zPlus, q0zMinus);

			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 3, dx, q1xPlus, q1xMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 4, dx, q1yPlus, q1yMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 5, dx, q1zPlus, q1zMinus);

			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 6, dx, q2xPlus, q2xMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 7, dx, q2yPlus, q2yMinus);
			stretchTriangleQuantitiesForDerivative(x, uv, v, bu, bv, 8, dx, q2zPlus, q2zMinus);

			// test first wu and wv derivatives:
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

			// test first C0 and C1 derivatives:
			Eigen::Vector3d dC0dP0;
			dC0dP0[0] = (q0xPlus.C0 - q0xMinus.C0) / (2 * dx);
			dC0dP0[1] = (q0yPlus.C0 - q0yMinus.C0) / (2 * dx);
			dC0dP0[2] = (q0zPlus.C0 - q0zMinus.C0) / (2 * dx);
			checkVectorEquality(dC0dP0, q.dC0dP0, tol);

			Eigen::Vector3d dC0dP1;
			dC0dP1[0] = (q1xPlus.C0 - q1xMinus.C0) / (2 * dx);
			dC0dP1[1] = (q1yPlus.C0 - q1yMinus.C0) / (2 * dx);
			dC0dP1[2] = (q1zPlus.C0 - q1zMinus.C0) / (2 * dx);
			checkVectorEquality(dC0dP1, q.dC0dP1, tol);

			Eigen::Vector3d dC0dP2;
			dC0dP2[0] = (q2xPlus.C0 - q2xMinus.C0) / (2 * dx);
			dC0dP2[1] = (q2yPlus.C0 - q2yMinus.C0) / (2 * dx);
			dC0dP2[2] = (q2zPlus.C0 - q2zMinus.C0) / (2 * dx);
			checkVectorEquality(dC0dP2, q.dC0dP2, tol);


			Eigen::Vector3d dC1dP0;
			dC1dP0[0] = (q0xPlus.C1 - q0xMinus.C1) / (2 * dx);
			dC1dP0[1] = (q0yPlus.C1 - q0yMinus.C1) / (2 * dx);
			dC1dP0[2] = (q0zPlus.C1 - q0zMinus.C1) / (2 * dx);
			checkVectorEquality(dC1dP0, q.dC1dP0, tol);

			Eigen::Vector3d dC1dP1;
			dC1dP1[0] = (q1xPlus.C1 - q1xMinus.C1) / (2 * dx);
			dC1dP1[1] = (q1yPlus.C1 - q1yMinus.C1) / (2 * dx);
			dC1dP1[2] = (q1zPlus.C1 - q1zMinus.C1) / (2 * dx);
			checkVectorEquality(dC1dP1, q.dC1dP1, tol);

			Eigen::Vector3d dC1dP2;
			dC1dP2[0] = (q2xPlus.C1 - q2xMinus.C1) / (2 * dx);
			dC1dP2[1] = (q2yPlus.C1 - q2yMinus.C1) / (2 * dx);
			dC1dP2[2] = (q2zPlus.C1 - q2zMinus.C1) / (2 * dx);
			checkVectorEquality(dC1dP2, q.dC1dP2, tol);


			// test second C0 and C1 derivatives:
			Eigen::Matrix3d d2C0dP0dP0;
			d2C0dP0dP0.col(0) = (q0xPlus.dC0dP0 - q0xMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP0.col(1) = (q0yPlus.dC0dP0 - q0yMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP0.col(2) = (q0zPlus.dC0dP0 - q0zMinus.dC0dP0) / (2 * dx);
			checkMatrixEquality(d2C0dP0dP0, q.d2C0dP0dP0, tol);

			Eigen::Matrix3d d2C0dP0dP1;
			d2C0dP0dP1.col(0) = (q1xPlus.dC0dP0 - q1xMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP1.col(1) = (q1yPlus.dC0dP0 - q1yMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP1.col(2) = (q1zPlus.dC0dP0 - q1zMinus.dC0dP0) / (2 * dx);
			checkMatrixEquality(d2C0dP0dP1, q.d2C0dP0dP1, tol);

			Eigen::Matrix3d d2C0dP0dP2;
			d2C0dP0dP2.col(0) = (q2xPlus.dC0dP0 - q2xMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP2.col(1) = (q2yPlus.dC0dP0 - q2yMinus.dC0dP0) / (2 * dx);
			d2C0dP0dP2.col(2) = (q2zPlus.dC0dP0 - q2zMinus.dC0dP0) / (2 * dx);
			checkMatrixEquality(d2C0dP0dP2, q.d2C0dP0dP2, tol);


			Eigen::Matrix3d d2C0dP1dP0;
			d2C0dP1dP0.col(0) = (q0xPlus.dC0dP1 - q0xMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP0.col(1) = (q0yPlus.dC0dP1 - q0yMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP0.col(2) = (q0zPlus.dC0dP1 - q0zMinus.dC0dP1) / (2 * dx);
			checkMatrixEquality(d2C0dP1dP0, q.d2C0dP1dP0, tol);

			Eigen::Matrix3d d2C0dP1dP1;
			d2C0dP1dP1.col(0) = (q1xPlus.dC0dP1 - q1xMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP1.col(1) = (q1yPlus.dC0dP1 - q1yMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP1.col(2) = (q1zPlus.dC0dP1 - q1zMinus.dC0dP1) / (2 * dx);
			checkMatrixEquality(d2C0dP1dP1, q.d2C0dP1dP1, tol);

			Eigen::Matrix3d d2C0dP1dP2;
			d2C0dP1dP2.col(0) = (q2xPlus.dC0dP1 - q2xMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP2.col(1) = (q2yPlus.dC0dP1 - q2yMinus.dC0dP1) / (2 * dx);
			d2C0dP1dP2.col(2) = (q2zPlus.dC0dP1 - q2zMinus.dC0dP1) / (2 * dx);
			checkMatrixEquality(d2C0dP1dP2, q.d2C0dP1dP2, tol);


			Eigen::Matrix3d d2C0dP2dP0;
			d2C0dP2dP0.col(0) = (q0xPlus.dC0dP2 - q0xMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP0.col(1) = (q0yPlus.dC0dP2 - q0yMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP0.col(2) = (q0zPlus.dC0dP2 - q0zMinus.dC0dP2) / (2 * dx);
			checkMatrixEquality(d2C0dP2dP0, q.d2C0dP2dP0, tol);

			Eigen::Matrix3d d2C0dP2dP1;
			d2C0dP2dP1.col(0) = (q1xPlus.dC0dP2 - q1xMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP1.col(1) = (q1yPlus.dC0dP2 - q1yMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP1.col(2) = (q1zPlus.dC0dP2 - q1zMinus.dC0dP2) / (2 * dx);
			checkMatrixEquality(d2C0dP2dP1, q.d2C0dP2dP1, tol);

			Eigen::Matrix3d d2C0dP2dP2;
			d2C0dP2dP2.col(0) = (q2xPlus.dC0dP2 - q2xMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP2.col(1) = (q2yPlus.dC0dP2 - q2yMinus.dC0dP2) / (2 * dx);
			d2C0dP2dP2.col(2) = (q2zPlus.dC0dP2 - q2zMinus.dC0dP2) / (2 * dx);
			checkMatrixEquality(d2C0dP2dP2, q.d2C0dP2dP2, tol);




			Eigen::Matrix3d d2C1dP0dP0;
			d2C1dP0dP0.col(0) = (q0xPlus.dC1dP0 - q0xMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP0.col(1) = (q0yPlus.dC1dP0 - q0yMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP0.col(2) = (q0zPlus.dC1dP0 - q0zMinus.dC1dP0) / (2 * dx);
			checkMatrixEquality(d2C1dP0dP0, q.d2C1dP0dP0, tol);

			Eigen::Matrix3d d2C1dP0dP1;
			d2C1dP0dP1.col(0) = (q1xPlus.dC1dP0 - q1xMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP1.col(1) = (q1yPlus.dC1dP0 - q1yMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP1.col(2) = (q1zPlus.dC1dP0 - q1zMinus.dC1dP0) / (2 * dx);
			checkMatrixEquality(d2C1dP0dP1, q.d2C1dP0dP1, tol);

			Eigen::Matrix3d d2C1dP0dP2;
			d2C1dP0dP2.col(0) = (q2xPlus.dC1dP0 - q2xMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP2.col(1) = (q2yPlus.dC1dP0 - q2yMinus.dC1dP0) / (2 * dx);
			d2C1dP0dP2.col(2) = (q2zPlus.dC1dP0 - q2zMinus.dC1dP0) / (2 * dx);
			checkMatrixEquality(d2C1dP0dP2, q.d2C1dP0dP2, tol);


			Eigen::Matrix3d d2C1dP1dP0;
			d2C1dP1dP0.col(0) = (q0xPlus.dC1dP1 - q0xMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP0.col(1) = (q0yPlus.dC1dP1 - q0yMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP0.col(2) = (q0zPlus.dC1dP1 - q0zMinus.dC1dP1) / (2 * dx);
			checkMatrixEquality(d2C1dP1dP0, q.d2C1dP1dP0, tol);

			Eigen::Matrix3d d2C1dP1dP1;
			d2C1dP1dP1.col(0) = (q1xPlus.dC1dP1 - q1xMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP1.col(1) = (q1yPlus.dC1dP1 - q1yMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP1.col(2) = (q1zPlus.dC1dP1 - q1zMinus.dC1dP1) / (2 * dx);
			checkMatrixEquality(d2C1dP1dP1, q.d2C1dP1dP1, tol);

			Eigen::Matrix3d d2C1dP1dP2;
			d2C1dP1dP2.col(0) = (q2xPlus.dC1dP1 - q2xMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP2.col(1) = (q2yPlus.dC1dP1 - q2yMinus.dC1dP1) / (2 * dx);
			d2C1dP1dP2.col(2) = (q2zPlus.dC1dP1 - q2zMinus.dC1dP1) / (2 * dx);
			checkMatrixEquality(d2C1dP1dP2, q.d2C1dP1dP2, tol);


			Eigen::Matrix3d d2C1dP2dP0;
			d2C1dP2dP0.col(0) = (q0xPlus.dC1dP2 - q0xMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP0.col(1) = (q0yPlus.dC1dP2 - q0yMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP0.col(2) = (q0zPlus.dC1dP2 - q0zMinus.dC1dP2) / (2 * dx);
			checkMatrixEquality(d2C1dP2dP0, q.d2C1dP2dP0, tol);

			Eigen::Matrix3d d2C1dP2dP1;
			d2C1dP2dP1.col(0) = (q1xPlus.dC1dP2 - q1xMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP1.col(1) = (q1yPlus.dC1dP2 - q1yMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP1.col(2) = (q1zPlus.dC1dP2 - q1zMinus.dC1dP2) / (2 * dx);
			checkMatrixEquality(d2C1dP2dP1, q.d2C1dP2dP1, tol);

			Eigen::Matrix3d d2C1dP2dP2;
			d2C1dP2dP2.col(0) = (q2xPlus.dC1dP2 - q2xMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP2.col(1) = (q2yPlus.dC1dP2 - q2yMinus.dC1dP2) / (2 * dx);
			d2C1dP2dP2.col(2) = (q2zPlus.dC1dP2 - q2zMinus.dC1dP2) / (2 * dx);
			checkMatrixEquality(d2C1dP2dP2, q.d2C1dP2dP2, tol);


		}


		TEST_METHOD(TestStretchConditionTriangleQuantities)
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

			double bu = 0.9;
			double bv = 1.1;

			StretchCondition<double>::TriangleQuantities q(
				x.segment<3>(3 * 0), x.segment<3>(3 * 1), x.segment<3>(3 * 2),
				uv.segment<2>(2 * 0), uv.segment<2>(2 * 1), uv.segment<2>(2 * 2),
				v.segment<3>(3 * 0), v.segment<3>(3 * 1), v.segment<3>(3 * 2),
				bu, bv
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

			checkStretchTriangleQuantities(x, uv, v, bu, bv, dx, tol);

			// test for 10 random configurations:
			for (int i = 0; i < 10; ++i)
			{
				// randomize the positions:
				x = Eigen::VectorXd::Random(3 * 3);
				uv = Eigen::VectorXd::Random(3 * 2);
				v = Eigen::VectorXd::Random(3 * 3);
				checkStretchTriangleQuantities(x, uv, v, bu, bv, dx, tol);
			}

		}


		TEST_METHOD(TestStretchCondition)
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
			StretchCondition<double> sc(0, 1, 2, bu, bv);

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
