#include "CppUnitTest.h"

#include "EqualityTests.h"

#include "ClothMesh.h"

#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTest1
{
	TEST_CLASS(ClothMeshTest)
	{
	public:

		TEST_METHOD(TestSmallClothMesh)
		{
			Eigen::VectorXd v(4 * 3);

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

			Eigen::VectorXd x(4 * 3);

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

			Eigen::VectorXd uv(4 * 2);
			uv[2 * 0 + 0] = -2.0;
			uv[2 * 0 + 1] = 0.0;

			uv[2 * 1 + 0] = 0.0;
			uv[2 * 1 + 1] = 1.0;

			uv[2 * 2 + 0] = 0.0;
			uv[2 * 2 + 1] = -1.0;

			uv[2 * 3 + 0] = 1.0;
			uv[2 * 3 + 1] = 0.0;

			std::vector<int> triangleIndices;
			triangleIndices.push_back(0);
			triangleIndices.push_back(1);
			triangleIndices.push_back(2);

			triangleIndices.push_back(2);
			triangleIndices.push_back(1);
			triangleIndices.push_back(3);

			ClothMesh<double> m(
				x, v, uv, triangleIndices,
				1.0, 1.0, 1.0,
				1.0, 1.0, 1.0,
				1.0
			);

			Assert::IsTrue(m.x() == x);
			Assert::IsTrue(m.v() == v);

			Assert::AreEqual((int)m.m().size(), 4);

			Assert::AreEqual(m.m()[0], 2.0 / 3.0);
			Assert::AreEqual(m.m()[1], 1.0);
			Assert::AreEqual(m.m()[2], 1.0);
			Assert::AreEqual(m.m()[3], 1.0 / 3.0);

			Assert::AreEqual((int)m.bendConditions().size(), 1);
			Assert::AreEqual((int)m.shearConditions().size(), 2);
			Assert::AreEqual((int)m.stretchConditions().size(), 2);

			Assert::AreEqual(m.shearConditions()[0].inds()[0], 0);
			Assert::AreEqual(m.shearConditions()[0].inds()[1], 1);
			Assert::AreEqual(m.shearConditions()[0].inds()[2], 2);

			Assert::AreEqual(m.shearConditions()[1].inds()[0], 2);
			Assert::AreEqual(m.shearConditions()[1].inds()[1], 1);
			Assert::AreEqual(m.shearConditions()[1].inds()[2], 3);

			Assert::AreEqual(m.stretchConditions()[0].inds()[0], 0);
			Assert::AreEqual(m.stretchConditions()[0].inds()[1], 1);
			Assert::AreEqual(m.stretchConditions()[0].inds()[2], 2);

			Assert::AreEqual(m.stretchConditions()[1].inds()[0], 2);
			Assert::AreEqual(m.stretchConditions()[1].inds()[1], 1);
			Assert::AreEqual(m.stretchConditions()[1].inds()[2], 3);

			Assert::AreEqual(m.bendConditions()[0].inds()[0], 0);
			Assert::AreEqual(m.bendConditions()[0].inds()[1], 1);
			Assert::AreEqual(m.bendConditions()[0].inds()[2], 2);
			Assert::AreEqual(m.bendConditions()[0].inds()[3], 3);

		}

		TEST_METHOD(TestLargerClothMesh)
		{
			Eigen::VectorXd v(16 * 3);
			Eigen::VectorXd x(16 * 3);
			Eigen::VectorXd uv(16 * 2);

			for (int i = 0; i <= 3; ++i)
			{
				for (int j = 0; j <= 3; ++j)
				{
					x[3 * i + 0] = uv[2 * i + 0] = i;
					x[3 * i + 1] = uv[2 * i + 1] = j;
					x[3 * i + 2] = 0;

					v[3 * i + 0] = v[3 * i + 1] = v[3 * i + 2] = 0;
				}
			}

			std::vector<int> triangleInds;
			int shuffledI[] = { 2, 0, 1 };
			int shuffledJ[] = { 1, 2, 0 };
			for (int i = 0; i < 3; ++i)
			{
				int si = shuffledI[i];
				for (int j = 0; j < 3; ++j)
				{
					int base = shuffledI[i] + 4 * shuffledJ[j];

					triangleInds.push_back(base + 0);
					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + 4);

					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + 5);
					triangleInds.push_back(base + 4);

				}
			}

			ClothMesh<double> m(
				x, v, uv, triangleInds,
				1.0, 1.0, 1.0,
				1.0, 1.0, 1.0,
				1.0
			);

			Assert::AreEqual((int)m.stretchConditions().size(), 18);
			Assert::AreEqual((int)m.shearConditions().size(), 18);
			Assert::AreEqual((int)m.bendConditions().size(), 9 + 6 + 6);

		}
	};
}
