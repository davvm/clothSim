#include "CppUnitTest.h"

#include "test\EqualityTests.h"

#include "simLib\ClothMesh.h"
#include "simLib\DirectSolver.h"

#include <set>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace ClothSim;

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

			for (size_t i = 0; i < m.stretchConditions().size(); ++i)
			{
				Assert::AreEqual(m.stretchConditions()[i].inds()[0], triangleInds[3 * i + 0]);
				Assert::AreEqual(m.stretchConditions()[i].inds()[1], triangleInds[3 * i + 1]);
				Assert::AreEqual(m.stretchConditions()[i].inds()[2], triangleInds[3 * i + 2]);

				Assert::AreEqual(m.shearConditions()[i].inds()[0], triangleInds[3 * i + 0]);
				Assert::AreEqual(m.shearConditions()[i].inds()[1], triangleInds[3 * i + 1]);
				Assert::AreEqual(m.shearConditions()[i].inds()[2], triangleInds[3 * i + 2]);
			}
			
			std::set<int> bendConditionsFound;

			// check all the shared edges and make sure there's a corresponding bend condition:
			for (size_t i = 0; i < triangleInds.size(); i += 3)
			{
				// find triangles sharing an edge with triangle i:
				for (size_t j = i + 3; j < triangleInds.size(); j += 3)
				{
					int iverts[] = { triangleInds[i + 0], triangleInds[i + 1], triangleInds[i + 2] };
					int jverts[] = { triangleInds[j + 0], triangleInds[j + 1], triangleInds[j + 2] };

					std::vector<int> sharedVerts;
					if (triangleInds[i + 0] == triangleInds[j + 0]) sharedVerts.push_back(triangleInds[i + 0]);
					if (triangleInds[i + 1] == triangleInds[j + 0]) sharedVerts.push_back(triangleInds[i + 1]);
					if (triangleInds[i + 2] == triangleInds[j + 0]) sharedVerts.push_back(triangleInds[i + 2]);
					if (triangleInds[i + 0] == triangleInds[j + 1]) sharedVerts.push_back(triangleInds[i + 0]);
					if (triangleInds[i + 1] == triangleInds[j + 1]) sharedVerts.push_back(triangleInds[i + 1]);
					if (triangleInds[i + 2] == triangleInds[j + 1]) sharedVerts.push_back(triangleInds[i + 2]);
					if (triangleInds[i + 0] == triangleInds[j + 2]) sharedVerts.push_back(triangleInds[i + 0]);
					if (triangleInds[i + 1] == triangleInds[j + 2]) sharedVerts.push_back(triangleInds[i + 1]);
					if (triangleInds[i + 2] == triangleInds[j + 2]) sharedVerts.push_back(triangleInds[i + 2]);

					// found a shared edge!
					if (sharedVerts.size() == 2)
					{
						// search for it in the bend conditions:
						bool found(false);
						for (size_t k = 0; k < m.bendConditions().size(); ++k)
						{
							if (m.bendConditions()[k].inds()[1] == sharedVerts[0] &&
								m.bendConditions()[k].inds()[2] == sharedVerts[1])
							{
								found = true;
							}
							else if (m.bendConditions()[k].inds()[2] == sharedVerts[0] &&
								m.bendConditions()[k].inds()[1] == sharedVerts[1])
							{
								found = true;
							}

							if (found)
							{
								bendConditionsFound.insert((int)k);
								std::set<int> trivertset;
								trivertset.insert(triangleInds[i + 0]);
								trivertset.insert(triangleInds[i + 1]);
								trivertset.insert(triangleInds[i + 2]);
								trivertset.insert(triangleInds[j + 0]);
								trivertset.insert(triangleInds[j + 1]);
								trivertset.insert(triangleInds[j + 2]);

								std::set<int> bcvertset;
								bcvertset.insert(m.bendConditions()[k].inds()[0]);
								bcvertset.insert(m.bendConditions()[k].inds()[1]);
								bcvertset.insert(m.bendConditions()[k].inds()[2]);
								bcvertset.insert(m.bendConditions()[k].inds()[3]);

								Assert::IsTrue( trivertset == bcvertset );

								break;
							}
						}
						Assert::IsTrue( found );
					}
				}
			}

			Assert::AreEqual(bendConditionsFound.size(), m.bendConditions().size());
		}

		double numericalForce(const ClothMesh<double> &m, const Eigen::VectorXd x, const Eigen::VectorXd &uv, int index, double dx)
		{
			Eigen::VectorXd xTest(x.size());
			xTest.setConstant(0);
			xTest += x;

			xTest[index] = x[index] - dx;
			double eMinus = m.energy(xTest, uv);

			xTest[index] = x[index] + dx;
			double ePlus = m.energy(xTest, uv);

			// f = -dE/dx
			return -(ePlus - eMinus) / (2 * dx);
		}

		Eigen::VectorXd numericalForceDerivative(const ClothMesh<double> &m, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, int i, double dx)
		{
			Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddv((int)x.size(), (int)x.size());
			Eigen::VectorXd dampingForces((int)x.size());
			Eigen::SparseMatrix<double> dampingPseudoDerivatives((int)x.size(), (int)x.size());

			double xOrig = x[i];

			x[i] = xOrig + dx;

			Eigen::VectorXd fPlus(x.size());
			fPlus.setConstant(0);
			m.forcesAndDerivatives(x, uv, v, fPlus, dampingForces, dfdx, dddx, dddv);

			x[i] = xOrig - dx;

			Eigen::VectorXd fMinus(x.size());
			fMinus.setConstant(0);
			m.forcesAndDerivatives(x, uv, v, fMinus, dampingForces, dfdx, dddx, dddv);

			x[i] = xOrig;

			// df/dx
			return (fPlus - fMinus) / (2 * dx);
		}

		Eigen::VectorXd numericalDampingDerivative(const ClothMesh<double> &m, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, int i, double dx)
		{
			Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddv((int)x.size(), (int)x.size());
			Eigen::VectorXd forces((int)x.size());
			Eigen::SparseMatrix<double> dampingPseudoDerivatives((int)x.size(), (int)x.size());

			double vOrig = v[i];

			v[i] = vOrig + dx;

			Eigen::VectorXd fPlus(v.size());
			fPlus.setConstant(0);
			m.forcesAndDerivatives(x, uv, v, forces, fPlus, dfdx, dddx, dddv);

			v[i] = vOrig - dx;

			Eigen::VectorXd fMinus(x.size());
			fMinus.setConstant(0);
			m.forcesAndDerivatives(x, uv, v, forces, fMinus, dfdx, dddx, dddv);

			v[i] = vOrig;

			// df/dx
			return (fPlus - fMinus) / (2 * dx);
		}

		TEST_METHOD(TestClothMeshForcesAndDerivatives)
		{
			int nx = 4;
			int ny = 4;

			Eigen::VectorXd v((nx + 1) * (ny + 1) * 3);
			Eigen::VectorXd x((nx + 1) * (ny + 1) * 3);
			Eigen::VectorXd uv((nx + 1) * (ny + 1) * 2);

			for (int i = 0; i <= nx; ++i)
			{
				for (int j = 0; j <= nx; ++j)
				{
					int base = i + (nx + 1) * j;
					x[3 * base + 0] = uv[2 * base + 0] = (float)i / nx;
					x[3 * base + 1] = uv[2 * base + 1] = (float)j / ny;
					x[3 * base + 2] = 0.1f * i * j / (nx * ny);

					v[3 * base + 0] = v[3 * base + 1] = 0;
					v[3 * base + 2] = -x[3 * base + 2];
				}
			}

			std::vector<int> triangleInds;
			for (int i = 0; i < nx; ++i)
			{
				for (int j = 0; j < ny; ++j)
				{
					int base = i + (nx + 1) * j;

					triangleInds.push_back(base + 0);
					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + (nx + 1));

					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + (nx + 2));
					triangleInds.push_back(base + (nx + 1));

				}
			}

			double kBend(1.5), kStretch(200.0), kShear(100.0);
			double dBend(2.0), dStretch(20.0), dShear(10.0);

			ClothMesh<double> m(
				x, v, uv, triangleInds,
				kBend, kStretch, kShear,
				dBend, dStretch, dShear,
				1.0
			);

			Eigen::VectorXd forces((nx + 1) * (ny + 1) * 3);
			Eigen::VectorXd dampingforces((nx + 1) * (ny + 1) * 3);
			Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dddv((int)x.size(), (int)x.size());
			m.forcesAndDerivatives(x, uv, v, forces, dampingforces, dfdx, dddx, dddv);

			double dx = 0.00001;
			double tol = 1.e-4;

			// check forces are correct:
			for( int i=0; i < forces.size(); ++i )
			{
				Assert::AreEqual( forces[i], numericalForce( m, x, uv, i, dx ), tol );
			}

			// check dfdx is correct:
			for (int i = 0; i < forces.size(); ++i)
			{
				Eigen::VectorXd fdN = dfdx.row(i);
				checkVectorEquality(fdN, numericalForceDerivative(m, uv, x, v, i, dx), tol, true);
			}

			// check damping forces are correct:

			// work out time derivative of c
			Eigen::VectorXd xTest;
			xTest = x + dx * v;
			Eigen::VectorXd cPlus;
			m.C(x + dx * v, uv, cPlus);
			Eigen::VectorXd cMinus;
			m.C(x - dx * v, uv, cMinus);

			Eigen::VectorXd dCdt = (cPlus - cMinus) / (2 * dx);

			Eigen::VectorXd dCdxi;
			for (int i = 0; i < dampingforces.size(); ++i)
			{

				// work out dCdxi:
				double xOrig = x[i];
				x[i] = xOrig + dx;
				m.C(x, uv, cPlus);

				x[i] = xOrig - dx;
				m.C(x, uv, cMinus);

				x[i] = xOrig;
				dCdxi = (cPlus - cMinus) / ( 2 * dx );

				double dampingForce(0);
				size_t n = 0;
				for (size_t j = 0; j < m.bendConditions().size(); ++j, ++n)
				{
					dampingForce -= dBend * dCdxi[n] * dCdt[n];
				}
				for (size_t j = 0; j < m.shearConditions().size(); ++j, ++n)
				{
					dampingForce -= dShear * dCdxi[n] * dCdt[n];
				}
				for (size_t j = 0; j < m.stretchConditions().size(); ++j, n += 2)
				{
					dampingForce -= dStretch * dCdxi[n] * dCdt[n];
					dampingForce -= dStretch * dCdxi[n+1] * dCdt[n+1];
				}

				Assert::AreEqual(dampingforces[i], dampingForce, tol);
			}

			// check dddv is correct:
			for (int i = 0; i < forces.size(); ++i)
			{
				Eigen::VectorXd fdN = dddv.row(i);
				checkVectorEquality(fdN, numericalDampingDerivative(m, uv, x, v, i, dx), tol, true);
			}




			// \todo: check dddx works... even though it doesn't really
			




			// check doubling up the arguments adds the results:
			Eigen::VectorXd f2((nx + 1) * (ny + 1) * 3);
			Eigen::SparseMatrix<double> dfdx2((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dfdv((int)x.size(), (int)x.size());

			m.forcesAndDerivatives(x, uv, v, f2, f2, dfdx2, dfdx2, dfdv);
			checkVectorEquality(f2, forces + dampingforces, 1.e-6, true);

			for (int i = 0; i < dfdx2.outerSize(); ++i)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(dfdx2, i); it; ++it)
				{
					Assert::AreEqual( it.value(), dfdx.coeff(it.row(),it.col()) + dddx.coeff( it.row(), it.col() ), 1.e-6 );
				}
			}

			for (int i = 0; i < dfdv.outerSize(); ++i)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(dfdv, i); it; ++it)
				{
					Assert::AreEqual(it.value(), dddv.coeff(it.row(), it.col()), 1.e-6);
				}
			}

		}

		TEST_METHOD(TestClothMeshImplicitUpdate)
		{
			int nx = 4;
			int ny = 4;

			Eigen::VectorXd v((nx + 1) * (ny + 1) * 3);
			Eigen::VectorXd x((nx + 1) * (ny + 1) * 3);
			Eigen::VectorXd uv((nx + 1) * (ny + 1) * 2);

			for (int i = 0; i <= nx; ++i)
			{
				for (int j = 0; j <= nx; ++j)
				{
					int base = i + (nx + 1) * j;
					x[3 * base + 0] = uv[2 * base + 0] = (float)i / nx;
					x[3 * base + 1] = uv[2 * base + 1] = (float)j / ny;
					x[3 * base + 2] = 0.1f * i * j / (nx * ny);

					v[3 * base + 0] = v[3 * base + 1] = 0;
					v[3 * base + 2] = -x[3 * base + 2];
				}
			}

			std::vector<int> triangleInds;
			for (int i = 0; i < nx; ++i)
			{
				for (int j = 0; j < ny; ++j)
				{
					int base = i + (nx + 1) * j;

					triangleInds.push_back(base + 0);
					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + (nx + 1));

					triangleInds.push_back(base + 1);
					triangleInds.push_back(base + (nx + 2));
					triangleInds.push_back(base + (nx + 1));

				}
			}

			double kBend(1.5), kStretch(200.0), kShear(100.0);
			double dBend(2.0), dStretch(20.0), dShear(10.0);

			ClothMesh<double> m(
				x, v, uv, triangleInds,
				kBend, kStretch, kShear,
				dBend, dStretch, dShear,
				1.0
				);



			Eigen::VectorXd forces((nx + 1) * (ny + 1) * 3);
			Eigen::SparseMatrix<double> dfdx((int)x.size(), (int)x.size());
			Eigen::SparseMatrix<double> dfdv((int)x.size(), (int)x.size());

			Eigen::SparseMatrix<double> A((int)x.size(), (int)x.size());
			Eigen::VectorXd rhs((int)x.size());

			for (int i = 0; i < 20; ++i)
			{
				m.forcesAndDerivatives(m.x(), uv, m.v(), forces, forces, dfdx, dfdx, dfdv);

				// assemble equation:
				double dt(0.1);

				m.assembleImplicitUpdateEquations(dt, forces, m.v(), dfdx, dfdv, A, rhs);

				// solve to find changes in velocity and position:
				Eigen::SparseLU< Eigen::SparseMatrix<double> > s;
				s.analyzePattern(A);
				s.factorize(A);

				Eigen::VectorXd dv = s.solve(rhs);
				Eigen::VectorXd dx = (m.v() + dv) * dt;

				// now check the following:

				// df = dfdx * dx + dfdv * dv
				// fNext = forces + df
				// aNext = M^-1 fNext
				// dv = aNext * dt

				Eigen::VectorXd df = dfdx * dx + dfdv * dv;
				Eigen::VectorXd fNext = forces + df;
				Eigen::VectorXd aNext(fNext.size());
				for (int i = 0; i < m.m().size(); ++i)
				{
					aNext.segment<3>(3 * i) = fNext.segment<3>(3 * i) / m.m()[i];
				}

				checkVectorEquality(dv, aNext * dt, 1.e-6, true);

				Eigen::VectorXd vOrig = m.v();
				Eigen::VectorXd xOrig = m.x();

				DirectSolver<double> solver;
				std::vector< ForceField<double>* > forceFields;
				m.advance(forceFields, dt, solver);

				checkVectorEquality(m.v() - vOrig, dv, 1.e-6, true);
				checkVectorEquality(m.x() - xOrig, dx, 1.e-6, true);
			}

		}
	};
}
