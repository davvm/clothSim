#ifndef EQUALITYTESTS_H
#define EQUALITYTESTS_H

#include "EnergyCondition.h"

void checkVectorEquality(Eigen::Vector3d v0, Eigen::Vector3d v1, double tol);
void checkVectorEquality(Eigen::VectorXd v0, Eigen::VectorXd v1, double tol, bool relative = false);
void checkMatrixEquality(Eigen::Matrix3d m0, Eigen::Matrix3d m1, double tol);
double numericalForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, const Eigen::VectorXd &x, double k, int i, double dx);
double numericalDampingForce(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, const Eigen::VectorXd &x, const Eigen::VectorXd &v, double d, int i, double dx);
Eigen::VectorXd numericalForceDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, double d, int i, double dx);
Eigen::VectorXd numericalDampingForceDerivative(const EnergyCondition<double> &c, const Eigen::VectorXd &uv, Eigen::VectorXd &x, Eigen::VectorXd &v, double k, double d, int i, double dx);

#endif // EQUALITYTESTS_H