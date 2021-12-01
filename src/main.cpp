
#include <iostream>
#include <Eigen/Dense>
#include "GeneralEigenSolver.h"
#include "PowerMethod.h"

int main() {
    constexpr unsigned int n = 100;
    Eigen::Matrix<double, n, n> A;
    A.setZero();
    A.diagonal() = Eigen::Vector<double, n>::Constant(1.);
    A(0,0) = 2.;

    //Eigen::Vector<double, n> V;
    //V.setRandom();

    double threshold = 1e-8;
    int maxIter = 100;
    auto solver = new PowerMethod<double>();
    solver->setThreshold(threshold);
    solver->setMaxIter(maxIter);
    solver->setMatrix(A);
    //solver->setEigenVector(V);

    std::cout << "threshold: " << solver->getThreshold() << ", max iters: " << solver->getMaxIter() << std::endl;

    auto lambda = solver->solve();
    std::cout << "Solution: " << lambda;

    return 0;
}
