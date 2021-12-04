
#include <iostream>
#include <Eigen/Dense>
#include "GeneralEigenSolver.h"
#include "PowerMethod.h"

int main() {
    constexpr unsigned int n = 10;
    Eigen::Matrix<std::complex<double>, n, n> A;

    A.setZero();
    A.diagonal() = Eigen::Vector<std::complex<double>, n>::Constant(1.);
    A(0,0) = std::complex<double>(0., 2.);

    Eigen::Vector<double, n> V;
    V.setOnes();

    double threshold = 1e-15;
    int maxIter = 100;

    auto solver = new PowerMethod<std::complex<double>>();
    solver->setThreshold(threshold);
    solver->setMaxIter(maxIter);
    solver->setMatrix(A);
    solver->setEigenVector(V);

    std::cout << "threshold: " << solver->getThreshold() << ", max iters: " << solver->getMaxIter() << std::endl;

    try {
        auto lambda = solver->solve();
        std::cout << "Solution: " << lambda;
    } catch(std::invalid_argument& e) {
        std::cerr << e.what();
    }

    return 0;
}
