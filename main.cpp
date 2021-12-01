#include <iostream>
#include <Eigen/Dense>

//#include "GeneralEigenSolver.h"
#include "src/GeneralEigenSolver.h"

int main() {
    double threshold = 1e-15;
    int maxIter = 1000;
    auto solver = new GeneralEigenSolver<double>();
    solver->setThreshold(threshold);
    solver->setMaxIter(maxIter);

    std::cout << "threshold: " << solver->getThreshold() << ", max iters: " << solver->getMaxIter();


    //auto solver = new GeneralEigenSolver<double>();
    //auto m = Eigen::Matrix<double, 3, 3>::Ones();
    //solver->setMatrix(m);
    //auto m_c = solver->getMatrix();
    //std::cout << m_c << std::endl;

    return 0;
}
