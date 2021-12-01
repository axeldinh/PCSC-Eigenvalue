#include <iostream>
#include <Eigen/Dense>

//#include "GeneralEigenSolver.h"
#include "src/foo.h"

int main() {
    foo A;
    A.print();
    //auto solver = new GeneralEigenSolver<double>();
    //auto m = Eigen::Matrix<double, 3, 3>::Ones();
    //solver->setMatrix(m);
    //auto m_c = solver->getMatrix();
    //std::cout << m_c << std::endl;
    return 0;
}
