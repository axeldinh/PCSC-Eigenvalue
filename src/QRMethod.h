//
// Created by axeld on 08/12/2021.
//

// TODO Documentation

#ifndef EIGENVALUE_PROJECT_QRMETHOD_H
#define EIGENVALUE_PROJECT_QRMETHOD_H

#include "GeneralEigenSolver.h"

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template <typename ScalarType>
class QRMethod: public GeneralEigenSolver<ScalarType> {

public:
    ScalarType solve();
    ScalarType solve(int n);
    VectorType<ScalarType> solve_all();
};

template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve() {
    return this->solve(1);
}

// TODO need to sort diagonal and return n largest eigenvalues
template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve(int n) {

    VectorType<ScalarType> eigenValues = this->solve_all();

    std::sort(eigenValues.begin(), eigenValues.end());

    return eigenValues(n-1);
}

template<typename ScalarType>
VectorType<ScalarType> QRMethod<ScalarType>::solve_all() {

    auto A = GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 1; i < this->getMaxIter(); i++) {
        MatrixType<ScalarType> Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    return A.diagonal();
}

#endif //EIGENVALUE_PROJECT_QRMETHOD_H
