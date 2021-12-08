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

    auto A = GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 0; i < this->getMaxIter(); i++) {
        auto Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    ScalarType eigenValue = 0.;

    for (int i = 0; i < A.cols(); i++) {
        if (std::abs(A(i,i)) > std::abs(eigenValue)) {
            eigenValue = A(i,i);
        }
    }

    return eigenValue;
}

// TODO need to sort diagonal and return n largest eigenvalues
template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve(int n) {

    auto A = GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 0; i < this->getMaxIter(); i++) {
        auto Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    ScalarType eigenValue = 0.;

    for (int i = 0; i < A.cols(); i++) {
        if (std::abs(A(i,i)) > std::abs(eigenValue)) {
            eigenValue = A(i,i);
        }
    }

    return eigenValue;
}

template<typename ScalarType>
VectorType<ScalarType> QRMethod<ScalarType>::solve_all() {
    auto A = GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 0; i < this->getMaxIter(); i++) {
        auto Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    return A.diagonal();
}

#endif //EIGENVALUE_PROJECT_QRMETHOD_H
