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

template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve(int n) {

    if (n > GeneralEigenSolver<ScalarType>::mMatrix->rows()) {
        throw std::invalid_argument("INVALID EIGENVALUES NUMBER: n is larger than the maximum number of eigenvalues.\n");
    }

    VectorType<ScalarType> eigenValues;

    try {
        eigenValues = this->solve_all();
    } catch(UninitializedSolver& e) {
        throw e;
    }

    std::sort(eigenValues.begin(), eigenValues.end(), std::greater<ScalarType>());

    return eigenValues(n-1);
}

template<typename ScalarType>
VectorType<ScalarType> QRMethod<ScalarType>::solve_all() {

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // A is a copy of mMatrix, so the original matrix is not modified
    auto A = *GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 1; i < this->getMaxIter(); i++) {
        MatrixType<ScalarType> Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    return A.diagonal();
}

#endif //EIGENVALUE_PROJECT_QRMETHOD_H
