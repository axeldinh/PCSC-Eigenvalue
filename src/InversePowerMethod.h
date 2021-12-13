//
// Created by axeld on 08/12/2021.
//

// TODO Documentation

#ifndef EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H
#define EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H

#include <iostream>

#include "GeneralEigenSolver.h"

template <typename ScalarType>
class InversePowerMethod: public GeneralEigenSolver<ScalarType> {
public:
    InversePowerMethod();
    ~InversePowerMethod();

    ScalarType solve();
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

template <typename ScalarType>
InversePowerMethod<ScalarType>::InversePowerMethod()
    : GeneralEigenSolver<ScalarType>() {}

template <typename ScalarType>
InversePowerMethod<ScalarType>::~InversePowerMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

template <typename ScalarType>
ScalarType InversePowerMethod<ScalarType>::solve() {

    if (!this->getIsVectorInit()) {
        this->initRandomEigenVector();
    }

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // Computing the inverse of A
    auto I = Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>(GeneralEigenSolver<ScalarType>::mMatrix->rows(), GeneralEigenSolver<ScalarType>::mMatrix->cols());
    I.setZero();
    I.diagonal() = Eigen::Vector<ScalarType,Eigen::Dynamic>(GeneralEigenSolver<ScalarType>::mMatrix->rows()).setOnes();
    auto B = GeneralEigenSolver<ScalarType>::mMatrix->partialPivLu().inverse();

    double threshold = this->getThreshold();
    int maxIter = this->getMaxIter();
    double error;
    int iter = 0;
    ScalarType lambda;

    while (iter < maxIter) {
        // Do one iteration of the Power Method. The matrix is fetched from the mother class
        iter++;
        auto temp = B * GeneralEigenSolver<ScalarType>::mEigenVector;
        auto temp_norm = temp.norm();
        if (temp_norm < 1e-15) {
            throw std::invalid_argument("INVALID STARTING VECTOR: The guessed eigenvector "
                                        "is currently in the null "
                                        "space of the matrix,\ntry to initialize "
                                        "the vector with a different value "
                                        "(e.g using GeneralEigenSolver<ScalarType>::setEigenVector(...))\n");
        }
        GeneralEigenSolver<ScalarType>::mEigenVector = temp / temp_norm; // Update the

        // Compute the corresponding eigenvalue
        auto Bv = B * GeneralEigenSolver<ScalarType>::mEigenVector;
        lambda = GeneralEigenSolver<ScalarType>::mEigenVector.transpose() * Bv;
        error = (Bv - lambda * GeneralEigenSolver<ScalarType>::mEigenVector).norm();

        if (error < threshold) {
            return 1. / lambda;
        }
    }
    std::cout << "The Shifted Inverse Power Method did not converge after " << iter << " iterations\n";
    return 1. / lambda;
}

#endif //EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H
