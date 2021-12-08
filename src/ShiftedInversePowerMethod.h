//
// Created by axeld on 08/12/2021.
//

#ifndef EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H
#define EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H

#include "GeneralEigenSolver.h"
#include <iostream>

template<typename ScalarType>
class ShiftedInversePowerMethod: public GeneralEigenSolver<ScalarType> {

private:
    ScalarType mShift;
    bool isShiftInit;

public:
    ShiftedInversePowerMethod();
    ~ShiftedInversePowerMethod();

    /** @name Setters
     *
     */
    ///@{
    void setShift(ScalarType shift);
    ///@}

    /** @name Getters
     *
     */
    ///@{
    ScalarType getShift();
    bool getIsShiftInit();
    ///@}

    ScalarType solve();
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

template<typename ScalarType>
ShiftedInversePowerMethod<ScalarType>::ShiftedInversePowerMethod()
    : GeneralEigenSolver<ScalarType>() {
        isShiftInit = false;
    }

template<typename ScalarType>
ShiftedInversePowerMethod<ScalarType>::~ShiftedInversePowerMethod() {}

/*==============================================*//*=
 =  Setters
 ================================================*/

template <typename ScalarType>
void ShiftedInversePowerMethod<ScalarType>::setShift(ScalarType shift) {
    mShift = shift;
    isShiftInit = true;
}

/*==========================================*//*=
=  Getters
==============================================*/

template <typename ScalarType>
ScalarType ShiftedInversePowerMethod<ScalarType>::getShift() {
    return mShift;
}

template <typename ScalarType>
bool ShiftedInversePowerMethod<ScalarType>::getIsShiftInit() {
    return isShiftInit;
}

/*==========================================*//*=
=  Solver
==============================================*/

template <typename ScalarType>
ScalarType ShiftedInversePowerMethod<ScalarType>::solve() {

    if (!this->getIsVectorInit()) {
        this->initRandomEigenVector();
    }

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    if (!this->getIsShiftInit()) {
        throw UninitializedSolver("shift", "please initialize with ShiftedInversePowerMethod<typename ScalarType>::setShift");
    }

    // Computing the inverse of A - shift * I
    auto I = Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>(GeneralEigenSolver<ScalarType>::mMatrix.rows(), GeneralEigenSolver<ScalarType>::mMatrix.cols());
    I.setZero();
    I.diagonal() = Eigen::Vector<ScalarType,Eigen::Dynamic>(GeneralEigenSolver<ScalarType>::mMatrix.rows()).setOnes();
    auto B = (GeneralEigenSolver<ScalarType>::mMatrix - mShift * I).partialPivLu().inverse();

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
            return 1. / lambda + mShift;
        }
    }
    std::cout << "The Shifted Inverse Power Method did not converge after " << iter << " iterations\n";
    return 1. / lambda + mShift;


}


#endif //EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H
