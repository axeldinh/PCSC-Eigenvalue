//
// Created by axeld on 08/12/2021.
//

// TODO Documentation

#ifndef EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H
#define EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H

#include "GeneralPowerMethod.h"
#include <iostream>

template<typename ScalarType>
class ShiftedInversePowerMethod: public GeneralPowerMethod<ScalarType> {

protected:
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
    ScalarType getShift() const;
    bool getIsShiftInit() const;
    ///@}

    ScalarType solve();
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

template<typename ScalarType>
ShiftedInversePowerMethod<ScalarType>::ShiftedInversePowerMethod()
    : GeneralPowerMethod<ScalarType>() {
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
ScalarType ShiftedInversePowerMethod<ScalarType>::getShift() const {
    return mShift;
}

template <typename ScalarType>
bool ShiftedInversePowerMethod<ScalarType>::getIsShiftInit() const {
    return isShiftInit;
}

/*==========================================*//*=
=  Solver
==============================================*/

template <typename ScalarType>
ScalarType ShiftedInversePowerMethod<ScalarType>::solve() {

    if (!this->getIsShiftInit()) {
        throw UninitializedSolver("shift", "please initialize with ShiftedInversePowerMethod<typename ScalarType>::setShift");
    }

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // Computing the inverse of A - shift * I
    auto I = Eigen::Matrix<ScalarType,Eigen::Dynamic,Eigen::Dynamic>(this->mMatrix->rows(), this->mMatrix->cols());
    I.setIdentity();
    MatrixType<ScalarType> matrixInv = (*this->mMatrix - mShift*I).partialPivLu().inverse();
    ScalarType lambda;

    try {
        lambda = GeneralPowerMethod<ScalarType>::solve(matrixInv);
    } catch (std::invalid_argument& e) {
        throw;
    } catch (std::exception& e) {
        throw;
    }

    return 1. / lambda + mShift;
}


#endif //EIGENVALUE_PROJECT_SHIFTEDINVERSEPOWERMETHOD_H
