
// TODO Documentation

#ifndef EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H
#define EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H

#include <iostream>

#include "GeneralPowerMethod.h"

template <typename ScalarType>
class InversePowerMethod: public GeneralPowerMethod<ScalarType> {
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
    : GeneralPowerMethod<ScalarType>() {}

template <typename ScalarType>
InversePowerMethod<ScalarType>::~InversePowerMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

template <typename ScalarType>
ScalarType InversePowerMethod<ScalarType>::solve() {

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // Computing the inverse of A
    MatrixType<ScalarType> matrixInv = this->mMatrix->partialPivLu().inverse();

    ScalarType lambda;

    try {
        lambda = GeneralPowerMethod<ScalarType>::solve(matrixInv);
    } catch (std::invalid_argument& e) {
        throw;
    } catch (std::exception& e) {
        throw;
    }

    return 1. / lambda;
}

#endif //EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H
