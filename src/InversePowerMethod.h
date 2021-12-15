
// TODO Documentation

#ifndef EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H
#define EIGENVALUE_PROJECT_INVERSEPOWERMETHOD_H

#include <iostream>
#include "GeneralPowerMethod.h"

/**
 * Class to solve an eigenvalue problem using the inverse power method.
 * This class aims at solving \f$Ax = \lambda x\f$, see GeneralPowerMethod for more information about the algorithm.
 *
 * For the inverse power method, the power method algorithm is used on \f$A^{-1}\f$.
 *
 * If the eigenvalues are such that \f$|\lambda_1| \ge |\lambda_2| \ge |\lambda_n|\f$, then \f$\lambda_n\f$ should be returned,
 * unless the starting vector is in the null space of \f$A\f$ or the starting vector is the eigenvector corresponding to another eigenvalue.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */

template <typename ScalarType>
class InversePowerMethod: public GeneralPowerMethod<ScalarType> {
public:
    InversePowerMethod();
    ~InversePowerMethod();

    /** @name Solve method
     *
     */
    ///@{
    ScalarType solve();
    ///@}
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

/**
 * Basic constructor.
 * Uses the constructor of the GeneralPowerMethod class.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template <typename ScalarType>
InversePowerMethod<ScalarType>::InversePowerMethod()
    : GeneralPowerMethod<ScalarType>() {}

/**
* Destructor.
* Uses the destructor of the GeneralPowerMethod class.
* @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
*/
template <typename ScalarType>
InversePowerMethod<ScalarType>::~InversePowerMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

/**
 * Solver for the eigenvalue problem.
 * Solves the eigenvalue problem \f$Ax=\lambda x\f$ by returning the lowest eigenvalue with respect to its absolute value.
 * The corresponding eigenvector is stored in the #mEigenVector attribute.
 *
 * Throws an UninitializedSolver exception if #mMatrix has not been initialized.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return ScalarType, \f$\lambda\f$ the eigenvalue.
 * @throws UninitializedSolver
 */
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
