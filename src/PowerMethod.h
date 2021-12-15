
#ifndef EIGENVALUE_PROJECT_POWERMETHOD_H
#define EIGENVALUE_PROJECT_POWERMETHOD_H

#include <iostream>
#include "GeneralPowerMethod.h"
#include "Exceptions/UninitializedSolver.h"

/**
 * Class to solve eigenvalues problems using the power method.
 * Aims at solving the eigenvalue problem \f$Ax=\lambda x\f$, see GeneralPowerMethod for more information about the algorithm.
 *
 * If the eigenvalues are such that \f$|\lambda_1| \ge |\lambda_2| \ge |\lambda_n|\f$, then \f$\lambda_1\f$ should be returned,
 * unless the starting vector is in the null space of \f$A\f$ or the starting vector is the eigenvector corresponding to another eigenvalue.
 *  @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */

template<typename ScalarType>
class PowerMethod: public GeneralPowerMethod<ScalarType> {
public:
    // Constructors and Destructors
    PowerMethod();
    ~PowerMethod();

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
template<typename ScalarType>
PowerMethod<ScalarType>::PowerMethod()
    : GeneralPowerMethod<ScalarType>() {}

/**
 * Destructor.
 * Uses the destructor of the GeneralPowerMethod class.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
PowerMethod<ScalarType>::~PowerMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

/**
 * Solver for the eigenvalue problem.
 * Solves the eigenvalue problem \f$Ax=\lambda x\f$ by returning the highest eigenvalue with respect to its absolute value.
 * The corresponding eigenvector is stored in the #mEigenVector attribute.
 *
 * Throws an UninitializedSolver exception if #mMatrix has not been initialized.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return ScalarType, \f$\lambda\f$ the eigenvalue.
 * @throws UninitializedSolver
 */
template<typename ScalarType>
ScalarType PowerMethod<ScalarType>::solve() {

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    ScalarType lambda;

    try {
        lambda = GeneralPowerMethod<ScalarType>::solve(*this->mMatrix);
    } catch (std::invalid_argument& e) {
        throw;
    } catch (std::exception& e) {
        throw;
    }

}

#endif //EIGENVALUE_PROJECT_POWERMETHOD_H
