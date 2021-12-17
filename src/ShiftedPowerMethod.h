#include "GeneralPowerMethod.h"
#include <iostream>

#ifndef EIGENVALUE_PROJECT_SHIFTEDPOWERMETHOD_H
#define EIGENVALUE_PROJECT_SHIFTEDPOWERMETHOD_H

/**
 * Class to solve an eigenvalue problem using the shifted power method.
 * This class aims at solving \f$Ax = \lambda x\f$, see GeneralPowerMethod for more information about the algorithm.
 *
 * For the shifted  power method, the power method algorithm is used on \f$A - \sigma I\f$, where \f$\sigma\f$
 * is some shift (#mShift) given by the user. #mShift should be given has a close approximate of the desired eigenvalue.
 *
 * If the eigenvalues are such that \f$|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|\f$, then the largest \f$|\lambda_i - \sigma|\f$ should be returned,
 * unless the starting vector is in the null space of \f$A\f$ or the starting vector is the eigenvector corresponding to another eigenvalue.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */

template<typename ScalarType>
class ShiftedPowerMethod: public GeneralPowerMethod<ScalarType> {

    using MatrixType = Eigen::Matrix<ScalarType, -1, -1>;

protected:
    ScalarType mShift; /**< Shift for the shifted power method algorithm. */
    bool isShiftInit; /**< Bool, true if the shift has been initialized. */

public:
    ShiftedPowerMethod();
    ~ShiftedPowerMethod();

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
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
ShiftedPowerMethod<ScalarType>::ShiftedPowerMethod()
        : GeneralPowerMethod<ScalarType>() {
    isShiftInit = false;
}

/**
* Destructor.
* Uses the destructor of the GeneralPowerMethod class.
* @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
*/
template<typename ScalarType>
ShiftedPowerMethod<ScalarType>::~ShiftedPowerMethod() {}

/*==============================================*//*=
 =  Setters
 ================================================*/

/**
 * Setter for #mShift.
 * Modifies the value of #mShift.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param shift, ScalarType, the desired shift.
 */
template <typename ScalarType>
void ShiftedPowerMethod<ScalarType>::setShift(ScalarType shift) {
    mShift = shift;
    isShiftInit = true;
}

/*==========================================*//*=
=  Getters
==============================================*/

/**
 * Getter for #mShift.
 * Returns #mShift.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return ScalarType #mShift, the shift for the shifted power method (see the detailed description ShiftedInversePowerMethod).
 */
template <typename ScalarType>
ScalarType ShiftedPowerMethod<ScalarType>::getShift() const {
    return mShift;
}

/**
 * Getter for #isShiftInit.
 * Returns #isShiftInit.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return bool #isShiftInit, true if #mShift has been initialized.
 */
template <typename ScalarType>
bool ShiftedPowerMethod<ScalarType>::getIsShiftInit() const {
    return isShiftInit;
}

/*==========================================*//*=
=  Solver
==============================================*/

/**
 * Solver for the eigenvalue problem.
 * Solves the eigenvalue problem \f$Ax=\lambda x\f$ by returning the closest eigenvalue to #mShift.
 * The corresponding eigenvector is stored in the #mEigenVector attribute.
 *
 * Throws an UninitializedSolver exception if #mMatrix or #mShift have not been initialized.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return ScalarType, \f$\lambda\f$ the eigenvalue.
 * @throws UninitializedSolver
 */
template <typename ScalarType>
ScalarType ShiftedPowerMethod<ScalarType>::solve() {

    if (!this->getIsShiftInit()) {
        throw UninitializedSolver("shift", "please initialize with ShiftedPowerMethod<ScalarType>::setShift");
    }

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<ScalarType>::setMatrix");
    }

    // Computing A-mShiftI
    auto I = MatrixType(this->mMatrix->rows(), this->mMatrix->cols());
    I.setIdentity();
    MatrixType matrixShifted = *this->mMatrix - mShift*I;
    ScalarType lambda;

    try {
        lambda = GeneralPowerMethod<ScalarType>::solve(matrixShifted);
    } catch (std::invalid_argument& e) {
        throw;
    } catch (std::exception& e) {
        throw;
    }

    return lambda + mShift;
}


#endif //EIGENVALUE_PROJECT_SHIFTEDPOWERMETHOD_H
