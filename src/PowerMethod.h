/**
 * Class to solve eigenvalues problems using the power method.
 * The power method is given by:
 *      Given \f$ A\in \mathbb{C}^{n \times n} \f$:
 *          Choose a starting vector \f$ x_0 \in \mathbb{C}^n \f$
 *          k = 0
 *          repeat:
 *              k = k+1
 *              Compute y_k = Ax_{k-1}
 *              Normalize x_k = \frac{y_k}{||y_k||_2}
 *          End at convergence or when the number of iterations has exceeded maxIter.
 *
 */

#ifndef EIGENVALUE_PROJECT_POWERMETHOD_H
#define EIGENVALUE_PROJECT_POWERMETHOD_H

#include <iostream>
#include "GeneralPowerMethod.h"
#include "Exceptions/UninitializedSolver.h"

/**
 * Class to solve eigenvalues problems using the power method.
 * The power method is given by:
 *      Given \f$ A\in \mathbb{C}^{n \times n} \f$:
 *          Choose a starting vector \f$ x_0 \in \mathbb{C}^n \f$
 *          k = 0
 *          repeat:
 *              k = k+1
 *              Compute y_k = Ax_{k-1}
 *              Normalize x_k = \frac{y_k}{||y_k||_2}
 *          End at convergence or when the number of iterations has exceeded maxIter.
 *  @tparam ScalarType
 */

template<typename ScalarType>
class PowerMethod: public GeneralPowerMethod<ScalarType> {
public:
    // Constructors and Destructors
    PowerMethod();
    ~PowerMethod();

    // Solver
    ScalarType solve();
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

template<typename ScalarType>
PowerMethod<ScalarType>::PowerMethod()
    : GeneralPowerMethod<ScalarType>() {}

template<typename ScalarType>
PowerMethod<ScalarType>::~PowerMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

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

    return lambda;

    /*
    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // If the user did not initialize the EigenVector
    // We do it now

    if (!this->getIsVectorInit()) {
        this->initRandomEigenVector();
    }

    ScalarType lambda;

    double threshold = this->getThreshold();
    int maxIter = this->getMaxIter();
    double error;
    int iter = 0;

    while (iter < maxIter) {
        // Do one iteration of the Power Method. The matrix is fetched from the mother class
        iter++;
        auto temp = *GeneralEigenSolver<ScalarType>::mMatrix * GeneralEigenSolver<ScalarType>::mEigenVector;
        auto temp_norm = temp.norm();
        if (temp_norm < 1e-15) {
            throw std::invalid_argument("INVALID STARTING VECTOR: The guessed eigenvector "
                                        "is currently in the null "
                                        "space of the matrix,\ntry to initialize "
                                        "the vector with a different value "
                                        "(e.g using PowerMethod<ScalarType>::setEigenVector(...))\n");
        }
        GeneralEigenSolver<ScalarType>::mEigenVector = temp / temp_norm; // Update the

        // Compute the corresponding eigenvalue
        auto Av = *GeneralEigenSolver<ScalarType>::mMatrix * GeneralEigenSolver<ScalarType>::mEigenVector;
        lambda = GeneralEigenSolver<ScalarType>::mEigenVector.transpose() * Av;
        error = (Av - lambda * GeneralEigenSolver<ScalarType>::mEigenVector).norm();

        if (error < threshold) {
            return lambda;
        }
    }
    std::cout << "The method did not converge after " << iter << " iterations\n";
    return lambda;
     */
}

#endif //EIGENVALUE_PROJECT_POWERMETHOD_H
