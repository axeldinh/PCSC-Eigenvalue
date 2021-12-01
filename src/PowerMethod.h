//
// Created by axeld on 01/12/2021.
//

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
#include "GeneralEigenSolver.h"

template<typename ScalarType>
class PowerMethod: public GeneralEigenSolver<ScalarType> {
public:
    // Constructors and Destructors
    PowerMethod();
    ~PowerMethod();

    // Solver
    ScalarType solve();
};

/********************************************//**
 *  Constructors and Destructors
 ***********************************************/

template<typename ScalarType>
PowerMethod<ScalarType>::PowerMethod()
    : GeneralEigenSolver<ScalarType>() {}

template<typename ScalarType>
PowerMethod<ScalarType>::~PowerMethod() {}

/********************************************//**
 *  Solver
 ***********************************************/

template<typename ScalarType>
ScalarType PowerMethod<ScalarType>::solve() {

    // If the user did not initialize the EigenVector
    // We do it now

    if (!this->getIsVectorInit()) {
        this->initRandomEigenVector();
    }

    assert(this->getIsMatrixInit());

    ScalarType lambda;

    double threshold = this->getThreshold();
    int maxIter = this->getMaxIter();
    double error;
    int iter = 0;

    while (iter < maxIter) {
        // Do one iteration of the Power Method. The matrix is fetched from the mother class
        iter++;
        auto temp = GeneralEigenSolver<ScalarType>::mMatrix * GeneralEigenSolver<ScalarType>::mEigenVector;
        auto temp_norm = temp.norm();
        assert (temp_norm > 1e-15); // TODO provide exceptions (here it means that init vector is in the nullspace of A)
        GeneralEigenSolver<ScalarType>::mEigenVector = temp / temp_norm; // Update the

        // Compute the corresponding eigenvalue
        lambda = (GeneralEigenSolver<ScalarType>::mEigenVector.transpose() * GeneralEigenSolver<ScalarType>::mMatrix) * GeneralEigenSolver<ScalarType>::mEigenVector;
        error = (GeneralEigenSolver<ScalarType>::mMatrix * GeneralEigenSolver<ScalarType>::mEigenVector - lambda * GeneralEigenSolver<ScalarType>::mEigenVector).norm();

        if (error < threshold) {
            return lambda;
        }
    }
    std::cout << "The Power Method did not converge after " << iter << " iterations\n";
    return lambda;
}




#endif //EIGENVALUE_PROJECT_POWERMETHOD_H
