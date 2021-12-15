
#ifndef EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
#define EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include "Exceptions/UninitializedSolver.h"

// TODO Change the place of the usings, here they propagate to all files
template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename ScalarType>
using VectorType = Eigen::Vector<ScalarType, Eigen::Dynamic>;

/**
 * Abstract class for EigenValue solvers.
 * Mother class for classes used to solve eigenvalues
 * problem of the form \f$Ax=\lambda b\f$, where \f$A \in \mathbb{C}^{n\times n}\f$ is known and \f$b \in \mathbb{C}^n\f$.
 * This abstract class is based on the assumption
 * that \f$A\f$ and \f$b\f$ are declared using the Eigen library
 * (<a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">link</a>).
 *
 * The maximum number of iterations can be defined with the #mMaxIter attribute.
 * This class also allows the initialization of the #mMatrix attribute, which is a pointer to the matrix operator \f$A\f$-
 *
 * A more general class could have been used (with a fully tempated vector space), but this approach seemed too abstract
 * for the scope of the project.
 *
 * All methods are implemented in the header file as the compiler needs the code when calling on
 * a particular template (e.g GeneralEigenSolver<double>)
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
class GeneralEigenSolver {

protected:
    MatrixType<ScalarType> *mMatrix;         /**< Pointer to the matrix operator */
    unsigned int mMaxIter;                  /**< Maximum number of iterations for the solve() method */

    bool isMatrixInit;                      /**< Boolean, true if the matrix has been initialized */

public:

    // Constructors and Destructors
    GeneralEigenSolver();
    ~GeneralEigenSolver();

    /** @name Setters
     *
     */
    ///@{
    void setMatrix(MatrixType<ScalarType>& A);
    void setMaxIter(int maxIter);
    ///@}

    /** @name Getters
     *
     */
    ///@{
    int getMaxIter() const;
    bool getIsMatrixInit() const;
    ///@}


    /** @name Other functions
     *
     */
     ///@{
    virtual ScalarType solve() = 0;
    ///@}
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

/**
 * Initialization constructor.
 * Basic constructor, sets #mMaxIter to 1000 and #isMatrixInit to false.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
GeneralEigenSolver<ScalarType>::GeneralEigenSolver() {
    mMaxIter = 1000;
    isMatrixInit = false;
}

/**
 * Destructor.
 * Destructs the pointer to the matrix operator, #mMatrix.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
GeneralEigenSolver<ScalarType>::~GeneralEigenSolver() {
    // TODO investigate the destructor
    delete mMatrix;
}

/*==============================================*//*=
 =  Setters
 ================================================*/

/**
 * Setter for #mMatrix.
 * Sets #mMatrix as a pointer to A, also sets #isMatrixInit to true.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param A Matrix the solver should point to.
 */

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMatrix(MatrixType<ScalarType>& A) {
    mMatrix = &A;
    isMatrixInit = true;
}

/**
 * Setter for #mMaxIter
 * Modifies the value of #mMaxIter. In case maxIter is negative, throws an
 * an <a href="https://en.cppreference.com/w/cpp/error/invalid_argument"> std::invalid_argument </a> exception.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param maxIter Maximum number of iteration for the algorithm.
 */
template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMaxIter(const int maxIter) {
    if (maxIter <= 0) {
        throw std::invalid_argument("INVALID MAXIMUM ITERATIONs: The maximum number of iterations need to be positive\n");
    }
    mMaxIter = maxIter;
}

/*==============================================*//*=
 =  Getters
 ================================================*/

/**
 * Getter for #mMaxIter.
 * Returns #mMaxIter.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return int #mMaxIter, Maximum number of iteration for the algorithm.
 */
template<typename ScalarType>
int GeneralEigenSolver<ScalarType>::getMaxIter() const{
    return mMaxIter;
}
/**
 * Getter for #isMatrixInit.
 * Returns #isMatrixInit.
 * If false, #mMatrix has not been initialized. In this case, the solve() method should not be called.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return bool #isMatrixInit, true is #mMatrix has been initialized.
 */
template<typename ScalarType>
bool GeneralEigenSolver<ScalarType>::getIsMatrixInit() const {
    return isMatrixInit;
}

#endif //EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
