//
// Created by axeld on 01/12/2021.
//



#ifndef EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
#define EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H

#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include "Exceptions/UninitializedSolver.h"

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename ScalarType>
using VectorType = Eigen::Vector<ScalarType, Eigen::Dynamic>;

/**
 * Abstract class for EigenValue solvers.
 * Mother class for classes used to solve eigenvalues
 * problem of the form Ax=\lambda b, where A \in C^nn is known.
 * This abstract class is based on the assumption
 * that A and b are declared using the Eigen library
 * (<a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">link</a>).
 *
 * A more general class could have been used (with templated
 * operators and vectors), but this approach seemed too abstract
 * for the scope of the project.
 *
 * All methods are implemented in the header file as the compiler needs the code when calling on
 * a particular template (e.g GeneralEigenSolver<double>)
 *
 */
template<typename ScalarType>
class GeneralEigenSolver {

protected:
    // TODO make mMatrix a pointer to a given matrix
    MatrixType<ScalarType> mMatrix;         /**< Matrix operator */
    VectorType<ScalarType> mEigenVector;    /**< Eigenvector */
    double mThreshold;                      /**< Threshold for the solve() method */
    unsigned int mMaxIter;                  /**< Maximum number of iterations for the solve() method */

    bool isMatrixInit;                      /**< To know if the matrix has been initialized */
    bool isVectorInit;                      /**< To know if the vector has been initialized */

public:

    // Constructors and Destructors
    GeneralEigenSolver();
    ~GeneralEigenSolver();

    /** @name Setters
     *
     */
    ///@{
    void setMatrix(MatrixType<ScalarType> A);
    void setEigenVector(VectorType<ScalarType> V);
    void setThreshold(double threshold);
    void setMaxIter(int maxIter);
    ///@}

    /** @name Getters
     *
     */
    ///@{
    double getThreshold() const;
    int getMaxIter() const;
    bool getIsMatrixInit() const;
    bool getIsVectorInit() const;
    ///@}


    /** @name Other functions
     *
     */
     ///@{
    virtual ScalarType solve() = 0;
    void initRandomEigenVector();
    ///@}
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

/**
 * Initialization constructor.
 * Basic constructor, sets #mThreshold to 1e-15, #mMaxIter to 1000, #isMatrixInit and #isVectorInit to false
 * @tparam ScalarType
 */
template<typename ScalarType>
GeneralEigenSolver<ScalarType>::GeneralEigenSolver() {
    mThreshold = 1e-15;
    mMaxIter = 1000;
    isMatrixInit = false;
    isVectorInit = false;
}

// TODO Doc
template<typename ScalarType>
GeneralEigenSolver<ScalarType>::~GeneralEigenSolver() {
    // TODO investigate the destructor
}

/*==============================================*//*=
 =  Setters
 ================================================*/

// TODO test this method
/**
 * Setter for #mMatrix.
 * Sets #mMatrix as a copy of A, also sets #isMatrixInit to true.
 *
 * @tparam ScalarType
 * @param A
 */

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMatrix(const MatrixType<ScalarType> A) {
    mMatrix.resize(A.rows(), A.cols());
    mMatrix = A;
    isMatrixInit = true;
}

// TODO test this method
/**
 * Setter for #mEigenVector.
 * Sets #mEigenVector as a copy of V, also sets isVectorInit to true.
 *
 * @tparam ScalarType
 * @param V
 */
template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setEigenVector(const VectorType<ScalarType> V) {
    mEigenVector.resize(V.rows(), V.cols());
    mEigenVector = V;
    isVectorInit = true;
}

/**
 * Setter for #mThreshold.
 * Modifies the value of #mThreshold. In case threshold is negative, throws
 * an <a href="https://en.cppreference.com/w/cpp/error/invalid_argument"> std::invalid_argument </a> exception.
 * In case threshold is less than 1e-20, prints a warning.
 *
 * @tparam ScalarType
 * @param threshold double
 */
template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setThreshold(const double threshold) {
    if (threshold < 0) {
        throw std::invalid_argument("INVALID THRESHOLD: please choose a positive threshold "
                                         "(current value " + std::to_string(threshold) + ").\n");
    }
    if (threshold < 1e-20) {
        std::cerr << "WARNING: Threshold < 1e-16, the computation might take a long time.\n";
    }
    mThreshold = threshold;
}

/**
 * Setter for #mMaxIter
 * Modifies the value of #mMaxIter. In case maxIter is negative, throws an
 * an <a href="https://en.cppreference.com/w/cpp/error/invalid_argument"> std::invalid_argument </a> exception.
 *
 * @tparam ScalarType
 * @param maxIter
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
 * Getter for #mThreshold
 * Returns #mThreshold
 *
 * @tparam ScalarType
 * @return
 */
template<typename ScalarType>
double GeneralEigenSolver<ScalarType>::getThreshold() const {
    return mThreshold;
}

/**
 * Getter for #mMaxIter
 * Returns #mMaxIter
 * @tparam ScalarType
 * @return
 */
template<typename ScalarType>
int GeneralEigenSolver<ScalarType>::getMaxIter() const{
    return mMaxIter;
}
/**
 * Getter for #isMatrixInit.
 * Returns #isMatrixInit.
 * If true, #mMatrix has not been initialized. In this case, the solve() method should not be called.
 * @tparam ScalarType
 * @return
 */
template<typename ScalarType>
bool GeneralEigenSolver<ScalarType>::getIsMatrixInit() const {
    return isMatrixInit;
}

/**
 * Getter for #isVectorInit.
 * Returns #isVectorInit.
 * If true, #mEigenVector has not been initialized. In this case, the solve() method should not be called.
 * @tparam ScalarType
 * @return
 */
template<typename ScalarType>
bool GeneralEigenSolver<ScalarType>::getIsVectorInit() const {
    return isVectorInit;
}

/*==========================================*//*=
 = Other functions
 ==============================================*/

/**
 * Initialize #mEigenVector.
 * Initialize #mEigenVector with a uniform distribution in [-1,1].
 * The number of rows of #mEigenVector is changed to the number of columns of #mMatrix.
 * If #isMatrixInit is false, throws an UninitailizedSolver exception.
 * @tparam ScalarType
 */
template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::initRandomEigenVector() {
    if (isMatrixInit) {
        mEigenVector.resize(mMatrix.cols(), 1);
        mEigenVector.setRandom();
        isVectorInit = true;
    }
    else {
        throw UninitializedSolver(
                "matrix",
                "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }
}


#endif //EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
