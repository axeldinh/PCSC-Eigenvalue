//
// Created by axeld on 01/12/2021.
//

/**
 * Implementation of an abstract class for the EigenValue Solver
 * All methods are implemented in the header file as the compiler needs the code when calling on
 * a particular template (e.g GeneralEigenSolver<double>)
 *
 */

#ifndef EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
#define EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H

#include <Eigen/Dense>
#include <cassert>
#include <iostream>

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename ScalarType>
using VectorType = Eigen::Vector<ScalarType, Eigen::Dynamic>;

template<typename ScalarType>
class GeneralEigenSolver {

protected:
    // TODO make mMatrix a pointer to a given matrix
    MatrixType<ScalarType> mMatrix;
    VectorType<ScalarType> mEigenVector;
    double mThreshold;
    unsigned int mMaxIter;

    bool isMatrixInit;
    bool isVectorInit;

public:

    // Constructors and Destructors
    GeneralEigenSolver();
    ~GeneralEigenSolver();

    // Setters
    void setMatrix(MatrixType<ScalarType> A);
    void setEigenVector(VectorType<ScalarType> V);
    void setThreshold(double threshold);
    void setMaxIter(int maxIter);

    // Getters
    double getThreshold() const;
    int getMaxIter() const;
    bool getIsMatrixInit() const;
    bool getIsVectorInit() const;

    // Other functions
    virtual ScalarType solve() = 0;
    void initRandomEigenVector();
};

/********************************************//**
 *  Constructors and Destructors
 ***********************************************/

template<typename ScalarType>
GeneralEigenSolver<ScalarType>::GeneralEigenSolver() {
    mThreshold = 1e-15;
    mMaxIter = 1000;
    isMatrixInit = false;
    isVectorInit = false;
}

template<typename ScalarType>
GeneralEigenSolver<ScalarType>::~GeneralEigenSolver() {
    // TODO investigate the destructor
}

/********************************************//**
 *  Setters
 ***********************************************/

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMatrix(const MatrixType<ScalarType> A) {
    mMatrix.resize(A.rows(), A.cols());
    mMatrix = A;
    isMatrixInit = true;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setEigenVector(const VectorType<ScalarType> V) {
    mEigenVector.resize(V.rows(), V.cols());
    mEigenVector = V;
    isVectorInit = true;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setThreshold(const double threshold) {
    if (threshold < 0) {
        throw std::invalid_argument("INVALID THRESHOLD: please choose a positive threshold "
                                         "(current value " + std::to_string(threshold) + ").\n");
    }
    if (threshold < 1e-16) {
        std::cerr << "WARNING: Threshold < 1e-16, the computation might take a long time.\n";
    }
    mThreshold = threshold;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMaxIter(const int maxIter) {
    if (maxIter <= 0) {
        throw std::invalid_argument("INVALID MAXIMUM ITERATIONs: The maximum number of iterations need to be positive\n");
    }
    mMaxIter = maxIter;
}

/********************************************//**
 *  Getters
 ***********************************************/

template<typename ScalarType>
double GeneralEigenSolver<ScalarType>::getThreshold() const {
    return mThreshold;
}

template<typename ScalarType>
int GeneralEigenSolver<ScalarType>::getMaxIter() const{
    return mMaxIter;
}

template<typename ScalarType>
bool GeneralEigenSolver<ScalarType>::getIsMatrixInit() const {
    return isMatrixInit;
}

template<typename ScalarType>
bool GeneralEigenSolver<ScalarType>::getIsVectorInit() const {
    return isVectorInit;
}

/********************************************//**
 *  Other functions
 ***********************************************/

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::initRandomEigenVector() {
    if (isMatrixInit) {
        mEigenVector.resize(mMatrix.cols(), 1);
        mEigenVector.setRandom();
        isVectorInit = true;
    }
    else {
        // TODO transform it in an exception
        std::cerr << "Matrix not initialized, please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix";
    }
}


#endif //EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
