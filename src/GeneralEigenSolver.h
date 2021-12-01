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

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename ScalarType>
using VectorType = Eigen::Vector<ScalarType, Eigen::Dynamic>;

template<typename ScalarType>
class GeneralEigenSolver {
private:
    // TODO make mMatrix a pointer to a given matrix
    MatrixType<ScalarType> mMatrix;
    VectorType<ScalarType> mEigenVector;
    double mThreshold;
    unsigned int mMaxIter;

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

    // TODO add a pure virtual solve method
    virtual ScalarType solve(bool logging = false) = 0;
};

/********************************************//**
 *  Constructors and Destructors
 ***********************************************/

template<typename ScalarType>
GeneralEigenSolver<ScalarType>::GeneralEigenSolver() {
    mEigenVector.setRandom();
    mThreshold = 1e-15;
    mMaxIter = 1000;
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
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setEigenVector(const VectorType<ScalarType> V) {
    mEigenVector.resize(V.rows(), V.cols());
    mEigenVector = V;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setThreshold(const double threshold) {
    // TODO add exceptions (>0 and not too small) ValueErrorException
    assert(threshold > 0);
    mThreshold = threshold;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMaxIter(const int maxIter) {
    // TODO add exceptions
    assert(maxIter > 0);
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


#endif //EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
