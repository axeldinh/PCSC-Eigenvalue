//
// Created by axeld on 01/12/2021.
//

/**
 * Implementation of an abstract class for the EigenValue Solver
 * All methods are implemented in the header file as the compiler needs the code when calling on
 * a particular template (e.g GeneralEigenSolver<double>)
 *
 * TODO Make the class pure virtual
 */

#ifndef EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
#define EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H

#include <Eigen/Dense>

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

template<typename ScalarType>
using VectorType = Eigen::Vector<ScalarType, Eigen::Dynamic>;

template<typename ScalarType>
class GeneralEigenSolver {
private:
    MatrixType<ScalarType> mMatrix;
    VectorType<ScalarType> mEigenVector;
    double mThreshold;
    unsigned int mMaxIter;

public:

    // Constructors and Destructors
    GeneralEigenSolver();
    ~GeneralEigenSolver();

    // Setters
    void setThreshold(double threshold);
    void setMaxIter(int maxIter);

    // Getters
    double getThreshold();
    int getMaxIter();
};

template<typename ScalarType>
GeneralEigenSolver<ScalarType>::GeneralEigenSolver() {
    mMatrix.fill(1.);
}

template<typename ScalarType>
GeneralEigenSolver<ScalarType>::~GeneralEigenSolver() {}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setThreshold(double threshold) {
    mThreshold = threshold;
}

template<typename ScalarType>
void GeneralEigenSolver<ScalarType>::setMaxIter(int maxIter) {
    mMaxIter = maxIter;
}

template<typename ScalarType>
double GeneralEigenSolver<ScalarType>::getThreshold() {
    return mThreshold;
}

template<typename ScalarType>
int GeneralEigenSolver<ScalarType>::getMaxIter() {
    return mMaxIter;
}


#endif //EIGENVALUE_PROJECT_GENERALEIGENSOLVER_H
