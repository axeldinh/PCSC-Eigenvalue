
#ifndef EIGENVALUE_PROJECT_QRMETHOD_H
#define EIGENVALUE_PROJECT_QRMETHOD_H

#include "GeneralEigenSolver.h"

template<typename ScalarType>
using MatrixType = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * Solve an eigenvalue problem using the QR method.
 *
 * Aims at solving the eigenvalue problem \f$Ax = \lambda x\f$ using the QR method.
 * The QR method is given as:
 * > Given \f$ A\in \mathbb{C}^{n \times n} \f$: <br>
 *
 * > Set \f$A_0 = A\f$ <br>
 * > \f$k = 0\f$ <br>
 * > Repeat while \f$k <\f$#mMaxIter: <br>
 * > - Compute the QR decomposition of \f$A_k = QR\f$ <br>
 * > - \f$k = k+1\f$ <br>
 * > - \f$A_k = Q^T A_{k-1} Q\f$ <br>
 *
 * > The eigenvalues are given in the diagonal of \f$A_k\f$.
 *
 * The QR method can be used only on real matrices, not complex.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */

template <typename ScalarType>
class QRMethod: public GeneralEigenSolver<ScalarType> {
public:

    // Constructors and destructors
    QRMethod();
    ~QRMethod();

    /** @name Solve methods
     *
     */
    ///@{
    ScalarType solve();
    ScalarType solve(int n);
    VectorType<ScalarType> solve_all();
    ///@}
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

/**
 * Basic constructor.
 * Uses the constructor of the GeneralEigenSolver class.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
QRMethod<ScalarType>::QRMethod()
        : GeneralEigenSolver<ScalarType>() {}

/**
 * Destructor.
 * Uses the destructor of the GeneralEigenSolver class.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
QRMethod<ScalarType>::~QRMethod() {}

/*==========================================*//*=
=  Solver
==============================================*/

/**
 * Call to solve(int n) with n = 1.
 * Returns the highest eigenvalue.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return ScalarType, the highest eigenvalue.
 */
template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve() {
    return this->solve(1);
}

/**
 * Returns the n-th highest eigenvalue.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param n int, the rank of the desired eigenvalue.
 * @return ScalarType, the n-th highest eigenvalue.
 */

template<typename ScalarType>
ScalarType QRMethod<ScalarType>::solve(int n) {

    if (n > GeneralEigenSolver<ScalarType>::mMatrix->rows()) {
        throw std::invalid_argument("INVALID EIGENVALUES NUMBER: n is larger than the maximum number of eigenvalues.\n");
    }

    // TODO test this exception
    if (n<1) {
        throw std::invalid_argument("INVALID EIGENVALUES NUMBER: n is smaller than 1.\n");
    }

    VectorType<ScalarType> eigenValues;

    try {
        eigenValues = this->solve_all();
    } catch(UninitializedSolver& e) {
        throw e;
    }

    return eigenValues(n-1);
}

/**
 * Solves the eigenvalue problem \f$Ax = \lambda\f$ using the QR method.
 *
 * The eigenvalues are returned into a vector, in descending order.
 *
 * Throws an UninitializedSolver exception if #mMatrix has not been initialized.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return VectorType<ScalarType> name for Eigen::Vector<ScalarType,Eigen::Dynamic,Eigen::Dynamic>, vector containing the eigenvalues of #mMatrix in descending order.
 * @throws UninitializedSolver
 */

template<typename ScalarType>
VectorType<ScalarType> QRMethod<ScalarType>::solve_all() {

    if (!this->getIsMatrixInit()) {
        throw UninitializedSolver("matrix", "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }

    // A is a copy of mMatrix, so the original matrix is not modified
    auto A = *GeneralEigenSolver<ScalarType>::mMatrix;

    for (int i = 1; i < this->getMaxIter(); i++) {
        MatrixType<ScalarType> Q = Eigen::HouseholderQR<MatrixType<ScalarType>>(A).householderQ();
        A = Q.transpose() * A * Q;
    }

    VectorType<ScalarType> eigenValues = A.diagonal();

    std::sort(eigenValues.begin(), eigenValues.end(), std::greater<ScalarType>());

    return eigenValues;
}

#endif //EIGENVALUE_PROJECT_QRMETHOD_H
