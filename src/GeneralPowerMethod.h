
#ifndef EIGENVALUE_PROJECT_GENERALPOWERMETHOD_H
#define EIGENVALUE_PROJECT_GENERALPOWERMETHOD_H

#include "GeneralEigenSolver.h"
#include "Exceptions/UninitializedSolver.h"

/**
 * Abstract class for Power Methods.
 *
 * Mother class for eigenvalue solvers which are variants of the power method.
 * The power method is given by:
 * > Given \f$ A\in \mathbb{C}^{n \times n} \f$: <br>
 *
 * > Choose a starting vector \f$ x_0 \in \mathbb{C}^n \f$ <br>
 * > \f$k = 0\f$ <br>
 * > Repeat while \f$k <\f$#mMaxIter: <br>
 * > - \f$k = k+1\f$ <br>
 * > - Compute \f$y_k = Ax_{k-1}\f$ <br>
 * > - Normalize \f$x_k = \frac{y_k}{||y_k||_2}\f$ <br>
 * > - Approximate the eigenvalue as \f$ \lambda = x_k^T A x_k \f$
 * > - End if \f$ ||Ax_k - \lambda x_k|| <\f$#mThreshold
 *
 * > End at convergence or when the number of iterations has exceeded #mMaxIter.
 *
 * The matrix \f$A\f$ can vary depending on the solver and the eigenvalue wished.
 *
 * The resulting \f$\lambda\f$ is the largest eigenvalue of \f$A\f$, with respect to their absolute values.
 *
 * This class implements the methods to set and get #mThreshold, and the starting vector #mEigenVector.
 *
 * It also implements the method to initialize #mEigenVector randomly, initRandomEigenVector(),
 *  and to use the power method, solve(MatrixType<ScalarType>& A), which takes a pointer to the matrix \f$A\f$ as parameter.
 * Note that when calling solve(MatrixType<ScalarType>& A), where #mEigenVector is randomly initialized if it has not been done before-hand,
 * an <a href="https://en.cppreference.com/w/cpp/error/invalid_argument">std::invalid_argument</a> exception is thrown if #mMatrix has not been initialized.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */

template <typename ScalarType>
class GeneralPowerMethod: public GeneralEigenSolver<ScalarType>{
public:
    // Constructors and Destructors
    GeneralPowerMethod();
    ~GeneralPowerMethod();

    /** @name Setters
     *
     */
    ///@{
    void setEigenVector(VectorType<ScalarType> V);
    void setThreshold(double threshold);
    ///@}

    /** @name Getters
     *
     */
    ///@{
    double getThreshold() const;
    bool getIsVectorInit() const;
    VectorType<ScalarType> getEigenVector() const;
    ///@}


    /** @name Other functions
     *
     */
    ///@{
    //virtual ScalarType solve() = 0;
    void initRandomEigenVector();
    ///@}

protected:

    VectorType<ScalarType> mEigenVector;    /**< Starting vector for the algorithm. In the end becomes the eigenvector. */
    bool isVectorInit;                      /**< Boolean, true if the vector has been initialized */
    double mThreshold;                      /**< Threshold for the solve() method */

    /** @name Other functions
     *
     */
    ///@{
    ScalarType solve(MatrixType<ScalarType>& A);
    ///@}
};

/*==========================================*//*=
 =  Constructors and Destructors
 ==============================================*/

/**
 * Initialization constructor.
 * Basic constructor, sets #mThreshold to 1e-15 and #isVectorInit to false
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
GeneralPowerMethod<ScalarType>::GeneralPowerMethod() {
    mThreshold = 1e-15;
    isVectorInit = false;
}

/**
 * Destructor.
 * Calls the destructor of GeneralEigenSolver.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
GeneralPowerMethod<ScalarType>::~GeneralPowerMethod() {}

/*==============================================*//*=
 =  Setters
 ================================================*/

/**
 * Setter for #mEigenVector.
 * Sets #mEigenVector as a copy of V, also sets #isVectorInit to true.
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param V VectorType<ScalarType> (Eigen::Vector<ScalarType, -1,-1>), Starting vector for the power method.
 */
template<typename ScalarType>
void GeneralPowerMethod<ScalarType>::setEigenVector(const VectorType<ScalarType> V) {
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
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param threshold double, Threshold for the power method algorithm (see the detailed description GeneralPowerMethod).
 */
template<typename ScalarType>
void GeneralPowerMethod<ScalarType>::setThreshold(const double threshold) {
    if (threshold < 0) {
        throw std::invalid_argument("INVALID THRESHOLD: please choose a positive threshold "
                                         "(current value " + std::to_string(threshold) + ").\n");
    }
    if (threshold < 1e-20) {
        std::cerr << "WARNING: Threshold < 1e-20, the computation might take a long time.\n";
    }
    mThreshold = threshold;
}

/*==============================================*//*=
 =  Getters
 ================================================*/

/**
 * Getter for #mThreshold.
 * Returns #mThreshold
 *
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return double, Threshold for the power method algorithm (see the detailed description GeneralPowerMethod).
 */
template<typename ScalarType>
double GeneralPowerMethod<ScalarType>::getThreshold() const {
    return mThreshold;
}

 /**
 * Getter for #isVectorInit.
 * Returns #isVectorInit.
 * If false, #mEigenVector has not been initialized. If solve() is called, #mEigenVector will be initialized randomly.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @return bool #isVectorInit, true if #mEigenVector has been initialized.
 */
template<typename ScalarType>
bool GeneralPowerMethod<ScalarType>::getIsVectorInit() const {
    return isVectorInit;
}

/**
 * Getter for #mEigenVector.
 * Returns a copy of #mEigenVector.
 * @tparam VectorType naming for Eigen::Vector<ScalarType, -1, -1>
 * @return VectorType, the starting vector, or the eigenvector after the power method algorithm.
 */
template<typename ScalarType>
VectorType<ScalarType> GeneralPowerMethod<ScalarType>::getEigenVector() const {
    return this->mEigenVector;
}

/*==========================================*//*=
 = Other functions
 ==============================================*/

/**
 * Solves the eigenvalue problem \f$Ax=\lambda x\f$.
 * For more information about the algorithm, see the detailed description, GeneralPowerMethod.
 * If #mMatrix has not been initialized (see GeneralEigenSolver), throws an UninitializedSolver exception.
 * In case #mEigenVector ends up in the null space of \f$A\f$, throws an <a href="https://en.cppreference.com/w/cpp/error/invalid_argument">std::invalid_argument</a> exception.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 * @param A pointer to matrix operator \f$A\f$
 * @return ScalarType The resulting eigenvalue.
 * @throws UninitializedSolver, invalid_argument
 */
template <typename ScalarType>
ScalarType GeneralPowerMethod<ScalarType>::solve(MatrixType<ScalarType>& A) {

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
        auto temp = A * this->mEigenVector;
        auto temp_norm = temp.norm();
        if (temp_norm < 1e-15) {
            throw std::invalid_argument("INVALID STARTING VECTOR: The guessed eigenvector "
                                        "is currently in the null "
                                        "space of the matrix,\ntry to initialize "
                                        "the vector with a different value "
                                        "(e.g using PowerMethod<ScalarType>::setEigenVector(...))\n");
        }
        this->mEigenVector = temp / temp_norm; // Update the

        // Compute the corresponding eigenvalue
        auto Av = A * this->mEigenVector;
        lambda = this->mEigenVector.transpose().conjugate() * Av;
        error = (Av - lambda *  this->mEigenVector).norm();

        if (error < threshold) {
            return lambda;
        }
    }
    std::cout << "The method did not converge after " << iter << " iterations\n";
    return lambda;
}

/**
 * Initialize #mEigenVector.
 * Initialize #mEigenVector with a uniform distribution in [-1,1].
 * The number of rows of #mEigenVector is changed to the number of columns of #mMatrix.
 * If #isMatrixInit is false, throws an UninitializedSolver exception.
 * @tparam ScalarType The type of the scalars used in the eigenvalue problem (usually of type int, double or <a href="https://en.cppreference.com/w/cpp/numeric/complex">std::complex</a>)
 */
template<typename ScalarType>
void GeneralPowerMethod<ScalarType>::initRandomEigenVector() {
    if (this->isMatrixInit) {
        mEigenVector.resize(this->mMatrix->cols(), 1);
        mEigenVector.setRandom();
        isVectorInit = true;
    }
        // If the matrix has not been initialized, the dimensions of #mEigenVector cannot be known
    else {
        throw UninitializedSolver(
                "matrix",
                "please initialize with GeneralEigenSolver<typename ScalarType>::setMatrix");
    }
}

#endif //EIGENVALUE_PROJECT_GENERALPOWERMETHOD_H
