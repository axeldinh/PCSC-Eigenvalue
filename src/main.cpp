
/**
 * \mainpage EigenValue Solver
 *
 * \tableofcontents
 *
 * \section improvements To Be Improved
 * - Could have a warning to prevent the user from using QRMethod with complex matrices.
 * - More options to QRMethod::solve_all() could allow the user to choose the ordering of the eigenvalues
 * - A logging option could be added to monitor the convergence of the methods.
 * - The current code does not handle sparse matrices, as a lot of matrices used in practice, e.g. in physics,
 * the usage of sparse matrices can be preferred and should be implemented in the future.
 *
 * ===
 *
 * \section intro Introduction
 *
 * This code aims at solving eigenvalue problems of the form \f$Ax=\lambda x\f$,
 * where \f$A\in \mathbb{C}^{n\times n}\f$, \f$x \in \mathbb{C}^n\f$ and \f$\lambda \in \mathbb{C}\f$
 * (The vector space can also be \f$\mathbb{Z}\f$ or \f$\mathbb{R}\f$).
 *
 * The code uses the following libraries:
 * - Eigen 3.4.0 (<a href="https://eigen.tuxfamily.org/index.php?title=Main_Page">Link</a>).
 * - GoogleTest 1.11.0 (<a href="https://github.com/google/googletest">Link</a>).
 *
 * 4 classes have been implemented to solve this problem in 4 different algorithms:
 * - PowerMethod: Solves the eigenvalue problem using the standard power method (see section \ref power_method).
 * - InversePowerMethod: Solves the eigenvalue problem using the inverse power method (see section \ref inverse_power_method).
 * - ShiftedInversePowerMethod: Solves the eigenvalue problem using the shifted inverse power method (see section \ref shifted_inverse_power_method).
 * - QRMethod: Solves the eigenvalue problem using the QR method (see section \ref qr_method).
 *
 * More details can be found on the algorithms used in the detailed description of the classes.
 *
 * ===
 *
 * \section solvers Solvers
 *
 * An abstraction of an eigenvalue solver is implemented under the GeneralEigenSolver class.
 * All solvers inherit (directly or indirectly) from this class. This class is templated such that any type of scalar can be used,
 * as long as Eigen::Matrix<ScalarType,-1,-1> can be called.
 *
 * Two attributes are initialized with this class:
 * - The maximum number of iterations.
 * - The matrix \f$A\f$ from which we wish to extract the eigenvalues.
 *
 * Then 2 types of daughter classes inherit from GeneralEigenSolver:
 * - GeneralPowerMethod for the power methods and its variants (see the next section \ref general_power_method).
 * - QRMethod (see section \ref qr_method).
 *
 * ===
 *
 * \subsection general_power_method Power Method and Variants
 *
 * An abstraction of a solver using the power method is implemented under the GeneralPowerMethod class.
 * All variants of the power method inherit from this class.
 *
 * Two attributes are initialized with this class:
 * - The threshold for the power method. (GeneralPowerMethod#mThreshold)
 * - The starting vector for the power method. After solving, this variable becomes the corresponding eigenvector.  (GeneralPowerMethod#mEigenVector)
 *
 * Two methods are also implemented in this class:
 * - GeneralPowerMethod#solve(MatrixType<ScalarType>& A) which solves an eigenvalue problem using the power method.
 * - GeneralPowerMethod#initRandomEigenVector() which randomly initialize the starting vector.
 *
 * The variants then use solve(MatrixType<ScalarType>& A) with different matrices depending on the desired algorithm.
 *
 * Then 3 daughter classes inherit from GeneralPowerMethod:
 * - PowerMethod, the standard power method. Computes the highest eigenvalue, in absolute value (see section \ref power_method).
 * - InversePowerMethod. Computes the smallest eigenvalue, in absolute value (see section \ref inverse_power_method).
 * - ShiftedInversePowerMethod. Computes the closest eigenvalue to some shift \f$\sigma\f$ (see section \ref shifted_inverse_power_method).
 *
 * ---
 *
 * \subsubsection power_method Power Method
 *
 * The PowerMethod class implements the standard power method.
 * Given a matrix \f$A \in \mathbb{C}^{n\times n}\f$ with eigenvalues
 * \f$|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|\f$, then the returned eigenvalue is \f$\lambda_1\f$.
 * This will not be the case if the starting vector is in the null space of \f$A\f$ or an eigenvector of \f$A\f$.
 *
 * ---
 *
 * \subsubsection inverse_power_method Inverse Power Method
 *
 * The InversePowerMethod class implements the inverse power method algorithm.
 * Given a matrix \f$A \in \mathbb{C}^{n\times n}\f$ with eigenvalues
 * \f$|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|\f$, then the returned eigenvalue is \f$\lambda_n\f$.
 * This will not be the case if the starting vector is in the null space of \f$A\f$ or an eigenvector of \f$A\f$.
 *
 * This class assumes that \f$A\f$ is invertible.
 *
 * ---
 *
 * \subsubsection shifted_inverse_power_method Shifted Inverse Power Method
 *
 * The ShiftedInversePowerMethod class implements the inverse power method algorithm.
 * Given a matrix \f$A \in \mathbb{C}^{n\times n}\f$ with eigenvalues
 * \f$|\lambda_1| \ge |\lambda_2| \ge \cdots \ge |\lambda_n|\f$,
 * then the returned eigenvalue is the one closest to some shift (ShiftedInversePowerMethod#mShift).
 * This will not be the case if the starting vector is in the null space of \f$A\f$ or an eigenvector of \f$A\f$.
 *
 * Here, the shift parameter, ShiftedInversePowerMethod#mShift, must be initialized before calling ShiftedInversePowerMethod#solve().
 *
 * ===
 *
 *  \subsection qr_method QR Method
 *
 *  The QRMethod class implements a solver using the QR method.
 *  In order for the algorithm to retrieve the correct eigenvalues,
 *  the matrix \f$A\f$ must real and not complex.
 *
 *  The QRMethod#solve_all() method returns a Eigen::Vector<ScalarType> containing the eigenvalues in descending order.
 *  It is also possible to retrieve specific eigenvalues using QRMethod#solve(int n) and QRMethod#solve().
 *  The former will retrieve the n-th largest eigenvalue while the latter the largest one.
 *
 *  Here, special care should be taken in choosing the GeneralEigenSolver#mMaxIter attribute,
 *  as the algorithm \b will iterate this number of times.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include "GeneralEigenSolver.h"
#include "PowerMethod.h"
#include "ShiftedInversePowerMethod.h"
#include "InversePowerMethod.h"
#include "QRMethod.h"

//https://stackoverflow.com/questions/865668/parsing-command-line-arguments-in-c

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return nullptr;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

Eigen::MatrixXd createDoubleMatrix(int n) {
    Eigen::MatrixXd Q(n,n);
    Q.setRandom();
    Q = Eigen::HouseholderQR<Eigen::MatrixXd>(Q).householderQ();

    Eigen::MatrixXd A(n,n);
    A.setZero();
    A.diagonal() = Eigen::VectorXd::LinSpaced(n, 1, n);
    A = Q*A*Q.transpose();

    return A;
}

Eigen::MatrixXcd createComplexMatrix(int n) {

    Eigen::MatrixXcd Q(n,n);
    Q.setRandom();
    Q = Eigen::HouseholderQR<Eigen::MatrixXcd>(Q).householderQ();

    Eigen::MatrixXcd A(n,n);
    A.setZero();
    for (int i = 0; i < n; i++) {
        A(i,i) = std::complex<double>(0., (double) i+1.);
    }
    A = Q*A*Q.transpose();

    return A;
}

template <typename ScalarType>
GeneralEigenSolver<ScalarType>* createSolver(const std::string& algorithm, int maxIter, double threshold, double shift) {

    GeneralEigenSolver<ScalarType>* solver;

    if (algorithm == "pw") {
        solver =  new PowerMethod<ScalarType>();
        if (threshold!=-1.) {
            ((PowerMethod<ScalarType> *) solver)->setThreshold(threshold);
        }
    }
    else if (algorithm == "ipw") {
        solver = new InversePowerMethod<ScalarType>();
        if (threshold!=-1.) {
            ((InversePowerMethod<ScalarType> *) solver)->setThreshold(threshold);
        }

    }
    else if (algorithm == "sipw") {
        solver = new ShiftedInversePowerMethod<ScalarType>();
        if (threshold!=-1.) {
            ((ShiftedInversePowerMethod<ScalarType> *) solver)->setThreshold(threshold);
        }
        if (shift!=-1.) {
            ((ShiftedInversePowerMethod<ScalarType> *) solver)->setShift(shift);
        }
        else {
            throw std::invalid_argument("Shift not given, use -shift [option] with option being a double \n");
        }
    }

    // Set directly the ScalarType to double so that the compiler do not try to compile QRMethod<std::complex<>>
    else if (algorithm == "qr") {
        return solver;
        //solver = new QRMethod<ScalarType>();
    }
    else {
        throw std::invalid_argument(
            "Wrong algorithm, use -algorithm [option] with option being"
            " one of \"qr\", \"pw\", \"ipw\", \"sipw\" \n"
        );
    }

    if (maxIter!=-1.) {
        solver->setMaxIter(maxIter);
    }

    return solver;
}

template <typename ScalarType>
void printSolverConfig(const GeneralEigenSolver<ScalarType>* pSolver, const std::string& algorithm) {

    // Print the configuration onscreen
    std::cout << "Solving the eigenvalue problem using the ";
    if (algorithm=="pw") {
        std::cout << "PowerMethod";
    }
    if (algorithm=="ipw") {
        std::cout << "InversePowerMethod";
    }
    if (algorithm=="sipw") {
        std::cout << "ShiftedInversePowerMethod";
    }
    if (algorithm=="qr") {
        std::cout << "QRMethod";
    }
    std::cout << " class with maxIter=" << pSolver->getMaxIter();
    if (algorithm=="pw") {
        std::cout << ", threshold=" << ((PowerMethod<ScalarType>*) pSolver)->getThreshold();
    }
    if (algorithm=="ipw") {
        std::cout << ", threshold=" << ((InversePowerMethod<ScalarType>*) pSolver)->getThreshold();
    }
    if (algorithm=="sipw") {
        std::cout << ", threshold=" << ((ShiftedInversePowerMethod<ScalarType>*) pSolver)->getThreshold()
                  << ", shift=" << ((ShiftedInversePowerMethod<ScalarType>*) pSolver)->getShift();
    }
    std::cout << "\n";
}

void solveDouble(const std::string& algorithm, int maxIter, double threshold, double shift, int n, int nEigenValue) {

    // First create a double matrix
    auto A = createDoubleMatrix(n);

    // Create the solver and initialize them with specific attributes
    auto solver = createSolver<double>(algorithm, maxIter, threshold, shift);

    // If algorithm is "qr" the pointer has not been initialized,
    // this is to prevent the compiler from compiling a QRMethod<std::complex<>>
    // (std::greater is not implemented for this ScalarType).
    if (algorithm == "qr") {
        solver = new QRMethod<double>();
        if (maxIter!=-1) {
            solver->setMaxIter(maxIter);
        }
    }

    // Initialize the solver's matrix pointer
    solver->setMatrix(A);

    // Print the configuration onscreen
    printSolverConfig(solver, algorithm);

    double lambda;

    try {
        // If algorithm is QR, return the n-th eigenvalue
        if (algorithm=="qr") {
            lambda = ((QRMethod<double>*)solver)->solve(nEigenValue);
        }
        else {
            lambda = solver->solve();
        }
    } catch(UninitializedSolver& e) {
        throw;
    } catch(std::invalid_argument& e) {
        throw;
    } catch(std::exception& e) {
        throw;
    }

    std::cout << "Solution: " << std::scientific << lambda << "\n";
}

void solveComplex(const std::string& algorithm, int maxIter, double threshold, double shift, int n) {

    // First create a double matrix
    auto A = createComplexMatrix(n);

    // Create the solver and initialize them with specific attributes
    auto solver = createSolver<std::complex<double>>(algorithm, maxIter, threshold, shift);

    // Initialize the solver's matrix pointer
    solver->setMatrix(A);

    // Print the configuration onscreen
    printSolverConfig(solver, algorithm);

    std::complex<double> lambda;

    try {
        lambda = solver->solve();
    } catch(UninitializedSolver& e) {
        throw;
    } catch(std::invalid_argument& e) {
        throw;
    } catch(std::exception& e) {
        throw;
    }

    std::cout << "Solution: " << std::scientific << lambda << "\n";
}

int main(int argc, char* argv[]) {

    // Options
    std::string algorithm;
    std::string scalarType = "double";
    int maxIter = -1;
    double threshold = -1.;
    double shift = -1.;
    int n = 100;
    int nEigenValue = 1;

    // Parse mandatory options
    if (cmdOptionExists(argv, argv+argc, "-algorithm")) {
        algorithm = getCmdOption(argv, argv + argc, "-algorithm");
    }
    else {
        throw std::invalid_argument(
                "Algorithm not given, use -algorithm [option] with option being"
                " one of \"qr\", \"pw\", \"ipw\", \"sipw\" \n"
        );
    }

    // Parse optional options
    if (cmdOptionExists(argv, argv+argc, "-scalartype")) {
        scalarType = getCmdOption(argv, argv+argc, "-scalartype");
    }

    if (cmdOptionExists(argv, argv+argc, "-maxiter")) {
        maxIter = std::stoi(getCmdOption(argv, argv+argc, "-maxiter"));
    }

    if (cmdOptionExists(argv, argv+argc, "-threshold")) {
        threshold = std::stod(getCmdOption(argv, argv+argc, "-threshold"));
    }

    if (cmdOptionExists(argv, argv+argc, "-shift")) {
        shift = std::stod(getCmdOption(argv, argv+argc, "-shift"));
    }

    if (cmdOptionExists(argv, argv+argc, "-matrixsize")) {
        n = std::stoi(getCmdOption(argv, argv+argc, "-matrixsize"));
    }

    if (cmdOptionExists(argv, argv+argc, "-neigenvalue") && algorithm=="qr") {
        nEigenValue = std::stoi(getCmdOption(argv, argv+argc, "-neigenvalue"));
        if (nEigenValue <= 0) {
            throw std::invalid_argument(
                    "-neigenvalue should be positive."
                    );
        }
    }

    // If the scalarType is complex and algorithm is QRMethod, throw an exception
    if (scalarType=="complex" && algorithm=="qr") {
        throw std::invalid_argument(
                "-scalartype cannot be \"complex\" if -algorithm is \"qr\""
                );
    }

    std::cout << "Matrix Size = " << n << "\n";

    // Solve the eigenvalue problem.
    if (scalarType == "double") {
        try {
            solveDouble(algorithm, maxIter, threshold, shift, n, nEigenValue);
        } catch(UninitializedSolver& e) {
            // TODO .what() does not print on screen if they come from deeper than solveDouble / createSolver.
            std::cerr<<e.what();
        } catch(std::invalid_argument& e) {
            std::cerr<<e.what();
        } catch(std::exception& e) {
            std::cerr<<e.what();
        }
    }
    else if (scalarType == "complex") {
        try {
            solveComplex(algorithm, maxIter, threshold, shift, n);
        } catch(UninitializedSolver& e) {
            // TODO .what() does not print on screen if they come from deeper than solveDouble / createSolver.
            std::cerr<<e.what();
        } catch(std::invalid_argument& e) {
            std::cerr<<e.what();
        } catch(std::exception& e) {
            std::cerr<<e.what();
        }
    }
    else {
        std::cerr << "Wrong -scalartype [option], option should be one of: "
                  << "\"double\", \"complex\"\n";
    }

    return 0;
}
