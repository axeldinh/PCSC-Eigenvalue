/**
 * We implement here the tests for the PowerMethod class.
 */

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"

class TestPowerMethod: public ::testing::Test {
    /**
     * Test suite for the PowerMethod class
     * At initialization, a PowerMethod is instantiated
     * with a diagonal Eigen::Matrix<double,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new PowerMethod<double>();
        M.resize(n,n);
        M << 2., -1.,
            -1., 2.;
        solver->setMatrix(M);
    }

    void TearDown() override {}

    PowerMethod<double>* solver;
    Eigen::MatrixXd M;
    static constexpr int n {2};
};

TEST_F(TestPowerMethod, solveReturnsBiggestEigenvalue) {
    /**
     * Tests that if A = [[2,-1,[-1,2]], then the solve method returns 3.
     */
     Eigen::VectorXd V(n);
     V.setRandom();
     solver->setEigenVector(V);
     auto lambda = solver->solve();
     ASSERT_NEAR(lambda, 3., 1e-15);
}

TEST_F(TestPowerMethod, nullEigenVector) {
    /**
     * Checks if a null vector sends an std::invalid_argument exception
     */
     Eigen::VectorXd V(n);
     V.setZero();
     solver->setEigenVector(V);
     //auto lambda = solver->solve();
     ASSERT_THROW(solver->solve(), std::invalid_argument);
}

TEST_F(TestPowerMethod, uninitMatrix) {
    /**
     * Checks that trying to call .solve()
     * without matrix initialization returns an
     * UninitializedSolver Exception
     */
     solver = new PowerMethod<double>(); // Need new solver
     Eigen::VectorXd V(n);
     V.setRandom();
     solver->setEigenVector(V);
     ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestPowerMethod, uninitVector) {
    /**
     * Checks if the code still runs without
     * initializing the EigenVector
     * i.e the eigenvector is automatically initialized
     */
    auto lambda = solver->solve();
    ASSERT_NEAR(lambda, 3., 1e-15);
}

TEST_F(TestPowerMethod, noConvergencePrintsToScreen) {
    /**
     * Checks if the absence of convergence prints a warning to screen
     */
     solver = new PowerMethod<double>;
     Eigen::MatrixXd M(n,n);
     M.setRandom();
     solver->setMaxIter(1);
     solver->setMatrix(M);
     testing::internal::CaptureStdout();
     auto lambda = solver->solve();
     std::string output = testing::internal::GetCapturedStdout();
     std::string expected = "The method did not converge after ";
     expected += std::to_string(solver->getMaxIter()) + " iterations\n";
     EXPECT_STREQ(output.c_str(), expected.c_str());
}

class TestPowerMethodComplex: public ::testing::Test {
    /**
     * Test suite for the PowerMethod class with std::complex<double> values.
     * At initialization, a PowerMethod is instantiated
     * with a diagonal Eigen::Matrix<std::complex<double>,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new PowerMethod<std::complex<double>>();
        M.resize(n,n);
        M << std::complex<double>(0.,2.), std::complex<double>(0.,-1.),
                std::complex<double>(0.,-1.), std::complex<double>(0.,2.);
        solver->setMatrix(M);
    }

    void TearDown() override {}

    PowerMethod<std::complex<double>>* solver;
    Eigen::MatrixXcd M;
    static constexpr int n {2};
};

TEST_F(TestPowerMethodComplex, solveReturnsBiggestEigenvalue) {
    /**
     * Tests that if A = [[2i,-1i,[-1i,2i]], then the solve method returns 3i.
     */
    Eigen::VectorXcd V(n);
    V.setRandom();
    solver->setEigenVector(V);
    auto lambda = solver->solve();
    ASSERT_NEAR(std::abs(lambda-std::complex<double>(0.,3.)), 0.0, 1e-15);
}
