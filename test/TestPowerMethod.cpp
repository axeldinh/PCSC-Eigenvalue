//
// Created by axeld on 01/12/2021.
//

/**
 * We implement here the tests for the PowerMethod class.
 */

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"
#include "../src/Exceptions/UninitializedSolver.h"

class TestPowerMethod: public ::testing::Test {
    /**
     * Test suite for the Power method class
     * At initialization, a PowerMethod is instantiated
     * with a diagonal Eigen::Matrix<double,5,5>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new PowerMethod<double>;
        Eigen::MatrixXd M(n,n);
        M.setZero();
        M.diagonal() = Eigen::Vector<double,n>::Constant(3.);
        solver->setMatrix(M);
    }

    void TearDown() override {
        // TODO inspect why cannot delete solver
        //delete solver;
    }

    PowerMethod<double>* solver;
    static constexpr int n {5};
};

TEST_F(TestPowerMethod, diagonalMatrix) {
    /**
     * Tests if give A = 3*I, then the solve method return 3.
     */
     Eigen::VectorXd V(n);
     V.setOnes();
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
     ASSERT_THROW(solver->solve(), std::invalid_argument);
}

TEST_F(TestPowerMethod, uninitMatrix) {
    /**
     * Checks that trying to call .solve()
     * without matrix initialization returns an
     * UninitializedSolver Exception
     */
     solver = new PowerMethod<double>(); // Need new solver
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
     std::string expected = "The Power Method did not converge after ";
     expected += std::to_string(solver->getMaxIter()) + " iterations\n";
     EXPECT_STREQ(output.c_str(), expected.c_str());
}
