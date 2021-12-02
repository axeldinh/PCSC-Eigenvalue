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
        Eigen::Matrix<double, n, n> M;
        M.setZero();
        M.diagonal() = Eigen::Vector<double,n>::Constant(3.);
        solver.setMatrix(M);
    }

    void TearDown() override {}

    PowerMethod<double> solver;
    static constexpr int n {5};
};

TEST_F(TestPowerMethod, diagonalMatrix) {
    /**
     * Tests if give A = 3*I, then the solve method return 3.
     */
     Eigen::Vector<double,5> V;
     V.setOnes();
     auto lambda = solver.solve();
     ASSERT_NEAR(lambda, 3., 1e-15);
}

TEST_F(TestPowerMethod, nullEigenVector) {
    /**
     * Checks if a null vector sends an std::invalid_argument exception
     */
     Eigen::Vector<double, 5> V;
     V.setZero();
     solver.setEigenVector(V);
     ASSERT_THROW(solver.solve(), std::invalid_argument);
}

TEST_F(TestPowerMethod, uninitMatrix) {
    /**
     * Checks that trying to call .solve()
     * without matrix initialization returns an
     * UninitializedSolver Exception
     */
     auto solver2 = new PowerMethod<double>(); // Need new solver
     ASSERT_THROW(solver2->solve(), UninitializedSolver);
     delete solver2;
}

TEST_F(TestPowerMethod, uninitVector) {
    /**
     * Checks if the code still runs without
     * initializing the EigenVector
     * i.e the eigenvector is automatically initialized
     */
    auto lambda = solver.solve();
    ASSERT_NEAR(lambda, 3., 1e-15);
}

TEST_F(TestPowerMethod, noConvergencePrintsToScreen) {
    /**
     * Checks if the absence of convergence prints a warning to screen
     */
     auto solver2 = new PowerMethod<double>;
     Eigen::Matrix<double,n,n> M;
     M.setRandom();
     solver2->setMaxIter(1);
     solver2->setMatrix(M);
     testing::internal::CaptureStdout();
     auto lambda = solver2->solve();
     std::string output = testing::internal::GetCapturedStdout();
     std::string expected = "The Power Method did not converge after ";
     expected += std::to_string(solver2->getMaxIter()) + " iterations\n";
     EXPECT_STREQ(output.c_str(), expected.c_str());
}
