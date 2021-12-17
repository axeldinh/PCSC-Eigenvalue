#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/QRMethod.h"

class TestQRMethod: public ::testing::Test {
    /**
     * Test suite for the QRMethod class
     * At initialization, a QRMethod is instantiated
     * with a Eigen::Matrix<double,2,2>
     */
protected:
    void SetUp() override {
        solver = new QRMethod<double>();
        M.resize(n,n);
        M << 2, -1,
                -1, 2;
        solver->setMatrix(M);
    }

    void TearDown() override {}

    QRMethod<double>* solver;
    Eigen::MatrixXd M;
    static constexpr int n {2};
};

TEST_F(TestQRMethod, solveReturnsBiggestEigenvalue) {
    /**
     * Checks that the solve() method returns the biggest eigenvalue
     */
     auto lambda = solver->solve();
     ASSERT_NEAR(lambda, 3., 1e-8);
}

TEST_F(TestQRMethod, solve2ReturnsSecondBiggestEigenvalue) {
    /**
     * Checks that the solve(2) method returns the second biggest eigenvalue
     */
    auto lambda = solver->solve(2);
    ASSERT_NEAR(lambda, 1., 1e-8);
}

TEST_F(TestQRMethod, solveAllReturnsAllEigenvalue) {
    /**
     * Checks that the solveAll() method returns the right eigenvalues
     * And that they are sorted in a decreasing manner
     */
    auto lambdas = solver->solveAll();
    ASSERT_EQ(lambdas.rows(), 2);
    ASSERT_NEAR(lambdas(0), 3., 1e-8);
    ASSERT_NEAR(lambdas(1), 1., 1e-8);
}

TEST_F(TestQRMethod, negativeRankEigenvalueThrowsException) {
    /**
     * Checks that calling solve(-1) throws an std::invalid_argument exception
     */
    ASSERT_THROW(solver->solve(-1), std::invalid_argument);
}

TEST_F(TestQRMethod, tooBigRankEigenvalueThrowsException) {
    /**
     * Checks that calling solve(3) throws an std::invalid_argument exception
     * as 3 is bigger than the matrix size
     */
    ASSERT_THROW(solver->solve(3), std::invalid_argument);
}

TEST_F(TestQRMethod, uninitSolverThrowsException) {
    /**
     * Checks that calling solve(..) or solveAll()
     * without matrix initialization throws an UninitializedSolver exception
     */
    ASSERT_THROW((new QRMethod<double>())->solve(), UninitializedSolver);
    ASSERT_THROW((new QRMethod<double>())->solve(2), UninitializedSolver);
    ASSERT_THROW((new QRMethod<double>())->solveAll(), UninitializedSolver);
}



