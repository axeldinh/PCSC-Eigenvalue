#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/InversePowerMethod.h"

class TestInversePowerMethod: public ::testing::Test {
    /**
     * Test suite for the InversePowerMethod class
     * At initialization, a ShiftedInversePowerMethod is instantiated
     * with a diagonal Eigen::Matrix<double,5,5>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new InversePowerMethod<double>;
        M.resize(n,n);
        M << 2, -1,
                -1, 2;
        solver->setMatrix(M);
    }

    void TearDown() override {
    }

    InversePowerMethod<double>* solver;
    Eigen::MatrixXd M;
    static constexpr int n {2};
};

TEST_F(TestInversePowerMethod, solveReturnsSmallestEigenvalue) {
    /**
     * Checks if the InversePowerMethod returns the smallest eigenvalue
     */

    auto lambda = solver->solve();
    ASSERT_NEAR(lambda, 1., 1e-15);
}

TEST_F(TestInversePowerMethod, noInitMakesRandomInitOfVector) {
    /**
     * Checks that calling solve() without initializing EigenVector
     * makes it being initialized through initRandomEigenVector()
     */

    solver->setMaxIter(1);
    auto lambda = solver->solve();
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestInversePowerMethod, noMatrixInitReturnsException) {
    /**
     * Checks that calling solve without matrix initialization
     * returns an UninitializedSolver Exception
     */

    solver = new InversePowerMethod<double>();
    ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestInversePowerMethod, nullSpaceEigenVectorReturnsException) {
    /**
     * Checks that if the eigenvector is initialized in the null space
     * of (A-shift*I)^-1 then solve() returns an std::invalid_argument exception
     */

    Eigen::VectorXd V(n);
    V.setZero();
    solver->setEigenVector(V);
    ASSERT_THROW(solver->solve(), std::invalid_argument);
}

TEST_F(TestInversePowerMethod, noConvergencePrintsToScreen) {
    /**
     * Checks if the absence of convergence prints a warning to screen
     */
    solver = new InversePowerMethod<double>;
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
