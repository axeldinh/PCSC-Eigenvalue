#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/ShiftedPowerMethod.h"

class TestShiftedPowerMethod: public ::testing::Test {
    /**
     * Test suite for the ShiftedInversePowerMethod class
     * At initialization, a ShiftedInversePowerMethod is instantiated
     * with a Eigen::Matrix<double,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new ShiftedPowerMethod<double>();
        M.resize(n,n);
        M << 2, -1,
                -1, 2;
        solver->setMatrix(M);
    }

    void TearDown() override {}

    ShiftedPowerMethod<double>* solver;
    Eigen::MatrixXd M;
    static constexpr int n {2};
};

TEST_F(TestShiftedPowerMethod, setShiftMakesIsShiftInitToTrue) {
    /**
     * Checks if initializing the shift sets isShiftInit to true
     */

    solver->setShift(1.);
    ASSERT_TRUE(solver->getIsShiftInit());
}

TEST_F(TestShiftedPowerMethod, testGetShift) {
    /**
     * Tests if getShift returns the right shift
     */

    solver->setShift(2.5);
    ASSERT_EQ(solver->getShift(), 2.5);
}

TEST_F(TestShiftedPowerMethod, solveReturnsClosestEigenvalueToShift) {
    /**
     * Checks if the ShiftedPowerMethod returns the right eigenvalue to the shift
     */

    solver->setShift(2.5);
    auto lambda = solver->solve();
    ASSERT_NEAR(lambda, 1., 1e-15);
}

TEST_F(TestShiftedPowerMethod, noInitMakesRandomInitOfVector) {
    /**
     * Checks that calling solve() without initializing EigenVector
     * makes it being initialized through initRandomEigenVector()
     */

    solver->setShift(2.5);
    solver->setMaxIter(1);
    auto lambda = solver->solve();
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestShiftedPowerMethod, noMatrixInitReturnsException) {
    /**
     * Checks that calling solve without matrix initialization
     * returns an UninitializedSolver Exception
     */

    solver = new ShiftedPowerMethod<double>();
    solver->setShift(1.);
    ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestShiftedPowerMethod, noShiftInitReturnsException) {
    /**
     * Checks that calling solve without shift initialization
     * returns an UninitializedSolver Exception
     */

    ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestShiftedPowerMethod, nullSpaceEigenVectorReturnsException) {
    /**
     * Checks that if the eigenvector is initialized in the null space
     * of (A-shift*I) then solve() returns an std::invalid_argument exception
     */

    Eigen::VectorXd V(n);
    V.setZero();
    ASSERT_THROW(solver->solve(), std::invalid_argument);
}


TEST_F(TestShiftedPowerMethod, noConvergencePrintsToScreen) {
    /**
     * Checks if the absence of convergence prints a warning to screen
     */
    solver = new ShiftedPowerMethod<double>();
    Eigen::MatrixXd M(n,n);
    M.setRandom();
    solver->setMaxIter(1);
    solver->setMatrix(M);
    solver->setShift(1.);
    testing::internal::CaptureStdout();
    auto lambda = solver->solve();
    std::string output = testing::internal::GetCapturedStdout();
    std::string expected = "The method did not converge after ";
    expected += std::to_string(solver->getMaxIter()) + " iterations\n";
    EXPECT_STREQ(output.c_str(), expected.c_str());
}

class TestShiftedPowerMethodComplex: public ::testing::Test {
    /**
     * Test suite for the ShiftedPowerMethod class with std::complex<double> values
     * At initialization, a ShiftedPowerMethod is instantiated
     * with a diagonal Eigen::Matrix<std::complex<double>,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new ShiftedPowerMethod<std::complex<double>>();
        M.resize(n,n);
        M << std::complex<double>(0.,2.), std::complex<double>(0.,-1.),
                std::complex<double>(0.,-1.), std::complex<double>(0.,2.);
        solver->setMatrix(M);
    }

    void TearDown() override {}

    ShiftedPowerMethod<std::complex<double>>* solver;
    Eigen::MatrixXcd M;
    static constexpr int n {2};
};

TEST_F(TestShiftedPowerMethodComplex, solveWorksCorrectly) {
    /**
     * Checks if the ShiftedPowerMethod returns the right eigenvalue to the shift
     */

    solver->setShift(std::complex<double>(0.,2.5));
    auto lambda = solver->solve();
    ASSERT_NEAR(std::abs(lambda-std::complex<double>(0.,1.)), 0.0, 1e-15);
}
