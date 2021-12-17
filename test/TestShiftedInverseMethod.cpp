#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/ShiftedInversePowerMethod.h"

class TestShiftedInversePowerMethod: public ::testing::Test {
    /**
     * Test suite for the ShiftedInversePowerMethod class
     * At initialization, a ShiftedInversePowerMethod is instantiated
     * with a diagonal Eigen::Matrix<double,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new ShiftedInversePowerMethod<double>();
        M.resize(n,n);
        M << 2, -1,
            -1, 2;
        solver->setMatrix(M);
    }

    void TearDown() override {}

    ShiftedInversePowerMethod<double>* solver;
    Eigen::MatrixXd M;
    static constexpr int n {2};
};

TEST_F(TestShiftedInversePowerMethod, setShiftMakesIsShiftInitToTrue) {
    /**
     * Checks if initializing the shift sets isShiftInit to true
     */

    solver->setShift(1.);
    ASSERT_TRUE(solver->getIsShiftInit());
}

TEST_F(TestShiftedInversePowerMethod, testGetShift) {
    /**
     * Tests if getShift returns the right shift
     */

    solver->setShift(2.5);
    ASSERT_EQ(solver->getShift(), 2.5);
}

TEST_F(TestShiftedInversePowerMethod, solveReturnsClosestEigenvalueToShift) {
    /**
     * Checks if the ShiftedInversePowerMethod returns the closest eigenvalue to the shift
     */

    solver->setShift(2.5);
    auto lambda = solver->solve();
    ASSERT_NEAR(lambda, 3., 1e-15);
}

TEST_F(TestShiftedInversePowerMethod, noInitMakesRandomInitOfVector) {
    /**
     * Checks that calling solve() without initializing EigenVector
     * makes it being initialized through initRandomEigenVector()
     */

    solver->setShift(2.5);
    solver->setMaxIter(1);
    auto lambda = solver->solve();
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestShiftedInversePowerMethod, noMatrixInitReturnsException) {
    /**
     * Checks that calling solve without matrix initialization
     * returns an UninitializedSolver Exception
     */

    solver = new ShiftedInversePowerMethod<double>();
    solver->setShift(1.);
    ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestShiftedInversePowerMethod, noShiftInitReturnsException) {
    /**
     * Checks that calling solve without shift initialization
     * returns an UninitializedSolver Exception
     */

    ASSERT_THROW(solver->solve(), UninitializedSolver);
}

TEST_F(TestShiftedInversePowerMethod, nullSpaceEigenVectorReturnsException) {
    /**
     * Checks that if the eigenvector is initialized in the null space
     * of (A-shift*I)^-1 then solve() returns an std::invalid_argument exception
     */

    Eigen::VectorXd V(n);
    V.setZero();
    ASSERT_THROW(solver->solve(), std::invalid_argument);
}

TEST_F(TestShiftedInversePowerMethod, returnShiftIfItIsEigenValue) {
    /**
     * Checks that the shift is directly returned if it is an eigenvalue
     */
     solver->setShift(3.);
     solver->setMaxIter(1);
     ASSERT_EQ(solver->solve(), 3.);
}

TEST_F(TestShiftedInversePowerMethod, noConvergencePrintsToScreen) {
    /**
     * Checks if the absence of convergence prints a warning to screen
     */
    solver = new ShiftedInversePowerMethod<double>;
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

class TestShiftedInversePowerMethodComplex: public ::testing::Test {
    /**
     * Test suite for the ShiftedInversePowerMethod class with std::complex<double> values
     * At initialization, a ShiftedInversePowerMethod is instantiated
     * with a diagonal Eigen::Matrix<std::complex<double>,2,2>
     * (No vector initialization)
     */
protected:
    void SetUp() override {
        solver = new ShiftedInversePowerMethod<std::complex<double>>();
        M.resize(n,n);
        M << std::complex<double>(0.,2.), std::complex<double>(0.,-1.),
                std::complex<double>(0.,-1.), std::complex<double>(0.,2.);
        solver->setMatrix(M);
    }

    void TearDown() override {}

    ShiftedInversePowerMethod<std::complex<double>>* solver;
    Eigen::MatrixXcd M;
    static constexpr int n {2};
};

TEST_F(TestShiftedInversePowerMethodComplex, solveReturnsClosestEigenvalueToShift) {
    /**
     * Checks if the ShiftedInversePowerMethod returns the closest eigenvalue to the shift
     */

    solver->setShift(std::complex<double>(0.,2.5));
    auto lambda = solver->solve();
    ASSERT_NEAR(std::abs(lambda-std::complex<double>(0.,3.)), 0.0, 1e-15);
}

TEST_F(TestShiftedInversePowerMethodComplex, returnShiftIfItIsEigenValue) {
    /**
     * Checks that the shift is directly returned if it is an eigenvalue
     */
    solver->setShift(std::complex<double>(0.,3.));
    solver->setMaxIter(1);
    auto lambda = solver->solve();
    ASSERT_NEAR(std::abs(lambda-std::complex<double>(0.,3.)), 0.0, 1e-15);
}
