
#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"

class TestGeneralPowerMethod: public ::testing::Test {
protected:
    void SetUp() override {
        solver = new PowerMethod<double>();
    }

    void TearDown() override {
    }
    PowerMethod<double>* solver;
};

TEST_F(TestGeneralPowerMethod, setEigenVectorActuallyUpdates) {
    /**
     * Checks if calling setEigenVector actually updates mEigenVector
     */

    Eigen::VectorXd V(1);
    V(0) = 2.;
    solver->setEigenVector(V);
    ASSERT_EQ(solver->getEigenVector()(0), 2.);
}

TEST_F(TestGeneralPowerMethod, goodInit) {
    /**
     * Checks if the solver is well initialized
     */
    ASSERT_EQ(solver->getThreshold(), 1e-15);
    ASSERT_FALSE(solver->getIsVectorInit());
}

TEST_F(TestGeneralPowerMethod, changePositiveThreshold) {
    /**
     * Check that the threshold actually updates
     */
    solver->setThreshold(1);
    ASSERT_EQ(solver->getThreshold(), 1);
}

TEST_F(TestGeneralPowerMethod, changeNegativeThreshold) {
    /**
     * Check that changing to a <0 threshold throws an std::invalid_argument
     */
    ASSERT_THROW(solver
    ->setThreshold(-1), std::invalid_argument);
}

TEST_F(TestGeneralPowerMethod, changeSmallThreshold) {
    /**
     * Check that a warning displays if the desired threshold is small
     */
    testing::internal::CaptureStderr();
    solver->setThreshold(1e-21);
    std::string output = testing::internal::GetCapturedStderr();
    std::string expected = "WARNING: Threshold < 1e-20, the computation might take a long time.\n";
    EXPECT_STREQ(output.c_str(), expected.c_str());
}

TEST_F(TestGeneralPowerMethod, changeVectorUpdatesIsInitVector) {
    /**
     * Tests if changing the vector updates the boolean isVectorInit
     */
    Eigen::VectorXd V(3);
    V.setZero();
    solver->setEigenVector(V);
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestGeneralPowerMethod, randomInitVectorSetBoolToTrue) {
    /**
     * Checks if Randomly initialize the EigenVector
     * sets isVectorInit to true
     */
    Eigen::MatrixXd M(3,3);
    M.setZero();
    solver->setMatrix(M);
    solver->initRandomEigenVector();
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestGeneralPowerMethod, cannotRandomInitVectorWithoutMatrix) {
    /**
     * Checks that trying to randomly initialize the EigenVector without matrix
     * throws an UninitializedSolver exception
     */
    ASSERT_THROW(solver->initRandomEigenVector(), UninitializedSolver);
}