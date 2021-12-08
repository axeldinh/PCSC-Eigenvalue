//
// Created by axeld on 02/12/2021.
//

/**
 * We test the GeneralEigenSolver class
 * As it is a pure virtual class, we test through the usage of a PowerMethod class
 */

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"
#include "../src/Exceptions/UninitializedSolver.h"

class TestGeneralEigenSolver: public ::testing::Test {
protected:
    void SetUp() override {
        solver = new PowerMethod<double>;
    }

    void TearDown() override {
        delete solver;
    }
    PowerMethod<double>* solver;
};

TEST_F(TestGeneralEigenSolver, goodInit) {
    /**
     * Checks if the solver is well initialized
     */
     ASSERT_EQ(solver->getMaxIter(), 1000);
     ASSERT_EQ(solver->getThreshold(), 1e-15);
     ASSERT_FALSE(solver->getIsVectorInit());
     ASSERT_FALSE(solver->getIsMatrixInit());
}

TEST_F(TestGeneralEigenSolver, changePositiveThreshold) {
    /**
     * Check that the threshold actually updates
     */
    solver->setThreshold(1);
    ASSERT_EQ(solver->getThreshold(), 1);
}

TEST_F(TestGeneralEigenSolver, changeNegativeThreshold) {
    /**
     * Check that changing to a <0 threshold throws an std::invalid_argument
     */
    ASSERT_THROW(solver->setThreshold(-1), std::invalid_argument);
}

TEST_F(TestGeneralEigenSolver, changeSmallThreshold) {
    /**
     * Check that a warning displays if the desired threshold is small
     */
     testing::internal::CaptureStderr();
     solver->setThreshold(1e-21);
     std::string output = testing::internal::GetCapturedStderr();
     std::string expected = "WARNING: Threshold < 1e-20, the computation might take a long time.\n";
     EXPECT_STREQ(output.c_str(), expected.c_str());
}


TEST_F(TestGeneralEigenSolver, changePositiveMaxIter) {
    /**
     * Check if the maximum number of iterations
     * actually updates
     */
    solver->setMaxIter(2);
    ASSERT_EQ(solver->getMaxIter(), 2);
}

TEST_F(TestGeneralEigenSolver, changeNegativeMaxIter) {
    /**
     * Checks thats setting the maximum numer of iterations to a
     * negative number throws an std::invalid_argument
     */
    ASSERT_THROW(solver->setMaxIter(-1), std::invalid_argument);
}

TEST_F(TestGeneralEigenSolver, changeMatrixUpdatesIsInitMatrix) {
    /**
     * Tests if changing the matrix updates the boolean isMatrixInit
     */
     Eigen::Matrix<double,3,3> M;
     M.setZero();
     solver->setMatrix(M);
     ASSERT_TRUE(solver->getIsMatrixInit());
}

TEST_F(TestGeneralEigenSolver, changeVectorUpdatesIsInitMatrix) {
    /**
     * Tests if changing the matrix updates the boolean isMatrixInit
     */
    Eigen::Vector<double,3> V;
    V.setZero();
    solver->setEigenVector(V);
    ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestGeneralEigenSolver, randomInitVectorSetBoolToTrue) {
    /**
     * Checks if Randomly initialize the EigenVector
     * sets isVectorInit to true
     */
     Eigen::Matrix<double,3,3> M;
     M.setZero();
     solver->setMatrix(M);
     solver->initRandomEigenVector();
     ASSERT_TRUE(solver->getIsVectorInit());
}

TEST_F(TestGeneralEigenSolver, cannotRandomInitVectorWithoutMatrix) {
    /**
     * Checks that trying to randomly initialize the EigenVector without matrix
     * throws an UninitializedSolver exception
     */
     ASSERT_THROW(solver->initRandomEigenVector(), UninitializedSolver);
}

