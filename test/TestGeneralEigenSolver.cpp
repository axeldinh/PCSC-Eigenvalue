/**
 * We test the GeneralEigenSolver class
 * As it is a pure virtual class, we test through the usage of a PowerMethod class
 */

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"

class TestGeneralEigenSolver: public ::testing::Test {
    /**
     * Test suite for the GeneralEigenSolver class.
     * The test are done using the PowerMethod subclass.
     */
protected:
    void SetUp() override {
        solver = new PowerMethod<double>();
    }

    void TearDown() override {}
    PowerMethod<double>* solver;
};

TEST_F(TestGeneralEigenSolver, goodInit) {
    /**
     * Checks if the solver is well initialized
     */
     ASSERT_EQ(solver->getMaxIter(), 1000);
     ASSERT_FALSE(solver->getIsMatrixInit());
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
     * Checks that setting the maximum number of iterations to a
     * negative number throws an std::invalid_argument
     */
    ASSERT_THROW(solver->setMaxIter(-1), std::invalid_argument);
}

TEST_F(TestGeneralEigenSolver, changeMatrixUpdatesIsInitMatrix) {
    /**
     * Tests if changing the matrix updates the boolean isMatrixInit
     */
     Eigen::MatrixXd M(3,3);
     M.setZero();
     solver->setMatrix(M);
     ASSERT_TRUE(solver->getIsMatrixInit());
}

TEST_F(TestGeneralEigenSolver, nonSquareMatrixThrowsException) {
    /**
     * Checks that the user cannot instantiate a solver with a non-square matrix
     */
    Eigen::MatrixXd M(2,3);
    M.setZero();
    ASSERT_THROW(solver->setMatrix(M), std::invalid_argument);
}
