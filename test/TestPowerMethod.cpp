//
// Created by axeld on 01/12/2021.
//

/**
 * We implement here the tests for the PowerMethod class.
 */

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include "../src/PowerMethod.h"

TEST(TestPowerMethod, diagonalMatrix) {
    /**
     * Tests if give A = c*I, then the solve method return c.
     */
    Eigen::Matrix<double, 5, 5> M;
    M.setZero();
    M.diagonal() = Eigen::Vector<double,5>::Constant(3.);
    auto solver = new PowerMethod<double>();
    solver->setMatrix(M);
    auto lambda = solver->solve();
    ASSERT_NEAR(lambda, 3., 1e-15);
}