cmake_minimum_required(VERSION 3.20)
project(eigenvalue_project)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(eigenvalue_solver src/PowerMethod.cpp src/GeneralEigenSolver.cpp src/ShiftedInversePowerMethod.cpp)
target_link_libraries(eigenvalue_solver Eigen3::Eigen)

add_executable(main src/main.cpp)
target_link_libraries(main eigenvalue_solver Eigen3::Eigen)

add_subdirectory(googletest)

add_executable(test_eigenvalue_solver test/TestPowerMethod.cpp test/TestGeneralEigenSolver.cpp)
target_link_libraries(test_eigenvalue_solver gtest_main gtest pthread eigenvalue_solver Eigen3::Eigen)

add_custom_target(test ./test_eigenvalue_solver DEPENDS test_eigenvalue_solver)
