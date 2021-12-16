cmake_minimum_required(VERSION 3.0)
project(eigenvalue_project)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(googletest)
add_subdirectory(eigen)

add_library(eigensolver src/PowerMethod.cpp src/GeneralEigenSolver.cpp src/ShiftedInversePowerMethod.cpp
        src/InversePowerMethod.cpp src/QRMethod.cpp)
target_link_libraries(eigensolver Eigen3::Eigen)

add_executable(main src/main.cpp)
target_link_libraries(main eigensolver Eigen3::Eigen)

add_executable(test_eigensolver test/TestPowerMethod.cpp test/TestInverseMethod.cpp test/TestShiftedInverseMethod.cpp
        test/TestGeneralPowerMethod.cpp test/TestGeneralEigenSolver.cpp)
target_link_libraries(test_eigensolver gtest_main gtest pthread eigensolver Eigen3::Eigen)

add_custom_target(test ./test_eigensolver DEPENDS test_eigensolver)
