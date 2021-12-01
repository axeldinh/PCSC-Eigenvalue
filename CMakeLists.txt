cmake_minimum_required(VERSION 3.20)
project(eigenvalue_project)

set(CMAKE_CXX_STANDARD 20)

include_directories(eigen)
include_directories(googletest)

add_subdirectory(src)

add_executable(eigenvalue_project main.cpp src/GeneralEigenSolver.cpp)