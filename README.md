PCSC-Eigenvalue
===
---

> C++ implementation of the eigenvalue problem from the EPFL Programming Concepts in Scientific Computing course (MATH-458)

Requirements
===
---

In order to compile the code you will need:
* C++ compiler (compatible with C++20)
* CMake (version >= 3.16.3)
* GoogleTest (version >= 1.11.0) (Downloaded via git submodules)
* Eigen (version >= 3.4.0) (Downloaded via git submodules)
* doxygen (version >= 1.9.2.) (For the documentation)

Compilation
===
---

In order to compile the code, you should:

* Clone the repository using in a terminal, then enter the directory:
```
git clone https://github.com/axeldinh/PCSC-Eigenvalue.git
cd PCSC-Eigenvalue
```

* Clone the "googletest" and "Eigen" libraries using:
```
git submodule update --init
```

* The code can now be built using CLion or the command prompt with:
```
mkdir build
cd build
cmake ..
make
```

Once this is done, two executables are created:
* *main.exe*, which can be executed using:
```
./main -option [value]
```
* *test_eigensolver.exe*, which can be executed using:
```
./test_eigensolver
```

The former allows the usage of the code via different inputs and the latter performs the tests, 
more information can be found later in the README.


Documentation
===
---

The documentation was made using doxygen 1.9.2. 
The style adopted is the JavaDoc style ([link](https://www.doxygen.nl/manual/docblocks.html) to documentation). <br>
<br>
In order to generate the documentation, first create a "doc" directory in the project directory using:
```
mkdir doc
```

Then the documentation can be generated using:
```
doxygen doc_config
```

The documentation is generated through a html file in *"doc/html/index.html"* that can be opened in a browser.

Main
===
---

Structure of the code
---

The *main.exe* file aims at being able to use the library through user inputs. The classes' inheritance diagram is as such:

![Inheritance Diagram](class_general_eigen_solver.png)

There are two mother classes, one for all eigenvalue solvers, and one dedicated for the power methods and its variants.
Doing so, we were able to avoid repeating the code for each `GeneralPowerMethod` subclass. <br>

You can notice the presence of a templated parameter `ScalarType`, this allows the solver to be used on any type as long as `Eigen::Matrix<ScalarType,-1,-1>` can be called.

The *main* executable can handle two matrices random and non-diagonal matrices which are instantiated in the code, even though it would have been better to allow the usage of .mtx files.

Each matrix correspond to a different `ScalarType`, each of them is obtained by multiplying a diagonal matrix with a random orthogonal matrix and its adjoint of the right type,
we give here the eigenvalues of the matrices:
* `double` matrix, the eigenvalues are *1*, *2*, ..., *n* where *n* is the size of the matrix.
* `std::complex<double>` matrix, the eigenvalues are *i*, *2i*, ..., *ni*  where *n* is the size of the matrix.

Usage
---

As previously stated, the *main* executable can take different kind of inputs by using:
```
./main -option [value]
```
We will enumerate the options here:
* `-option=-algorithm`, this option is mandatory, defines which algorithm to use, `value` can be:
    - `value=pw`, to use the `PowerMethod` class.
    - `value=spw`, to use the `ShiftedPowerMethod` class.
    - `value=ipw`, to use the `InversePowerMethod` class.
    - `value=sipw`, to use the `ShiftedInversePowerMethod` class.
    - `value=qr`, to use the `QRMethod` class. This can not be used if `-scalartype=complex`.
  

* `-option=-matrixsize`, this option is optional, this is an `int` option. Gives the size of the matrix on which to solve the eigenvalue problem.
If not given, the default value (100) is used.


* `-option=-scalartype`, this option is optional, defines the type of scalar to be used, value can be:
    - `value=double`, to use `double` inputs.
    - `value=complex`, to use `std::complex<double>` inputs. This is not compatible with the `QRMethod` class.


* `-option=-maxiter`, this option is optional, this is an `int` option. Defines the maximum number of iterations of the solver.
If not given, the default value (1000) is used. This option can be used on every `-algorithm`.


* `-option=-threshold`, this option is optional, this is a `double` option. Defines the threshold at which the solver should stop.
If not given, the default value (1e-15) is used. This option can be used for `-algorithm=pw,spw,ipw,sipw`, or in other words, all subclasses of `GeneralPowerMethod`.


* `-option=-shiftdouble`, this option is optional (except if `ShiftedPowerMethod` is used), this is a `double` option.
Defines the real part of the shift for the `ShiftedPowerMethod` and `ShiftedInversePowerMethod` class.


* `-option=-shiftcomplex`, this option is optional, this is a `double` option.
  Defines the complex part of the shift for the `ShiftedInversePowerMethod` class.


* `-option=-neigenvalue`, this option is optional, this is an `int` option. When using the `QRMethod`, allows the user
to choose the rank of the eigenvalue. If it is 1, the biggest eigenvalue is returned.


Notice that, if the user wants to enter a complex shift, two options must be entered (unless one of them is 0.), and the shift
will in the end defined as `shift = std::complex<double>(shiftdouble, shiftcomplex)`.

Examples
---

We give here examples of usage of the *main*:

* Use the `ShiftedInversePowerMethod` on a `std::complex<double>` matrix of size `200`, with a `threshold` of `1e-6`, a `shift` of `1. + 1.1i` 
and a maximum number of iterations of `200`:
```
./main -algorithm sipw -scalartype complex -matrixsize 200 -threshold 1e-6 -shiftdouble 1. -shiftcomplex 1.1 -maxiter 200
```
This should return:
```
Matrix Size = 200, Type = complex
Solving the eigenvalue problem using the ShiftedInversePowerMethod class with maxIter=200, threshold=1e-06, shift=(1.,1.1)
Solution: (-7.376322e-13,1.000000e+00)
```

* Use the `QRMethod` on a `double` matrix of size `50`, with a maximum number of iterations of `500`, and aiming to extract the second biggest eigenvalue:
```
./main -algorithm qr -matrixsize 50 -maxiter 500 -neigenvalue 2
```
This should return:
```
Matrix Size = 50, Type = double
Solving the eigenvalue problem using the QRMethod class with maxIter=500
Solution: 4.900000e+01 
```

Test
===
---

You can run all tests by using:
```
./test_eigensolver
```
or
```
make test
```
in the *build* folder. Doing so, 51 tests should run and pass.

The test suite is separated in 6 files, each one corresponds to a different source file.

* `TestGeneralEigenSolver.cpp`
* `TestGeneralPowerMethod.cpp`
* `TestPowerMethod.cpp`
* `TestShiftedPowerMethod.cpp`
* `TestInversePowerMethod.cpp`
* `TestShiftedInversePowerMethod.cpp`
* `TestQRMethod.cpp`

In each test suite, we test:
* If the initialization are well performed, e.g. are the default parameters well assigned
* If the exception are indeed thrown
* If the setters and getters work correctly
* If the different `solve()` methods work correctly


For `TestPowerMethod.cpp`, `TestInversePowerMethod.cpp`, `TestShiftedInversePowerMethod.cpp` and `TestQRMethod.cpp`
the `solve()` methods are tested on a non-diagonal `Eigen::Matrix<double,2,2>`, with eigenvalues `1` and `3`.

Additionally, for `TestPowerMethod.cpp`, `TestShiftedPowerMethod.cpp`, `TestInversePowerMethod.cpp` and `TestShiftedInversePowerMethod.cpp`,
the `solve()` are tested on a non-diagonal `Eigen::Matrix<std::complex<double>,2,2>`, with eigenvalues `i` and `3i`
to assert that the code works on complex values.

Library Features
===
---

In this project, we implemented the following features:
* Implementation of an abstract templated class `GeneralEigenSolver` which allows all kind of scalar inputs, as 
long as `Eigen::Matrix<ScalarType-1,-1>` can be instantiated.
More algorithms can be based on this class. Also, it only takes a pointer to the matrix, so it does not have to 
be copied (unless necessary).


* Implementation of an abstract class `GeneralPowerSolver`, which allows the implementation of variants of the power method
by simply defining the matrix to be used by the power method (e.g. for the inverse method, `InversePowerMethod` only needs to compute the inverse of the matrix).


* Implementation of the QR method, `QRMethod`, which allows the user to recover all eigenvalues of a matrix (which is limited to real eigenvalues though, and it has heavy computations).


* The solver attributes can be handled on the command line. Those options are:
  - The maximum number of iterations (all solvers)
  - The threshold (useful for power method family)
  - The shift (useful for shifted power methods)


Issues
===
---

* The user cannot use its own matrices through the program. 
A matrix reader class should have been created to read .mtx files (not done due to time constraints).


* The program does solve eigenvalue problems using sparse matrices, can be a problem for big matrices.
Given the current code, implementing such a feature would ask for a huge refactoring (to our opinion).


* Code could have been refactored more, by grouping shifted (inverse) power methods under a common superclass.


