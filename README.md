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

* The code can now be built using CLion or the command line with:
```
mkdir build
cd build
cmake ..
make
```

Once this is done, two executables are created:
* *main.exe*, which can executed using:
```
./main -option [value]
```
* *test_eigensolver.exe*, which can executed using:
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

![Inheritance Diagram](classes_diagram.png)

There are two mother classes, one for all eigenvalue solvers, and one dedicated for the power methods and its variants.
Doing so, we were able to avoid repeating the code for each `GeneralPowerMethod` subclass. <br>

You can notice the presence of a templated parameter `ScalarType`, this allows the solver to be used on any type as long as `Eigen::Matrix<ScalarType,-1,-1>` can be called.

The *main* executable can handle two matrices which instantiated in the code, even though it would have been better to allow the usage of .mtx files.

Each matrix correspond to a different `ScalarType`, each of them is obtained by multiplying a diagonal matrix with an orthogonal matrix of the right type,
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
If not given, the default value (1e-15) is used. This option can be used for `-algorithm=pw,ipw,sipw`, mainly all subclasses of `GeneralPowerMethod`.


* `-option=-neigenvalue`, this option is optional, this is an `int` option. When using the `QRMethod`, allows the user
to choose the rank of the eigenvalue. If it is 1, the biggest eigenvalue is returned.

Examples
---

We give here examples of usage of the *main*:

* 

TODOs
===
---

Check how to handle libraries (make cmake download them, let the user clone the gits, etc...)


---

Examples of markdown writing:

Write code:
`int main (int argc, char* argv[]){}`

Make a line

---

Link: [power method](https://en.wikipedia.org/wiki/Power_iteration)

Insert Image ![alt test](imag.jpg)