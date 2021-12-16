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

Compilation
===
---

In order to compile the code, you should:

* Clone the repository using in a terminal, then enter the directory:
> git clone https://github.com/axeldinh/PCSC-Eigenvalue.git <br>
> cd PCSC-Eigenvalue

* Clone the "googletest" and "Eigen" libraries using:
> git submodule update --init

* The code can now be built using CLion or the command line with:
> mkdir build <br>
> cd build <br>
> cmake .. <br>
> make <br>

Once this is done, two executables are created:
* *main.exe*, which can executed using:
> ./main -option [value]
* *test_eigensolver.exe*, which can executed using:
> ./test_eigensolver

The former allows the usage of the code via different inputs and the latter performs the tests, 
more information can be found later in the README.


Documentation
===
---

The documentation was made using doxygen 1.9.2. The style adopted is the JavaDoc style ([link](https://www.doxygen.nl/manual/docblocks.html) to documentation).

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