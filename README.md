## Introduction

This code is an implementation of a discontinuous Galerkin finite element solver for acoustic waves on complex and unstructured 3D meshes. The code relies on the open-source Gmsh library for building the initial domain, and on the Eigen library for performing algebraic operations.

## Installation

First, make sure that the Eigen library is installed and that your GCC supports OpenMP. Then move to the source folder and compile the project with CMake. An example of compiler command is given in compile.sh.
```sh
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```
The executable requires a Mesh and a Config files containing the simulation parameters such as the initial conditions, some examples are given in the example folder. The output can then be opened in your Gmsh application to visualize the solution.
```sh
./Acoustic ${path-to-mesh.msh} ${path-to-config.conf}
```