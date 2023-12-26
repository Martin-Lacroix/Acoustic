## Introduction

Algorithm developed for the **Multiphysics Integrated Computational Project**. The project consists in the implementation of a discontinuous Galerkin finite element solver for acoustic waves on 3-dimensional complex and unstructured meshes using the Gmsh library.

## Installation

First, make sure that the Eigen library is installed and that your GCC supports OpenMP. Then move to the source folder and compile the project with CMake. An example of compiler command is given in compile.sh.
```css
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
```
The executable requires a Mesh and a Config files containing the simulation parameters such as the initial conditions, some examples are given in the example folder. The output can then be opened in your Gmsh application to visualize the solution.
```css
.\Acoustic ${path-to-mesh.msh} ${path-to-config.conf}
```

## Authors

* Martin Lacroix
* Pierre-Olivier Vanberg
* Tom Servais
