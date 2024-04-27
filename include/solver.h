#include "configParser.h"
#include "Mesh.h"

#ifndef DGALERKIN_SOLVER_H
#define DGALERKIN_SOLVER_H
using namespace std;

namespace solver{
    
    // Solve using forward explicit scheme
    // @param u initial nodal solution vector
    // @param mesh
    // @param config
     
    void forwardEuler(vector<vector<double>> &u,Mesh &mesh,Config config);

    // Solve using explicit Runge-Kutta integration method
    // @param u initial nodal solution vector
    // @param mesh
    // @param config

    void rungeKutta(vector<vector<double>> &u,Mesh &mesh,Config config);
}

#endif