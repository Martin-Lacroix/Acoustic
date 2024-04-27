#include <string>
#include <map>

#ifndef DGALERKIN_CONFIG_H
#define DGALERKIN_CONFIG_H
using namespace std;

struct Config{

    // Initial, final time and time step(t>0)

    double timeEnd =  1;
    double timeStart = 0;
    double timeStep = 0.1;
    double timeRate = 0.1;

    // Element Type

    string elementType = "Lagrange";

    // Time integration method

    string timeIntMethod = "Euler";

    // Boundary condition
    // key: physical group Tag
    // value: tuple<BCType,BCValue>

    map<int,pair<string,double>> physBCs;

    // Mean flow parameters, number of threads

    vector<double> v0 = {0,0,0};
    int numThreads = 1;
    double rho0 = 1;
    double c0 = 1;

    // Sources, initial conditions, save file

    vector<vector<double>> sources;
    vector<vector<double>> initConditions;
    string saveFile = "results.msh";
};

namespace config{Config parseConfig(string name);}
#endif