#include <iostream>
#include <errno.h>
#include <gmsh.h>
#include <cstdio>
#include <omp.h>

#include "Mesh.h"
#include "solver.h"
#include "configParser.h"

using namespace std;

int main(int argc,char**argv){

    if(argc!=3){return E2BIG;}
    string msh_name = argv[1];
    string config_name = argv[2];

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal",1);
    gmsh::open(msh_name);

    Config config = config::parseConfig(config_name);
    gmsh::logger::write("Config loaded : "+config_name);
    Mesh mesh(msh_name,config);

    // Initializes the solution

    vector<vector<double>> u(4,vector<double>(mesh.getNumNodes(),0));

    for(int i=0; i<config.initConditions.size(); ++i){

	double x = config.initConditions[i][1];
	double y = config.initConditions[i][2];
	double z = config.initConditions[i][3];
	double size = config.initConditions[i][4];
	double amp = config.initConditions[i][5];

        for(int n=0; n<mesh.getNumNodes(); n++){

            vector<double> coord,paramCoord;
            int dim,tag;
            gmsh::model::mesh::getNode(mesh.getElNodeTags()[n],coord,paramCoord,dim,tag);
	        u[0][n] += amp*exp(-((coord[0]-x)*(coord[0]-x)+(coord[1]-y)*(coord[1]-y)+(coord[2]-z)*(coord[2]-z))/size);
        }
    }
   
    // Starts solver

    if(config.timeIntMethod == "Euler"){solver::forwardEuler(u,mesh,config);}
    else if(config.timeIntMethod == "Runge-Kutta"){solver::rungeKutta(u,mesh,config);}

    gmsh::finalize();
    return EXIT_SUCCESS;
}
