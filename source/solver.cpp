#include <iostream>
#include <sstream>
#include <gmsh.h>
#include <vector>
#include <chrono>
#include <omp.h>

#include "utils.h"
#include "Mesh.h"
#include "configParser.h"

using namespace std;

namespace solver{

    // Common variables to all solver
    
    int numNodes;
    int elNumNodes;
    vector<std::size_t> elTags;
    vector<double> elFlux;
    vector<string> g_names;
    vector<double> elStiffvector;
    vector<vector<vector<double>>> Flux;

    // Perform a numerical step for all elements in mesh object
    // @param mesh Mesh object
    // @param config Configuration file
    // @param u Nodal solution vector
    // @param Flux Nodal physical Flux
    // @param beta double coefficient
    
    void numStep(Mesh &mesh,Config config,vector<vector<double>> &u,vector<vector<vector<double>>> &Flux,double beta){

        for(int eq=0; eq<4; ++eq){
            mesh.precomputeFlux(u[eq],Flux[eq],eq);

            #pragma omp parallel for schedule(static) firstprivate(elFlux,elStiffvector) num_threads(config.numThreads)
            for(int el=0; el<mesh.getElNum(); ++el){

                mesh.getElFlux(el,elFlux.data());
                mesh.getElStiffVector(el,Flux[eq],u[eq],elStiffvector.data());
                eigen::minus(elStiffvector.data(),elFlux.data(),elNumNodes);
                eigen::linEq(&mesh.elMassMatrix(el),&elStiffvector[0],&u[eq][el*elNumNodes],config.timeStep,beta,elNumNodes);
            }
        }
    }

    // Solve using forward explicit scheme
    // @param u initial nodal solution vector
    // @param mesh
    // @param config

    void forwardEuler(vector<vector<double>> &u,Mesh &mesh,Config config){

        // Memory allocation

        elNumNodes = mesh.getElNumNodes();
        numNodes = mesh.getNumNodes();
        elTags = vector<std::size_t>(&mesh.elTag(0),&mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(elNumNodes);
        elStiffvector.resize(elNumNodes);
        Flux = vector<vector<vector<double>>>(4,vector<vector<double>>(mesh.getNumNodes(),vector<double>(3)));

        // Gmsh save init

        gmsh::model::list(g_names);
        int gp_viewTag = gmsh::view::add("Pressure");
        int gv_viewTag = gmsh::view::add("Velocity");
        int grho_viewTag = gmsh::view::add("Density");
        vector<vector<double>> g_p(mesh.getElNum(),vector<double>(elNumNodes));
        vector<vector<double>> g_rho(mesh.getElNum(),vector<double>(elNumNodes));
        vector<vector<double>> g_v(mesh.getElNum(),vector<double>(3*elNumNodes));

        // Precomputation

        mesh.precomputeMassMatrix();

        // Source

        vector<vector<int>> srcIndices;
        for(int i=0;i<config.sources.size(); ++i){

            vector<int> indice;
            for(int n=0; n<mesh.getNumNodes(); ++n){

                vector<double> coord,paramCoord;
                int dim,tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n],coord,paramCoord,dim,tag);
                if(pow(coord[0]-config.sources[i][1],2)+pow(coord[1]-config.sources[i][2],2)+pow(coord[2]-config.sources[i][3],2)<pow(config.sources[i][4],2)){
                    indice.push_back(n);
                }
            }
            srcIndices.push_back(indice);
        }

        // Main Loop: time iteration
         
        auto start = chrono::system_clock::now();
        for(double t=config.timeStart,step=0,tDisplay=0; t<=config.timeEnd; t+=config.timeStep,tDisplay+=config.timeStep,++step){
            
            // Savings and prints

            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;

                // [1] Copy solution to match GMSH format

                for(int el = 0; el < mesh.getElNum(); ++el){
                    for(int n = 0; n < mesh.getElNumNodes(); ++n){

                        int elN = el*elNumNodes+n;
                        g_p[el][n] = u[0][elN];
                        g_rho[el][n] = u[0][elN]/(config.c0*config.c0);
                        g_v[el][3*n+0] = u[1][elN];
                        g_v[el][3*n+1] = u[2][elN];
                        g_v[el][3*n+2] = u[3][elN];
                    }
                }

                gmsh::view::addModelData(gp_viewTag,step,g_names[0],"ElementNodeData",elTags,g_p,t,1);
                gmsh::view::addModelData(grho_viewTag,step,g_names[0],"ElementNodeData",elTags,g_rho,t,1);
                gmsh::view::addModelData(gv_viewTag,step,g_names[0],"ElementNodeData",elTags,g_v,t,3);

                // [2] Print and compute iteration time
                
                auto end = chrono::system_clock::now();
                auto elapsed = chrono::duration_cast<chrono::seconds>(end-start);
                gmsh::logger::write("["+to_string(t)+"/"+to_string(config.timeEnd)+"s] Step number : "+ to_string((int) step)+",Elapsed time: "+to_string(elapsed.count())+"s");
            }
        
            // Updates Source
            
            for(int src=0; src<config.sources.size(); ++src){

                double amp = config.sources[src][5];
                double freq = config.sources[src][6];
                double phase = config.sources[src][7];
                double duration = config.sources[src][8];

                if(t<duration){
                    for(int n=0; n<srcIndices[src].size(); ++n){
                        u[0][srcIndices[src][n]] = amp*sin(2*M_PI*freq*t+phase);
                    }
                }
            }

            // First order Euler
            
            mesh.updateFlux(u,Flux,config.v0,config.c0,config.rho0);
            numStep(mesh,config,u,Flux,1);
        }

        // Saves to file

        gmsh::view::write(gp_viewTag,config.saveFile,true);
        gmsh::view::write(grho_viewTag,config.saveFile,true);
        gmsh::view::write(gv_viewTag,config.saveFile,true);
    }

    // Solves using explicit Runge-Kutta integration method
    // @param u initial nodal solution vector
    // @param mesh
    // @param config

    void rungeKutta(vector<vector<double>> &u,Mesh &mesh,Config config){

        // Memory allocation
         
        elNumNodes = mesh.getElNumNodes();
        numNodes = mesh.getNumNodes();
        elTags = vector<std::size_t>(&mesh.elTag(0),&mesh.elTag(0)+mesh.getElNum());
        elFlux.resize(elNumNodes);  
        elStiffvector.resize(elNumNodes);
        vector<vector<double>> k1,k2,k3,k4;
        Flux = vector<vector<vector<double>>>(4,vector<vector<double>>(mesh.getNumNodes(),vector<double>(3)));

        // Gmsh save init

        gmsh::model::list(g_names);
        int gp_viewTag = gmsh::view::add("Pressure");
        int gv_viewTag = gmsh::view::add("Velocity");
        int grho_viewTag = gmsh::view::add("Density");
        vector<vector<double>> g_p(mesh.getElNum(),vector<double>(elNumNodes));
        vector<vector<double>> g_rho(mesh.getElNum(),vector<double>(elNumNodes));
        vector<vector<double>> g_v(mesh.getElNum(),vector<double>(3*elNumNodes));

        // Precomputation
         
        mesh.precomputeMassMatrix();

        /// Source

        vector<vector<int>> srcIndices;
        for(int i=0;i<config.sources.size(); ++i){

            vector<int> indice;
            for(int n=0; n<mesh.getNumNodes(); n++){

                vector<double> coord,paramCoord;
                int dim,tag;
                gmsh::model::mesh::getNode(mesh.getElNodeTags()[n],coord,paramCoord,dim,tag);
                if(pow(coord[0]-config.sources[i][1],2)+pow(coord[1]-config.sources[i][2],2)+pow(coord[2]-config.sources[i][3],2)<pow(config.sources[i][4],2)){
                    indice.push_back(n);
                }
            }
            srcIndices.push_back(indice);
        }
        
        // Main Loop: time iterations

        auto start = chrono::system_clock::now();
        for(double t=config.timeStart,step=0,tDisplay=0; t<=config.timeEnd; t+=config.timeStep,tDisplay+=config.timeStep,++step){
            
            //  Savings and prints
            
            if(tDisplay>=config.timeRate || step==0){
                tDisplay = 0;

                // [1] Copy solution to match GMSH format

                for(int el=0; el<mesh.getElNum(); ++el){
                    for(int n=0; n<mesh.getElNumNodes(); ++n){

                        int elN = el*elNumNodes+n;
                        g_p[el][n] = u[0][elN];
                        g_rho[el][n] = u[0][elN]/(config.c0*config.c0);
                        g_v[el][3*n+0] = u[1][elN];
                        g_v[el][3*n+1] = u[2][elN];
                        g_v[el][3*n+2] = u[3][elN];
                    }
                }

                gmsh::view::addModelData(gp_viewTag,step,g_names[0],"ElementNodeData",elTags,g_p,t,1);
                gmsh::view::addModelData(grho_viewTag,step,g_names[0],"ElementNodeData",elTags,g_rho,t,1);
                gmsh::view::addModelData(gv_viewTag,step,g_names[0],"ElementNodeData",elTags,g_v,t,3);

                // [2] Print and compute iteration time

                auto end = chrono::system_clock::now();
                auto elapsed = chrono::duration_cast<chrono::seconds>(end-start);
                gmsh::logger::write("["+to_string(t)+"/"+to_string(config.timeEnd)+"s] Step number : "+ to_string((int) step)+",Elapsed time: "+to_string(elapsed.count()) +"s");
            }

            // Source

            for(int src=0; src<config.sources.size(); ++src){

                double amp = config.sources[src][5];
                double freq = config.sources[src][6];
                double phase = config.sources[src][7];
                double duration = config.sources[src][8];

                if(t<duration){
                    for(int n=0; n<srcIndices[src].size(); ++n){
                        u[0][srcIndices[src][n]] = amp*sin(2*M_PI*freq*t+phase);
                    }
                }
            }

            // Fourth order Runge-Kutta algorithm
            
            k1 = k2 = k3 = k4 = u;

            mesh.updateFlux(k1,Flux,config.v0,config.c0,config.rho0);
            numStep(mesh,config,k1,Flux,0);
            for(int eq=0; eq<u.size(); ++eq){eigen::plusTimes(k2[eq].data(),k1[eq].data(),0.5,numNodes);}

            mesh.updateFlux(k2,Flux,config.v0,config.c0,config.rho0);
            numStep(mesh,config,k2,Flux,0);
            for(int eq=0; eq<u.size(); ++eq){eigen::plusTimes(k3[eq].data(),k2[eq].data(),0.5,numNodes);}

            mesh.updateFlux(k3,Flux,config.v0,config.c0,config.rho0);
            numStep(mesh,config,k3,Flux,0);
            for(int eq=0; eq<u.size(); ++eq){eigen::plusTimes(k4[eq].data(),k3[eq].data(),1,numNodes);}

            mesh.updateFlux(k4,Flux,config.v0,config.c0,config.rho0);
            numStep(mesh,config,k4,Flux,0);

            // Concat results of R-K iterations

            for(int eq=0; eq<u.size(); ++eq){
                for(int i=0; i<numNodes; ++i){
                    u[eq][i] += (k1[eq][i]+2*k2[eq][i]+2*k3[eq][i]+k4[eq][i])/6.0;
                }
            }
        }

        // Save to file

        gmsh::view::write(gp_viewTag,config.saveFile,true);
        gmsh::view::write(grho_viewTag,config.saveFile,true);
        gmsh::view::write(gv_viewTag,config.saveFile,true);
    }
}