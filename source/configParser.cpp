#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdio>
#include <gmsh.h>
#include <map>

#include "utils.h"
#include "configParser.h"

using namespace std;

namespace config{

    vector<string> split(string str,char delimiter){
        vector<string> internal;
        stringstream ss(str);
        string tok;

        while(getline(ss,tok,delimiter)){internal.push_back(tok);}
        return internal;
    }

    Config parseConfig(string name){

        Config config;
        ifstream cFile(name);

        if(cFile.is_open()){
            
            string line;
            map<string,string> configMap;

            while(getline(cFile,line)){

                auto is_space = [](unsigned char x){return isspace(x);};
                line.erase(remove_if(line.begin(),line.end(),is_space),line.end());
                if(line[0] == '#' || line.empty()){continue;}

                auto delimiterPos = line.find("=");
                auto name = line.substr(0,delimiterPos);
                auto value = line.substr(delimiterPos+1);
                configMap[name] = value;
            }

            config.numThreads = stoi(configMap["numThreads"]);
            config.numThreads = config.numThreads==1? 0 : config.numThreads;
            config.timeStart = stod(configMap["timeStart"]);
            config.timeStep = stod(configMap["timeStep"]);
            config.timeRate = stod(configMap["timeRate"]);
            config.timeIntMethod = configMap["timeIntMethod"];
            config.timeEnd = stod(configMap["timeEnd"]);
            config.elementType = configMap["elementType"];
            config.v0[0] = stod(configMap["v0_x"]);
            config.v0[1] = stod(configMap["v0_y"]);
            config.v0[2] = stod(configMap["v0_z"]);
            config.rho0 = stod(configMap["rho0"]);
            config.saveFile = configMap["saveFile"];
            config.c0 = stod(configMap["c0"]);

            for(map<string,string>::iterator iter=configMap.begin(); iter!=configMap.end(); ++iter){

                string key = iter->first;

                if(key.find("source") == 0){

                    vector<string> sep = split(iter->second,',');
                    double x = stod(sep[1]);
                    double y = stod(sep[2]);
                    double z = stod(sep[3]);
                    double size = stod(sep[4]);
                    double amp = stod(sep[5]);
                    double freq = stod(sep[6]);
                    double phase = stod(sep[7]);
                    double duration = stod(sep[8]);
                    double pole;

                    if(sep[0] == "dipole"){

                        pole = 1;
                        vector<double> source1 ={pole,x-size,y,z,size/2.,amp,freq,phase,duration};
                        vector<double> source2 ={pole,x+size,y,z,size/2.,amp,freq,phase+M_PI,duration};
                        config.sources.push_back(source1);
                        config.sources.push_back(source2);
                    }

                    else if(sep[0] == "quadrupole"){

                        pole = 2;
                        vector<double> source1 ={pole,x-size,y,z,size/2.,amp,freq,phase,duration};
                        vector<double> source2 ={pole,x+size,y,z,size/2.,amp,freq,phase,duration};
                        vector<double> source3 ={pole,x,y-size,z,size/2.,amp,freq,phase+M_PI,duration};
                        vector<double> source4 ={pole,x,y+size,z,size/2.,amp,freq,phase+M_PI,duration};
                        config.sources.push_back(source1);
                        config.sources.push_back(source2);
                        config.sources.push_back(source3);
                        config.sources.push_back(source4);
                    }

                    else{

                        pole = 0;
                        vector<double> source ={pole,x,y,z,size,amp,freq,phase,duration};
                        config.sources.push_back(source);
                    }
                }

                else if(key.find("initialCondtition") == 0){

                    vector<string> sep = split(iter->second,',');
                    double x = stod(sep[1]);
                    double y = stod(sep[2]);
                    double z = stod(sep[3]);
                    double size = stod(sep[4]);
                    double amp = stod(sep[5]);
                    vector<double> init1 ={0,x,y,z,size,amp};
                    config.initConditions.push_back(init1);
                }
            }

            string physName;
            remove(&config.saveFile[0]);
            gmsh::vectorpair m_physicalDimTags;
            int bcDim = gmsh::model::getDimension()-1;
            gmsh::model::getPhysicalGroups(m_physicalDimTags,bcDim);

            for(int p=0; p<m_physicalDimTags.size(); ++p){

                gmsh::model::getPhysicalName(m_physicalDimTags[p].first,m_physicalDimTags[p].second,physName);
                if(configMap[physName].find("Absorbing") == 0){config.physBCs[m_physicalDimTags[p].second] = make_pair("Absorbing",0);}
                else if(configMap[physName].find("Reflecting") == 0){config.physBCs[m_physicalDimTags[p].second] = make_pair("Reflecting",0);}
                else{gmsh::logger::write("Not specified or supported boundary conditions.");}
            }

            gmsh::logger::write("==================================================");
            gmsh::logger::write("Simulation parameters : ");
            gmsh::logger::write("Time step : "+to_string(config.timeStep));
            gmsh::logger::write("Final time : "+to_string(config.timeEnd));
            gmsh::logger::write("Mean flow velocity : ("+to_string(config.v0[0])+","+ to_string(config.v0[1])+","+to_string(config.v0[2])+")");
            gmsh::logger::write("Mean density : "+to_string(config.rho0));
            gmsh::logger::write("Speed of sound : "+to_string(config.c0));
            gmsh::logger::write("Solver : "+config.timeIntMethod);
        }

        else{throw;}
        return config;
    }
}