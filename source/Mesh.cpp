#include <string>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <gmsh.h>
#include <omp.h>

#include "utils.h"
#include "Mesh.h"
#include "configParser.h"

using namespace std;

Mesh::Mesh(string name,Config config) : name(name),config(config){

    // Elements

    int numPrimaryNodes = 0;

    m_elDim = gmsh::model::getDimension();
    gmsh::model::mesh::getElementTypes(m_elType,m_elDim);
    gmsh::model::mesh::getElementProperties(m_elType[0],m_elName,m_elDim,m_elOrder,m_elNumNodes,m_elParamCoord,numPrimaryNodes);
    gmsh::model::mesh::getElementsByType(m_elType[0],m_elTags,m_elNodeTags);
    m_elNum = (int) m_elTags.size();
    m_elIntType = "Gauss"+to_string(2*m_elOrder);

    int numOrientations;
    int numComponents;

    gmsh::model::mesh::getIntegrationPoints(m_elType[0], m_elIntType, m_elParamCoord, m_elWeight);
    gmsh::model::mesh::getBasisFunctions(m_elType[0],m_elParamCoord,config.elementType,numComponents,m_elBasisFcts,numOrientations);
    gmsh::model::mesh::getBasisFunctions(m_elType[0],m_elParamCoord,"Grad"+config.elementType,numComponents,m_elUGradBasisFcts,numOrientations);
    gmsh::model::mesh::getJacobians(m_elType[0],m_elParamCoord,m_elJacobians,m_elJacobianDets, m_elIntPtCoords);

    m_elNumIntPts = (int) m_elJacobianDets.size() / m_elNum;

    // Gmsh provides the derivative of the shape functions along the parametric directions.
    // We therefore compute their derivative along the physical directions thanks to composed derivative.
    // The system can be expressed as J^T*df/dx = df/du

    // |dx/du dx/dv dx/dw|^T   |df/dx|   |df/du|
    // |dy/du dy/dv dy/dw|  *  |df/dy| = |df/dv|
    // |dz/du dz/dv dz/dw|     |df/dz|   |df/dw|

    // (x,y,z) are the physical coordinates
    // (u,v,w) are the parametric coordinates

    vector<double> jacobian(m_elDim*m_elDim);
    m_elGradBasisFcts.resize(m_elNum*m_elNumNodes*m_elNumIntPts*3);

    for (int el=0; el<m_elNum; ++el){
        for (int g=0; g<m_elNumIntPts; ++g){
            for (int f=0; f<m_elNumNodes; ++f){

                // The copy are not required but are enforced to ensure that the inputs (jacobian,grad) remains unchanged

                for(int i=0; i<m_elDim; ++i){
                    for(int j=0; j<m_elDim; ++j){
                        jacobian[i*m_elDim+j] = elJacobian(el,g,i,j);
                    }
                }
                copy(&elUGradBasisFct(g,f),&elUGradBasisFct(g,f)+m_elDim,&elGradBasisFct(el,g,f));
                eigen::solve(jacobian.data(),&elGradBasisFct(el,g,f),m_elDim);
            }
        }
    }

    assert(m_elType.size() == 1);
    assert(m_elNodeTags.size() == m_elNum*m_elNumNodes);
    assert(m_elJacobianDets.size() == m_elNum*m_elNumIntPts);
    assert(m_elBasisFcts.size() == m_elNumNodes*m_elNumIntPts);
    assert(m_elGradBasisFcts.size() == m_elNum*m_elNumIntPts*m_elNumNodes*3);

    gmsh::logger::write("==================================================");
    gmsh::logger::write("Number of Elements : "+to_string(m_elNum));
    gmsh::logger::write("Element dimension : "+to_string(m_elDim));
    gmsh::logger::write("Element Type : "+m_elName);
    gmsh::logger::write("Element Order : "+to_string(m_elOrder));
    gmsh::logger::write("Element Nbr Nodes : "+to_string(m_elNumNodes));
    gmsh::logger::write("Integration type : "+m_elIntType);
    gmsh::logger::write("Integration Nbr points : "+to_string(m_elNumIntPts));

    // Faces

    m_fDim = m_elDim-1;
    m_fName = m_fDim==0? "point": m_fDim==1? "line": m_fDim==2? "triangle": "None";
    m_fNumNodes = m_fDim==0? 1: m_fDim==1? 1+m_elOrder: m_fDim==2? (m_elOrder+1)*(m_elOrder+2)/2: 0;
    m_fType = gmsh::model::mesh::getElementType(m_fName,m_elOrder);

    // [1] Get Faces for all elements

    if(m_fDim<2){gmsh::model::mesh::getElementEdgeNodes(m_elType[0],m_elFNodeTags,-1);}
    else{gmsh::model::mesh::getElementFaceNodes(m_elType[0],3,m_elFNodeTags,-1);}
    m_fNumPerEl = m_elFNodeTags.size()/(m_elNum*m_fNumNodes);

    // [2] Remove the faces counted two times

    getUniqueFaceNodeTags();

    // [3] Create a single entity containing all the unique faces

    m_fEntity = gmsh::model::addDiscreteEntity(m_fDim);
    gmsh::model::mesh::addElementsByType(m_fEntity, m_fType, {}, m_fNodeTags);
    m_fNodeTags.clear();
    gmsh::model::mesh::getElementsByType(m_fType,m_fTags,m_fNodeTags,m_fEntity);
    m_fNum = m_fTags.size();


    // The same integration type and order is applied to the surface and to the volume integrals

    m_fIntType = m_elIntType;

    gmsh::model::mesh::getIntegrationPoints(m_fType, m_fIntType, m_fIntParamCoords, m_fWeight);
    gmsh::model::mesh::getBasisFunctions(m_fType,m_fIntParamCoords,config.elementType,*new int,m_fBasisFcts,numOrientations);
    gmsh::model::mesh::getBasisFunctions(m_fType,m_fIntParamCoords,"Grad"+config.elementType, *new int, m_fUGradBasisFcts, numOrientations);
    gmsh::model::mesh::getJacobians(m_fType,m_fIntParamCoords,m_fJacobians,m_fJacobianDets,m_fIntPtCoords,m_fEntity);

    m_fNumIntPts = (int)m_fJacobianDets.size() / m_fNum;
    
    // See element part for explanations

    m_fGradBasisFcts.resize(m_fNum*m_fNumNodes*m_fNumIntPts*3);

    for (int f=0; f<m_fNum; ++f){
        for (int g=0; g<m_fNumIntPts; ++g){
            for (int n=0; n<m_fNumNodes; ++n){
                for(int i=0; i<m_elDim; ++i){
                    for(int j=0; j<m_elDim; ++j){
                        jacobian[i*m_elDim+j] = fJacobian(f,g,i,j);
                    }
                }
                copy(&fUGradBasisFct(g,n),&fUGradBasisFct(g,n)+m_elDim,&fGradBasisFct(f,g,n));
                eigen::solve(jacobian.data(),&fGradBasisFct(f,g,n),m_elDim);
            }
        }
    }

    // Define a normal associated to each surface.

    vector<double> normal(m_Dim);

    for(int f=0; f<m_fNum; ++f){
        for(int g=0; g<m_fNumIntPts; ++g){

            switch (m_fDim){

                case 0:{
                    normal = {1,0,0};
                    break;
                }

                case 1:{
                    vector<double> normalPlane = {0,0,-1};
                    eigen::cross(&fGradBasisFct(f,g,0),normalPlane.data(),normal.data());
                    if(eigen::dot(&fGradBasisFct(f,g),&fGradBasisFct(f,0),m_Dim)<0){
                        for (int x=0; x<m_Dim; ++x){normal[x] = -normal[x];}
                    }
                    break;
                }

                case 2:{
                    eigen::cross(&fGradBasisFct(f,g,0),&fGradBasisFct(f,g,1),normal.data());
                    if(g!=0 && eigen::dot(&fNormal(f,0),normal.data(),m_Dim)<0){
                        for (int x=0; x<m_Dim; ++x){normal[x] = -normal[x];}
                    }
                    break;
                }
            }
            eigen::normalize(normal.data(),m_Dim);
            m_fNormals.insert(m_fNormals.end(),normal.begin(),normal.end());
        }
    }

    if(m_elDim==3 && m_elOrder!=1){fc = -1;}
    assert(m_elFNodeTags.size() == m_elNum*m_fNumPerEl*m_fNumNodes);
    assert(m_fJacobianDets.size() == m_fNum*m_fNumIntPts);
    assert(m_fBasisFcts.size() == m_fNumNodes*m_fNumIntPts);
    assert(m_fNormals.size() == m_Dim*m_fNum*m_fNumIntPts);

    gmsh::logger::write("==================================================");
    gmsh::logger::write("Number of Faces : "+to_string(m_fNum));
    gmsh::logger::write("Faces per Element : "+to_string(m_fNumPerEl));
    gmsh::logger::write("Face dimension : "+to_string(m_fDim));
    gmsh::logger::write("Face Type : "+m_fName);
    gmsh::logger::write("Face Nbr Nodes : "+to_string(m_fNumNodes));
    gmsh::logger::write("Integration type : "+m_fIntType);
    gmsh::logger::write("Integration Nbr points : "+to_string(m_fNumIntPts));

    // Connectivity
    // Assign corresponding faces to each element, we use the fact that the node tags per face has already been ordered
    
    m_fNbrElIds.resize(m_fNum);

    for(int el=0; el<m_elNum; ++el){
        for(int elF=0; elF<m_fNumPerEl; ++elF){
            for(int f=0; f<m_fNum; ++f){

                if(equal(&fNodeTagOrdered(f),&fNodeTagOrdered(f)+m_fNumNodes,&elFNodeTagOrdered(el,elF))){
                    m_elFIds.push_back(f);
                    m_fNbrElIds[f].push_back(el);
                }
            }
        }
    }

    // For efficiency purposes we also directly store the mapping between face node id and element node id
    // For example, the 3rd node of the face correspond to the 7th of the element

    m_fNToElNIds.resize(m_fNum);

    for(int f=0; f<m_fNum; ++f){
        for(int nf=0; nf<m_fNumNodes; ++nf){
            for(int el : m_fNbrElIds[f]){
                for(int nel=0; nel<m_elNumNodes; ++nel){
                    if(fNodeTag(f,nf) == elNodeTag(el,nel)){m_fNToElNIds[f].push_back(nel);}
                }
            }
        }
    }
    
    // Up to now, the normals are associated to the faces
    // We still need to know how the normal is oriented with respect to its neighbouring elements
    // For instance, element1 normal has the same orientation, we therefore assign the orientation +1
    // Reciprocally we set the orientation to -1 for e2.

    //  ______      f1       ______
    // |      |     |       |      |
    // |  e1  |->   |->   <-|  e2  |
    // |______|     |       |______|

    double dotProduct;
    vector<double> m_elBarycenters,fNodeCoord(3),elOuterDir(3),paramCoords;
    gmsh::model::mesh::getBarycenters(m_elType[0],-1,false,true,m_elBarycenters);

    for(int el=0; el<m_elNum; ++el){
        for(int f=0; f<m_fNumPerEl; ++f){

            dotProduct = 0.0;
            int dim,tag;
            gmsh::model::mesh::getNode(elFNodeTag(el,f),fNodeCoord,paramCoords,dim,tag);

            for(int x=0; x<m_Dim; x++){
                elOuterDir[x] = fNodeCoord[x] - m_elBarycenters[el*3+x];
                dotProduct += elOuterDir[x]*fNormal(elFId(el,f),0,x);
            }
            if(dotProduct>= 0){m_elFOrientation.push_back(1);}
            else{m_elFOrientation.push_back(-1);}
        }
    }

    // Once the orientation known, we reclassify the neighbouring elements by imposing
    // the first one to be oriented in the same direction than the corresponding face

    int elf;

    for(int f=0; f<m_fNum; ++f){
        for(int lf=0; lf<m_fNumPerEl; ++lf){
            if(elFId(fNbrElId(f,0),lf) == f){elf = lf;}
        }

        if(m_fNbrElIds.size() == 2){
            if(elFOrientation(fNbrElId(f,0),elf)<= 0){
                swap(m_fNbrElIds[f][0],m_fNbrElIds[f][1]);
                for(int nf=0; nf<m_fNumNodes; ++nf){swap(fNToElNId(f,nf,0),fNToElNId(f,nf,1));}
            }
        }
    }

    assert(m_elFIds.size() == m_elNum*m_fNumPerEl);
    assert(m_elFOrientation.size() == m_elNum*m_fNumPerEl);
    gmsh::logger::write("==================================================");
    gmsh::logger::write("Element-Face connectivity retrieved");

    // Boundary conditions
    // Check if a face is a boundary or not and orientate the normal at boundaries in the outward direction
    // This convention is particularly useful for boundary conditions

    for(int f=0; f<m_fNum; ++f){
        if(m_fNbrElIds[f].size()<2){

            m_fIsBoundary.push_back(true);
            for(int lf=0; lf<m_fNumPerEl; ++lf){
                if(elFId(fNbrElId(f,0),lf) == f){
                    for(int g=0; g<m_fNumIntPts; ++g){

                        fNormal(f,g,0) *= elFOrientation(fNbrElId(f,0),lf);
                        fNormal(f,g,1) *= elFOrientation(fNbrElId(f,0),lf);
                        fNormal(f,g,2) *= elFOrientation(fNbrElId(f,0),lf);
                    }
                    elFOrientation(fNbrElId(f,0),lf) = 1;
                }
            }
        }
        else{m_fIsBoundary.push_back(false);}
    }


    // Iterate over the physical boundaries and over each nodes belonging to that boundary
    // Retrieve the associated face and assign it an unique integer representing the BC type
    // Default = Absorbing, 1 = Reflecting, 2 = Absorbing
     
    m_fBC.resize(m_fNum);
    vector<std::size_t> nodeTags;
    vector<double> coord;

    for (auto const& physBC : config.physBCs){

        auto physTag = physBC.first;
        auto BCtype = physBC.second.first;
        auto BCvalue = physBC.second.second;
        gmsh::model::mesh::getNodesForPhysicalGroup(m_fDim,physTag,nodeTags,coord);

        if(BCtype == "Reflecting"){
            for(int f=0; f<m_fNum; ++f){
                if(m_fIsBoundary[f] && find(nodeTags.begin(),nodeTags.end(),fNodeTag(f)) != nodeTags.end()){m_fBC[f] = 1;}
            }
        }
        else{
            for(int f=0; f<m_fNum; ++f){
                if(m_fIsBoundary[f] && find(nodeTags.begin(),nodeTags.end(),fNodeTag(f)) != nodeTags.end()){m_fBC[f] = 0;}
            }
        }
    }

    // Computes the R*K*R matrix product
    // This matrix product is related to the absorbing boundary conditions in the specific context of acoustic waves
    // It is mainly used to suppress the outgoing solution along characteristics lines
    
    RKR.resize(m_fNum*m_fNumIntPts);

    for(int f=0; f<m_fNum; ++f){
        if(m_fBC[f] == 0){
            for(int g=0; g<m_fNumIntPts; ++g){
                
                int i = f*m_fNumIntPts+g;
                RKR[i].resize(16);

                RKR[i][0] = 0.25*config.c0;
                RKR[i][1] = 0.25*config.c0*config.c0*config.rho0*fNormal(f,g,0);
                RKR[i][2] = 0.25*config.c0*config.c0*config.rho0*fNormal(f,g,1);
                RKR[i][3] = 0.25*config.c0*config.c0*config.rho0*fNormal(f,g,2);

                RKR[i][4] = 0.25*fNormal(f,g,0)/config.rho0;
                RKR[i][5] = 0.25*config.c0*fNormal(f,g,0)*fNormal(f,g,0);
                RKR[i][6] = 0.25*config.c0*fNormal(f,g,0)*fNormal(f,g,1);
                RKR[i][7] = 0.25*config.c0*fNormal(f,g,0)*fNormal(f,g,2);

                RKR[i][8] = 0.25*fNormal(f,g,1)/config.rho0;
                RKR[i][9] = 0.25*config.c0*fNormal(f,g,1)*fNormal(f,g,0);
                RKR[i][10] = 0.25*config.c0*fNormal(f,g,1)*fNormal(f,g,1);
                RKR[i][11] = 0.25*config.c0*fNormal(f,g,1)*fNormal(f,g,2);

                RKR[i][12] = 0.25*fNormal(f,g,2)/config.rho0;
                RKR[i][13] = 0.25*config.c0*fNormal(f,g,2)*fNormal(f,g,0);
                RKR[i][14] = 0.25*config.c0*fNormal(f,g,2)*fNormal(f,g,1);
                RKR[i][15] = 0.25*config.c0*fNormal(f,g,2)*fNormal(f,g,2);
            }
        }
    }

    assert(m_fIsBoundary.size() == m_fNum);

    // Extra Memory allocation
    // Instantiate Ghost Elements and numerical flux storage
    
    m_fFlux.resize(m_fNum*m_fNumNodes);
    uGhost = vector<vector<double>>(4,vector<double>(m_fNum*m_fNumIntPts));
    FluxGhost = vector<vector<vector<double>>>(4,vector<vector<double>>(m_fNum*m_fNumIntPts,vector<double>(3)));
    gmsh::logger::write("Boundary conditions successfuly loaded");
    gmsh::logger::write("==================================================");
}

// Precompute and store the mass matris for all elements in m_elMassMatrix

void Mesh::precomputeMassMatrix(){

    m_elMassMatrices.resize(m_elNum*m_elNumNodes*m_elNumNodes);
    for(int el=0; el<m_elNum; ++el){getElMassMatrix(el,true,&elMassMatrix(el));}
}

// Compute the element mass matrix
// @param el integer: element id (!= gmsh tag,it is the location in memory storage)
// @param inverse boolean: whether or not the mass matrix must be inverted before returned
// @param elMassMatrix double array: output storage of the element mass matrix

void Mesh::getElMassMatrix(const int el,const bool inverse,double *elMassMatrix){

    for(int i=0; i<m_elNumNodes; ++i){
        for (int j=0; j<m_elNumNodes; ++j){
            elMassMatrix[i*m_elNumNodes+j] = 0.0;
            for(int g=0; g<m_elNumIntPts; g++){
                elMassMatrix[i*m_elNumNodes+j] += elBasisFct(g,i)*elBasisFct(g,j)*m_elWeight[g]*elJacobianDet(el,g);
            }
        }
    }
    if(inverse){eigen::inverse(elMassMatrix,m_elNumNodes);}
}

// Compute the element stiffness/convection matrix
// @param el integer: element id
// @param Flux double array: physical flux
// @param u double array: solution at element node
// @param elStiffVector double array: output storage of the element stiffness vector

void Mesh::getElStiffVector(const int el,vector<vector<double>> &Flux,vector<double> &u,double *elStiffVector){

    int jId;

    for(int i=0; i<m_elNumNodes; ++i){
        elStiffVector[i] = 0.0;
        for (int j=0; j<m_elNumNodes; ++j){
            jId = el*m_elNumNodes+j;
            for(int g=0; g<m_elNumIntPts; g++){
                elStiffVector[i] += eigen::dot(Flux[jId].data(),&elGradBasisFct(el,g,i),m_Dim)*elBasisFct(g,j)*m_elWeight[g]*elJacobianDet(el,g);
            }
        }
    }
}

// Precompute the numerical flux through all the faces
// The flux implemented is the Rusanov Flux
// Also note that the following code is paralelized using openMP
// @param Flux double array: physical flux
// @param u double array: solution at the node
// @param eq: equation id (0 = pressure, 1 = Vx, 2= Vy, 3= Vz)

void Mesh::precomputeFlux(vector<double> &u,vector<vector<double>> &Flux,int eq){

    #pragma omp parallel num_threads(config.numThreads)
    {
        // Memory allocation (Cross-plateform compatibility)

        int elUp,elDn;
        vector<double> FIntPts(m_fNumIntPts,0);
        vector<double> Fnum(m_Dim,0);

        #pragma omp parallel for schedule(static)
        for(int f=0; f<m_fNum; ++f){

            fill(FIntPts.begin(),FIntPts.end(),0);

            // Numerical Flux at Integration points

            if (m_fIsBoundary[f]){
                for (int g=0; g<m_fNumIntPts; ++g){FIntPts[g] = FluxGhost[eq][f*m_fNumIntPts+g][0];}
            }
            else{
                for (int i=0; i<m_fNumNodes; ++i){
                    elUp = fNbrElId(f,0)*m_elNumNodes+fNToElNId(f,i,0);
                    elDn = fNbrElId(f,1)*m_elNumNodes+fNToElNId(f,i,1);

                    for (int g=0; g<m_fNumIntPts; ++g){
                        for (int x=0; x<m_Dim; ++x){Fnum[x] = 0.5*((Flux[elUp][x]+Flux[elDn][x])+fc*config.c0*fNormal(f,g,x)*(u[elUp]-u[elDn]));}
                        FIntPts[g] += eigen::dot(&fNormal(f,g),Fnum.data(),m_Dim)*fBasisFct(g,i);
                    }
                }
            }
            
            // Surface integral

            for (int n=0; n<m_fNumNodes; ++n){
                fFlux(f,n) = 0;
                for (int g=0; g<m_fNumIntPts; ++g){fFlux(f,n) += m_fWeight[g]*fBasisFct(g,n)*FIntPts[g]*fJacobianDet(f,g);}
            }
        }
    }
}

// Compute flux through a given element from the value of the flux at the face
// @param el integer: element id
// @param F double array: Output element flux

void Mesh::getElFlux(const int el,double* F){

    int i;
    fill(F,F+m_elNumNodes,0);

    for(int f=0; f<m_fNumPerEl; ++f){
        el == fNbrElId(elFId(el,f),0)? i = 0 : i = 1;
        for(int nf=0; nf<m_fNumNodes; ++nf){F[fNToElNId(elFId(el,f),nf,i)] += elFOrientation(el,f)*fFlux(elFId(el,f),nf);}
    }
}

// Compute physical flux from the nodal solution. Also update the ghost element and numerical flux
// @param u: nodal solution vector
// @param Flux: Physical flux
// @param v0: mean flow speed (V0x,V0y,V0z)
// @param c0: speed of sound
// @param rho0: mean flow density

void Mesh::updateFlux(vector<vector<double>> &u,vector<vector<vector<double>>> &Flux,vector<double> &v0,double c0,double rho0){

    for(int el=0; el<m_elNum; ++el){
        for(int n=0; n<m_elNumNodes; ++n){
            int i = el*m_elNumNodes+n;

            // Pressure flux, Vx, Vy, Vy

            Flux[0][i] ={v0[0]*u[0][i]+rho0*c0*c0*u[1][i],v0[1]*u[0][i]+rho0*c0*c0*u[2][i],v0[2]*u[0][i]+rho0*c0*c0*u[3][i]};
            Flux[1][i] ={v0[0]*u[1][i]+u[0][i]/rho0,v0[1]*u[1][i],v0[2]*u[1][i]};
            Flux[2][i] ={v0[0]*u[2][i],v0[1]*u[2][i]+u[0][i]/rho0,v0[2]*u[2][i]};
            Flux[3][i] ={v0[0]*u[3][i],v0[1]*u[3][i],v0[2]*u[3][i]+u[0][i]/rho0};
        }

        // Ghost elements

        for(int f=0; f<m_fNumPerEl; ++f){
            int fId = elFId(el,f);
            if(m_fIsBoundary[fId]){
                for(int g=0; g<m_fNumIntPts; ++g){
                    int gId = fId*m_fNumIntPts+g;

                    // Interpolate solution at integration points

                    uGhost[0][gId] = 0;
                    uGhost[1][gId] = 0;
                    uGhost[2][gId] = 0;
                    uGhost[3][gId] = 0;

                    for(int n=0; n<m_fNumNodes; ++n){

                        int nId = el*m_elNumNodes+fNToElNId(fId,n,0);
                        uGhost[0][gId] += u[0][nId]*fBasisFct(g,n);
                        uGhost[1][gId] += u[1][nId]*fBasisFct(g,n);
                        uGhost[2][gId] += u[2][nId]*fBasisFct(g,n);
                        uGhost[3][gId] += u[3][nId]*fBasisFct(g,n);
                    }

                    if(m_fBC[fId] == 1){

                        // Normal component of velocity

                        double dot = fNormal(fId,g,0)*uGhost[1][gId] +fNormal(fId,g,1)*uGhost[2][gId] +fNormal(fId,g,2)*uGhost[3][gId];

                        // Removes normal component (rigid wall BC)

                        uGhost[1][gId] -= dot*fNormal(fId,g,0);
                        uGhost[2][gId] -= dot*fNormal(fId,g,1);
                        uGhost[3][gId] -= dot*fNormal(fId,g,2);

                        // Flux at integration points (Pressure flux, Vx, Vy, Vz)

                        FluxGhost[0][gId] ={v0[0]*uGhost[0][gId]+rho0*c0*c0*uGhost[1][gId],v0[1]*uGhost[0][gId]+rho0*c0*c0*uGhost[2][gId],v0[2]*uGhost[0][gId]+rho0*c0*c0*uGhost[3][gId]};
                        FluxGhost[1][gId] ={v0[0]*uGhost[1][gId]+uGhost[0][gId]/rho0,v0[1]*uGhost[1][gId],v0[2]*uGhost[1][gId]};
                        FluxGhost[2][gId] ={v0[0]*uGhost[2][gId],v0[1]*uGhost[2][gId]+uGhost[0][gId]/rho0,v0[2]*uGhost[2][gId]};
                        FluxGhost[3][gId] ={v0[0]*uGhost[3][gId],v0[1]*uGhost[3][gId],v0[2]*uGhost[3][gId]+uGhost[0][gId]/rho0};

                        // Project Flux on the normal

                        for(int eq=0; eq<4; ++eq){FluxGhost[eq][gId][0] = eigen::dot(&fNormal(fId,g),&FluxGhost[eq][gId][0],m_Dim);}
                    }

                    else{

                        // Absorbing boundary conditions, the flux already projected on the normal

                        FluxGhost[0][gId][0] = RKR[gId][0]*uGhost[0][gId] +RKR[gId][1]*uGhost[1][gId] +RKR[gId][2]*uGhost[2][gId] +RKR[gId][3]*uGhost[3][gId];
                        FluxGhost[1][gId][0] = RKR[gId][4]*uGhost[0][gId] +RKR[gId][5]*uGhost[1][gId] +RKR[gId][6]*uGhost[2][gId] +RKR[gId][7]*uGhost[3][gId];
                        FluxGhost[2][gId][0] = RKR[gId][8]*uGhost[0][gId] +RKR[gId][9]*uGhost[1][gId] +RKR[gId][10]*uGhost[2][gId] +RKR[gId][11]*uGhost[3][gId];
                        FluxGhost[3][gId][0] = RKR[gId][12]*uGhost[0][gId] +RKR[gId][13]*uGhost[1][gId] +RKR[gId][14]*uGhost[2][gId] +RKR[gId][15]*uGhost[3][gId];
                    }
                }
            }
        }
    }
}

// List of nodes for each unique face given a list of node per face and per elements

void Mesh::getUniqueFaceNodeTags(){

    // Ordering per face for efficient comparison

    m_elFNodeTagsOrdered = m_elFNodeTags;
    for(int i=0; i<m_elFNodeTagsOrdered.size(); i+=m_fNumNodes){
        sort(m_elFNodeTagsOrdered.begin()+i,m_elFNodeTagsOrdered.begin()+(i+m_fNumNodes));
    }

    // Unordered keep gmsh order while ordered array are used for comparison

    m_fNodeTags = m_elFNodeTags;
    m_fNodeTagsOrdered = m_elFNodeTagsOrdered;

    // Remove identical faces by comparing ordered arrays

    vector<std::size_t>::iterator it_delete;
    vector<std::size_t>::iterator it_deleteUnordered;
    vector<std::size_t>::iterator it_unordered = m_fNodeTags.begin();

    for(vector<std::size_t>::iterator it_ordered = m_fNodeTagsOrdered.begin(); it_ordered != m_fNodeTagsOrdered.end();){
        it_deleteUnordered = it_unordered+m_fNumNodes;

        for(it_delete=it_ordered+m_fNumNodes; it_delete != m_fNodeTagsOrdered.end(); it_delete+=m_fNumNodes){
            if(equal(it_ordered,it_ordered+m_fNumNodes,it_delete)){break;}
            it_deleteUnordered+=m_fNumNodes;
        }

        if(it_delete != m_fNodeTagsOrdered.end()){
            m_fNodeTagsOrdered.erase(it_delete,it_delete+m_fNumNodes);
            m_fNodeTags.erase(it_deleteUnordered,it_deleteUnordered+m_fNumNodes);
        }

        else{
            it_ordered+=m_fNumNodes;
            it_unordered+=m_fNumNodes;
        }
    }
}
