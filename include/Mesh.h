#include "configParser.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <string>
#include <map>

#ifndef DGALERKIN_MESH_H
#define DGALERKIN_MESH_H
using namespace std;

class Mesh{

public:

    Mesh(string name,Config config);

    // List of getters used for vector access and to improve readability
    // The inline optimization (equivalent of c macro) is used to remove the associated function overhead
    // el = element indices in memory storage (!= element tag)
    // n,i = element or face node (0,..,elNumNode)
    // u = parametric coordinates (0=u,1=v,2=w)
    // x = physical coordinates (0=x,1=y,2=z)
    // f = face indices in memory storage
    // g = integration point indices
    
    inline std::size_t& elTag(int el){
        return m_elTags[el];
    };
    inline std::size_t& elNodeTag(int el,int n=0){
        return m_elNodeTags[el*m_elNumNodes+n];
    };
    inline double& elJacobian(int el,int g=0,int x=0,int u=0){
        return m_elJacobians[el*m_elNumIntPts*9+g*9+u*3+x];
    };
    inline double& elJacobianDet(int el,int g=0){
        return m_elJacobianDets[el*m_elNumIntPts+g];
    };
    inline double& elWeight(int g){
        return m_elIntParamCoords[g*4+3];
    };
    inline double& elBasisFct(int g,int i=0){
        return m_elBasisFcts[g*m_elNumNodes+i];
    };
    inline double& elUGradBasisFct(int g,int i=0,int u=0){
        return m_elUGradBasisFcts[g*m_elNumNodes*3+i*3+u];
    };
    inline double& elGradBasisFct(int el,int g=0,int i=0,int x=0){
        return m_elGradBasisFcts[el*m_elNumIntPts*m_elNumNodes*3+g*m_elNumNodes*3+i*3+x];
    };
    inline std::size_t& elFNodeTag(int el,int f=0,int i=0){
        return m_elFNodeTags[el*m_fNumPerEl*m_fNumNodes+f*m_fNumNodes+i];
    };
    inline std::size_t& fNodeTag(int f,int i=0){
        return m_fNodeTags[f*m_fNumNodes+i];
    };
    inline std::size_t& elFNodeTagOrdered(int el,int f=0,int i=0){
        return m_elFNodeTagsOrdered[el*m_fNumPerEl*m_fNumNodes+f*m_fNumNodes+i];
    };
    inline std::size_t& fNodeTagOrdered(int f,int i=0){
        return m_fNodeTagsOrdered[f*m_fNumNodes+i];
    };
    inline double& fJacobian(int f,int g=0,int x=0,int u=0){
        return m_fJacobians[f*m_fNumIntPts*9+g*9+u*3+x];
    };
    inline double& fJacobianDet(int f,int g=0){
        return m_fJacobianDets[f*m_fNumIntPts+g];
    };
    inline double& fUGradBasisFct(int g,int i=0,int u=0){
        return m_fUGradBasisFcts[g*m_fNumNodes*3+i*3+u];
    };
    inline double& fGradBasisFct(int f,int g=0,int i=0,int x=0){
        return m_fGradBasisFcts[f*m_fNumIntPts*m_fNumNodes*3+g*m_fNumNodes*3+i*3+x];
    };
    inline double& fIntPtCoord(int f,int g=0,int x=0){
        return m_fIntPtCoords[f*m_fNumIntPts*3+g*3+x];
    };
    inline double& fWeight(int g){
        return m_fIntParamCoords[g*4+3];
    };
    inline double& fBasisFct(int g,int i=0){
        return m_fBasisFcts[g*m_fNumNodes+i];
    };
    inline double& fNormal(int f,int g=0,int x=0){
        return m_fNormals[f*3*m_fNumIntPts+g*3+x];
    };
    inline int &elFId(int el,int f=0){
        return m_elFIds[el*m_fNumPerEl+f];
    };
    inline int &fNbrElId(int f,int el=0){
        return m_fNbrElIds[f][el];
    };
    inline int &fNToElNId(int f,int nf=0,int el=0){
        return m_fNToElNIds[f][nf*m_fNbrElIds[f].size()+el];
    };
    inline int &elFOrientation( int el,int f){
        return m_elFOrientation[el*m_fNumPerEl+f];
    };
    inline double &elMassMatrix(int el,int i=0,int j=0){
        return m_elMassMatrices[el*m_elNumNodes*m_elNumNodes+i*m_elNumNodes+j];
    };
    inline double &fFlux(int f,int n=0){
        return m_fFlux[f*m_fNumNodes+n];
    };

    // Extra getters
    
    vector<std::size_t> const &getElNodeTags(){return m_elNodeTags;}
    int getNumNodes(){return m_elNodeTags.size();}
    int getElNumNodes(){return m_elNumNodes;}
    int getElNum(){return m_elNum;}

    // Matrices and vectors assembly
    
    void updateFlux(vector<vector<double>> &u,vector<vector<vector<double>>> &Flux,vector<double> &v0,double c0,double rho0);
    void getElStiffVector(int el,vector<vector<double>> &Flux,vector<double> &u,double *elStiffVector);
    void precomputeFlux(vector<double> &u,vector<vector<double>> &Flux,int eq);
    void getElMassMatrix(int el,bool inverse,double *elMassMatrix);
    void getElFlux(int el,double* F);
    void getUniqueFaceNodeTags();
    void precomputeMassMatrix();

private:

    string name;                                // Associated mesh file string
    Config config;                              // Configuration object

    int fc=1;					                // Numerical flux coefficient	
    int m_elNum;                                // Number of elements in dim
    int m_elDim;                                // Dimension of the element (and the domain)
    int m_Dim = 3;                              // Physical space dimension
    int m_elOrder;                              // Element Order
    int m_elNumNodes;                           // Number of nodes per element
    int m_elNumIntPts;                          // Number of integration points

    string m_elName;                            // Element Type name
    string m_elIntType;                         // Integration type name
    vector<int> m_elType;                       // Element Types (integer)
    vector<std::size_t> m_elTags;                       // Tags of the elements
    vector<int> m_elFIds;                       // Faces ids for each element [e1f1,e1f2,..,e2f1,e2f2,..]
    vector<std::size_t> m_elNodeTags;                   // Tags of the nodes associated to each element [e1n1,e1n2,..,e2n1,e2n2,..]
    vector<std::size_t> m_elFNodeTags;                  // Node tags for each face and each element
    vector<int> m_elFOrientation;               // Contains 1 or -1, if the outward element face is in the same direction as the face normal or -1 if not [e1f1,e1f2,..,e2f1,e2f2]
    vector<double> m_elBasisFcts;               // Evaluation of the basis functions at the integration points [g1f1,g1f2,..,g2f1,g2f2,..]
    vector<double> m_elJacobians;               // Jacobian evaluated at each integration points (dx/du) [e1g1Jxx,e1g1Jxy,e1g1Jxz,..,e1gGJzz,e2g1Jxx,..]
    vector<double> m_elParamCoord;              // Parametric coordinates of the element
    vector<double> m_elIntPtCoords;             // x,y,z coordinates of the integration points element by element [e1g1x,e1g1y,e1g1z,.. ,e1gGz,e2g1x,..]
    vector<double> m_elJacobianDets;            // Determinants of the jacobian evaluated at each integration points [e1g1DetJ,e1g2DetJ,.. e2g1DetJ,e2g2DetJ,..]
    vector<double> m_elGradBasisFcts;           // Derivatives of the basis functions at the integration points [e1g1df1/dx,e1g1df1/dy,..,e1g2df1/dx,e1g2df1/dy,..,e1g1df2/dx,e1g1df2/dy,..]
    vector<std::size_t> m_elFNodeTagsOrdered;           // Ordered used for comparison, Unordered preserve GMSH ordering and locality [e1f1n1,e1f1n2,..,e1f2n1,e1f2n2,..,e2f1n1,e2f1n2,..]
    vector<double> m_elIntParamCoords;          // u,v,w coordinates and the weight q for each integration point [g1u,g1v,g1w,g1q,g2u,..]
    vector<double> m_elUGradBasisFcts;          // Derivatives of the basis functions at the integration points [g1df1/du,g1df1/dv,..,g2df1/du,g2df1/dv,..,g1df2/du,g1df2/dv,..]

    std::vector<double> m_elWeight;

    int m_fDim;                                 // Face dimension
    string m_fName;                             // Face type name
    int m_fType;                                // Face type
    int m_fNumNodes;                            // Number of nodes per face
    int m_fNumPerEl;                            // Number of faces per element
    int m_fEntity;                              // Entity containing all the faces
    int m_fNum;                                 // Number of unique faces
    int m_fNumIntPts;                           // Number of integration points on each face

    string m_fIntType;                          // Face integration type (same as element)
    vector<std::size_t> m_fNodeTags;                    // Node tags for each unique face
    vector<std::size_t> m_fNodeTagsOrdered;             // [f1n1,f1n2,..,f2n1,f2n2,..]
    vector<std::size_t> m_fTags;                        // Tag for each unique face [f1,f2,f3,..]
    vector<double> m_fJacobians;                // Jacobian evaluated at each integration points (dx/du) [f1g1Jxx,f1g1Jxy,f1g1Jxz,.. f1g1Jzz,f1g2Jxx,..,f1gGJzz,f2g1Jxx,..]
    vector<double> m_fJacobianDets;             // Determinants of the jacobian evaluated at each integration points [f1g1DetJ,f1g2DetJ,.. f2g1DetJ,f2g2DetJ,..]
    vector<double> m_fIntPtCoords;              // x,y,z coordinates of the integration points for each faces [f1g1x,f1g1y,f1g1z,.. ,f1gGz,f2g1x,..]
    vector<double> m_fIntParamCoords;           // u,v,w coordinates and the weight q for each integration point [g1u,g1v,g1w,g1q,g2u,..]
    vector<double> m_fBasisFcts;                // Basis functions at the integration points [g1f1,g1f2,..,g2f1,g2f2,..]
    vector<double> m_fUGradBasisFcts;           // Derivatives of the basis functions at the integration points [g1df1/du,g1df1/dv,..,g2df1/du,g2df1/dv,..,g1df2/du,g1df2/dv,..]
    vector<double> m_fGradBasisFcts;            // Derivatives of the basis functions at the integration points [e1g1df1/dx,e1g1df1/dy,..,e1g2df1/dx,e1g2df1/dy,..,e1g1df2/dx,e1g1df2/dy,..]
    vector<double> m_fNormals;                  // Normal for each face at each int point [f1g1Nx,f1g1Ny,f1g1Nz,f1g2Nx,..,f2g1Nx,f2g1Ny,f2g1Nz,..]

    std::vector<double> m_fWeight;

    vector<vector<int>> m_fNbrElIds;            // Id of element of each side of the face
    vector<vector<int>> m_fNToElNIds;           // Map face node Ids to element node Ids
    vector<double> m_elMassMatrices;            // Element mass matrix stored contiguously (row major) [e1m11,e1m12,..,e1m21,e1m22,..,e2m11,..]
    vector<double> m_fFlux;                     // Flux through all faces [f1n1,f1n2,..,f2n1,f2n2,..]
    vector<bool> m_fIsBoundary;                 // Is Face a boundary
    vector<int> m_fBC;                          // Boundary type

    vector<vector<double>> RKR;                 // R*K*R^-1 matrix product,absorbing boundary
    vector<vector<double>> uGhost;              // Ghost element nodal solution
    vector<vector<vector<double>>> FluxGhost;   // Ghost flux
};

#endif