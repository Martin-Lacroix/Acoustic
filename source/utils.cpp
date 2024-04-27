#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

namespace eigen{

    void solve(double *A,double *B,int &N){
        Eigen::Map<Eigen::MatrixXd> A_eigen(A,N,N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B,N);
        B_eigen = A_eigen.lu().solve(B_eigen);
    }

    void inverse(double *A,int &N){
        Eigen::Map<Eigen::MatrixXd> A_eigen(A,N,N);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_eigen);
        A_eigen = A_eigen.inverse();
    }

    void normalize(double *A,int &N){
        Eigen::Map<Eigen::VectorXd> A_eigen(A,N);
        A_eigen.normalize();
    }

    void linEq(double *A,double *X,double *Y,double &alpha,double beta,int &N){
        Eigen::Map<Eigen::VectorXd> X_eigen(X,N);
        Eigen::Map<Eigen::VectorXd> Y_eigen(Y,N);
        Eigen::Map<Eigen::MatrixXd> A_eigen(A,N,N);
        Y_eigen = beta*Y_eigen +  alpha*A_eigen*X_eigen;
    }

    double dot(double *A,double *B,int N){
        Eigen::Map<Eigen::VectorXd> A_eigen(A,N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B,N);
        return A_eigen.dot(B_eigen);
    }

    void minus(double *A,double *B,int N){
        Eigen::Map<Eigen::VectorXd> A_eigen(A,N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B,N);
        A_eigen -= B_eigen;
    }

    void plus(double *A,double *B,int N){
        Eigen::Map<Eigen::VectorXd> A_eigen(A,N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B,N);
        A_eigen += B_eigen;
    }

    void plusTimes(double *A,double *B,double c,int N){
        Eigen::Map<Eigen::VectorXd> A_eigen(A,N);
        Eigen::Map<Eigen::VectorXd> B_eigen(B,N);
        A_eigen += c*B_eigen;
    }

    void cross(double *A,double *B,double *OUT){
        Eigen::Map<Eigen::Vector3d> A_eigen(A);
        Eigen::Map<Eigen::Vector3d> B_eigen(B);
        Eigen::Map<Eigen::Vector3d> OUT_eigen(OUT);
        OUT_eigen = A_eigen.cross(B_eigen);
    }
}