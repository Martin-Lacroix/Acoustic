#include <iostream>
#include <iomanip>

#ifndef DGALERKIN_UTILS_H
#define DGALERKIN_UTILS_H
using namespace std;

namespace eigen{

    void inverse(double *A,int &N);
    void solve(double *A,double *B,int &N);
    void normalize(double *A,int &N);
    double dot(double *A,double *b,int N);
    void linEq(double *A,double *X,double *Y,double &alpha,double beta,int &N);
    void minus(double *A,double *B,int N);
    void plus(double *A,double *B,int N);
    void plusTimes(double *A,double *B,double c,int N);
    void cross(double *A,double *B,double *OUT);
}

namespace display{

    template<typename Container>
    void print(const Container& cont,int row=1,bool colMajor=false){

        if(colMajor){
            for(int rowIt=0; rowIt<row; ++rowIt){

                int colIt = 0;
                for (auto const& x : cont){

                    if(colIt%row == rowIt){cout << setprecision(4) << left << setw(10) << x << " ";}
                    colIt++;
                }
                cout << endl;
            }
        }
        else{
            int colIt = 0;
            for (auto const& x : cont){

                colIt++;
                cout << setprecision(4) << left << setw(10) << x << " ";
                if(colIt%row == 0){cout << endl;}
            }
        }
    }
}

#endif