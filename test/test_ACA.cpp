#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "lrmat.hpp"
using namespace std;

int main(){
    
    // Build matrix A with property for ACA
    int nr = 100;
    // p1: random points in a unit disk, plane z=z1
    srand (time(NULL));
    double z1 = 1;
    vectR3 p1(nr);
    for(int j=0; j<nr; j++){
        double rho = ((double) rand() / (double)(RAND_MAX)); // (double) otherwise integer division!
        double theta = ((double) rand() / (double)(RAND_MAX));
        p1[j][0] = rho*cos(2*M_PI*theta); p1[j][1] = rho*sin(2*M_PI*theta); p1[j][2] = z1;
    }
    // p2: random points in a unit disk, plane z=z2
    srand (time(NULL));
    double z2 = 10;
    vectR3 p2(nr);
    for(int j=0; j<nr; j++){
        double rho = ((double) rand() / (RAND_MAX)); // (double) otherwise integer division!
        double theta = ((double) rand() / (RAND_MAX));
        p2[j][0] = rho*cos(2*M_PI*theta); p2[j][1] = rho*sin(2*M_PI*theta); p2[j][2] = z2;
    }
    Matrix A(nr,nr);
    for(int j=0; j<nr; j++){
        for(int k=0; k<nr; k++){
            A(j,k) = 1./(4*M_PI*norm(p1[j]-p2[k]));
        }
    }
    
    LowRankMatrix B(A, 1e-4); //ACA
    
    // Vecteur (pseudo-)aleatoire
    vectCplx u(nr);
    int NbSpl = 1000;
    double du = 5./double(NbSpl);
    srand (time(NULL));
    for(int j=0; j<nr; j++){
        int n = rand()%(NbSpl+1);
        u[j] = n*du;}
    
    /* // Vector of 1
     vectCplx u(nr);
     for(int j=0; j<nr; j++){
     u[j] = 1.;}*/
    
    vectCplx ua(nr),ub(nr);
    MvProd(ua,A,u);
    MvProd(ub,B,u);
    Real err = norm(ua-ub)/norm(ua);
    cout << "Erreur:\t" << err << endl;
    
    //cout << "Taux de compression:\t";
    //cout << CompressionRate(B) << endl;

}
