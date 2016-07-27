#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "lrmat.hpp"
using namespace std;

/**************************************************************************//**
 * Authomatic check of the LowRankMatrix structure
 * (for the test of the ACA algorithm precision see the src file test_ACA.cpp)
 *****************************************************************************/

int main(){
    
    srand (1);
    // we set a constant seed for rand because we want always the same result if we run the check many times
    // (two different initializations with the same seed will generate the same succession of results in the subsequent calls to rand)
    
    // Build matrix A with property for ACA
    int nr = 100;
    // p1: points in a unit disk of the plane z=z1
    double z1 = 1;
    vectR3 p1(nr);
    for(int j=0; j<nr; j++){
        double rho = ((double) rand() / (double)(RAND_MAX)); // (double) otherwise integer division!
        double theta = ((double) rand() / (double)(RAND_MAX));
        p1[j][0] = sqrt(rho)*cos(2*M_PI*theta); p1[j][1] = sqrt(rho)*sin(2*M_PI*theta); p1[j][2] = z1;
        // sqrt(rho) otherwise the points would be concentrated in the center of the disk
    }
    // p2: points in a unit disk of the plane z=z2
    double z2 = 10;
    vectR3 p2(nr);
    for(int j=0; j<nr; j++){
        double rho = ((double) rand() / (RAND_MAX)); // (double) otherwise integer division!
        double theta = ((double) rand() / (RAND_MAX));
        p2[j][0] = sqrt(rho)*cos(2*M_PI*theta); p2[j][1] = sqrt(rho)*sin(2*M_PI*theta); p2[j][2] = z2;
    }
    Matrix A(nr,nr);
    for(int j=0; j<nr; j++){
        for(int k=0; k<nr; k++){
            A(j,k) = 1./(4*M_PI*norm(p1[j]-p2[k]));
        }
    }
    
    LowRankMatrix B(A); // construct a low rank matrix B applying ACA to matrix A
    
    // Vecteur
    vectCplx u(nr);
    int NbSpl = 1000;
    double du = 5./double(NbSpl);
    for(int j=0; j<nr; j++){
        int n = rand()%(NbSpl+1);
        u[j] = n*du;}
    
    vectCplx ua(nr),ub(nr);
    MvProd(ua,A,u);
    MvProd(ub,B,u);
    Real err = norm(ua-ub)/norm(ua);
    //cout << "Erreur:\t" << err << endl;
    
    assert(abs(err-1.10797e-07)<1e-12);
}
