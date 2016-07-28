#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "cluster.hpp"
#include "hmatrix.hpp"
#include "loading.hpp"


using namespace std;

int main(){
    
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    Build a matrix A 	////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

	int nr = 100;
	// p1: random points in a unit disk, plane z=z1
	srand (1);
	double z1 = 1;
	vectR3 p1(nr);
	vectReal r1(nr);
	for(int j=0; j<nr; j++){
		double rho = ((double) rand() / (double)(RAND_MAX)); // (double) otherwise integer division!
		double theta = ((double) rand() / (double)(RAND_MAX));
		p1[j][0] = sqrt(rho)*cos(2*M_PI*theta); p1[j][1] = sqrt(rho)*sin(2*M_PI*theta); p1[j][2] = z1;
		r1[j]=1.e-16;
	}
	// p2: random points in a unit disk, plane z=z2
	double z2 = 10;
	vectR3 p2(nr);
	vectReal r2(nr);
	for(int j=0; j<nr; j++){
		double rho = ((double) rand() / (RAND_MAX)); // (double) otherwise integer division!
		double theta = ((double) rand() / (RAND_MAX));
		p2[j][0] = sqrt(rho)*cos(2*M_PI*theta); p2[j][1] = sqrt(rho)*sin(2*M_PI*theta); p2[j][2] = z2;
		r2[j]=1.e-16;
	}
	
	Matrix A(nr,nr);
	for(int j=0; j<nr; j++){
	for(int k=0; k<nr; k++){
	   	A(j,k) = 1./(4*M_PI*norm(p1[j]-p2[k]));
	}
	}
    
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    Build Hmatrix 	////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
	// Parametres
	Param Parametre(-1.,1.e-1);
	
//	HMatrix B(A,p1,r1,p2,r2);


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    Test MvProd 	////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
	// Vecteur (pseudo-)aleatoire
//    	vectCplx u(nr);
//   	int NbSpl = 1000;
//    	double du = 5./double(NbSpl);
//
//	for(int j=0; j<nr; j++){
//		int n = rand()%(NbSpl+1);
//		u[j] = n*du;
//	}
//	
//
//    
//	/* // Vector of 1
//	vectCplx u(nr);
//	for(int j=0; j<nr; j++){
//	u[j] = 1.;}*/
//    
//	vectCplx ua(nr),ub(nr);
//	MvProd(ua,A,u);
//	MvProd(ub,B,u);
//	Real err = norm(ua-ub)/norm(ua);
//	cout << "Erreur:\t" << err << endl;
	
    //cout << "Taux de compression:\t";
    //cout << CompressionRate(B) << endl;

}
