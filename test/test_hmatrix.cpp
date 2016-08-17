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
	SetNdofPerElt(1); // ndofperelt
	int nr = 100;
	// p1: random points in a unit disk, plane z=z1
	srand (1);
	double z1 = 1;
	vectR3 p1(nr);
	vectReal r1(nr);
	vectInt tab1(nr);
	for(int j=0; j<nr; j++){
		double rho = ((double) rand() / (double)(RAND_MAX)); // (double) otherwise integer division!
		double theta = ((double) rand() / (double)(RAND_MAX));
		p1[j][0] = sqrt(rho)*cos(2*M_PI*theta); p1[j][1] = sqrt(rho)*sin(2*M_PI*theta); p1[j][2] = z1;
		r1[j]=1.e-16;
		tab1[j]=j;
	}
	// p2: random points in a unit disk, plane z=z2
	double z2 = 10;
	vectR3 p2(nr);
	vectReal r2(nr);
	vectInt tab2(nr);
	for(int j=0; j<nr; j++){
		double rho = ((double) rand() / (RAND_MAX)); // (double) otherwise integer division!
		double theta = ((double) rand() / (RAND_MAX));
		p2[j][0] = sqrt(rho)*cos(2*M_PI*theta); p2[j][1] = sqrt(rho)*sin(2*M_PI*theta); p2[j][2] = z2;
		r2[j]=1.e-16;
		tab2[j]=j;
	}
	
	Matrix A1(nr,nr);
	Matrix A2(nr,nr);
	
	for(int j=0; j<nr; j++){
	for(int k=0; k<nr; k++){
	   	A1(j,k) = 1./(4*M_PI*norm(p1[j]-p2[k]));
		if (j!=k){
			A2(j,k) = 1./(4*M_PI*norm(p1[j]-p1[k]));
		}
		else{
			A2(j,k)=0;
		}
	}
	}
    
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    Build Hmatrix 	////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
	// Parametres
	SetEta(-1); // Pas de low rank matrices pour eta=-1
	SetEpsilon(1e-1);

	HMatrix B1(A1,p1,r1,tab1,p2,r2,tab2);
	HMatrix B2(A2,p1,r1,tab1,p1,r1,tab2);
	HMatrix B3(A2,p1,r1,tab1);


////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////    Test MvProd 	////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
	// Vecteur (pseudo-)aleatoire
    	vectCplx u(nr);
   	int NbSpl = 1000;
    	double du = 5./double(NbSpl);

	for(int j=0; j<nr; j++){
		int n = rand()%(NbSpl+1);
		u[j] = n*du;
	}
    
	vectCplx ua(nr),ub(nr);
	MvProd(ua,A1,u);
	MvProd(ub,B1,u);
	Real err = norm(ua-ub)/norm(ua);
	
	assert(abs(CompressionRate(B1))<1e-10);
	assert(err<1e-16);
	
	MvProd(ua,B2,u);
	MvProd(ub,B3,u);
	err = norm(ua-ub)/norm(ua);
	assert(abs(CompressionRate(B2)-CompressionRate(B3))<1e-16);
	assert(err<1e-16);
				      

}
