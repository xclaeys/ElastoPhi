#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "cluster.hpp"
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
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Build cluster 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	Cluster t(p1,r1); Cluster s(p2,r2);
	
	
	

}
