#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "lrmat.hpp"
#include "cluster.hpp"
#include "hmatrix.hpp"
#include "loading.hpp"
#include "export.hpp"
#include "user.hpp"
#include "parametres.hpp"


#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace std;

/**************************************************************************//**
* It builds the hierarchical matrix with compressed and dense blocks,
* computes the consistency error for a matrix vector product and
* the relative error in Frobenius norm with respect to the dense matrix.
* It also produces an output file to visualize the compression of the matrix
* (use graphes_output_local_compression.py in postprocessing folder).
*
* (To be run it requires the input file with the desidered parameters)
*****************************************************************************/

int main(int argc, char* argv[]){
	
	
	////////////////========================================================////////////////
	////////////////////////////////========  Input ========////////////////////////////////
	
	// Check the number of parameters
	if (argc < 2) {
		// Tell the user how to run the program
		cerr << "Usage: " << argv[0] << " input name" << endl;
		/* "Usage messages" are a conventional way of telling the user
		 * how to run a program if they enter the command incorrectly.
		 */
		return 1;
	}
	
	// Load the inputs
	string inputname = argv[1];
	LoadParamIO(inputname);
	LoadParam(inputname);
	
	cout<<"############# Inputs #############"<<endl;
	cout<<"Eta : "+NbrToStr(GetEta())<<endl;
	cout<<"Epsilon : "+NbrToStr(GetEpsilon())<<endl;
	cout<<"Output path : "+GetOutputPath()<<endl;
	cout<<"Mesh path : "+GetMeshPath()<<endl;
	cout<<"Matrix path : "+GetMatrixPath()<<endl;
	cout<<"##################################"<<endl;
 
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////    Build Hmatrix 	////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	//vector<double> times;
	vectReal r;
	vectR3   x;
	Matrix   A;
	tic();
	LoadMatrix((GetMatrixPath()).c_str(),A);
	LoadPoints((GetMeshPath()).c_str(),x,r);
	vectInt tab(nb_rows(A));
	for (int j=0;j<x.size();j++){
		tab[3*j]  = j;
		tab[3*j+1]= j;
		tab[3*j+2]= j;
	}
	toc();
	tic();
	HMatrix B(A,x,r,tab);
	toc();
	
	
	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////     Errors 	//////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	
	// Vecteur
	int nr  = nb_rows(A);
	vectCplx u(nr);
	int NbSpl = 1000;
	double du = 5./double(NbSpl);
	srand (1);
	for(int j=0; j<nr; j++){
		int n = rand()%(NbSpl+1);
		u[j] = n*du;}
	
	vectCplx ua(nr),ub(nr);
	MvProd(ua,A,u);
	MvProd(ub,B,u);
	Real err = norm(ua-ub)/norm(ua);
	Real compression=CompressionRate(B);
	
	cout<<"Matrix-vector product relative error : "<<err<<endl;
	cout<<"Compression rate: "<<compression<<endl;
	
	Real normA = NormFrob(A);
	//cout << "Frobenius norm of the dense matrix: " << normA << endl;
	
	Real froberrH = sqrt(squared_absolute_error(B,A))/normA;
	cout << "Relative error in Frobenius norm: " << froberrH << endl;

	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////      Output 	//////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	
	tic();
	Output(B, "output_local_comp_"+NbrToStr(GetEta())+"_"+NbrToStr(GetEpsilon())+"_"+GetMatrixName()); // to visualize the compression of the matrix
	toc();
	
}

