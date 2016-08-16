//#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "lrmat.hpp"
#include "cluster.hpp"
#include "hmatrix.hpp"
#include "sparsematrix.hpp"
#include "loading.hpp"
#include "export.hpp"
#include "user.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace std;

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
	
	cout<<"############# Inputs #############"<<endl;
	cout<<"Output path : "+GetOutputPath()<<endl;
	cout<<"Mesh path : "+GetMeshPath()<<endl;
	cout<<"Matrix path : "+GetMatrixPath()<<endl;
	cout<<"##################################"<<endl;
	
	
	vector<double> times2;
	
	vectReal r;
	vectR3   x;
	Matrix   A;
	SpMatrix spA;
	
	tic();
    LoadSpMatrix((GetDataPath()+'/'+split(GetMatrixName(),'.').at(0)+"Creuse.txt").c_str(),spA);
	toc();
	
	tic();
	LoadPoints(GetMeshPath().c_str(),x,r);
	toc();
	tic();
	LoadMatrix(GetMatrixPath().c_str(),A);
	toc();
	
	// Vecteur pour le produit matrice vecteur
	int nr  = nb_rows(A);
	vectCplx u(nr);
	int NbSpl = 1000;
	double du = 5./double(NbSpl);
	srand (1); // !! pour reproducibilite'
	for(int j=0; j<nr; j++){
		int n = rand()%(NbSpl+1);
		u[j] = n*du;}
	
	vectCplx ua(nr);
	MvProd(ua,A,u);
		
	tic();
	Real normA = NormFrob(A);
	toc();
	cout << "Frobenius norm of the dense matrix: " << normA << endl;
	
	// Vecteur renvoyant pour chaque dof l'indice de l'entite geometrique correspondante dans x
	vectInt tab(nb_rows(A));
	for (int j=0;j<x.size();j++){
		tab[3*j]  = j;
		tab[3*j+1]= j;
		tab[3*j+2]= j;
	}
	
	// Values of eta and epsilon
	const int neta = 1;
	double eta [neta] = {5e-1};
	const int nepsilon = 3;
	double epsilon[nepsilon] = {-1, 1e0, 1e-1};
	
	// for output file
	string filename=GetOutputPath()+"/output_compression_16_08_2016"+GetMatrixName();
	ofstream output(filename.c_str());
	if (!output){
		cerr<<"Output file cannot be created"<<endl;
		exit(1);
	}
	output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur_MvProd "<<"Erreur_Frob"<<endl;
	
	for(int iepsilon=0; iepsilon<nepsilon; iepsilon++)
	{
		cout << "iepsilon: " << iepsilon << endl;
		
		for(int ieta=0; ieta<neta; ieta++)
		{
			cout << "ieta: " << ieta << endl;
			SetEta(eta[ieta]);
			SetEpsilon(epsilon[iepsilon]);
			cout<<"Eta : "+NbrToStr(GetEta())<<endl;
			cout<<"Epsilon : "+NbrToStr(GetEpsilon())<<endl;
			
			vector<double> times;
			
			tic();
			// Build the hierarchical matrix with compressed and dense blocks
			HMatrix B(A,x,r,tab,GetEpsilon()==-1 ? 0 : -1);
			// if epsilon=-1 rank 0 blocks, otherwise aca compression with the given precision
			toc(times);
			
			vectCplx ub(nr);
			tic();
			MvProd(ub,B,u); // Do the matrix vector product
			toc(times);
			
			Real errH = norm(ua-ub)/norm(ua);
			cout << "Matrix-vector product relative error with HMatrix: " << errH << endl;
			
			Real compressionH=CompressionRate(B);
			cout << "Compression rate with HMatrix:" << compressionH << endl;
			
			Real froberrH = sqrt(squared_absolute_error(B,A))/normA;
			cout << "Relative error in Frobenius norm with HMatrix: " << froberrH << endl;
			
			// write in output file
			output<<GetEta()<<" "<<GetEpsilon()<<" "<<compressionH<<" "<<errH<<" "<<froberrH<<endl;
			//cout<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<errH<<" "<<froberrH<<endl;
		}
	}
	output.close();
    
    
    ///// With the sparse matrix obtained with the heuristic strategy:
    
    vectCplx uasp(nr);
    tic();
    MvProd(uasp,spA,u);
    toc();
    
    Real errSp = norm(ua-uasp)/norm(ua);
    cout << "Matrix-vector product relative error with Sparse matrix: " << errSp << endl;
    
    Real compressionSp=CompressionRate(spA);
    cout << "Compression rate with Sparse matrix:" << compressionSp << endl;
	
	for(int i=0; i<nb_coeff(spA); i++){
		A(spA.I_(i),spA.J_(i)) = 0; // now A = A-spA !! (to save memory)
	}
	
	Real froberrSp = NormFrob(A)/normA;
	cout << "Relative error in Frobenius norm with Sparse Matrix: " << froberrSp << endl;
    
    
    // for output file of sparse matrix
    string filename2=GetOutputPath()+"/output_sparseMatrix_16_08_2016"+GetMatrixName();
    ofstream output2(filename2.c_str());
    if (!output2){
        cerr<<"Output file cannot be created"<<endl;
        exit(1);
    }
    output2<<"Compression "<<"Erreur_MvProd "<<"Erreur_Frob"<<endl;
    output<<compressionSp<<" "<<errSp<<" "<<froberrSp<<endl;

	
}
