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
	Param Parametres(inputname);
 
	cout<<"############# Inputs #############"<<endl;
	//cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
	//cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
	cout<<"Data path : "+Parametres.datapath<<endl;
	cout<<"Output path : "+Parametres.outputpath<<endl;
	cout<<"Mesh name : "+Parametres.meshname<<endl;
	cout<<"Matrix name : "+Parametres.matrixname<<endl;
	cout<<"##################################"<<endl;
	
	vector<double> times2;
	
	vectReal r;
	vectR3   x;
	Matrix   A;
	SpMatrix spA;
	
	tic();
	LoadSpMatrix((Parametres.datapath+"/"+(split(Parametres.matrixname,'.')).at(0)+"Creuse.txt").c_str(),spA);
	toc();
	
	tic();
	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
	toc();
	tic();
	LoadMatrix((Parametres.datapath+"/"+Parametres.matrixname).c_str(),A);
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
	
	vectCplx uasp(nr);
	tic();
	MvProd(uasp,spA,u);
	toc();
	
	Real errSp = norm(ua-uasp)/norm(ua);
	cout << "Matrix-vector product relative error with Sparse matrix: " << errSp << endl;
	
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
	string filename=Parametres.outputpath+"/output_compression_12_08_2016"+Parametres.matrixname;
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
			Param Parametre(eta[ieta],epsilon[iepsilon]);
			cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
			cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
			
			vector<double> times;
			
			tic();
            // Build the hierarchical matrix with compressed and dense blocks
			HMatrix B(A,x,r,tab,Parametres.epsilon==-1 ? 0 : -1);
            // if epsilon=-1 rank 0 blocks, otherwise aca compression with the given precision
			toc(times);			
			
			vectCplx ub(nr);
			tic();
			MvProd(ub,B,u); // Do the matrix vector product
			toc(times);
			
			Real errH = norm(ua-ub)/norm(ua);
			cout << "Matrix-vector product relative error with HMatrix: " << errH << endl;
			
			Real compression=CompressionRate(B);
			cout << "Compression rate :" << compression << endl;
			
			Real froberrH = sqrt(squared_absolute_error(B,A))/normA;
			cout << "Relative error in Frobenius norm with HMatrix: " << froberrH << endl;
			
			// write in output file
			output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<errH<<" "<<froberrH<<endl;
			//cout<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<errH<<" "<<froberrH<<endl;
		}
	}
	output.close();
	
	
	for(int i=0; i<nb_coeff(spA); i++){
		A(spA.I_(i),spA.J_(i)) = 0; // now A = A-spA !! (to save memory)
	}
	
	Real froberrSp = NormFrob(A)/normA;
	cout << "Relative error in Frobenius norm with Sparse Matrix: " << froberrSp << endl;
	
}
