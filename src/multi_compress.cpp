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
	cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
	cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
	cout<<"Data path : "+Parametres.datapath<<endl;
	cout<<"Output path : "+Parametres.outputpath<<endl;
	cout<<"Mesh name : "+Parametres.meshname<<endl;
	cout<<"Matrix name : "+Parametres.matrixname<<endl;
	cout<<"##################################"<<endl;
 
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Build Hmatrix 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	
	vectReal r;
	vectR3   x;
	Matrix   A;
	
	LoadMatrix((Parametres.datapath+"/"+Parametres.matrixname).c_str(),A);
	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
	
	HMatrix B(A,x,r,x,r);
 
	
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Test MvProd 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	
	// Vecteur (pseudo-)aleatoire
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
	
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Fichier de sortie 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	string filename=Parametres.outputpath+"/output_compression_"+Parametres.matrixname;
	ifstream infile(filename);
	ofstream output;
	output.open(filename,ios::app);
	if (!output){
		cerr<<"Output file cannot be created"<<endl;
		exit(1);
	}
	else{
		if (!infile.good()){
			output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur"<<endl;
		}
		else{
			
		}
		output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<endl;
	}
	output.close();
	
}
