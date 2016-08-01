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
	//cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
	//cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
	cout<<"Data path : "+Parametres.datapath<<endl;
	cout<<"Output path : "+Parametres.outputpath<<endl;
	cout<<"Mesh name : "+Parametres.meshname<<endl;
	cout<<"Matrix name : "+Parametres.matrixname<<endl;
	cout<<"##################################"<<endl;
    
	
	vectReal r;
	vectR3   x;
	Matrix   A;
	LoadMatrix((Parametres.datapath+"/"+Parametres.matrixname).c_str(),A);
	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
    
    // Vecteur pour le produit matrice vecteur
    int nr  = nb_rows(A);
    vectCplx u(nr);
    int NbSpl = 1000;
    double du = 5./double(NbSpl);
    srand (1); //!!
    for(int j=0; j<nr; j++){
        int n = rand()%(NbSpl+1);
        u[j] = n*du;}
    
    vectCplx ua(nr);
    MvProd(ua,A,u);
    
    // Values of eta and epsilon
    int neta = 4;
    double eta[neta];
    eta[0] = 1e+1; eta[1] = 1e+0; eta[2] = 1e-1; eta[3] = 1e-2;
    int nepsilon = 4;
    double epsilon[nepsilon];
    epsilon[0] = 1e+1; epsilon[1] = 1e+0; epsilon[2] = 1e-1; epsilon[3] = 1e-2;
    
    // for output file
    string filename=Parametres.outputpath+"/output_compression_"+Parametres.matrixname;
    ofstream output(filename.c_str());
    if (!output){
        cerr<<"Output file cannot be created"<<endl;
        exit(1);
    }
    //output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur "<<"TimeHMatrix "<<"TimeMvProd"<<endl;
    
    for(int iepsilon=0; iepsilon<nepsilon; iepsilon++)
    {
        cout << "iepsilon: " << iepsilon << endl;
        
        for(int ieta=0; ieta<neta; ieta++)
        {
            cout << "ieta: " << ieta << endl;
            Param Parametre(eta[ieta],epsilon[iepsilon]);
            
            vector<double> times;
            
            tic();
            HMatrix B(A,x,r,x,r); // Build the hierarchical matrix with compressed and dense blocks
            toc(times);
            
            vectCplx ub(nr);
            tic();
            MvProd(ub,B,u); // Do the matrix vector product
            toc(times);
            
            Real err = norm(ua-ub)/norm(ua);
            Real compression=CompressionRate(B);
            
            // write in output file
            output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<" "<<times[0]<<" "<<times[1]<<endl;
        }
    }
	output.close();
	
}
