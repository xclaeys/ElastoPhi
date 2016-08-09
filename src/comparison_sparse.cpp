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
#include "sparsematrix.hpp"
#include "loading.hpp"
#include "export.hpp"
#include "user.hpp"


#include <stdio.h>
#include <stdlib.h>
#include <time.h>


using namespace std;

/**************************************************************************//**
 * 
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
	
//	vectReal r;
//	vectR3   x;
	Matrix   A;
    SpMatrix spA;
    
    tic();
    LoadSpMatrix((Parametres.datapath+"/"+"matrice5364FracsTrianglesCreuse.txt").c_str(),spA);
    toc();
    
    tic();
	LoadMatrix((Parametres.datapath+"/"+"matrice5364FracsTriangles.txt").c_str(),A);
    toc();
//    tic();
//	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
//    toc(times2);
    
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
    MvProd(uasp,spA,u);
    
    Real errsp = norm(ua-uasp)/norm(ua);
    cout << "Erreur Sparse: " << errsp << endl;
    
    tic();
    cout << "Norme de Frobenius de A: " << NormFrob(A) << endl;
    toc();
    
    for(int i=0; i<nb_coeff(spA); i++){
        A(spA.I_(i),spA.J_(i)) = 0; // now A = A-spA !!
    }
    
    cout << "Norme de Frobenius de A-spA: " << NormFrob(A) << endl;
    
    
//    // Values of eta and epsilon
//    int neta = 4;
//    double eta[neta];
//    eta[0] = 1e+1; eta[1] = 1e+0; eta[2] = 1e-1; eta[3] = 1e-2;
//    //eta[0] = 1e+1; eta[1] = 8e-1; eta[2] = 1e-1; eta[3] = 1e-2;
//    int nepsilon = 4;
//    double epsilon[nepsilon];
//    epsilon[0] = 1e+0; epsilon[1] = 5e-1; epsilon[2] = 1e-1; epsilon[3] = 1e-2;
//    
//    // for output file
//    string filename=Parametres.outputpath+"/output_compression_"+Parametres.matrixname;
//    //string filename=Parametres.outputpath+"/output_compression_CondAdmMin08_"+Parametres.matrixname;
//    ofstream output(filename.c_str());
//    if (!output){
//        cerr<<"Output file cannot be created"<<endl;
//        exit(1);
//    }
//    //output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur "<<"TimeHMatrix "<<"TimeMvProd"<<endl;
//    
//    for(int iepsilon=0; iepsilon<nepsilon; iepsilon++)
//    {
//        cout << "iepsilon: " << iepsilon << endl;
//        
//        for(int ieta=0; ieta<neta; ieta++)
//        {
//            cout << "ieta: " << ieta << endl;
//            Param Parametre(eta[ieta],epsilon[iepsilon]);
//            
//            vector<double> times;
//            
//            tic();
//            HMatrix B(A,x,r); // Build the hierarchical matrix with compressed and dense blocks
//            toc(times);
//            
//            vectCplx ub(nr);
//            tic();
//            MvProd(ub,B,u); // Do the matrix vector product
//            toc(times);
//            
//            Real err = norm(ua-ub)/norm(ua);
//            Real compression=CompressionRate(B);
//            
//            // write in output file
//            output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<" "<<times[0]<<" "<<times[1]<<endl;
//            cout<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<" "<<times[0]<<" "<<times[1]<<endl;
//        }
//    }
//	output.close();
	
}
