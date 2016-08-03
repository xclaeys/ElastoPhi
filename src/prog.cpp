#include <iostream>
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
 cout<<"##################################"<<endl;
 
 ////////////////========================================================////////////////

 vectReal r;
 vectR3   x;
 Matrix   A;
 
 LoadMatrix( "../data/matrice7.txt",A);
 LoadPoints( "../data/maillage7.txt",x,r);
 ExportMEDIT("../data/maillage7.txt");
 
 HMatrix B(A,x,r,x,r); 

 // Vecteur (pseudo-)aleatoire
 int nr  = nb_rows(A);
 vectCplx u(nr);
 int NbSpl = 1000; 
 double du = 5./double(NbSpl);
 srand (time(NULL));  
 for(int j=0; j<nr; j++){
   int n = rand()%(NbSpl+1);
   u[j] = n*du;}
 
 vectCplx ua(nr),ub(nr);
 MvProd(ua,A,u);
 MvProd(ub,B,u);  
 Real err = norm(ua-ub)/norm(ua);
 cout << "Erreur:\t" << err << endl;
 
 cout << "Taux de compression:\t";
 cout << CompressionRate(B) << endl;  
  
}
