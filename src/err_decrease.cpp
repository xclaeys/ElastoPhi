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
	
	
	//	////////////////========================================================////////////////
	//	////////////////////////////////========  Input ========////////////////////////////////
	//
	//	// Check the number of parameters
	//	if (argc < 2) {
	//		// Tell the user how to run the program
	//		cerr << "Usage: " << argv[0] << " input name" << endl;
	//		/* "Usage messages" are a conventional way of telling the user
	//		 * how to run a program if they enter the command incorrectly.
	//		 */
	//		return 1;
	//	}
	//
	//	// Load the inputs
	//	string inputname = argv[1];
	//	Param Parametres(inputname);
	//
	//	cout<<"############# Inputs #############"<<endl;
	//	cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
	//	cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
	//	cout<<"Data path : "+Parametres.datapath<<endl;
	//	cout<<"Output path : "+Parametres.outputpath<<endl;
	//	cout<<"Mesh name : "+Parametres.meshname<<endl;
	//	cout<<"Matrix name : "+Parametres.matrixname<<endl;
	//	cout<<"##################################"<<endl;
	
	
	
	
	srand (1);
	// we set a constant seed for rand because we want always the same result if we run the check many times
	// (two different initializations with the same seed will generate the same succession of results in the subsequent calls to rand)
	
	// Build matrix A with property for ACA
	int nr = 100;
	vectInt Ir(nr); // row indices for the lrmatrix
	vectInt Ic(nr); // column indices for the lrmatrix
	// p1: points in a unit disk of the plane z=z1
	double z1 = 1;
	vectR3 p1(nr);
	for(int j=0; j<nr; j++){
		Ir[j] = j;
		Ic[j] = j;
		double rho = ((double) rand() / (double)(RAND_MAX)); // (double) otherwise integer division!
		double theta = ((double) rand() / (double)(RAND_MAX));
		p1[j][0] = sqrt(rho)*cos(2*M_PI*theta); p1[j][1] = sqrt(rho)*sin(2*M_PI*theta); p1[j][2] = z1;
		// sqrt(rho) otherwise the points would be concentrated in the center of the disk
	}
	
	
	// Parametres
	// Load the inputs
	string inputname = argv[1];
	Param Parametre(inputname);
	Param Parametres(1e-1,1e-1); // eta (only for hmatrix), ACA epsilon
	string filename=Parametres.outputpath+"/output_err_decrease.txt";
	ofstream output(filename);
//	output.open(filename,ios::app);
	if (!output){
		cerr<<"Output file cannot be created"<<endl;
		exit(1);
	}

	
	
	for (int j=1;j<6;j++){
		// p2: points in a unit disk of the plane z=z2
		
		double z2 = 1+double(j)*0.2;
		vectR3 p2(nr);
		for(int j=0; j<nr; j++){
			double rho = ((double) rand() / (RAND_MAX)); // (double) otherwise integer division!
			double theta = ((double) rand() / (RAND_MAX));
			p2[j][0] = sqrt(rho)*cos(2*M_PI*theta); p2[j][1] = sqrt(rho)*sin(2*M_PI*theta); p2[j][2] = z2;
		}
		Matrix A(nr,nr);
		for(int j=0; j<nr; j++){
			for(int k=0; k<nr; k++){
				A(j,k) = 1./(4*M_PI*norm(p1[j]-p2[k]));
			}
		}
		
		SubMatrix subm(A,Ir,Ic); // A viewed as a SubMatrix
		
		for (int i=1;i<50;i++){
			
			tic();
			LowRankMatrix lrm(subm,Ir,Ic,i);
			toc();
			tic();
			LowRankMatrixSVD lrmsvd(subm,Ir,Ic,i); // construct a low rank matrix B applying ACA to matrix A
			toc();
			tic();
			Real err=0;
			Real err_SVD=0;
			err=squared_relative_error(lrm,subm);
			err_SVD=squared_relative_error(lrmsvd,subm);
//			cout<<z2-z1<<" "<<i<<" "<<err<<endl;
			output <<z2-z1<<" "<<i<<" "<<err<<" "<<err_SVD<<endl;
			
			toc();


			//    ofstream output;
			//    output.open(filename,ios::app);
			//    if (!output){
			//        cerr<<"Output file cannot be created"<<endl;
			//        exit(1);
			//    }
			//    else{
			//        if (!infile.good()){
			//            output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur"<<endl;
			//        }
			//        else{
			//
			//        }
			//        output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<endl;
			//    }
			//    output.close();
		}
	}
	output.close();
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Build Hmatrix 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
//	vector<double> times;
	//	vectReal r;
	//	vectR3   x;
	//	Matrix   A;
	//	tic();
	//	LoadMatrix((Parametres.datapath+"/"+Parametres.matrixname).c_str(),A);
	//	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
	//	toc(times);
	//	tic();
	//	HMatrix B(A,x,r);
	//	toc(times);
	
	
	////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////    Test MvProd 	////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	//	tic();
	//	Real err=0;
	//	LowRankMatrix lrm = GetLowRankMatrix(B, 0);
	//	SubMatrix    subm = SubMatrix(A,ir_(lrm),ic_(lrm));
	//	squared_absolute_error(err,lrm,subm);
	//	toc();
	
	//	// Vecteur
	//	int nr  = nb_rows(A);
	//	vectCplx u(nr);
	//	int NbSpl = 1000;
	//	double du = 5./double(NbSpl);
	//	srand (1);
	//	for(int j=0; j<nr; j++){
	//		int n = rand()%(NbSpl+1);
	//		u[j] = n*du;}
	//
	//	vectCplx ua(nr),ub(nr);
	//	MvProd(ua,A,u);
	//	MvProd(ub,B,u);
	//	Real err = norm(ua-ub)/norm(ua);
	//	Real compression=CompressionRate(B);
	//
	//	cout<<"Erreur : "<<err<<endl;
	//	cout<<"Compression :"<<compression<<endl;
	
	//    // Ecriture dans un fichier avec append:
	//    ////////////////////////////////////////////////////////////////////////////////////////
	//    ////////////////////////////////    Fichier de sortie 	////////////////////////////////
	//    ////////////////////////////////////////////////////////////////////////////////////////
	//    string filename=Parametres.outputpath+"/output_compression_"+Parametres.matrixname;
	//    ifstream infile(filename);
	//    ofstream output;
	//    output.open(filename,ios::app);
	//    if (!output){
	//        cerr<<"Output file cannot be created"<<endl;
	//        exit(1);
	//    }
	//    else{
	//        if (!infile.good()){
	//            output<< "Eta "<<"Epsilon "<<"Compression "<<"Erreur"<<endl;
	//        }
	//        else{
	//
	//        }
	//        output<<Parametres.eta<<" "<<Parametres.epsilon<<" "<<compression<<" "<<err<<endl;
	//    }
	//    output.close();
	
}
