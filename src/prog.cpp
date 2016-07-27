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

	Real eps;
	Real eta;
	
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
	GetInput(inputname,eps,eta);
	
	cout<<"############# Inputs #############"<<endl;
	cout<<"Eta : "+NbrToStr(eta)<<endl;
	cout<<"Epsilon : "+NbrToStr(eps)<<endl;
	cout<<"##################################"<<endl;
	
////////////////========================================================////////////////
	
}
