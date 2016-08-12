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

/**************************************************************************//**
* It converts the mesh file from the format given by Ibtihel to gmsh format (.msh)
* and to medit format (.mesh) for visualization of the mesh
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
	Param Parametres(inputname);
 
	cout<<"############# Inputs #############"<<endl;
	cout<<"Eta : "+NbrToStr(Parametres.eta)<<endl;
	cout<<"Epsilon : "+NbrToStr(Parametres.epsilon)<<endl;
	cout<<"Data path : "+Parametres.datapath<<endl;
	cout<<"Output path : "+Parametres.outputpath<<endl;
	cout<<"Mesh name : "+Parametres.meshname<<endl;
	cout<<"Matrix name : "+Parametres.matrixname<<endl;
	cout<<"##################################"<<endl;
 
	////////////////========================================================////////////////
	
	vectReal r;
	vectR3   x;
	tic();
	LoadPoints((Parametres.datapath+"/"+Parametres.meshname).c_str(),x,r);
	vectInt tab(Parametres.ndofperelt*x.size());
	for (int j=0;j<x.size();j++){
		tab[3*j]  = j;
		tab[3*j+1]= j;
		tab[3*j+2]= j;
	}
	toc();
	tic();
	Cluster t(x,r,tab);
	toc();
	
	for(int idepth=1; idepth<4; idepth++){
		VisuPartitionedMesh(t, Parametres.datapath+"/"+Parametres.meshname, "VisuPart"+(split(Parametres.meshname,'.')).at(0)+"depth"+NbrToStr(idepth)+".msh", idepth);
	}
	
	
}
