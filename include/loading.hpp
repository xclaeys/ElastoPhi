#ifndef LOADING_HPP
#define LOADING_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "matrix.hpp"
#include "user.hpp"
#include "sparsematrix.hpp"

using namespace std;
//==================================================//
//
//  DESCRIPTION:
//
//
//
//  INPUT:
//  filename: nom du fichier d'input
//
//  OUTPUT:
//
//
//
//==================================================//

class Param{
public:
	static int  dim;
	static Real eta;
	static Real epsilon;
	static string datapath;
	static string outputpath;
	static string meshname;
	static string matrixname;
	
	Param (string inputname); // Lecture du fichier
	Param (); // Constructeur par defaut
	Param (Real , Real, string , string, string , string); // Valeurs données à la main
	Param (Real , Real); // Valeurs données à la main quand on construit la matrice
	Param (int); // Nombre d'inconnus par points géométriques
};

// Allocation de la mémoire pour les valeurs statiques (obligatoire)
Real Param::eta;
Real Param::epsilon;
int Param::dim;
string Param::datapath;
string Param::outputpath;
string Param::meshname;
string Param::matrixname;

Param::Param (){
	
}
Param::Param(string inputname){
	ifstream data(inputname.c_str());
	
	// Si le fichier n'existe pas
	if (!data){
		cerr << "Input file doesn't exist" << endl;
		exit(1);
	}
	// Lecture du fichier
	else {
		while (data){
			string strInput;
			getline(data,strInput);
			
			vector<string> line = split (strInput,' ');
			if (!line.empty()){
				if (line.at(0)=="Eta"){
					eta=StrToNbr<Real>(line.back());
				}
				else if (line.at(0)=="Epsilon"){
					epsilon=StrToNbr<Real>(line.back());
				}
				else if (line.at(0)=="Data_path"){
					datapath=line.back();
				}
				else if (line.at(0)=="Output_path"){
					outputpath=line.back();
				}
				else if (line.at(0)=="Mesh_name"){
					meshname=line.back();
				}
				else if (line.at(0)=="Matrix_name"){
					matrixname=line.back();
				}
			}
		}
		ifstream Matrixname(datapath+"/"+matrixname);
		if (!Matrixname){
			cerr << "Matrix file does not exist" << endl;
			exit(1);
		}
		ifstream Meshname(datapath+"/"+meshname);
		if (!Meshname){
			cerr << "Mesh file does not exist" << endl;
			exit(1);
		}
	}
}
Param::Param(Real eta0, Real epsilon0){
	eta=eta0;
	epsilon=epsilon0;
}
Param::Param(Real eta0, Real epsilon0, string datapath0, string outputpath0, string meshname0, string matrixname0){
	eta=eta0;
	epsilon=epsilon0;
	datapath=datapath0;
	outputpath=outputpath0;
	meshname=meshname0;
	matrixname=matrixname0;
	
}
Param::Param(int dim0){
	dim=dim0;
}

Param ParametreDefaut(1.0,1.e-2,"../data","../output","maillage1.txt","matrice1.txt");
Param ParametreDefautDim(3);

//==================================================//
//
//  DESCRIPTION:
//  Charge la matrice dense
//
//  INPUT:
//  filename: nom du fichier de la matrice de Ibtihel
//
//  OUTPUT:
//  m: matrice dense
//
//==================================================//

void LoadMatrix(const char* filename, Matrix& m){
	
	int NbRow, NbCol;
	string      line;
	int        j0,k0;
	Cplx         val;
	Param Parametres;
	
	// Ouverture fichier
	ifstream file; file.open(filename);
	if(!file.good()){
		cout << "LoadMatrix in loading.hpp: error opening the matrix file" << endl;
		abort();}
	
	// Lecture nombre de lignes et de colonnes
	file >> NbRow; file >> NbCol;
	m.resize(NbRow,NbCol);
	
	getline(file,line);
	getline(file,line);
	while(!file.eof()){
		
		// Lecture de la ligne
		istringstream iss(line);
		
		// Pour chaque ligne, stockage
		// du bloc d'interaction
		iss >> j0; j0 = Parametres.dim*(j0-1);
		iss >> k0; k0 = Parametres.dim*(k0-1);
		
		for(int j=0; j<Parametres.dim; j++){
			for(int k=0; k<Parametres.dim; k++){
				iss >> val;
				m(j0+j,k0+k) = val;
			}
		}
		getline(file,line);
	}
	
	file.close();
}

//==================================================//
//
//  DESCRIPTION:
//  Charge la matrice sparse
//
//  INPUT:
//  filename: nom du fichier de la matrice sparse de Ibtihel
//
//  OUTPUT:
//  m: matrice sparse
//
//==================================================//

void LoadSpMatrix(const char* filename, SpMatrix& m){
	
	int NbRow, NbCol;
	string      line;
	int        j0,k0;
	Cplx         val;
	Param Parametres;
	
	// Ouverture fichier
	ifstream file; file.open(filename);
	if(!file.good()){
		cout << "LoadMatrix in loading.hpp: error opening the matrix file" << endl;
		abort();}
	
	// Lecture nombre de lignes et de colonnes
	file >> NbRow; file >> NbCol;
	m.resize(NbRow,NbCol,NbRow*max(NbCol/10,10));
	
	int NbCoef=0;
	
	getline(file,line);
	getline(file,line);
	while(!file.eof()){
		
		// Lecture de la ligne
		istringstream iss(line);
		
		// Pour chaque ligne, stockage
		// du bloc d'interaction
		iss >> j0; j0 = Parametres.dim*(j0-1);
		iss >> k0; k0 = Parametres.dim*(k0-1);
		
		for(int j=0; j<Parametres.dim; j++){
			for(int k=0; k<Parametres.dim; k++){
				iss >> val;
				m.I_(NbCoef) = j0+j;
				m.J_(NbCoef) = k0+k;
				m.K_(NbCoef) = val;
				//m(j0+j,k0+k) = val;
				NbCoef++;
			}
		}
		getline(file,line);
	}
	
	file.close();
	
	m.resize(NbRow,NbCol,NbCoef);
}


//==================================================//
//
//  DESCRIPTION:
//  Charge les donnees geometriques
//  associees au nuage de points
//
//  INPUT:
//  filename: nom du fichier de maillage
//
//  OUTPUT:
//  x: nuage de points (centre des elements)
//  r: rayon de champ proche associe a chaque point
//
//==================================================//

void LoadPoints(const char* filename, vectR3& x, vectReal& r){
	
	x.clear(); r.clear();
	int NbElt, NbPt, poubelle;
	R3 Pt[4]; R3 Ctr; Real Rmax,Rad;
	Param Parametres;
	
	// Ouverture fichier
	ifstream file; file.open(filename);
	if(!file.good()){
		cout << "LoadPoints in loading.hpp: " ;
		cout << "error opening the geometry file\n";
		abort();}
	
	// Nombre d'elements
	file >> NbElt;
	
	// Lecture elements
	for(int e=0; e<NbElt; e++){
		Ctr=0.; file >> poubelle;
		file >> NbPt;
		
		// Calcul centre element
		for(int j=0; j<NbPt; j++){
			file >> poubelle; file>>Pt[j];
			Ctr+= (1./Real(NbPt))*Pt[j];}
		
		//Ajout de trois points identiques
		//dans le nuage de points...
		for(int j=0; j<Parametres.dim; j++){
			x.push_back(Ctr);}
		
		// Calcul du rayon champ
		// proche associe a l'element
		Rmax = norm(Ctr-Pt[0]);
		
		for(int j=1; j<NbPt; j++){
			Rad = norm(Ctr-Pt[j]);
			if(Rad>Rmax){Rmax=Rad;}}
		
		for(int j=0; j<Parametres.dim; j++){
			r.push_back(Rmax);}
		
		// Separateur inter-element
		if(e<NbElt-1){file >> poubelle;}
	}
	
	// Fermeture fichier
	file.close();
	
}



#endif
