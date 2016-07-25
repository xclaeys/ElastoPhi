#ifndef LOADING_HPP
#define LOADING_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "matrix.hpp"

using namespace std;

//==================================================//
//
//  DESCRIPTION:
//  Charge la matrice 
//
//  INPUT:
//  filename: nom du fichier de maillage
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
    iss >> j0; j0 = 3*(j0-1);
    iss >> k0; k0 = 3*(k0-1);
    
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
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
    //dans le nuage de points
    for(int j=0; j<3; j++){
      x.push_back(Ctr);}

    // Calcul du rayon champ 
    // proche associe a l'element
    Rmax = norm(Ctr-Pt[0]);
    
    for(int j=1; j<NbPt; j++){
      Rad = norm(Ctr-Pt[j]);
      if(Rad>Rmax){Rmax=Rad;}}
    
    for(int j=0; j<3; j++){
      r.push_back(Rmax);}    

    // Separateur inter-element
    if(e<NbElt-1){file >> poubelle;}
  }

  // Fermeture fichier
  file.close();

}


#endif
