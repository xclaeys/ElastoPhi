#ifndef LOADING_HPP
#define LOADING_HPP

#include <fstream>
#include "matrix.hpp"

using namespace std;

//==================================================//
//
//  INPUT:
//  filename: nom du fishier de maillage
//
//  OUTPUT:
//  x: nuage de points (centre des elements)
//  r: rayon de champ proche associe a chaque point
//
//==================================================//

void load(const char* filename, const vectR3& x, const vectReal& r){
  
  x.clear(); r.clear();
  ifstream file; file.open(filename);
  
  int NbElt, NbPt, poubelle;
  R3 Pt[4]; R3 Ctr; Real Rmax,Rad;
  
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
      if(Rad>Rmax){Rmax=Rad;}
    }
    
    for(int j=0; j<3; j++){
      r.push_back(Rmax);}    
    
    if(e<NbElt-1){file >> poubelle;}
  }
  
  
  for(int j=0; j<4; j++){
    cout << "Pt["<< j << "]:\t" << Pt[j] << endl;}
  


}





#endif
