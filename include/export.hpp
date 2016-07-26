#ifndef EXPORT_HPP
#define EXPORT_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include "matrix.hpp"

//==================================================//
//
//  DESCRIPTION:
//  Convertit le maillage au format gmsh
//  pour visualisation. Le fichier de sortie
//  s'appelle "visu.msh"
//
//  INPUT:
//  filename: nom du fichier de maillage
//
//  OUTPUT:
//  none
//
//==================================================//

void ExportGMSH(const char* filename){

  vector<R3x4> Elt;
  int NbElt, NbPt, poubelle;
  
  // Ouverture fichier
  ifstream infile;
  infile.open(filename);
  if(!infile.good()){
    cout << "LoadPoints in loading.hpp: error opening the geometry file" << endl;
    abort();}
  
  // Nombre d'elements
  infile >> NbElt; 
  Elt.resize(NbElt);
  
  // Lecture elements
  for(int e=0; e<NbElt; e++){
    infile >> poubelle;
    infile >> NbPt;
    
    // Calcul centre element
    for(int j=0; j<NbPt; j++){
      infile >> poubelle;
      infile >> Elt[e][j];
    }

    // Separateur inter-element
    if(e<NbElt-1){infile >> poubelle;}    
    
  }
  infile.close();  

  // Ecriture fichier de sortie
  ofstream outfile;
  outfile.open("visu.msh");
  outfile << "$MeshFormat\n"; 
  outfile << "2.2 0 8\n";
  outfile << "$EndMeshFormat\n";
  outfile << "$Nodes\n";  
  outfile << 4*NbElt << endl;
  for(int j=0; j<NbElt; j++){
    for(int k=0; k<4; k++){
      outfile << 1+4*j+k << "\t";
      outfile << Elt[j][k] << "\n";
    }
  }
  outfile << "$EndNodes\n";    
  outfile << "$Elements\n";      
  outfile << NbElt << endl;
  for(int j=0; j<NbElt; j++){
    outfile << j  << "\t";
    outfile << 3  << "\t";
    outfile << 2  << "\t";
    outfile << 99 << "\t";
    outfile << 1  << "\t";
    for(int k=0; k<4; k++){
      outfile << 1+4*j+k << "\t";}
    outfile << "\n";
  }
  outfile << "$EndElements\n";        
  
}

#endif
