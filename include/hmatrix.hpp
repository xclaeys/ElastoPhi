#ifndef HMATRIX_HPP
#define HMATRIX_HPP

#include <cassert>
#include <fstream>
#include "matrix.hpp"

//===============================//
//     MATRICE HIERARCHIQUE      // 
//===============================//

class HMatrix{
  
private:
  
  const Matrix& mat;
  const vectR3& xt;
  const vectR3& xs;
  
  vector<Block>         FarField;
  vector<LowRankMatrix> FarFieldMat;
    
  vector<Block>         NearField;
  vector<Matrix>        NearFieldMat;
  
  void BuildBlockTree(const Cluster&, const Cluster&);
  
public:
  
  HMatrix(const Matrix&, const vectR3&, const vectReal&, const vectR3&, const vectReal&); 
  friend void DisplayPartition(const HMatrix&, char const* const);
  friend void MvProd(vectCplx&, const HMatrix&, const vectCplx&);
  friend Real CompressionRate(const HMatrix&);
  
};

void HMatrix::BuildBlockTree(const Cluster& t, const Cluster& s){
  Block B(t,s);
  if( B.IsAdmissible() ){
    FarField.push_back(B);
    const vectInt& I = num_(t);
    const vectInt& J = num_(s);
    SubMatrix submat = SubMatrix(mat,I,J);
    FarFieldMat.push_back(LowRankMatrix(submat));
  }
  else if( s.IsLeaf() ){
    if( t.IsLeaf() ){
      NearField.push_back(B);
      const vectInt& I = num_(t);
      const vectInt& J = num_(s);
      NearFieldMat.push_back(SubMatrix(mat,I,J));      
    }    
    else{
      BuildBlockTree(son_(t,0),s);
      BuildBlockTree(son_(t,1),s);
    }
  }
  else{    
    if( t.IsLeaf() ){
      BuildBlockTree(t,son_(s,0));
      BuildBlockTree(t,son_(s,1));
    }
    else{
      BuildBlockTree(son_(t,0),son_(s,0));
      BuildBlockTree(son_(t,0),son_(s,1));
      BuildBlockTree(son_(t,1),son_(s,0));
      BuildBlockTree(son_(t,1),son_(s,1));    
    }
  }
  
}

HMatrix::HMatrix(const Matrix& mat0,
		 const vectR3& xt0, const vectReal& rt,
		 const vectR3& xs0, const vectReal& rs):
  
  mat(mat0), xt(xt0), xs(xs0) {
  assert( nb_rows(mat)==xt.size() && nb_cols(mat)==xs.size() );
  
  // Construction arbre des paquets
  Cluster t(xt,rt); Cluster s(xs,rs);
  
  // Construction arbre des blocs
  BuildBlockTree(t,s);  
  
}


Real CompressionRate(const HMatrix& hmat){
  
  Real comp = 0.;
  Real size = ( (hmat.xt).size() )*( (hmat.xs).size() );
  const vector<LowRankMatrix>& FarFieldMat  = hmat.FarFieldMat;
  const vector<Matrix>&        NearFieldMat = hmat.NearFieldMat;
  
  for(int j=0; j<FarFieldMat.size(); j++){
    Real nr   = nb_rows(FarFieldMat[j]);
    Real nc   = nb_cols(FarFieldMat[j]);
    Real rank = rank_of(FarFieldMat[j]);
    comp += rank*(nr + nc)/size;
  }
  
  for(int j=0; j<NearFieldMat.size(); j++){
    Real nr   = nb_rows(NearFieldMat[j]);
    Real nc   = nb_cols(NearFieldMat[j]);
    comp += nr*nc/size;
  }

  return comp;
  
}


// Representation graphique de la partition en bloc
void DisplayPartition(const HMatrix& hmat, char const * const name){

  const vector<Block>& FarField = hmat.FarField;
  const vectR3& xt = hmat.xt;
  const vectR3& xs = hmat.xs;;
  
  // Representation graphique
  const int  Ns = xs.size();
  const Real ds = 1./Real(Ns-1);
  const int  Nt = xt.size();
  const Real dt = 1./Real(Nt-1);
  
  ofstream file; file.open(name);
  for(int j=0; j<FarField.size(); j++){
    
    const Cluster& t = tgt_(FarField[j]);
    const vectInt& It = num_(t);
    Real at = (It[0]-0.5)*dt;
    Real bt = (It[It.size()-1]+0.5)*dt;
    
    const Cluster& s = src_(FarField[j]); 
    const vectInt& Is = num_(s);       
    Real as = (Is[0]-0.5)*ds;
    Real bs = (Is[Is.size()-1]+0.5)*ds;    

    file << as << "\t" << at << "\n";
    file << bs << "\t" << at << "\n";
    file << bs << "\t" << bt << "\n";    
    file << as << "\t" << bt << "\n";
    file << as << "\t" << at << "\n";
    file << endl;
    
  }
  file.close();
 
}


void MvProd(vectCplx& f, const HMatrix& A, const vectCplx& x){
  assert(size(f)==size(x)); fill(f,0.);
  
  const vector<Block>&         FarField     = A.FarField;
  const vector<Block>&         NearField    = A.NearField;
  const vector<LowRankMatrix>& FarFieldMat  = A.FarFieldMat;
  const vector<Matrix>&        NearFieldMat = A.NearFieldMat;

  // Contribution champ lointain
  for(int b=0; b<FarField.size(); b++){
    const Block&          B  = FarField[b];
    const LowRankMatrix&  M  = FarFieldMat[b];
    const vectInt&        It = num_(tgt_(B));
    const vectInt&        Is = num_(src_(B));      
    
    ConstSubVectCplx xx(x,Is); 
    SubVectCplx ff(f,It); 
    MvProd(ff,M,xx);
  }
  
  // Contribution champ proche
  for(int b=0; b<NearField.size(); b++){
    const Block&   B  = NearField[b];
    const Matrix&  M  = NearFieldMat[b];
    const vectInt& It = num_(tgt_(B));
    const vectInt& Is = num_(src_(B));      
    
    ConstSubVectCplx xx(x,Is); 
    SubVectCplx ff(f,It); 
    MvProd(ff,M,xx);
  }
  
}



#endif
