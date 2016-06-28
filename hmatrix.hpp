#ifndef HMATRIX_HPP
#define HMATRIX_HPP

#include <cassert>
#include <fstream>

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
  
  HMatrix(const Matrix&, const vectR3&, const vectR3&); 
  friend void DisplayPartition(const HMatrix&, char const* const name);
  
};



void HMatrix::BuildBlockTree(const Cluster& t, const Cluster& s){
  Block B(t,s);
  if( B.IsAdmissible() ){
    FarField.push_back(B);}
  else if( s.IsLeaf() ){
    if( t.IsLeaf() ){
      NearField.push_back(B);}    
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

HMatrix::HMatrix(const Matrix& mat0, const vectR3& xt0, const vectR3& xs0):
  mat(mat0), xt(xt0), xs(xs0) {
  assert( nb_rows(mat)==xt.size() && nb_cols(mat)==xs.size() );
  
  // Construction arbre des paquets
  Cluster t(xt); Cluster s(xs);
  
  // Construction arbre des blocs
  BuildBlockTree(t,s);  
  
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



#endif
