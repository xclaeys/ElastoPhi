#ifndef HMATRIX_HPP
#define HMATRIX_HPP

#include <cassert>


//===============================//
//     MATRICE HIERARCHIQUE      // 
//===============================//

class HMatrix{

private:

  const Matrix&         Mat;
  vector<Block>         FarField;
  vector<Block>         NearField;
  vector<Matrix>        NearFieldMat;
  vector<LowRankMatrix> FarFieldMat;  
  
  
  void BuildBlockTree(const Cluster&, const Cluster&);
  
public:
  
  HMatrix(const Matrix&,const Cluster&, const Cluster&); 
  
};



void HMatrix::BuildBlockTree(const Cluster& t, const Cluster& s){
  Block B(t,s);
  if( B.IsAdmissible() ){FarField.push_back(B);}
  else if( s.IsLeaf() ){
    if( t.IsLeaf() ){NearField.push_back(B);}    
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


HMatrix(const Matrix& M, const vectR3& xt, const vectR3& xs): Mat(M) {
  assert( nb_rows(M)==xt.size() && nb_cols(M)==xs.size() );

  // Construction arbre des paquets cibles
  vectInt numt(nb_rows(M));  
  for(int j=0; j<nb_rows(M); j++){numt[j]=j;}
  Cluster t(xt); t.push_back(numt); t.Build();
  
  // Construction arbre des paquets sources
  vectInt nums(nb_cols(M));  
  for(int j=0; j<nb_cols(M); j++){nums[j]=j;}
  Cluster s(xs); s.push_back(nums); s.Build();
  
  // Construction arbre des block
  BuildBlockTree(t,s);  
  
}



#endif
