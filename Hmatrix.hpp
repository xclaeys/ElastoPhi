#ifndef HMATRIX_HPP
#define HMATRIX_HPP

//===============================//
//     MATRICE HIERARCHIQUE      // 
//===============================//

class HMatrix{


  
private:
  vector<Block>         FarField;
  vector<Block>         NearField;
  vector<Matrix>        NearFieldMat;
  vector<LowRankMatrix> FarFieldMat;  
  
public:
  BlockTree(){};  
  void Assemble(const Cluster&, const Cluster&, const Matrix&);

};



#endif
