#ifndef CLUSTER_HPP
#define CLUSTER_HPP
#include <cassert>
#include "matrix.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

typedef Eigen::Matrix3d               MatR3;
typedef Eigen::EigenSolver<MatR3>     EigenSolver;
typedef EigenSolver::EigenvectorsType EigenVector;
typedef EigenSolver::EigenvalueType   EigenValue;

//===============================//
//           PAQUET              // 
//===============================//
class Cluster{
  
private:
  const vectR3&   x;       // Nuage complet de points
  vectInt         num;     // Indices des noeuds du nuage
  Cluster*        son[2];  // Paquets enfants
  R3              ctr;     // Centre du paquet
  Real            rad;     // Rayon du champ proche

  void Build();
  
public:
  Cluster(const vectR3& x0): x(x0), ctr(0.), rad(0.) {
    son[0]=0;son[1]=0;
    for(int j=0; j<x.size(); j++){num.push_back(j);}
    Build();
    
    //=============== TEST ===============//
    Real h = 1./5000.;
    //=============== TEST ===============//
    
    NearFieldBall(h);
  }
  
  Cluster(const vectR3& x0, const int& j0): x(x0), ctr(0.), rad(0.) {
    son[0]=0;son[1]=0; num.push_back(j0);}
  
  Cluster(const Cluster& ); // Pas de recopie 
  
  void NearFieldBall(const Real&);
  bool IsLeaf() const { if(son[0]==0){return true;} return false; }
  void push_back(const int& j){num.push_back(j);}
  void push_back(const vectInt& init){num = init;}  
  
  friend const Real&    rad_(const Cluster& t){return t.rad;}
  friend const R3&      ctr_(const Cluster& t){return t.ctr;}
  friend Cluster&       son_(const Cluster& t,const int& j){return *(t.son[j]);}
  friend const vectInt& num_(const Cluster& t){return t.num;}
  friend ostream& operator<<(ostream& os, const Cluster& cl){
    for(int j=0; j<(cl.num).size(); j++){os<<cl.num[j]<< "\t";} return os;}
  friend void DisplayTree(const Cluster&);
  
};

void DisplayTree(const Cluster& cl){
  cout << cl << endl;
  if(!cl.IsLeaf()){
    DisplayTree(son_(cl,0));
    DisplayTree(son_(cl,1));}
}


void Cluster::Build(){
  
  // Calcul centre du paquet
  int nb_pt = num.size();
  R3 xc = 0.;
  for(int j=0; j<nb_pt; j++){
    xc += x[num[j]];}
  xc = (1./Real(nb_pt))*xc;
  
  // Calcul matrice de covariance
  MatR3 cov; cov.setZero();
  for(int j=0; j<nb_pt; j++){
    R3 u = x[num[j]] - xc; 
    for(int p=0; p<3; p++){
      for(int q=0; q<3; q++){
	cov(p,q) += u[p]*u[q];
      }
    }
  }
  
  // Calcul direction principale
  EigenSolver eig(cov); 
  EigenValue  lambda = eig.eigenvalues();
  EigenVector ev = eig.eigenvectors();
  int l = 0; Real max=abs(lambda[0]);
  if( max<abs(lambda[1]) ){l=1; max=abs(lambda[1]);}
  if( max<abs(lambda[1]) ){l=2;}
  R3 w;
  w[0] = ev(0,l).real();
  w[1] = ev(1,l).real();
  w[2] = ev(2,l).real();  
  
  // Construction des paquets enfants
  if(nb_pt>1){
    for(int j=0; j<nb_pt; j++){
      R3 dx = x[num[j]] - xc;
      if( (w,dx)>0 ){
	if(son[0]==0){son[0] = new Cluster(x,num[j]);}
	else{son[0]->push_back(num[j]);}
      }
      else{
	if(son[1]==0){son[1] = new Cluster(x,num[j]);}
	else{son[1]->push_back(num[j]);}
      }
    }
  }
  
  // Recursivite
  if(nb_pt>1){
    assert( son[0]!=0 );
    assert( son[1]!=0 );
    son[0]->Build();
    son[1]->Build();
  }
  
}

void Cluster::NearFieldBall(const Real& h){

  // Si deja construit
  if(rad>0){return;}
  
  // Feuille de l'arbre
  int nb_pt = num.size();  
  if(nb_pt==1){ctr=x[num[0]]; rad=h; return;}
  
  // Recursivite
  son[0]->NearFieldBall(h);
  son[1]->NearFieldBall(h); 
  
  // Centre et rayon champ proche
  const Real& r0 = son[0]->rad;
  const Real& r1 = son[1]->rad;  
  const R3& c0   = son[0]->ctr;
  const R3& c1   = son[1]->ctr;  
  
  Real l = 0.5*( 1 + (r1-r0)/norm(c1-c0) );
  ctr = (1-l)*c0 + l*c1;
  rad = l*norm(c1-c0)+r0;
  
}

//===============================//
//           BLOCK               // 
//===============================//
class Block{
  
private:
  static const Real eta;  
  const Cluster* t; 
  const Cluster* s; 
  
public:
  Block(const Cluster& t0, const Cluster& s0):  t(&t0), s(&s0) {};
  Block(const Block& b): t(b.t), s(b.s) {};
  Block& operator=(const Block& b){t=b.t; s=b.s; return *this;}
  friend const Cluster& tgt_(const Block& b){return *(b.t);}
  friend const Cluster& src_(const Block& b){return *(b.s);}
  bool IsAdmissible() const{
    return max(rad_(*t),rad_(*s)) < 2*eta*( norm(ctr_(*t)-ctr_(*s))-rad_(*t)-rad_(*s) );}

  friend ostream& operator<<(ostream& os, const Block& b){
    os << "src:\t" << src_(b) << endl; os << "tgt:\t" << tgt_(b); return os;}
  
};
const Real Block::eta = 0.8;


//===============================//
//       ARBRE DES BLOCS         // 
//===============================//
class BlockTree{
  
private:
  vector<Block> FarField;
  vector<Block> NearField;
  
public:
  BlockTree(){};  
  void Build(const Cluster&, const Cluster&);
  friend ostream& operator<<(ostream& os, const BlockTree& bt)
  {
    os << "__________" << endl;
    os << "NearField:" << endl;
    for(int j=0; j<bt.NearField.size(); j++){
      os << bt.NearField[j] << endl << endl;}
    os << "__________" << endl;
    os << "FarField:" << endl;
    for(int j=0; j<bt.FarField.size(); j++){
      os << bt.FarField[j] << endl << endl;}
    return os;
  }  
  
};


void BlockTree::Build(const Cluster& t, const Cluster& s){
  Block B(t,s);
  if( B.IsAdmissible() ){FarField.push_back(B);}
  else if( s.IsLeaf() ){
    if( t.IsLeaf() ){NearField.push_back(B);}    
    else{
      Build(son_(t,0),s);
      Build(son_(t,1),s);
    }
  }
  else{
    if( t.IsLeaf() ){
      Build(t,son_(s,0));
      Build(t,son_(s,1));
    }
    else{
      Build(son_(t,0),son_(s,0));
      Build(son_(t,0),son_(s,1));
      Build(son_(t,1),son_(s,0));
      Build(son_(t,1),son_(s,1));    
    }
  }
  
}







#endif
