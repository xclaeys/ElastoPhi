#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include "point.hpp"

//================================//
//      DECLARATIONS DE TYPE      //
//================================//
using namespace std;
typedef pair<int,int>            Int2;

//================================//
//      VECTEUR DE COMPLEXES      //
//================================//
typedef vector<Cplx>    vectCplx;
typedef vector<Real>    vectReal;
typedef vector<int>     vectInt;
typedef vector<R3>      vectR3;

void operator+=(vectInt& J, const int& inc){
  for(int k=0; k<J.size(); k++){J[k]+=inc;} }

vectInt operator+(const int& inc, vectInt& J){
  vectInt I(J); for(int k=0; k<I.size(); k++){I[k]+=inc;}
  return I;}

vectInt operator+(vectInt& J, const int& inc){
  vectInt I(J); for(int k=0; k<I.size(); k++){I[k]+=inc;}
  return I;}

int size(const vectCplx& u){return u.size();}

vectCplx operator*(const Cplx& z, const vectCplx& u){
  vectCplx v=u; for(int j=0; j<v.size(); j++){v[j] = z*v[j];}
  return v;}

vectCplx operator+(const vectCplx& u, const vectCplx& v){
  assert(u.size()==v.size());
  vectCplx w=u; for(int j=0; j<v.size(); j++){w[j] = w[j]+v[j];}
  return w;}

vectCplx operator-(const vectCplx& u, const vectCplx& v){
  assert(u.size()==v.size());
  vectCplx w=u; for(int j=0; j<v.size(); j++){w[j] = w[j]-v[j];}
  return w;}

Cplx operator,(const vectCplx& u, const vectCplx& v){
  assert(u.size()==v.size());
  Cplx dot_prod = 0.;
  for(int j=0; j<u.size(); j++){dot_prod += u[j]*v[j];}
  return dot_prod;}

Cplx dprod(const vectCplx& u, const vectCplx& v){
  assert(u.size()==v.size());
  Cplx dot_prod = 0.;
  for(int j=0; j<u.size(); j++){dot_prod += u[j]*conj(v[j]);}
  return dot_prod;}

Real norm(const vectCplx& u){return abs(dprod(u,u));}

template <typename T>
ostream& operator<<(ostream& os, const vector<T>& u){
  for(int j=0; j<u.size(); j++){ os << u[j] << "\t";}
  return os;}

int argmax(const vectCplx& u){
  int k = 0;
  for(int j=0; j<u.size(); j++){
    if( abs(u[j]) > abs(u[k]) ){k=j;}}
  return k;}

//================================//
//         CLASSE MATRICE         //
//================================//
class Matrix{

private:
  
  static const int Dynamic = Eigen::Dynamic;
  typedef Eigen::Matrix<Cplx, Dynamic, Dynamic>  DenseMatrix;
  typedef Eigen::JacobiSVD<DenseMatrix>          SVDType;
  typedef SVDType::SingularValuesType            SgValType;

  
  
  DenseMatrix  mat;
  const int    nr;
  const int    nc;  
  
public:
  Matrix(const int& nbr, const int& nbc):
    mat(nbr,nbc), nr(nbr), nc(nbc){}
  
  Matrix(const Matrix& A):
    nr(A.nr), nc(A.nc), mat(A.mat){}
    
  Cplx& operator()(const int& j, const int& k){
    return mat(j,k);}
  
  const Cplx& operator()(const int& j, const int& k) const {
    return mat(j,k);}  
  
  void operator=(const Matrix& A){
    assert( nr==A.nr && nc==A.nc);
    mat = A.mat;} 
  
  void operator=(const Cplx& z){
    for(int j=0; j<nr; j++){
      for(int k=0; k<nc; k++){
	mat(j,k)=z;}}
  }
  
  vectCplx operator*(const vectCplx& u){
    vectCplx v(nr,0.);
    for(int j=0; j<nr; j++){
      for(int k=0; k<nc; k++){
	v[j]+= mat(j,k)*u[k];
      }
    }
    return v;}
  
  friend ostream& operator<<(ostream& os, const Matrix& m){
    return os << m.mat;} 
  
  friend const int& nb_rows(const Matrix& A){ return A.nr;} 
  
  friend const int& nb_cols(const Matrix& A){ return A.nc;} 
  
  friend vectCplx col(const Matrix& A, const int& k){    
    vectCplx u(A.nr,0.);
    for(int j=0; j<A.nr; j++){u[j]=A(j,k);}
    return u;}
  
  friend vectCplx row(const Matrix& A, const int& j){    
    vectCplx u(A.nc,0.);
    for(int k=0; k<A.nc; k++){u[k]=A(j,k);}
    return u;}
  
  friend Int2 argmax(const Matrix& A){
    int jj=0,kk=0; Real Amax=0.;
    for(int j=0; j<A.nr; j++){
      for(int k=0; k<A.nc; k++){
	if(abs(A(j,k))>Amax){
	  jj=j; kk=k; Amax=abs(A(j,k));
	}
      }
    }
    return Int2(jj,kk);
  }
  
  friend vectReal SVD(const Matrix& A){
  // friend void SVD(const Matrix& A){
    SVDType svd(A.mat);
    const SgValType& sv = svd.singularValues();
    vectReal s(sv.size());
    for(int j=0; j<sv.size(); j++){s[j]=sv[j];}
    return s;
    
//    cout << "sv.size():\t" << sv.size() << endl;
  }

  
};


//================================//
//      CLASSE SOUS-MATRICE       //
//================================//

class SubMatrix{
  
  const Matrix&  mat;
  const vectInt& ir;
  const vectInt& ic;  
  const int nr;
  const int nc;  
  
public:
  
  SubMatrix(const Matrix& mat0, const vectInt& ir0, const vectInt& ic0):
    mat(mat0), ir(ir0), ic(ic0), nr(ir0.size()), nc(ic0.size()) {}

  SubMatrix(const SubMatrix& submat):
    mat(submat.mat), ir(submat.ir), ic(submat.ic), nr(submat.nr), nc(submat.nc) {}
  
  const Cplx& operator()(const int& j, const int& k) const {
    return mat(ir[j],ic[k]);}

  vectCplx operator*(const vectCplx& u){
    vectCplx v(nr,0.);
    for(int j=0; j<nr; j++){
      for(int k=0; k<nc; k++){
	v[j]+= mat(ir[j],ic[k])*u[k];
      }
    }
    return v;}
  
  friend ostream& operator<<(ostream& os, const SubMatrix& m){
    for(int j=0; j<m.nr; j++){ for(int k=0; k<m.nc; k++){ 
	os << m(j,k) << "\t";} os << "\n";}
    return os;} 
  
  friend const int& nb_rows(const SubMatrix& A){ return A.nr;} 
  
  friend const int& nb_cols(const SubMatrix& A){ return A.nc;} 
  
  friend vectCplx col(const SubMatrix& A, const int& k){    
    vectCplx u(A.nr,0.);
    for(int j=0; j<A.nr; j++){u[j]=A(j,k);}
    return u;}
  
  friend vectCplx row(const SubMatrix& A, const int& j){    
    vectCplx u(A.nc,0.);
    for(int k=0; k<A.nc; k++){u[k]=A(j,k);}
    return u;}
  
};



#endif
