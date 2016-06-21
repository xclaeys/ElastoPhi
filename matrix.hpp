#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include <gmm.h>
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
typedef vector<int>     vectInt;
typedef vector<R3>      vectR3;

int size(const vectCplx& u){return u.size();}

vectCplx operator*(const Cplx& z, const vectCplx& u){
  vectCplx v=u; for(int j=0; j<v.size(); j++){v[j] = z*v[j];}
  return v;}

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

ostream& operator<<(ostream& os, const vectCplx& u){
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
  gmm::dense_matrix<Cplx>  mat;
  const int nr;
  const int nc;  
  
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
    vectCplx v(u.size());
    gmm::mult(mat,u,v);
    return v;}
    
  friend ostream& operator<<(ostream& os, const Matrix& m){
    return os << m.mat;} 
  
  friend const int& nb_rows(const Matrix& A){ return A.nr;} 
  
  friend const int& nb_cols(const Matrix& A){ return A.nc;} 

  friend vectCplx col(const Matrix& A, const int& k){    
    vectCplx u(nb_rows(A),0.);
    for(int j=0; j<nb_rows(A); j++){u[j]=A(j,k);}
    return u;}
  
  friend vectCplx row(const Matrix& A, const int& j){    
    vectCplx u(nb_cols(A),0.);
    for(int k=0; k<nb_cols(A); k++){u[k]=A(j,k);}
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
  
};





#endif
