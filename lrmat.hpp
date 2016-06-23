#ifndef LRMAT_HPP
#define LRMAT_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "matrix.hpp"


//================================//
//   CLASSE MATRICE RANG FAIBLE   //
//================================//
class LowRankMatrix{

private:  
  int rank, nr, nc;
  vector<vectCplx> u, v;

public:
  
  LowRankMatrix(const int& nbr, const int& nbc){
    nr=nbr; nc=nbc; rank=0;}
  
  LowRankMatrix(const LowRankMatrix& m){
    nr=m.nr; nc=m.nc; rank = m.rank;
    u.resize(rank); v.resize(rank);
    for(int j=0; j<rank; j++){
      u[j] = m.u[j]; v[j] = m.v[j];}
  }
  
  template <typename mat>
  LowRankMatrix(const int&, const mat&);
  
  void operator=(const LowRankMatrix& m){
    nr=m.nr; nc=m.nc; rank = m.rank;
    u.resize(rank); v.resize(rank);
    for(int j=0; j<rank; j++){
      u[j] = m.u[j]; v[j] = m.v[j];}
  }
  
  void Append(const vectCplx& new_u, const vectCplx& new_v){    
    assert(new_u.size()==nr); u.push_back(new_u);
    assert(new_v.size()==nc); v.push_back(new_v);    
    rank++;
  }

  friend Real Frob(const LowRankMatrix& m){
    const vector<vectCplx>& u = m.u;
    const vector<vectCplx>& v = m.v;
    const int& rank = m.rank;
    
    Cplx frob = 0.;
    for(int j=0; j<rank; j++){
      for(int k=0; k<rank; k++){
	frob += dprod(v[j],v[k])*dprod(u[k],u[j]) ;
      }
    }
    return sqrt(abs(frob));
  }
  
  vectCplx operator*(const vectCplx& w){
    assert(w.size()==nc);
    vectCplx res(nr,0.);
    for(int k=0; k<v.size(); k++){
      Cplx pk = (v[k],w);
      for(int j=0; j<nr; j++){
	res[j] += pk*u[k][j];
      }
    }
    return res;
  }
  
  friend ostream& operator<<(ostream& os, const LowRankMatrix& m){    
    os << "rank:\t" << m.rank << endl;
    os << "nr:\t"   << m.nr << endl;
    os << "nc:\t"   << m.nc << endl;
    os << "\nu:\n";    
    for(int j=0; j<m.nr; j++){
      for(int k=0; k<m.rank; k++){
	cout << m.u[k][j] << "\t";}
      cout << "\n";}
    os << "\nv:\n";    
    for(int j=0; j<m.nc; j++){
      for(int k=0; k<m.rank; k++){
	cout << m.v[k][j] << "\t";
      }
      cout << "\n";
    }
    
    return os;
  }
  
};


//=========================//
//    PARTIAL PIVOT ACA    //
//=========================//
template <typename mat>
LowRankMatrix::LowRankMatrix(const int& rk, const mat& A){
  
  nr = nb_rows(A);
  nc = nb_cols(A);  
  vector<bool> visited_row(nr,false); 
  vector<bool> visited_col(nc,false);   
  
  int I=0, J=0;
  for(int q=0; q<rk; q++){

    //==================//
    // Recherche colonne
    Real rmax = 0.;
    vectCplx r(nc);
    for(int k=0; k<nc; k++){
      r[k] = A(I,k);
      for(int j=0; j<u.size(); j++){
	r[k] += -u[j][I]*v[j][k];}
      if( abs(r[k])>rmax && !visited_col[k] ){
	J=k; rmax=abs(r[k]);}
    }
    visited_row[I] = true;
    
    //==================//
    // Recherche ligne
    if( abs(r[J]) ){
      Cplx gamma = Cplx(1.)/r[J];
      Real cmax = 0.;
      vectCplx c(nr);
      for(int j=0; j<nr; j++){
	c[j] = A(j,J);
	for(int k=0; k<u.size(); k++){
	  c[j] += -u[k][j]*v[k][J];}
	c[j] = gamma*c[j];
	if( abs(c[j])>cmax && !visited_row[j] ){
	  I=j; cmax=abs(c[j]);}
      }
      visited_col[J] = true;
      
      //==================//
      // Nouvelle croix
      u.push_back(c);
      v.push_back(r);
    }
  }

  rank = u.size();
  
}


//======================//
//    FULL PIVOT ACA    //
//======================//
/*
LowRankMatrix::LowRankMatrix(const int& rk, const Matrix& A){
  rank = 0;
  nr = nb_rows(A);
  nc = nb_cols(A);  
  Matrix R=A;

  Int2 ind = argmax(R);
  int jj=ind.first, kk=ind.second;
  Real rmax=0.;
  
  for(int q=0; q<rk; q++){
    
    if( abs(R(jj,kk))<1e-10 ){
      break;}
    else{
      Cplx gamma = 1./R(jj,kk);
      vectCplx c =  gamma*col(A,kk);          
      vectCplx r =  row(A,jj);         
      
      rmax = 0.;
      for(int j=0; j<nr; j++){
	for(int k=0; k<nc; k++){
	  R(j,k) = R(j,k) - c[j]*r[k];
	  if(abs(R(j,k))>rmax){
	    rmax=abs(R(j,k)); jj=j; kk=k;}
	}
      }
      
      u.push_back(c);
      v.push_back(r);
      rank++;
    }
    
  }

}
==========================*/
#endif
