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
//
// Refs biblio:
//
//  -> slides de Stéphanie Chaillat:
//           http://uma.ensta-paristech.fr/var/files/chaillat/seance2.pdf
//           et en particulier la slide 25                    
//
//  -> livre de M.Bebendorf:  
//           http://www.springer.com/kr/book/9783540771463
//           et en particulier le paragraphe 3.4
//
//  -> livre de Rjasanow-Steinbach:
//           http://www.ems-ph.org/books/book.php?proj_nr=125
//           et en particulier le paragraphe 3.2
//
//=================================//

class LowRankMatrix{

private:  

  Real epsilon;
  int rank, nr, nc;
  vector<vectCplx> u, v;
  
public:
  
  LowRankMatrix(const int& nbr, const int& nbc){
    nr=nbr; nc=nbc; rank=0; epsilon = 1e-6;}
  
  //=========================//
  //    PARTIAL PIVOT ACA    //
  //=========================//
  template <typename mat>
  LowRankMatrix(const mat& A, const Real& eps = 1e-4){
    
    epsilon = eps;
    nr = nb_rows(A);
    nc = nb_cols(A);  
    vector<bool> visited_row(nr,false); 
    vector<bool> visited_col(nc,false);   
    
    Real frob = 0.;
    Real aux;
    Cplx frob_aux;
    
    int I=0, J=0;
    int q = 0;
    while( (aux/frob)>epsilon && q < min(nr,nc) ){
      
      
      q++;
      //======= Historique estimateur ==============//
      //cout << "___________________________" << endl;
      //cout << "Iteration:\t"  << q++ << endl; 
      //cout << "Estimateur:\t" << aux/frob << endl;
      //============================================//

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
	
	//====================//
	// Estimateur d'erreur
	frob_aux = 0.;
	aux = abs(dprod(c,c)*dprod(r,r));      
	for(int j=0; j<u.size(); j++){
	  frob_aux += dprod(v[j],r)*dprod(c,u[j]);}
	frob += aux + 2*frob_aux.real();
	
	//==================//
	// Nouvelle croix
	u.push_back(c);
	v.push_back(r);
      }
      else{ break; }
    }
    
    rank = u.size();
  }
    
  LowRankMatrix(const LowRankMatrix& m){
    nr=m.nr; nc=m.nc; rank = m.rank; epsilon = m.epsilon;
    u.resize(rank); v.resize(rank);
    for(int j=0; j<rank; j++){
      u[j] = m.u[j]; v[j] = m.v[j];}
  }
  
  void operator=(const LowRankMatrix& m){
    nr=m.nr; nc=m.nc; rank = m.rank;
    u.resize(rank); v.resize(rank);
    for(int j=0; j<rank; j++){
      u[j] = m.u[j]; v[j] = m.v[j];}
  }
  
  friend Real CompressionRate(const LowRankMatrix& m){
    return m.rank*( 1./Real(m.nr) + 1./Real(m.nc) );
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
  
  vectCplx operator*(const vectCplx& rhs){
    assert(rhs.size()==nc);
    vectCplx lhs(nr,0.);
    for(int k=0; k<v.size(); k++){
      Cplx pk = (v[k],rhs);
      for(int j=0; j<nr; j++){
	lhs[j] += pk*u[k][j];
      }
    }
    return lhs;
  }
  
  template <typename LhsType, typename RhsType>
  friend void MvProd(LhsType& lhs, const LowRankMatrix& m, const RhsType& rhs){
    const vector<vectCplx>& u = m.u;   
    const vector<vectCplx>& v = m.v;
    for(int k=0; k<v.size(); k++){
      Cplx pk = (rhs,v[k]); 
      for(int j=0; j<m.nr; j++){
	lhs[j] += pk*u[k][j];
      }
    }
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
  
  friend const int& rank_of(const LowRankMatrix& m){ return m.rank;}
  friend const int& nb_rows(const LowRankMatrix& m){ return m.nr;}   
  friend const int& nb_cols(const LowRankMatrix& m){ return m.nc;}     
  
};




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
