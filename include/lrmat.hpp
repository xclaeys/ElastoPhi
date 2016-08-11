#ifndef LRMAT_HPP
#define LRMAT_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "matrix.hpp"
#include "loading.hpp"
#include <Eigen/Dense>


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
	int rank, nr, nc;
	vector<vectCplx> u, v;
	vectInt ir;
	vectInt ic;
	
public:
	
	LowRankMatrix(const vectInt& ir0, const vectInt& ic0){
		nr=ir0.size();
		nc=ic0.size();
		ir=ir0;
		ic=ic0;
		rank=0;
	}
	
	
	//=========================//
	//    PARTIAL PIVOT ACA    //
	//=========================//
	LowRankMatrix(const SubMatrix& A, const vectInt& ir0, const vectInt& ic0, int reqrank=-1){
		Param Parametres;
		nr = nb_rows(A);
		nc = nb_cols(A);
		ir=ir0;
		ic=ic0;
		
		vector<bool> visited_row(nr,false);
		vector<bool> visited_col(nc,false);
		
		Real frob = 0.;
		Real aux  = 0.;
		Cplx frob_aux;
		
		int I=0, J=0;
		int q = 0;
    
        if(reqrank == 0)
            rank = 0; // approximate with a zero matrix
        else{
            vectCplx r(nc),c(nr);
            
            // Compute the first cross
            // (don't modify the code because we want to really use the Bebendorf stopping criterion (3.58),
            // i.e. we don't want to accept the new cross if it is not satisfied because otherwise the approximation would be more precise than desired)
            //==================//
            // Recherche colonne
            Real rmax = 0.;
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
            if( abs(r[J])>1e-10 ){
                Cplx gamma = Cplx(1.)/r[J];
                Real cmax = 0.;
                for(int j=0; j<nr; j++){
                    c[j] = A(j,J);
                    for(int k=0; k<u.size(); k++){
                        c[j] += -u[k][j]*v[k][J];}
                    c[j] = gamma*c[j];
                    if( abs(c[j])>cmax && !visited_row[j] ){
                        I=j; cmax=abs(c[j]);}
                }
                visited_col[J] = true;
                
                aux = abs(dprod(c,c)*dprod(r,r));
            }
            else{cout << "There is a zero coefficient in the original full matrix and ACA didn't work" << endl;}
            
            // (see Bebendorf stopping criterion (3.58) pag 141)
            while ( (q == 0) ||
                      ( (reqrank > 0) && (q < reqrank) ) ||
                     ( (reqrank < 0) && ( sqrt(aux/frob)>Parametres.epsilon * (1 - Parametres.eta)/(1 + Parametres.epsilon) ) ) ) {
                
                // We accept the cross
                q++;
                //====================//
                // Estimateur d'erreur
                frob_aux = 0.;
                //aux = abs(dprod(c,c)*dprod(r,r)); // (already computed to evaluate the test)
                // aux: terme quadratiques du developpement du carre' de la norme de Frobenius de la matrice low rank
                for(int j=0; j<u.size(); j++){
                    frob_aux += dprod(r,v[j])*dprod(c,u[j]);}
                // frob_aux: termes croises du developpement du carre' de la norme de Frobenius de la matrice low rank
                frob += aux + 2*frob_aux.real(); // frob: Frobenius norm of the low rank matrix                
                //==================//
                // Nouvelle croix
                u.push_back(c);
                v.push_back(r);

                if (q >= min(nr,nc) )
                    break;
                // Compute another cross
                //==================//
                // Recherche colonne
                rmax = 0.;
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
                if( abs(r[J])){//>1e-10 ){
                    Cplx gamma = Cplx(1.)/r[J];
                    Real cmax = 0.;
                    for(int j=0; j<nr; j++){
                        c[j] = A(j,J);
                        for(int k=0; k<u.size(); k++){
                            c[j] += -u[k][j]*v[k][J];}
                        c[j] = gamma*c[j];
                        if( abs(c[j])>cmax && !visited_row[j] ){
                            I=j; cmax=abs(c[j]);}
                    }
                    visited_col[J] = true;
                    
                    aux = abs(dprod(c,c)*dprod(r,r));
                }
                else{ cout << "ACA's loop broke" << endl; break; }
            }
        
        // old stopping criterion:
        //}while(sqrt(aux/frob)>Parametres.epsilon && q < min(nr,nc) );
        // (see stopping criterion in slide 26 of Stephanie Chaillat and Rjasanow-Steinbach)
        // si epsilon >=1, always 1 iteration because aux=frob since frob_aux = 0!
        // indeed, it's a sort of relative error
        
            rank = u.size();
        }
            
            
            //======= Historique estimateur ==============//
            //cout << "___________________________" << endl;
            //cout << "Iteration:\t"  << q++ << endl;
            //cout << "Estimateur:\t" << aux/frob << endl;
            //============================================//
    }
	
	LowRankMatrix(const LowRankMatrix& m){
		Param Parametre;
		ir=m.ir;
		ic=m.ic;
		nr=m.nr; nc=m.nc; rank = m.rank;
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
	
//	void Append(const vectCplx& new_u, const vectCplx& new_v){
//		assert(new_u.size()==nr); u.push_back(new_u);
//		assert(new_v.size()==nc); v.push_back(new_v);
//		rank++;
//	}
	
	friend Real NormFrob(const LowRankMatrix& m){
		const vector<vectCplx>& u = m.u;
		const vector<vectCplx>& v = m.v;
		const int& rank = m.rank;
		
		Cplx frob = 0.;
		for(int j=0; j<rank; j++){
			for(int k=0; k<rank; k++){
				frob += dprod(v[k],v[j])*dprod(u[k],u[j]) ;
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
	friend const vectInt& ir_(const LowRankMatrix& m){ return m.ir;}
	friend const vectInt& ic_(const LowRankMatrix& m){ return m.ic;}
	
	friend Real squared_relative_error (const LowRankMatrix& m, const SubMatrix& subm){
		Real norm= 0;
        Real err = 0;
		for (int j=0;j<m.nr;j++){
			for (int k=0;k<m.nc;k++){
				Cplx aux=subm(j,k);
				norm+= pow(abs(aux),2);
				for (int l=0;l<m.u.size();l++){
					aux = aux - m.u[l][j] * m.v[l][k];
				}
				err+=pow(abs(aux),2);
			}
		}
		err =err/norm;
        return err;
	}
	friend Real squared_absolute_error (const LowRankMatrix& m, const SubMatrix& subm){
        Real err=0;
        for (int j=0;j<m.nr;j++){
			for (int k=0;k<m.nc;k++){
				Cplx aux=subm(j,k);
				for (int l=0;l<m.u.size();l++){
					aux = aux - m.u[l][j] * m.v[l][k];
				}
				err+=pow(abs(aux),2);
			}
		}
        return err;
	}
	
};


class LowRankMatrixSVD{
	
private:
	
	int rank, nr, nc;
	vector<vectCplx> u, v;
	vectInt ir;
	vectInt ic;
	
public:
	
	LowRankMatrixSVD(const vectInt& ir0, const vectInt& ic0){
		nr=ir0.size();
		nc=ic0.size();
		ir=ir0;
		ic=ic0;
		rank=0;
	}
	
	
	// Construit une matrix low rank SVD à nombre de matrice de rang 1 fixé
	LowRankMatrixSVD(const SubMatrix& A, const vectInt& ir0, const vectInt& ic0, int k){
		nr = nb_rows(A);
		nc = nb_cols(A);
		ir=ir0;
		ic=ic0;


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
