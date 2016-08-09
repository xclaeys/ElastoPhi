#ifndef SPMATRIX_HPP
#define SPMATRIX_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "point.hpp"

//================================//
//      DECLARATIONS DE TYPE      //
//================================//
using namespace std;
//typedef pair<int,int>            Int2;

typedef vector<Cplx>    vectCplx;

//=================================================================//
//                         CLASS SPARSE MATRIX
/******************************************************************//**
* This class is a wrapper for the SparseMatrix class of
*  the [Eigen3](http://eigen.tuxfamily.org/dox/) library....
*********************************************************************/

//The simplest way to create a sparse matrix while guaranteeing good performance is thus to first build a list of so-called triplets, and then convert it to a SparseMatrix.
//typedef Eigen::Triplet<double> T;
//std::vector<T> tripletList;
//tripletList.reserve(estimation_of_entries);
//for(...)
//{
//    // ...
//    tripletList.push_back(T(i,j,v_ij));
//}
//SparseMatrixType mat(rows,cols);
//mat.setFromTriplets(tripletList.begin(), tripletList.end());
//// mat is ready to go!

class SpMatrix{
	
private:
		
    vector<int> I,J;
    vector<Cplx> K;
	int  nr;
	int  nc;
	
public:
	
	//! ### Default constructor
	/*!
	 Initialise the matrix to the size 0*0
  */
	SpMatrix(): nr(0), nc(0){}
	
	
	//! ### Another constructor
	/*!
	 Initialise the matrix with _nrp_ rows and _ncp_ columns, Ip, Jp, Kp
  */
	SpMatrix(const vector<int>& Ip, const vector<int>& Jp, vector<Cplx>& Kp, const int& nrp, const int& ncp):
    I(Ip), J(Jp), K(Kp), nr(nrp), nc(ncp) {}
	
	
	//! ### Copy constructor
	/*!
  */
	SpMatrix(const SpMatrix& A):
	I(A.I), J(A.J), K(A.K), nr(A.nr), nc(A.nc){}
		
	
	//! ### Access operator
	/*!
	 If _A_ is the instance calling the operator
	 _A(j,k)_ returns the entry of _A_ located
	 jth row and kth column.
  */
    /*
	Cplx& operator()(const int& j, const int& k){
		return mat.coeffRef(j,k);}
	*/
    
	//! ### Assignement operator with matrix input argument
	/*!
	 Copies the value of the entries of the input _A_
	 (which is a matrix) argument into the entries of
	 calling instance.
  */
	void operator=(const SpMatrix& A){
		assert( nr==A.nr && nc==A.nc);
        I = A.I; J = A.J;  K = A.K;}
	
    
	//! ### Matrix-vector product
	/*!
	 The input parameter _u_ is the input vector
	 (i.e. the right operand).
  */
	vectCplx operator*(const vectCplx& u){
        int ncoef = I.size();
        vectCplx v(nr,0.);
        for(int j=0; j<ncoef; j++)
                v[I[j]]+= K[j]*u[J[j]];
        return v;}
	
	//! ### Matrix-vector product
	/*!
	 Another instanciation of the matrix-vector product
	 that avoids the generation of temporary instance for the
	 output vector. This routine achieves the operation
	 
	 lhs = m*rhs
	 
	 The left and right operands (_lhs_ and _rhs_) are templated
	 and can then be of any type (not necessarily of type vectCplx).
	 The only requirement is that an overload of the parentesis-based
	 access operator be available for the operands.
  */
	template <typename LhsType, typename RhsType>
	friend void MvProd(LhsType& lhs, const SpMatrix& m, const RhsType& rhs){
        int ncoef = m.I.size();
        for(int j=0; j<ncoef; j++)
            lhs[m.I[j]]+= m.K[j]*rhs[m.J[j]];
	}
    
    
    //! ### Modifies the size of the matrix
    /*!
     Changes the size of the matrix so that
     the number of rows is set to _nbr_ and
     the number of columns is set to _nbc_ and ...
     */
    void resize(const int nbr, const int nbc, const int nbcoef){
        assert(nbcoef<=nbr*nbc);
        nr = nbr; nc = nbc;
        I.resize(nbcoef); J.resize(nbcoef); K.resize(nbcoef);
    }
    
    int& I_(const int i){
        assert(i<I.size());
        return I[i];
    }
    
    int& J_(const int i){
        assert(i<J.size());
        return J[i];
    }
    
    Cplx& K_(const int i){
        assert(i<K.size());
        return K[i];
    }
	
    
//	friend ostream& operator<<(ostream& os, const SpMatrix& m){
//		return os << m...;}
	
	
	//! ### Access to number of rows
	/*!
	 Returns the number of rows of the input argument _A_
  */
	friend const int& nb_rows(const SpMatrix& A){ return A.nr;}
	
	//! ### Access to number of columnss
	/*!
	 Returns the number of columns of the input argument _A_
  */
	friend const int& nb_cols(const SpMatrix& A){ return A.nc;}
	
    
    //! ### Access to number of non zero coefficients
    /*!
     Returns the number of non zero coefficients of the input argument _A_
     */
	friend int nb_coeff(const SpMatrix& A){ return A.I.size();}
};





#endif
