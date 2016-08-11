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

template <typename T>
int size(const vector<T>& u){return u.size();}

void fill(vectCplx& u, const Cplx& v){
	for(int j=0; j<u.size(); j++){u[j]=v;}}

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
//      CLASSE SUBVECTOR          //
//================================//

template <typename VecType>
class SubVec{
	
private:
	VecType&       U;
	const vectInt& I;
	const int      size;
	
	typedef typename VecType::value_type ValType;
	
public:
	
	SubVec(VecType& U0, const vectInt& I0): U(U0), I(I0), size(I0.size()) {}
	SubVec(const SubVec&); // Pas de constructeur par recopie
	
	Cplx& operator[](const int& k) {return U[I[k]];}
	const Cplx& operator[](const int& k) const {return U[I[k]];}
	
	void operator=(const ValType& v){
		for(int k=0; k<size; k++){ U[I[k]]=v;}}
	
	template <typename RhsType>
	ValType operator,(const RhsType& rhs) const {
		ValType lhs = 0.;
		for(int k=0; k<size; k++){lhs += U[I[k]]*rhs[k];}
		return lhs;
	}
	
	friend int size(const SubVec& sv){ return sv.size;}
	
	friend ostream& operator<<(ostream& os, const SubVec& u){
		for(int j=0; j<u.size; j++){ os << u[j] << "\t";}
		return os;}
	
};

typedef SubVec<vectCplx> SubVectCplx;
typedef SubVec<const vectCplx> ConstSubVectCplx;


void fill(SubVectCplx& u, const Cplx& v){
	for(int j=0; j<size(u); j++){u[j]=v;}}



//=================================================================//
//                         CLASS MATRIX
/******************************************************************//**
* This class is a wrapper for the Matrix class with Dynamic size of
*  the [Eigen3](http://eigen.tuxfamily.org/dox/) library.
*  Elements of the class Matrix represent dense matrices.
*
*  The number of rows and columns can be changed after initialisation.
*
*  Elements of this class store
*     - mat: an instance of an Eigen3 matrix
*     - nr:  the number of rows
*     - nc:  the number of columns
*
*********************************************************************/

class Matrix{
	
private:
	
	static const int Dynamic = Eigen::Dynamic;
	typedef Eigen::Matrix<Cplx, Dynamic, Dynamic>  DenseMatrix;
	typedef Eigen::JacobiSVD<DenseMatrix>          SVDType;
	typedef SVDType::SingularValuesType            SgValType;
	typedef SVDType::MatrixUType		       UMatrixType;
	typedef SVDType::MatrixVType		       VMatrixType;
	
	
	DenseMatrix  mat;
	int  nr;
	int  nc;
	
public:
	
	//! ### Default constructor
	/*!
	 Initialise the matrix to the size 0*0
  */
	Matrix(): mat(0,0), nr(0), nc(0){}
	
	
	//! ### Another constructor
	/*!
	 Initialise the matrix with _nbr_ rows and _nbc_ columns,
	 and fills the matrix with zeros
  */
	Matrix(const int& nbr, const int& nbc):
	mat(nbr,nbc), nr(nbr), nc(nbc){}
	
	
	//! ### Copy constructor
	/*!
	 Initialise the matrix with _nbr_ rows and _nbc_ columns,
	 and fills the matrix with zeros
  */
	Matrix(const Matrix& A):
	mat(A.mat),nr(A.nr), nc(A.nc){}
	
	//! ### Templated copy constructor
	/*!
	 Does the same as the copy constructor
	 but the input argument can _A_ be of any type.
	 The only requirement is that the parenthesis operator
	 be overloaded so that the expression _A(j,k)_ provides
	 access to the elements of _A_.
  */
	template <typename MatType>
	Matrix(const MatType& A):
	mat(nb_rows(A),nb_cols(A)),nr(nb_rows(A)), nc(nb_cols(A)){
		for(int j=0; j<nr; j++){
			for(int k=0; k<nc; k++){
				mat(j,k) = A(j,k);
			}
		}
	}
	
	
	//! ### Access operator
	/*!
	 If _A_ is the instance calling the operator
	 _A(j,k)_ returns the entry of _A_ located
	 jth row and kth column. Modification of the entries
	 are allowed.
  */
	Cplx& operator()(const int& j, const int& k){
		return mat(j,k);}
	
	
	//! ### Access operator
	/*!
	 If _A_ is the instance calling the operator
	 _A(j,k)_ returns the entry of _A_ located
	 jth row and kth column. Modification of the
	 entries are forbidden.
  */
	const Cplx& operator()(const int& j, const int& k) const {
		return mat(j,k);}
	
	//! ### Assignement operator with matrix input argument
	/*!
	 Copies the value of the entries of the input _A_
	 (which is a matrix) argument into the entries of
	 calling instance.
  */
	void operator=(const Matrix& A){
		assert( nr==A.nr && nc==A.nc);
		mat = A.mat;}
	
	//! ### Assignement operator with scalar input argument
	/*!
	 Sets the values of the entries of the calling instance
	 to the input value _z_
  */
	void operator=(const Cplx& z){
		for(int j=0; j<nr; j++){
			for(int k=0; k<nc; k++){
				mat(j,k)=z;}}
	}
	
	//! ### Matrix-vector product
	/*!
	 Naive self-contained implementation of matrix-vector product
	 for dense matrices. The input parameter _u_ is the input vector
	 (i.e. the right operand). This operator does not rely on the
	 matrix-vector operator obtained via the library eigen3
  */
	vectCplx operator*(const vectCplx& u){
		vectCplx v(nr,0.);
		for(int j=0; j<nr; j++){
			for(int k=0; k<nc; k++){
				v[j]+= mat(j,k)*u[k];
			}
		}
		return v;}
	
	//! ### Modifies the size of the matrix
	/*!
	 Changes the size of the matrix so that
	 the number of rows is set to _nbr_ and
	 the number of columns is set to _nbc_
  */
	void resize(const int nbr, const int nbc){
		mat.resize(nbr,nbc); nr = nbr; nc = nbc;}
	
	
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
	friend void MvProd(LhsType& lhs, const Matrix& m, const RhsType& rhs){
		for(int j=0; j<m.nr; j++){
			for(int k=0; k<m.nc; k++){
				lhs[j]+= m.mat(j,k)*rhs[k];
			}
		}
	}
	
	friend ostream& operator<<(ostream& os, const Matrix& m){
		return os << m.mat;}
	
	
	//! ### Access to number of rows
	/*!
	 Returns the number of rows of the input argument _A_
  */
	friend const int& nb_rows(const Matrix& A){ return A.nr;}
	
	//! ### Access to number of columnss
	/*!
	 Returns the number of columns of the input argument _A_
  */
	friend const int& nb_cols(const Matrix& A){ return A.nc;}
	
	
	//! ### Extraction of a column
	/*!
	 Returns, as a vector, the column numbered _k_ ofthe matrix _A_
  */
	friend vectCplx col(const Matrix& A, const int& k){
		vectCplx u(A.nr,0.);
		for(int j=0; j<A.nr; j++){u[j]=A(j,k);}
		return u;}
	
	//! ### Extraction of a row
	/*!
	 Returns, as a vector, the row numbered _j_ ofthe matrix  _A_
  */
	friend vectCplx row(const Matrix& A, const int& j){
		vectCplx u(A.nc,0.);
		for(int k=0; k<A.nc; k++){u[k]=A(j,k);}
		return u;}
	
	//! ### Looking for the entry of maximal modulus
	/*!
	 Returns the number of row and column of the entry
	 of maximal modulus in the matrix _A_
  */
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
	
	
	//! ### Computation of singular values
	/*!
	 Returns a vector of Real containing the singular values
	 of the input matrix _A_ in decreasing order
  */
	friend vectReal SVD(const Matrix& A){
		SVDType svd(A.mat);
		const SgValType& sv = svd.singularValues();
		vectReal s(sv.size());
		for(int j=0; j<sv.size(); j++){s[j]=sv[j];}
		return s;
	}
	
	friend void PartialSVD(const Matrix& A, vector<vectCplx>& u ,vector<vectCplx>& v, int k){
		assert(k<=min(A.nr,A.nc));
		SVDType svd(A.mat,Eigen::ComputeThinU | Eigen::ComputeThinV );
		const SgValType& sv = svd.singularValues();
		
		const UMatrixType& uu = svd.matrixU();
		const VMatrixType& vv = svd.matrixV();
//		cout<<A.nr<<" "<<A.nc<<endl;
//		cout<<uu.size()<<" "<<vv.size()<<endl;
//		cout<<uu.rows()<<" "<<vv.cols()<<endl;
//		cout<<uu.rows()<<" "<<vv.cols()<<endl;
		
		
		for (int i=0;i<k;i++){
			vector<Cplx> uuu;
			vector<Cplx> vvv;
			for (int j=0;j<A.nr;j++){
				uuu.push_back(uu(j,i)*sv[i]);
			}
			for (int j=0;j<A.nc;j++){
				vvv.push_back(vv(j,i));
			}
			u.push_back(uuu);
			v.push_back(vvv);
			
		}
	}
    
    friend Real NormFrob (const Matrix& A){
        Real norm=0;
        for (int j=0;j<A.nr;j++){
            for (int k=0;k<A.nc;k++){
                norm = norm + pow(abs(A(j,k)),2);
            }
        }
        return sqrt(norm);
    }

	
	
};


//================================//
//      CLASSE SOUS-MATRICE       //
//================================//

class SubMatrix{
	
	const Matrix&  mat;
	vectInt ir;
	vectInt ic;
	int nr,nc;
	
public:
	
	SubMatrix(const Matrix& mat0, const vectInt& ir0, const vectInt& ic0):
	mat(mat0), ir(ir0), ic(ic0), nr(ir0.size()), nc(ic0.size()) {}
	
	SubMatrix(const SubMatrix& m): mat(m.mat){
		ir=m.ir;
		ic=m.ic;
		nr=m.nr;
		nc=m.nc;
	}
	
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
	
	friend const vectInt& ir_(const SubMatrix& A){ return A.ir;}
	
	friend const vectInt& ic_(const SubMatrix& A){ return A.ic;}
	
	friend const Matrix& mat_(const SubMatrix& A){ return A.mat;}
	
	friend vectCplx col(const SubMatrix& A, const int& k){
		vectCplx u(A.nr,0.);
		for(int j=0; j<A.nr; j++){u[j]=A(j,k);}
		return u;}
	
	friend vectCplx row(const SubMatrix& A, const int& j){
		vectCplx u(A.nc,0.);
		for(int k=0; k<A.nc; k++){u[k]=A(j,k);}
		return u;}
	
	template <typename LhsType, typename RhsType>
	friend void MvProd(LhsType& lhs, const SubMatrix& m, const RhsType& rhs){
		for(int j=0; j<m.nr; j++){
			for(int k=0; k<m.nc; k++){
				lhs[j]+= m.mat(m.ir[j],m.ic[k])*rhs[k];
			}
		}
	}
	friend Real NormFrob (const SubMatrix& A){
		Real norm=0;
		for (int j=0;j<A.ir.size();j++){
			for (int k=0;k<A.ic.size();k++){
				norm = norm + pow(abs(A(j,k)),2);
			}
		}
		return sqrt(norm);
	}
	
};



#endif
