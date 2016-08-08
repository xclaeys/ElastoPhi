#ifndef CLUSTER_HPP
#define CLUSTER_HPP
#include <cassert>
#include "matrix.hpp"
#include "loading.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix3d               MatR3;
typedef Eigen::EigenSolver<MatR3>     EigenSolver;
typedef EigenSolver::EigenvectorsType EigenVector;
typedef EigenSolver::EigenvalueType   EigenValue;

//===============================//
//           PAQUETS             //
//===============================//
//
// Refs biblio:
//
//  -> livre de Sauter-Schwab:
//           http://www.springer.com/kr/book/9783540680925
//           et en particulier le paragraphe 7.1.2
//
//  -> livre de Borm:
//           http://www.ems-ph.org/books/book.php?proj_nr=125
//           et en particulier less paragraphes 3.1, 3.2 et 3.3
//
//  -> livre de Rjasanow-Steinbach:
//           http://www.ems-ph.org/books/book.php?proj_nr=125
//           et en particulier le paragraphe 3.1
//
//=================================//


class Cluster{
	
private:
	const vectR3&   x;       // Nuage complet des points
	const vectReal& r;       // Rayon de champs proche pour chaque point
	
	vectInt         num;     // Indices des noeuds du nuage
	
	Cluster*        son[2];  // Paquets enfants
	R3              ctr;     // Centre du paquet
	Real            rad;     // Rayon du champ proche
	
	void Build();
	
public:
	Cluster(const vectR3& x0, const vectReal& r0): x(x0), r(r0), ctr(0.), rad(0.) {
		son[0]=0;son[1]=0;
		for(int j=0; j<x.size(); j++){num.push_back(j);}
		Build();
		NearFieldBall();
	}
	
	Cluster(const vectR3& x0, const vectReal& r0, const int& j0): x(x0), r(r0), ctr(0.), rad(0.) {
		son[0]=0;son[1]=0; num.push_back(j0);}
	Cluster(const Cluster& ); // Pas de recopie
	~Cluster(){if (son[0]!=0){ delete son[0];}if (son[1]!=0){ delete son[1];}};
	void NearFieldBall();
	bool IsLeaf() const { if(son[0]==0){return true;} return false; }
	void push_back(const int& j){num.push_back(j);}
	
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
	
//	ctr = xc;
	
	Real radmax=0.;
	for(int j=0; j<nb_pt; j++){
		if( radmax<norm(xc-x[num[j]]) ){
			radmax=norm(xc-x[num[j]]);} }
	
	// Calcul matrice de covariance
	MatR3 cov; cov.setZero();
//	rad=0.;
	for(int j=0; j<nb_pt; j++){
		R3 u = x[num[j]] - xc;
		//		cout<<rad<<" "<<norm(u)<<" "<<max(rad,norm(u)) <<endl;
//		rad=max(rad,norm(u));
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
	if( max<abs(lambda[2]) ){l=2; }
	R3 w;
	w[0] = ev(0,l).real();
	w[1] = ev(1,l).real();
	w[2] = ev(2,l).real();
	
	// Construction des paquets enfants
	if(radmax>1e-10){
		for(int j=0; j<nb_pt; j++){
			R3 dx = x[num[j]] - xc;
			if( (w,dx)>0 ){
				if(son[0]==0){son[0] = new Cluster(x,r,num[j]);}
				else{ son[0]->push_back(num[j]); }
			}
			else{
				if(son[1]==0){son[1] = new Cluster(x,r,num[j]);}
				else{ son[1]->push_back(num[j]); }
			}
		}
	}
	
	// Recursivite
	if(radmax>1e-10){
		assert( son[0]!=0 );
		assert( son[1]!=0 );
		son[0]->Build();
		son[1]->Build();
	}
//	else{
//		rad=r[num[0]];
//		ctr=x[num[0]];
//	}
//	cout<<ctr<<" "<<rad<<endl;
	
}

void Cluster::NearFieldBall(){
//	// Calcul centre du paquet
//	int nb_pt = num.size();
//	R3 xc = 0.;
//	for(int j=0; j<nb_pt; j++){
//		xc += x[num[j]];}
//	ctr = (1./Real(nb_pt))*xc;
//	
//	rad=0;
//	for(int j=0; j<nb_pt; j++){
//		R3 u =ctr-x[num[j]];
//		Real n = norm(u);
//		rad=max(rad,n);
//	}
//	// Feuille de l'arbre
//	if(son[0]==0){
//		ctr=x[num[0]];
//		rad=r[num[0]];
//		return;
//	}
//	// Recursivite
//	son[0]->NearFieldBall();
//	son[1]->NearFieldBall();
	
	
	// Si deja construit
	if(rad>0){return;}
	
	// Feuille de l'arbre
	if(son[0]==0){
		ctr=x[num[0]];
		rad=r[num[0]];
		return;}
	
	// Recursivite
	son[0]->NearFieldBall();
	son[1]->NearFieldBall();
	
	// Centre et rayon champ proche
	const Real& r0 = son[0]->rad;
	const Real& r1 = son[1]->rad;
	const R3&   c0 = son[0]->ctr;
	const R3&   c1 = son[1]->ctr;
	
	Real l = 0.5*( 1 + (r1-r0)/norm(c1-c0) );
	ctr = (1-l)*c0 + l*c1;
	rad = l*norm(c1-c0)+r0;
	
	cout<<ctr<<" "<<rad<<endl;
	

	
}

//===============================//
//           BLOCK               //
//===============================//
class Block{
	
private:
	
	const Cluster* t;
	const Cluster* s;
	
public:
	Block(const Cluster& t0, const Cluster& s0):  t(&t0), s(&s0) {};
	Block(const Block& b): t(b.t), s(b.s) {};
	Block& operator=(const Block& b){t=b.t; s=b.s; return *this;}
	friend const Cluster& tgt_(const Block& b){return *(b.t);}
	friend const Cluster& src_(const Block& b){return *(b.s);}
	bool IsAdmissible() const{
		Param Parametres;
		return 2*min(rad_(*t),rad_(*s)) < Parametres.eta*( norm(ctr_(*t)-ctr_(*s))-rad_(*t)-rad_(*s) );}
	friend ostream& operator<<(ostream& os, const Block& b){
		os << "src:\t" << src_(b) << endl; os << "tgt:\t" << tgt_(b); return os;}
	
};







#endif
