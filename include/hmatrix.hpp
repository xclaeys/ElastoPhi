#ifndef HMATRIX_HPP
#define HMATRIX_HPP

#include <cassert>
#include <fstream>
#include "matrix.hpp"
#include "loading.hpp"
#include "lrmat.hpp"

//===============================//
//     MATRICE HIERARCHIQUE      //
//===============================//

class HMatrix: public Parametres{
	
private:
	
	const Matrix& mat;
	const vectR3& xt;
	const vectR3& xs;
	const vectInt& tabt;
	const vectInt& tabs;
	
	
	vector<LowRankMatrix> FarFieldMat;
	vector<SubMatrix>     NearFieldMat;
	
	void BuildBlockTree(const Cluster&, const Cluster&, int reqrank=-1);
	
public:
	
	HMatrix(const Matrix&, const vectR3&, const vectReal&, const vectInt&, const vectR3&, const vectReal&, const vectInt&); // To be used with two different clusters
	HMatrix(const Matrix&, const vectR3&, const vectReal&, const vectInt&, int reqrank=-1); // To be used with one cluster
	//	friend void DisplayPartition(const HMatrix&, char const* const);
	friend void MvProd(vectCplx&, const HMatrix&, const vectCplx&);
	friend Real CompressionRate(const HMatrix&);
	friend void Output(const HMatrix&, string filename);
	friend const LowRankMatrix& GetLowRankMatrix(HMatrix m, int i){
		assert(i<m.FarFieldMat.size());
		return m.FarFieldMat[i];}
	friend Real squared_absolute_error(const HMatrix& B, const Matrix& A);
	
};

void HMatrix::BuildBlockTree(const Cluster& t, const Cluster& s, int reqrank){
	Block B(t,s);
	if( B.IsAdmissible() ){
		const vectInt& I = num_(t);
		const vectInt& J = num_(s);
		SubMatrix submat = SubMatrix(mat,I,J);
		LowRankMatrix lrm(submat,I,J,t,s,reqrank);
		if(rank_of(lrm)!=-5){
			FarFieldMat.push_back(lrm);
			return;
		}
	}
	if( s.IsLeaf() ){
		if( t.IsLeaf() ){
			const vectInt& I = num_(t);
			const vectInt& J = num_(s);
			NearFieldMat.push_back(SubMatrix(mat,I,J));
		}
		else{
			BuildBlockTree(son_(t,0),s, reqrank);
			BuildBlockTree(son_(t,1),s, reqrank);
		}
	}
	else{
		if( t.IsLeaf() ){
			BuildBlockTree(t,son_(s,0), reqrank);
			BuildBlockTree(t,son_(s,1), reqrank);
		}
		else{
			BuildBlockTree(son_(t,0),son_(s,0), reqrank);
			BuildBlockTree(son_(t,0),son_(s,1), reqrank);
			BuildBlockTree(son_(t,1),son_(s,0), reqrank);
			BuildBlockTree(son_(t,1),son_(s,1), reqrank);
		}
	}
	
}

HMatrix::HMatrix(const Matrix& mat0,
		 const vectR3& xt0, const vectReal& rt, const vectInt& tabt0,
		 const vectR3& xs0, const vectReal& rs, const vectInt& tabs0):

mat(mat0), xt(xt0), xs(xs0), tabt(tabt0), tabs(tabs0) {
	assert( nb_rows(mat)==tabt.size() && nb_cols(mat)==tabs.size() );
	
	// Construction arbre des paquets
	Cluster t(xt,rt,tabt); Cluster s(xs,rs,tabs);
	
	// Construction arbre des blocs
	BuildBlockTree(t,s);
	
}

HMatrix::HMatrix(const Matrix& mat0,
		 const vectR3& xt0, const vectReal& rt, const vectInt& tabt0, int reqrank):

mat(mat0), xt(xt0), xs(xt0), tabt(tabt0), tabs(tabt0) {
	assert( nb_rows(mat)==tabt.size() && nb_cols(mat)==tabs.size() );
	
	// Construction arbre des paquets
	Cluster t(xt,rt,tabt);
	
	// Construction arbre des blocs
	BuildBlockTree(t,t,reqrank);
	
}

// 1- !!!
Real CompressionRate(const HMatrix& hmat){
	
	Real comp = 0.;
	Real size = ( (hmat.tabt).size() )*( (hmat.tabs).size() );
	const vector<LowRankMatrix>& FarFieldMat  = hmat.FarFieldMat;
	const vector<SubMatrix>&     NearFieldMat = hmat.NearFieldMat;
	
	for(int j=0; j<FarFieldMat.size(); j++){
		Real nr   = nb_rows(FarFieldMat[j]);
		Real nc   = nb_cols(FarFieldMat[j]);
		Real rank = rank_of(FarFieldMat[j]);
		comp += rank*(nr + nc)/size;
	}
	
	for(int j=0; j<NearFieldMat.size(); j++){
		Real nr   = nb_rows(NearFieldMat[j]);
		Real nc   = nb_cols(NearFieldMat[j]);
		comp += nr*nc/size;
	}
	
	return 1-comp;
	
}

void Output(const HMatrix& hmat, string filename){
	string path=GetOutputPath()+"/"+filename;
	
	ofstream outputfile(path.c_str());
	
	if (!outputfile){
		cerr << "Output file cannot be created in "+path << endl;
		exit(1);
	}
	else{
		
		const vector<LowRankMatrix>& FarFieldMat  = hmat.FarFieldMat;
		//		const vector<SubMatrix>&     NearFieldMat = hmat.NearFieldMat;
		
		for(int i=0; i<FarFieldMat.size(); i++){
			
			vectInt ir = ir_(FarFieldMat[i]);
			vectInt ic = ic_(FarFieldMat[i]);
			Real local_compression = CompressionRate(FarFieldMat[i]);
			
			for (int j=0;j<ir.size();j++){
				for (int k=0;k<ic.size();k++){
					outputfile<<ir[j]<<" "<<ic[k]<<" "<<local_compression<<endl;
				}
			}
		}
		
		//		for(int i=0; i<NearFieldMat.size(); i++){
		//			vectInt ir = ir_(NearFieldMat[i]);
		//			vectInt ic = ic_(NearFieldMat[i]);
		//			for (int j=0;j<ir.size();j++){
		//				for (int k=0;k<ic.size();k++){
		//					outputfile<<ir[j]<<" "<<ic[k]<<" "<<1<<endl;
		//				}
		//			}
		//		}
	}
	
}




// Representation graphique de la partition en bloc
//void DisplayPartition(const HMatrix& hmat, char const * const name){
//
//	const vector<Block>& FarField = hmat.FarField;
//	const vectR3& xt = hmat.xt;
//	const vectR3& xs = hmat.xs;;
//
//	// Representation graphique
//	const int  Ns = xs.size();
//	const Real ds = 1./Real(Ns-1);
//	const int  Nt = xt.size();
//	const Real dt = 1./Real(Nt-1);
//
//	ofstream file; file.open(name);
//	for(int j=0; j<FarField.size(); j++){
//
//		const Cluster& t = tgt_(FarField[j]);
//		const vectInt& It = num_(t);
//		Real at = (It[0]-0.5)*dt;
//		Real bt = (It[It.size()-1]+0.5)*dt;
//
//		const Cluster& s = src_(FarField[j]);
//		const vectInt& Is = num_(s);
//		Real as = (Is[0]-0.5)*ds;
//		Real bs = (Is[Is.size()-1]+0.5)*ds;
//
//		file << as << "\t" << at << "\n";
//		file << bs << "\t" << at << "\n";
//		file << bs << "\t" << bt << "\n";
//		file << as << "\t" << bt << "\n";
//		file << as << "\t" << at << "\n";
//		file << endl;
//
//	}
//	file.close();
//
//}


void MvProd(vectCplx& f, const HMatrix& A, const vectCplx& x){
	assert(size(f)==size(x)); fill(f,0.);
	
	const vector<LowRankMatrix>&    FarFieldMat  = A.FarFieldMat;
	const vector<SubMatrix>&        NearFieldMat = A.NearFieldMat;
	
	// Contribution champ lointain
	for(int b=0; b<FarFieldMat.size(); b++){
		const LowRankMatrix&  M  = FarFieldMat[b];
		const vectInt&        It = ir_(M);
		const vectInt&        Is = ic_(M);
		
		ConstSubVectCplx xx(x,Is);
		SubVectCplx ff(f,It);
		MvProd(ff,M,xx);
	}
	
	// Contribution champ proche
	for(int b=0; b<NearFieldMat.size(); b++){
		const SubMatrix&  M  = NearFieldMat[b];
		
		const vectInt& It = ir_ (M);
		const vectInt& Is = ic_ (M);
		
		ConstSubVectCplx xx(x,Is);
		SubVectCplx ff(f,It);
		MvProd(ff,M,xx);
	}
	
}

Real squared_absolute_error(const HMatrix& B, const Matrix& A){
	Real err = 0;
	for(int j=0; j<B.FarFieldMat.size(); j++){
		SubMatrix subm(B.mat,ir_(B.FarFieldMat[j]),ic_(B.FarFieldMat[j]));
		err += squared_absolute_error(B.FarFieldMat[j], subm);
	}
	return err;
}



#endif
