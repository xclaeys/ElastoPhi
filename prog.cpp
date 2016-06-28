#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
#include "lrmat.hpp"
#include "cluster.hpp"
#include "hmatrix.hpp"

#include <stdlib.h> 
#include <time.h>   


using namespace std;

int main(){


  /* ==========================================
     const int N = 10;
     Matrix A(N,N);
     for(int j=0; j<N; j++){
     for(int k=0; k<N; k++){
     A(j,k) = 1./(1+abs(j+k));      
     }
     }
     cout << "A:\n" << A << endl << endl;
     
     vectInt I(3), J(3);
     I[0]=5; I[1]=6; I[2]=7;
     J[0]=5; J[1]=6; J[2]=7;
     LowRankMatrix B(SubMatrix(A,I,J));
     
     cout << "SubMatrix(A,I,J):\n" << SubMatrix(A,I,J) << endl << endl;
     
     vectCplx u(N); for(int j=0; j<N; j++){u[j]=j+1;}
     vectCplx v1(3); for(int k=0; k<3; k++){v1[k]=J[k]+1;}
     SubVectCplx v2(u,I);
     vectCplx g(N); 
     
     SubVectCplx f(g,I);
     //  vectCplx f(3);
     MvProd(f,B,v2);
     cout << "v2:\t" << v2 << endl;
     cout << "g:\t"  << g << endl;  
     
     fill(g,0.);
     
     cout << endl;
     MvProd(f,B,v1);
     cout << "v1:\t" << v1 << endl;
     cout << "g:\t"  << g << endl;  
     ========================================== */
  
  
  /*  
      const int N = 20;
      vectCplx U(N);
      for(int j=0; j<N; j++){
      U[j] = 1.;}
      
      vectInt I(4);
  I[0]=2; I[1]=3; I[2]=4; I[3]=5;
  
  SubVectCplx V(U,I);
  V = 2.;

  cout << "U:\t"     << U     << endl;
  cout << "(V,V):\t" << (V,V) << endl;  
  */  
  
  
  /*
  const int N = 10;
  Matrix A(N,N);
  for(int j=0; j<N; j++){
    for(int k=0; k<N; k++){
      A(j,k) = 1./(1+abs(j+k));      
      
    }
  }

  cout << "A:\n" << A << endl;
  
  vectInt I(4), J(4);
  I[0]=0; I[1]=1; I[2]=2; I[3]=3;
  J[0]=0; J[1]=1; J[2]=2; J[3]=3;
  
  SubMatrix B(A,I,J);
  cout << "B:\n" << B << endl;

  Matrix C = B;
  cout << "C:\n" << C << endl;
  */  
  

  
  /* =============================
     const int N = 100;
     Matrix A(N,N);
     for(int j=0; j<N; j++){
     for(int k=0; k<N; k++){
     A(j,k) = 1./(1+abs(j+k));      
     
     }
     }
     
     // vectReal SgVal = SVD(A);
     // cout << "SgVal:\t" << SgVal << endl;
     

     
     srand (time(NULL));
     //int I = rand()%N;
     vectCplx u(N);
     //u[I] = 1.;
     for(int j=0; j<N;j++){
     int I = rand()%N;
     Real val = Real(I)/Real(N);
     u[j] = val;
     }
     
     vectCplx ua=A*u, ub=B*u;
     cout << "rang de B:\t" << rank_of(B) << endl;
     cout << "norm(ua-ub)/norm(ua):\t" << norm(ua-ub)/norm(ua) << endl;
     ============================= */  
  
  // Cluster cl(x);
  // cout << "cl:\t" << cl << endl;
  // DisplayTree(cl);

  
  const int N = 4; const int N2 = N*N;
  const Real dx = 1./Real(N2-1);
  vectR3 x(N2);
  for(int j=0; j<N2; j++){
    x[j][0] = j*dx;}
  
  Matrix A(N2,N2);
  for(int j=0; j<N2; j++){
    for(int k=0; k<N2; k++){
      //      A(j,k) = 1./(1+abs(j-k));
      A(j,k) = 1./(1+norm(x[j]-x[k]));
    }
  }
  
  HMatrix B(A,x,x);
  DisplayPartition(B,"PlotPartition.txt");
  
  for(int j=0; j<N2; j++){
    vectCplx u(N2,0.), f(N2,0.);
    u[j] = 1.;
    MvProd(f,B,u);
    cout << f << endl;
  }
  
  cout << endl << endl;
  for(int j=0; j<N2; j++){
    vectCplx u(N2,0.), f(N2,0.);
    u[j] = 1.;
    f = A*u;
    cout << f << endl;
  }
  

  
  /* =======================
     int N2 = 6;
     vectInt I(N2),J(N2);
     for(int j=0; j<N2; j++){I[j]=j;J[j]=j;}
     SubMatrix B(A,I,J);
     cout << "B:\n" << B << endl;
     
     LowRankMatrix D(N2,B);
     
     cout << endl << endl;
     cout << "B*u:\t " << endl;
     for(int j=0; j<N2; j++){
     vectCplx u(N2,0.);    
     u[j] = 1;
     cout << B*u << endl;}
     
     cout << endl;
     
     cout << "D*u:\t " << endl;
     for(int j=0; j<N2; j++){
     vectCplx u(N2,0.);    
     u[j] = 1;
     cout << D*u << endl;}
     
     ======================= */
  

  
  
  
     /* =================================
	int N=10; 
	int Npt = N;
	cout << "Npt:\t" << Npt  << endl;
	
	vectR3 x(Npt);
	Real dx = 1./Real(N-1);
	for(int j0=0; j0<N; j0++){
	//    for(int j1=0; j1<N; j1++){
	//    int I = j0+N*j1; 
	x[j0][0] = j0*dx;
	//  x[I][1] = j1*dx;
	// }
	}
     
     vectInt num(Npt);  
     for(int j=0; j<Npt; j++){num[j]=j;}
     
     Cluster root(x);
     root.push_back(num);
     root.Build();
     root.NearFieldBall(0.5*dx);
     BlockTree BT; BT.Build(root,root);
     cout << BT << endl;
     ================================ */  
  
  //  root.NearFieldBall(dx);
  
  
  /* =================================
     vectR3 x(8);
     x[1][0] = 2.;
     x[2][1] = 2.;
     x[3][0] = 2.; x[3][1] = 2.;
     R3 ez; ez[2]=1.;
     x[4] = x[0]+2*ez;
     x[5] = x[1]+2*ez;
     x[6] = x[2]+2*ez;
     x[7] = x[3]+2*ez;
     for(int j=0; j<8; j++){num[j]=j;}
     ================================= */
  
  /* =================================
     vector<Box> B; B = SubBox(num,x);
     const int N = 3;
     cout << "B[" << N << "]:\n" << B[N] << endl;
     ================================= */  
  
  //Cluster  tree;
  //OcTree<int>* son = son_(tree,0);
  
  /* =================================
     Matrix B(4,4);
     for(int j=0; j<4; j++){
     for(int k=0; k<4; k++){
     B(j,k) = 1./( 1. + j+k );
     }
     }
     
     cout << endl;
     cout << B << endl;
     
     B(2,3) = -3;
     LowRankMatrix A(4,B);
     
     cout << endl << endl;
     cout << "A*u:\t " << endl;
     for(int j=0; j<4; j++){
     vectCplx u(4,0.);    
     u[j] = 1;
     cout << A*u << endl;}
     
     cout << endl;
     
     cout << "B*u:\t " << endl;
     for(int j=0; j<4; j++){
     vectCplx u(4,0.);    
     u[j] = 1;
     cout << B*u << endl;}
     ================================= */

}
