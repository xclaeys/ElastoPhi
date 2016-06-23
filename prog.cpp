#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "lrmat.hpp"
#include "cluster.hpp"
using namespace std;

int main(){

  const int N = 100;
  Matrix A(N,N);
  for(int j=0; j<N; j++){
    for(int k=0; k<N; k++){
      A(j,k) = 1./(1+abs(j-k));
    }
  }

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
