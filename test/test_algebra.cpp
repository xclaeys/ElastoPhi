#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include "point.hpp"
#include "matrix.hpp"
using namespace std;

int main(){
	R3 x,y,z;
	x[0]=0;
	x[1]=1;
	x[2]=2;
	y=2*x;
  
  	z=x+y;
	
  	assert(z[0]==3*x[0]);


}
