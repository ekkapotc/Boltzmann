#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "../../boltzmann.hpp"

using namespace boltzmann;

void f001();
void f002();

int main( int argc , char ** argv )
{
  f002();
  return 0;
}

void f001()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = -2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }
 
  y = fabs(2.0*x[0])+max(10.0*x[0],x[1]);
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f002()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = -2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }
 
  y = fabs(x[0]*x[0]+x[1]*x[1]);
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

