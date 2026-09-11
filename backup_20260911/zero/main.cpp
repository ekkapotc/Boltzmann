#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "../../boltzmann.hpp"

using namespace boltzmann;

void f001();
void f002();
void f003();
void f004();
void f005();

int main( int argc , char ** argv )
{
  f001();
  return 0;
}

void f001()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 1.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }

  active a = x[0];
  cout << "a.idx = " << a.idx << endl;
  a = sqrt(a-1.0);
  cout << "a.idx = " << a.idx << endl;
  y = a + x[1];  
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

  x[0] = 1.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }

  active a = cbrt(x[0]-1.0);
  cout << "a.idx = " << a.idx << endl;
  active b = cbrt(x[1]);
  cout << "b.idx = " << b.idx << endl;
  y = a+b; 
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f003()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }

  active a = hypot(x[0],x[1]); 
  cout << "a.idx = " << a.idx << endl;
  y = hypot(a,0.0);
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;

  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f004()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }

  active a = 3.0;
  cout << "a.idx = " << a.idx << endl;
  active b = 4.0;
  cout << "b.idx = " << b.idx << endl;
  active c = hypot(a,b);
  cout << "c.idx = " << c.idx << endl;
  y = c*x[0]*x[1];
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f005()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }

  active b = hypot(4.0,x[1]);
  cout << "b.idx = " << b.idx << endl;
  cout << "b = " << b.val << endl;
  y = b*x[0]*x[1];
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}



