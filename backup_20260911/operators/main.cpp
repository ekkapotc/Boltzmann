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
void f006();
void f007();
void f008();
void f009();
void f010();
void f011();

int main( int argc , char ** argv )
{
  f010();
  return 0;
}

void f001()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 2;
  a += x[0];
  cout << "a.idx = " << a.idx << endl;
  active b = 3; 
  b += x[1];
  cout << "b.idx = " << b.idx << endl;
  y = a+b;
  cout << "y.idx = " << y.idx << endl;

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

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  x[0] += x[0];
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] *= x[1];
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0]+x[1];
  cout << "y.idx = " << y.idx << endl;  

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

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
  
  active a = 2;
  a += a;
  cout << "a.idx = " << a.idx << endl;
  y = x[0]+x[1]*a;
  cout << "y.idx = " << y.idx << endl;  

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

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 2;
  cout << "a.idx = " << a.idx << endl;
  x[0] *= a;
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] *= a;
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0]+x[1];
  cout << "y.idx =  " << y.idx << endl;
  
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

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 2;
  cout << "a.idx = " << a.idx << endl;
  x[0] *= a;
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] /= a;
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0]+x[1];
  cout << "y.idx = " << y.idx << endl;

  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f006()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 2;
  cout << "a.idx = " << a.idx << endl;
  x[0] -= a;
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] -= 2*x[0]+0.5;
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0]+x[1];
  cout << "y.idx = " << y.idx << endl;

  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f007()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 2;
  cout << "a.idx = " << a.idx << endl;
  ++a;//a = 3
  cout << "a.idx = " << a.idx << endl;
  ++x[0];
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] *= a;
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0]+x[1];
  cout << "y.idx = " << y.idx << endl;
  
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f008()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 5;
  cout << "a.idx = " << a.idx << endl;
  a++;//a=6
  cout << "a.idx = " << a.idx << endl;
  y = a*x[0]*x[1];
  cout << "y.idx = " << y.idx << endl;
  
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f009()
{
  double ** A;
  initialize(2,1,1024);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = x[0]++;//a = 2 , x[0] = 3
  cout << "a.idx = " << a.idx << endl;
  cout << "x[0].idx = " <<x[0].idx << endl;
  active b = x[1]++;//b = 3 , x[1] = 4
  cout << "b.idx = " << b.idx << endl;
  cout << "x[1].idx = " <<x[1].idx << endl;
  y = a+b+x[0]+x[1];
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f010()
{
  double ** A;
  initialize(2,1,512);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = x[0]+x[1];
  std::cout << "a.idx = " << a.idx << std::endl;
  active b = x[0]-x[1];
  std::cout << "b.idx = " << b.idx << std::endl;
  a = b++;  
  std::cout << "a.val = " << passive_value(a) << std::endl;
  std::cout << "a.idx = " << a.idx << std::endl;
  y = a*b;
  std::cout << "y.val = " << passive_value(y) << std::endl;
  std::cout << "y.idx = " << y.idx << std::endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

void f011()
{
  double ** A;
  initialize(2,1,512);
  set_elim_mode(REVERSE_ELIM);
  
  active x[2];
  active y;

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = x[0]+x[1];
  std::cout << "a.idx = " << a.idx << std::endl;
  active b = x[0]-x[1];
  std::cout << "b.idx = " << b.idx << std::endl;
  a = a++;  
  std::cout << "a.val = " << passive_value(a) << std::endl;
  std::cout << "a.idx = " << a.idx << std::endl;
  y = a*b;
  std::cout << "y.val = " << passive_value(y) << std::endl;
  std::cout << "y.idx = " << y.idx << std::endl;
  dependent(y);  

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

