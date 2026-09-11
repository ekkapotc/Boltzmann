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

int main( int argc , char ** argv )
{
  f007();
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
 
  active a = min(x[0],x[1]);
  cout << "a.idx = " << a.idx << endl;
  active b = max(x[0],x[1]);
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
 
  active a = 5;
  cout << "a.idx = " << a.idx << endl;
  x[0] += 2.0*a;//x[0] = 12
  cout << "x[0].idx = " << x[0].idx << endl;
  x[1] = a+min(x[0],x[1]); //x[1] = 5+max(12,3)=5+3 = 8
  cout << "x[1].idx = " << x[1].idx << endl;
  y = sin(x[0]*x[1]);
  cout << "y.idx " << y.idx << endl;
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

  for( int i=0 ; i<2 ; i++ )
  {
    independent(x[i]);  
  }
 
  active a = 5;
  cout << "a.idx = " << a.idx << endl;
  active b = max(2,a);
  cout << "b.idx = " << b.idx << endl;
  active c = x[0]+b*x[1];
  cout << "c.idx = " << c.idx << endl;  
  y = c;
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
 
  active a = max(x[0],5);
  cout << "a.idx " << a.idx << endl;
  active b = a*x[0]+2;
  cout << "b.idx " << b.idx << endl;  
  y = x[1]*b;
  cout << "y.idx " << y.idx << endl;
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
 
  active a = x[0]+max(2,x[1]);//a = 2+max(2,3) = 2+3 = 5
  cout << "a.idx = " << a.idx << endl;
  cout << "a = " << a.val << endl;
  active b = max(a,max(x[0],x[1]+10));//b = max(5,max(2,3+10)) = max(5,13) = 13
  cout << "b.idx = " << b.idx << endl;
  cout << "b = " << b.val << endl;
  y = 2*a+3*b;
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
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
 
  active a = +x[0];//a = 2.0
  cout << "a.idx = " << a.idx << endl;
  active b = max(x[0],x[1]);//b = 3.0 
  cout << "b.idx = " << b.idx << endl;
  y = a*b;//2*3
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
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
 
  active a = 5;
  cout << "a.idx = " << a.idx << endl;
  active b = a+3;
  cout << "b.idx = " << b.idx << endl;
  active c = min(a,max(min(a,2),b));//both passive
  cout << "c.idx = " << c.idx << endl;
  y = c+3*x[0]+2*x[1];  
  cout << "y.idx = " << y.idx << endl;
  cout << "y = " << y.val << endl;
  dependent(y);

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

