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

  x[0] = 2.0;
  x[1] = 3.0;

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }

  active a,b;

  a = x[0]+5.2;
  b = sin(x[1]);
  y = a+b;

  dependent(y);
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

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }

  active a = x[0]+5.2;//point A
  active b = sin(x[1]);//point B
  y = a+b;//point C

  dependent(y);
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
 
  {
    active a = sin(x[0]);
    active b = x[0]+x[1];
    x[0] = a+b;
    x[1] = a+b;
  }
  
  y = x[0]+x[1];
 
  dependent(y);
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
 
  active b = x[0]+x[1];
  cout << "b.idx = " << b.idx << endl;

  cout << "Before entering block" << endl;
  {
    active d = b;//d is reachable from indep but not from dep
    cout << "d.idx = " << d.idx << endl;
    d = b+x[0];
    cout << "d.idx = " << d.idx << endl;
    b = sin(d);
    cout << "b.idx = " << b.idx << endl;
  }
  cout << "After leaving block" << endl;
  active c = b;
  cout << "c.idx = " << c.idx << endl;
 
  y = c;
  cout << "y.idx = " << y.idx << endl;

  dependent(y);
  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}

active f005_a( active x0 , active & x1 )
{
  cout << "Start of f005_a()" << endl;
  cout << "x0.idx = " << x0.idx << endl;
  cout << "x1.idx = " << x1.idx << endl;
  active a = x0+x1;
  cout << "a.idx = " << a.idx << endl;
  active b = 2*a;
  cout << "b.idx = " << b.idx << endl;
  cout << "End of f005_a()" << endl;
  return b;
}

active f005_b( active x0 )
{
  cout << "Start of f005_b()" << endl;
  cout << "x0.idx = " << x0.idx << endl;
  active a = x0;
  cout << "a.idx = " << a.idx << endl;
  active b = 2*a;
  cout << "b.idx = " << b.idx << endl;
  cout << "End of f005_b()" << endl;
  return b;
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
 
  active c;
  c = f005_a(x[0],x[1]);//assignment operator called
  cout << "c.idx = " << c.idx << endl;
  active d = f005_b(2*x[0]);//copy elision
  cout << "d.idx = " << d.idx << endl; 
  y = c+d;
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

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }
 
  active c = x[0]+x[1];
  cout << "c.idx = " << c.idx << endl;
  c*=2;
  cout << "c.idx = " << c.idx << endl;
  y = c;
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

  for( int i=0 ; i<2 ; i++ ){
    independent(x[i]);  
  }
 
  active a;
  cout << "a.idx = " << a.idx << endl;
  active b = a + x[0];//(a+x[0])=>b
  cout << "b.idx = " << b.idx << endl;
  active  c = a + x[1];//(a+x[1])=>c
  cout <<  "c.idx = " << c.idx << endl;
  y = b+c;//(b+c).idx = 5
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
 
  active a;
  cout << "a.idx = " << a.idx << endl;
  a = a + x[0];//(a+x[0])=>b
  cout << "a.idx = " << a.idx << endl;
  active  b = a + x[1];//(a+x[1])=>b
  cout <<  "b.idx = " << b.idx << endl;
  y = (a+b);//(a+b).idx = 5
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
 
  x[0] = x[0]+5.0;
  cout << "x[0].idx = " << x[0].idx << endl;
  active a = x[1] + 5;
  cout << "a.idx = " << a.idx << endl;
  x[1] = a+2;
  cout << "x[1].idx = " << x[1].idx << endl;
  y = x[0] + x[1]; 
  cout << "y.idx = " << y.idx << endl;
  dependent(y);

  cout << "reverse_elimination()"<< endl;
  reverse_elimination(0);
  harvest(2,1,A);
  finalize();
}




