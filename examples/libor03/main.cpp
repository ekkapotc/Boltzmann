#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include "../../boltzmann.hpp"


/*
 * The library header no longer says `using namespace boltzmann;` at global
 * scope for you -- a public header that dumps its namespace into every
 * translation unit that includes it is a name-collision generator.  Say it
 * here, or qualify: boltzmann::active, boltzmann::sin.
 */
using namespace boltzmann;
using namespace std;

void path_calc( const int N , const int Nmat , const double delta , active* LL , double * lambda , double * z  )
{
  double sqez,lam,con1;
  active v , vrat;  
  
  for( int n=0 ; n<Nmat ; n++ )
  {
    sqez = sqrt(delta)*z[n];
    
    v = 0.0;
    
    for( int i=n+1 ; i<N ; i++ ) 
    {
      lam = lambda[i-n-1];
      con1 = delta*lam;
      v += (con1*LL[i])/(1.0+delta*LL[i]);
      vrat = exp(con1*v + lam*(sqez-0.5*con1));
      LL[i] = LL[i]*vrat;
    }
  }
}

active portfolio( const int N , const int Nmat , const double delta , const int Nopt ,  int * maturities ,  double * swaprates , active * LL )
{
  active b,s,swapval;

  std::unique_ptr<active[]>  B (new active[N]);
  std::unique_ptr<active[]>  S (new active[N]);

  b = 1.0;
  s = 0.0;
  
  for( int n = Nmat ; n<N ; n++ )
  {
    b  = b / (1.0 + delta*LL[n] );
    s  = s + delta*b;
    B[n] = b;
    S[n] = s;
  }

  active v = 0.0;

  for ( int i=0 ; i<Nopt ; i++ )
  {
     int m = maturities[i] + Nmat - 1;
     swapval = B[m] + swaprates[i]*S[m] - 1.0;

     if( swapval < 0 )
     {
       v = v - 100.0*swapval;
     }

     for( int n=0 ; n<Nmat ; n++ )
     {
       v = v / ( 1.0+delta*LL[n] );
     }
  }

  return v;
}

int main(int argc, char* argv[])
{
  double delta = 0.25;
  int Nopt = 15;
  int maturities [] = {4,4,4,8,8,8,20,20,20,28,28,28,40,40,40};
  double swaprates [] = {.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055};

  int npath = atoi(argv[1]);
  int it = atoi(argv[2]);

  int Nmat = 40;
  int N = Nmat + 40;
  double result;

  double * lambda = new double[N];
  
  for(  int i=0 ; i<N ; i++ ) lambda[i] = 0.20;

  double * z = new double[N];
 
  for( int i=0 ; i<N ; i++ ) z[i] = 0.30;

  //inputs
  active* L = new active[N];
  //output
  active v;  
  //local inputs
  active* LL = new active[N];

  for( int i=0 ; i<N ; i++ ) L[i] = 0.05;

  initialize(N,1,4381540*it);
 
  for( int i=0 ; i<N ; i++ )  independent(L[i]);   
  
  while(checkpoint(L,v))
  {
    try
    {
       v = 0.0;   
       for( int i=0 ; i<npath ; i++ )
       {
         //double mem_before = get_memory();
         for( int j=0 ; j<N ; j++ ) LL[j] = L[j];
          path_calc(N,Nmat,delta,LL,lambda,z);
          v +=  portfolio(N,Nmat,delta,Nopt,maturities,swaprates,LL);
	 //double mem_after = get_memory();
         //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
       }
       v = v/npath;
    }catch(BreakException const & e )
    {

    }
  }

  dependent(v);

  double ** A;
	
  harvest( 1 , N , A );

  //std::cout << "result = " << v.val << std::endl;

  finalize();

  delete [] z;
  delete [] lambda;
  delete [] LL;
  delete [] L;

  return 0;
}

