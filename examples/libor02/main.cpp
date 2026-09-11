#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory>

#include "../../boltzmann.hpp"

using namespace boltzmann;

void path_calc( const int N , const int Nmat , const double delta , active * L , active * L1 , double * lambda , double * z , int path_i )
{
  double sqez,lam,con1;
  active v,vrat;

  std::unique_ptr<active []>  L2 (new active[N]);

  //double mem_before = get_memory();

  for( int i=0 ; i<N ; i++ ){
    L2[i] = L[i];
  }

  for( int n=0 ; n<Nmat ; n++ ){
    sqez = std::sqrt(delta)*z[n];
    
    v = 0.0;
    
    for( int i=n+1 ; i<N ; i++ ) {
      lam = lambda[i-n-1];
      con1 = delta*lam;
      v += (con1*L2[i])/(1.0+delta*L2[i]);
      vrat = exp(con1*v + lam*(sqez-0.5*con1));
      L2[i] = L2[i]*vrat;
    }
  }

  for( int i=0 ; i<N ; i++ ){
    if(path_i){
      L1[i] = L1[i] + L2[i] + 0.00;
    }else{ 
      L1[i] = L1[i] + L2[i]  - 0.05;
    }
  }

  //double mem_after = get_memory();
  //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
}

void portfolio( const int N , const int Nmat , const double delta , const int Nopt ,  int * maturities ,  double * swaprates , active * L , 	active & v )
{
  active b,s,swapval;

  std::unique_ptr<active []>  B (new active[N]);
  std::unique_ptr<active []>  S (new active[N]);

  b = 1.0;
  s = 0.0;
  
  for( int n = Nmat ; n<N ; n++ )
  {
    b  = b / (1.0 + delta*L[n] );
    s  = s + delta*b;
    B[n] = b;
    S[n] = s;
  }

  v = 0.0;

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
       v = v / ( 1.0+delta*L[n] );
     }
  }
}

int main( int argc , char ** argv )
{ 
  double ** A;
  
  int npath = atoi(argv[1]);
  int it = atoi(argv[2]); 

  double delta = 0.25;
  int Nopt = 15;
  int maturities [] = {4,4,4,8,8,8,20,20,20,28,28,28,40,40,40};
  double swaprates [] = {.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055};

  int Nmat = 40;
  int N = Nmat + 40;

  double * lambda = new double[N];
  
  for(  int i=0 ; i<N ; i++ ) lambda[i] = 0.20;

  double * z = new double[N];
 
  for( int i=0 ; i<N ; i++ ) z[i] = 0.30;

  //inputs
  active * L = new active[N];
  //output
  active v;  

  active * L1 = new active[N];
  active * L2 = new active[N];
	
  for( int i=0 ; i<N ; i++ ) L[i] = 0.05;

  initialize(N,1,3988320*it);

  set_probe_frequency(5);

  for( int i=0 ; i<N ; i++ ) independent(L[i]);
	
  while(checkpoint(L,v)){
    try{
      for( int i=0 ; i<N ; i++ )
        L1[i] = L[i];//L will stay the same

      for( int i=0 ; i<npath ; i++ )
        path_calc(N,Nmat,delta,L,L1,lambda,z,i);
  
      //find the averages for each one
      for( int i=0 ; i<N ; i++ )
        L1[i] = L1[i]/npath;

      portfolio(N,Nmat,delta,Nopt,maturities,swaprates,L1,v);
    }catch(BreakException const & e){
    }
  }

  dependent(v);

  harvest( 1 , N , A );//A is an 1 x N matrix //4->5
	
  finalize();
	
  delete [] lambda;
  delete [] z;
  delete [] L1;
  delete [] L2;
  delete [] L;

  return 0;
}
