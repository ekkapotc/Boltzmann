#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>

#include "../../boltzmann.hpp"

using namespace boltzmann;
using namespace std;

void path_calc( int Nmat , int N  , double delta , double * z , active * L , active * L1 , active * L2 ,  double * sigma , double * rho  , int path_i )
{
  double mat_i;
  double mat_j;
  double con1 = 0.0;
  double delta_W = 0.0;
  double mat_end;

  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  double d = 0.2;
  double beta = 0.0;
  double tau = 0.0;

  int indx_n_j = 0;

  std::unique_ptr<double []>  maturity_i (new double[N]);

  active drift = 0.0;
  active value = 0.00;

  for( int i=0 ; i<N ; i++ )
  {
    L2[i] = L[i];
  }

  mat_end = delta*N;
  
  for( int i=0 ; i<N ; i++ )
  {
     mat_i = delta*i;
     
     for( int j=0; j<N ; j++ )
     {
       mat_j = delta*j;
       
       if(mat_i>mat_j){
	  rho[i*N+j] = exp( (0-1) * beta * (mat_i - mat_j) );
       }else{
          rho[i*N+j] = exp( (0-1) * beta * (mat_j - mat_i) );
       }
     }   
  }

  for( int i=0 ; i<N ; i++ )
  {
    maturity_i[i] = delta*i;
  }

  for( int i=0 ; i<N ; i++ )
  {
    tau = maturity_i[i];
    sigma[i] = (a + b * (tau)) * exp((0-1) * c * (tau)) + d;
  }

  for( int i=0 ; i<Nmat ; i++ )
  {
    delta_W = sqrt(delta)*z[i];
    drift = 0.0;
    for( int j=i+1 ; j<N ; j++ )
    {
       mat_j = j*delta;
       indx_n_j = (N-1)*N + j;
       con1 = delta*sigma[j]*sigma[N-1]*rho[indx_n_j];
       drift = drift + (con1*L2[j]) / (1+delta*L2[j]);
       value = exp(drift*delta + ( delta_W-0.5*sigma[j]*delta)*sigma[j]);
       L2[j] = L2[j]*value;
    }
  }

   

  for( int i=0 ; i<N  ; i++ )
  {
    if(path_i)
    {
      L1[i] = L1[i] + L2[i] + 0.00;
    }else{
      L1[i] = L1[i] + L2[i] - 0.05;
    }
  }

  //double mem_after = get_memory();
  //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
}

void portfolio( int Nmat , int N , int Nopt , active * L , double delta , double * swaprates , active & v , int * maturity ,  active * portfolio_val , double * v_observed )
{
  int indx = 0;  

  active v1;
  active v2;
  active m1;
  active m2;
  active val_p;
  active value;
  active s;
  active max_val = 0.0;
  
  //active * portfolio_val = new active[Nmat];

  for( int i=0 ; i<Nmat ; i++ )
  {
    portfolio_val[i] = 0.0;
  }

  v = 0.0;
  m1 = 1.0;

  for( int i=0 ; i<Nmat ; i++ )
  {
    v1 = L[i+indx]*delta;
    m1 = m1 / (v1+1);
  }

  value = m1;
  indx = Nmat;
  
  for( int i=0 ; i<Nopt ; i++ )
  {
    for( int j=0 ; j<maturity[i] ; j++ ){
      int jp = j+1;
      m2 = 1;

      for( int k=0 ; k<jp ; k++ ){
        v2 = L[k+indx]*delta;
        m2 = m2 / (v2+1);
      }

      val_p = m2;

      s = val_p*delta*(L[Nmat+j]-swaprates[i]);

      if( s > 0 )
      {
        max_val = s;
      }else{
        max_val = 0.0;
      }

      //portfolio_val[i] = portfolio_val[i] + max_val*value*100;
      portfolio_val[i] += max_val*value*100;

      //std::cout << "portfolio_val[" << i << "].idx = "  << portfolio_val[i].idx << std::endl;
    }
  }

   for( int i=0 ; i<Nmat ; i++ )
   {
     v = v + (portfolio_val[i] - v_observed[i])*(portfolio_val[i] - v_observed[i]);
   }

   //delete [] portfolio_val;
}

int main( int argc , char ** argv )
{ 
  int npath = atoi(argv[1]);
  int it = atoi(argv[2]);

  
  double delta = 0.25;
  int Nopt = 15;
  int maturity [] = {4,4,4,8,8,8,20,20,20,28,28,28,40,40,40};
  double swaprates [] = {.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055};

  double ** A;
  int Nmat = 40;
  int N = Nmat + 40;
  double result;

  double * sigma = new double[N];

  double * lambda = new double[N];
  
  for(  int i=0 ; i<N ; i++ ) lambda[i] = 0.20;

  double * z = new double[N];
 
  for( int i=0 ; i<N ; i++ ) z[i] = 0.30;

  double * rho = new double[N*N];

  for( int i= 0 ; i<N*N ; i++ ) rho[i] = 0.0;

  double * v_observed = new double[Nmat];

  //declare inputs
  active * L = new active[N];
  //declare output
  active v;

  active * L1 = new active[N];
  active * L2 = new active[N];  

  active * portfolio_val = new active[Nmat];
	
  initialize(N,1,it*4331040);
  
  set_probe_frequency(5);

  for( int i=0 ; i<N ; i++ )  L[i] = 0.05;

  for( int i=0 ; i<N ; i++ )  independent(L[i]);

  //int t=0;

  while(checkpoint(L,v))
  {
     try{
       /*** function body ***/

      for( int i=0 ; i<N ; i++ ) L1[i] = L[i];//L will stay the same

      for( int i=0 ; i<npath ; i++ ) path_calc(Nmat,N,delta,z,L,L1,L2,sigma,rho,i);

      //double mem_before = get_memory();

      //find the averages for each one
      for( int i=0 ; i<N ; i++ ) L1[i] = L1[i]/npath;

       portfolio(Nmat,N,Nopt,L1,delta,swaprates,v,maturity,portfolio_val,v_observed);
  
       //double mem_after = get_memory();
       //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;

       /*** function body ***/
    }catch(BreakException)
    {
      //std::cout << t+1 << "th iteration..." << std::endl;
      //t++;
    }
 }
  
  dependent(v);

  //std::cout << "cost = " << reverse_elimination() << std::endl;
	
  harvest( 1 , N , A );

  finalize();

  //std::cout << "result = " << v.val << std::endl;

  delete [] rho;
  delete [] lambda;
  delete [] z;
  delete [] L1;
  delete [] L2;
  delete [] sigma;
  delete [] portfolio_val;
  delete [] v_observed;
  delete [] L;
  
  return 0;
}

