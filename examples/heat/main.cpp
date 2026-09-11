#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <memory>

#include "../../boltzmann.hpp"

using namespace boltzmann;

//single time step
void ts(
	int & nx, 
        double & delta_t, 
        double & c,
        active * temp) 
{
  std::unique_ptr<active[]> old_temp (new active[nx+1]);
  //active * old_temp = new active[nx+1];
  
  for( int j=0 ; j<=nx ; j++ ){
    old_temp[j] = temp[j];
  }

  for( int j=0 ; j<=nx ; j++ ){
    if(j>0 && j<nx){    
      int jp1=j+1;
      int jm1=j-1;
      temp[j] = old_temp[j]+(nx*nx)*c*delta_t*(old_temp[jp1]-2.0*old_temp[j]+old_temp[jm1]);
    }
  }

  //delete [] old_temp;
}

//time stepping scheme
void tss(
        int & nx,
        int & from,
        int & to,
        int & stride,
        double & delta_t,
        double & c,
        active * temp,
        double * temp_obs,
        active & cost)
{
  for( int i=from ; i<=to ; i++ ){
    //double mem_before = get_memory();
    ts(nx,delta_t,c,temp);
    //double mem_after = get_memory();
    //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
  }
}

void f(
       int nx, 
       int nt, 
       active * temp,
       active & cost
       ) 
{
  int stride=10;
  int zero=0;
  double delta_t=1./nt;
  double c=0.001;
  std::unique_ptr<double[]> temp_obs (new double[nx+1]);
  //double * temp_obs=new double[nx+1]; 

  for(int i=0;i<=nx;i++){
    temp_obs[i] = 2.-i/100.;
  }
  	
  tss(nx,zero,nt,stride,delta_t,c,temp,temp_obs.get(),cost);
  	
  for( int j=1 ; j<nx ; j++ ){
    temp[j] = (temp[j]-temp_obs[j])*(temp[j]-temp_obs[j]);
  }

  for( int j=1 ; j<nx ; j++ ){
    cost += temp[j];
  }

  //delete [] temp_obs;
}

int main( int argc , char ** argv )
{
  int nx = atoi(argv[1]);
  int nt = atoi(argv[2]);
  int it = atoi(argv[3]);
  int probe_freq = atoi(argv[4]);
  
  double ** A;
  
  //int k=0;

  //while(k<5){

  //std::cout << "k = " << k << std::endl;
  active * temp=new active[nx+1];
  active cost;

  temp[0] = 2.0;
 	
  for( int i=1 ; i<=nx ; i++ ){
    temp[i] = 0.0;
  }

  initialize(nx+1,1,(it*215208));//nx=200 //0->1
  //initialize(nx+1,1,(it*431208));//nx=400 //0->1
  //initialize(nx+1,1,(it*863208));//nx=800 //0->1

  set_probe_frequency(probe_freq);

  for( int i=0 ; i<=nx ; i++ ){
    independent(temp[i]);//1->2
  }
	
  while(checkpoint(temp,cost)){
    try{
      f(nx,nt,temp,cost);//2->3
    }catch(BreakException const & e){
    }
  }

  dependent(cost);//3->4

  harvest(1,nx+1,A);//A is an 1 x (nx+1) matrix //4->5
	
  finalize();//5->0
	
  delete [] temp;

  //k++;

  //}

  return 0;
}

