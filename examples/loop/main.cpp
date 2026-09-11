#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "../../boltzmann.hpp"

using namespace std;
using namespace boltzmann;

void foo( int n , int l , active * x )
{
  active * h = new active[n];

  for( int i = 0 ; i < l ; i++ ){

    //double mem_before = get_memory();

    if(i%2==0){

      h[0] = sin(x[n-1]*x[0]);
      //std::cout << "h[" << 0 << "].idx = " << h[0].idx << std::endl;

      for( int j = 1 ; j < n ; j++ ){
        h[j] = sin(x[j-1]*x[j]);
        //std::cout << "h[" << j << "].idx = " << h[j].idx << std::endl;
      }
    }else{

      x[0] = cos(h[n-1]*h[0]);
      //std::cout << "x[" << 0 << "].idx = " << x[0].idx << std::endl;

      for( int j = 1 ; j < n ; j++ ){
        x[j] = cos(h[j-1]*h[j]);
	//std::cout << "x[" << j << "].idx = " << x[j].idx << std::endl;
      }
    }

    //double mem_after = get_memory();
    //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
  }

  delete [] h;
}

int main( int argc , char ** argv )
{
  int n = atoi(argv[1]);
  int l = atoi(argv[2]);
  int it = atoi(argv[3]);

  //std::cout << "Runtime Settings : " << std::endl;
  //std::cout << "\tnx = " << n << std::endl;
  //std::cout << "\tnt = " << l << std::endl;
  //std::cout << "\tit = " << it << std::endl;

  double ** A;
  active * x = new active[n];
	
  for ( int i = 0 ; i < n ; i++ )
  { 
    x[i].val = (i+1)*0.25;
  }

  initialize(n,n,it*9120);//n=20

  for( int i=0 ; i<n ; i++ )
  {
    independent(x[i]);
  }

  while(checkpoint(x,x)){
    try{
      foo(n,l,x);
    }catch(BreakException const & e){

    }
  }

  for( int i=0 ; i<n ; i++ )
  {
    dependent(x[i]);
  }
	
  harvest(n,n,A);

  finalize();

  delete [] x;
  
  return 0;
}
	
