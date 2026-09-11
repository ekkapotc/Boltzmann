#include "../../boltzmann.hpp"


/*
 * The library header no longer says `using namespace boltzmann;` at global
 * scope for you -- a public header that dumps its namespace into every
 * translation unit that includes it is a name-collision generator.  Say it
 * here, or qualify: boltzmann::active, boltzmann::sin.
 */
using namespace boltzmann;
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>

using namespace std;

/**
 * Laplace 2D discretization. n is the dimension of x, 
 * m die number of time iterations
 */
void bratu( int  n , active ** x ) 
{
  const double lambda=1e-2;
  const double h = 1./(n-1);
 
  active ** r = new active*[n];

  for( int i=0 ; i<n ; i++ ){
    r[i] = new active[n];
  }

  //iterate over inner points
  for(int i = 1 ; i < (n-1) ; i++){
    for(int j = 1 ; j < (n-1) ; j++){
        r[i][j] = 0. - ( (x[i+1][j] - 2 * x[i][j] + x[i-1][j] ) / (h*h) )
              - ( (x[i][j+1] - 2 * x[i][j] + x[i][j-1] ) / (h*h) )
              - lambda * exp(x[i][j]);
    }
  }
 
  //updating the inner points
  for (int i = 1 ; i < n-1 ; i++){
    for (int j = 1 ; j < n-1 ; j++){
      x[i][j] = r[i][j];
    }
  }

  for( int i=0 ; i<n ; i++ ){
    delete [] r[i];
  }

  delete [] r;
}

int main(int argc,char** argv) 
{
  int nx = atoi(argv[1]);
  int blocksize = atof(argv[2]);
  double ** A;

  active** x = new active* [nx];

  for(int i = 0; i < nx; i++)
    x[i] = new active[nx];

  for (int i = 1; i < nx-1; i++)  
    for (int j = 1; j < nx-1; j++) 
      x[i][j].val = double(i+j);
    
  // enforce boundary condition
  for (int i = 0 ; i<nx ; i++){
    x[i][0].val = 0.; 
    x[i][nx-1].val = 0.; 
    x[0][i].val = 0.; 
  }

  for (int i = 0 ; i<nx ; i++){
    x[nx-1][i].val = 1.;
  }

  initialize(nx*nx,nx,nx,nx*nx,nx,nx,blocksize);

  for(int i = 0 ; i<nx ; i++){ 
    for(int j = 0 ; j<nx ; j++){
      independent(x[i][j]);
    } 
  }

  while(checkpoint(x,x)){
    try{
      bratu(nx,x);
    }catch(BreakException const & e){

    }
  }
  
  for(int i = 0 ; i<nx ; i++){
    for(int j = 0 ; j<nx ; j++){
      dependent(x[i][j]);
    }
  }
  
  harvest(nx*nx,nx*nx,A);
	
  finalize();

  for(int i = 0; i < nx; i++)
  {
    delete [] x[i];
  }
    
  delete [] x;

  return 0;
}//end of main
