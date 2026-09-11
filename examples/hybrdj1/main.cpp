#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../../boltzmann.hpp"

using namespace boltzmann;
using namespace std;

/**      
 *subroutine fcn for hybrj1 example. 
 */
int fcn(int n, const active *x, active *fvec)
{
  int k;
  double one=1, three=3, two=2, zero=0, four=4;
  active temp, temp1, temp2;

  for (k = 0; k < n; k++)
  {
    //double mem_before = get_memory();

    temp = (three - two*x[k])*x[k];
    temp1 = zero;

    if (k != 0) 
	temp1 = x[k-1];

    temp2 = zero;
    
    if (k != n-1) 
	temp2 = x[k+1];

    fvec[k] = temp - temp1 - two*temp2 + one;

    //double mem_after = get_memory();
    //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
  }

  return 0;//flowing off the end of a non-void function is undefined behaviour:
           //gcc emits ud2 here and the example died with SIGILL the moment
           //the suite actually RAN it.  Nothing ran it before.
}

int main(int argc,char** argv)
{
  double ** A;
  const int n  = atoi(argv[1]);
  const int it = atoi(argv[2]);

  //std::cout << "Runtime Settings : " << std::endl;
  //std::cout << "\tn = " << n << std::endl;
  //std::cout << "\tit = " << it << std::endl;

  initialize(n,n,1024/*it*1656*/);

  /*
   * This example is the living argument for both of the library's answers to
   * the try/catch problem, and it now uses both.
   *
   * Its checkpoint loop below has no try/catch -- the section is a bare call
   * to fcn() -- so under the old protocol the BreakException that
   * check_memory() throws from inside an operator escaped main and the job
   * died with `terminate called after throwing an instance of
   * boltzmann::BreakException`, on every rank, at a 1024-byte budget.  It did
   * exactly that the first time anything in this tree actually ran it.
   *
   * run_tape() below owns the try/catch, so the section cannot escape even
   * under the default BREAK_ON_TARGET.  RUN_TO_END goes further and removes
   * the throw altogether, so the section reaches its own cleanup rather than
   * being unwound out of the middle of an operator.  They are complements:
   * one makes the throw harmless, the other makes it not happen.  The
   * Jacobian is identical with either, both, or neither.
   */
  set_break_mode(RUN_TO_END);

  active* x = new active[n];
  active* fvec = new active[n];

  for( int j=0 ; j<n ; j++ )
  {
    x[j].val = -1;
  }

  for( int j=0; j<n; j++)
  {
    independent(x[j]);
  }

  //no try/catch anywhere in this file, deliberately: run_tape() has it
  run_tape( x , fvec , [&]{ fcn(n, x, fvec); } );

  for( int j=0 ; j<n; j++)
  {
    dependent(fvec[j]);
  }
 
  harvest( n , n , A );
 
  finalize(); 

  return 0;
}
