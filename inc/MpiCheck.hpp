#ifndef INCLUDE_MPI_CHECK_HPP
#define INCLUDE_MPI_CHECK_HPP

/*
 * CHECKING MPI RETURN CODES, AND WHY THAT NEEDED TWO CHANGES RATHER THAN ONE.
 *
 * No MPI return code was checked anywhere in this library.  The obvious fix --
 * wrap every call and test the result -- would have changed nothing on its
 * own, because a communicator's default error handler is
 * MPI_ERRORS_ARE_FATAL: the implementation aborts the job from inside the
 * call and the return value is MPI_SUCCESS or nothing at all.  Checking it is
 * dead code until somebody asks MPI to hand errors back.
 *
 * So the library now sets MPI_ERRORS_RETURN on the communicator it owns (see
 * open_tape()) and checks every call through this macro.  What that buys is
 * the message.  Instead of the implementation's own abort -- which names an
 * internal frame, no rank, and nothing about what the library was doing -- a
 * failure reads
 *
 *     boltzmann: MPI error on rank 3, Process.cpp:585
 *                MPI_Recv( ... )
 *                MPI_ERR_TRUNCATE: message truncated
 *
 * and then aborts the job deliberately.  It still aborts: a pipeline that has
 * lost a message has no correct way to continue, and a rank that returned
 * early would hang its neighbour in MPI_Probe rather than fail.  The gain is
 * that you are told which of the library's own calls failed, on which rank,
 * and what MPI called the problem.
 *
 * MPI_Init, MPI_Initialized, MPI_Finalized and MPI_Finalize run before any
 * communicator of ours exists, so those are covered by MPI_COMM_WORLD's
 * handler and checked here only for the diagnostic.
 */
#include <mpi.h>
#include <cstdio>
#include <cstdlib>

namespace boltzmann
{
namespace internals
{

inline void mpi_check( int code , const char * call , const char * file , int line )
{
  if( code==MPI_SUCCESS ) return;

  char text[MPI_MAX_ERROR_STRING];
  int  len = 0;

  if( MPI_Error_string(code,text,&len)!=MPI_SUCCESS || len<=0 ){
    len = 0;
  }
  text[len] = '\0';

  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);//best effort; -1 if even this fails

  std::fprintf(stderr,
               "boltzmann: MPI error on rank %d, %s:%d\n"
               "           %s\n"
               "           %s\n",
               rank,file,line,call,text);
  std::fflush(stderr);

  MPI_Abort(MPI_COMM_WORLD,code);
  std::abort();//not reached; keeps the compiler from assuming a return
}

}//end of namespace internals
}//end of namespace boltzmann

#define BZ_MPI(call) ::boltzmann::internals::mpi_check((call),#call,__FILE__,__LINE__)

#endif
