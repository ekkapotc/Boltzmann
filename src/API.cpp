#include "../inc/API.hpp"
#include "../inc/Vertex.hpp"
#include "../inc/Edge.hpp"
#include "../inc/Process.hpp"
#include "../inc/TapeState.hpp"
#include "../inc/MpiCheck.hpp"
#include "../inc/BreakException.hpp"

#include <mpi.h>
#include <sys/time.h>//gettimeofday, for get_wall_time(); kept out of the public header
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>//atexit
#include <cmath>
#include <ctime>

using namespace boltzmann;
using namespace boltzmann::internals;//in a .cpp, not a header

namespace boltzmann
{
  namespace internals
  {
    largeint vertex_counter = 0;
    largeint edge_counter   = 0;
  }
}


namespace boltzmann
{
  namespace internals
  {
    static inline void restore_values( active * x , active & y );
    static inline void restore_values( active ** x , active & y );
    static inline void restore_values( active * x , active * y );
    static inline void restore_values( active ** x , active ** y );
  }
}


/* -------------------------------------------------------------- the tape --
 *
 * SVEGP-26 : what used to be eleven file-statics sitting beside a Process
 * singleton -- two global states that had to agree with each other -- are
 * TapeState members now, reached through one thread_local pointer.
 */
namespace boltzmann
{
namespace internals
{

static thread_local TapeState * g_current_tape = NULL;

TapeState * current_tape()
{
  return g_current_tape;
}

void set_current_tape( TapeState * t )
{
  g_current_tape = t;
}

TapeState::TapeState( largeint n , largeint m , largeint mem ):
independent_size(n),
dependent_size(m),
memory_available(mem),
indep_rows(0),
indep_cols(0),
dep_rows(0),
dep_cols(0),
run_counter(0),
indep_shadow_copy(new double[size_t(n?n:1)]),
dep_shadow_copy(new double[size_t(m?m:1)])
{
  /*
   * One communicator per tape, duplicated from MPI_COMM_WORLD so the ring
   * traffic cannot be confused with the caller's.  The Process owns it and
   * this object frees it.  It used to be a file-global: a second
   * initialize() overwrote the first tape's handle, and a later finalize()
   * freed it out from under the first tape.
   */
  MPI_Comm dup = MPI_COMM_NULL;
  BZ_MPI( MPI_Comm_dup(MPI_COMM_WORLD,&dup) );

  /*
   * Hand errors back instead of aborting from inside the call.  Without this,
   * every BZ_MPI() in the library is dead code: the default handler is
   * MPI_ERRORS_ARE_FATAL and the return value never carries a failure.  Set
   * on our own duplicate only -- MPI_COMM_WORLD belongs to the caller, and a
   * library has no business changing how the caller's own communications
   * fail.
   */
  BZ_MPI( MPI_Comm_set_errhandler(dup,MPI_ERRORS_RETURN) );

  proc.initialize(n,m,mem,dup);
}

TapeState::~TapeState()
{
  delete [] indep_shadow_copy;
  delete [] dep_shadow_copy;

  if(proc.comm!=MPI_COMM_NULL){
    BZ_MPI( MPI_Comm_free(&proc.comm) );//one per tape; this used to leak one per call
    proc.comm = MPI_COMM_NULL;
  }
}

}//end of namespace internals
}//end of namespace boltzmann

namespace boltzmann//process-wide, and correctly so: MPI's lifecycle is not a tape's
{
  /*
   * THE LIBRARY DOES NOT OWN THE MPI LIFECYCLE.
   *
   * initialize() called MPI_Init unconditionally and finalize() called
   * MPI_Finalize, which had two consequences neither of them documented.
   *
   *   1  A PROCESS COULD OPEN EXACTLY ONE TAPE, EVER.  MPI_Init after
   *      MPI_Finalize is an error, so the second initialize() in a process
   *      aborted the job.  That is why no example in this tree differentiated
   *      more than one function, and why there was no assertion suite: a test
   *      that checks forty-eight elementary functions needs forty-eight
   *      tapes.  The missing tests and that line of code were the same fact.
   *
   *   2  A HOST APPLICATION THAT USES MPI ITSELF HAD ITS MPI TORN DOWN.  An
   *      optimiser that asks this library for a Jacobian and then wants to
   *      MPI_Allreduce the step found MPI already finalized.
   *
   * So MPI comes up on demand, only if nobody else has, and goes down at
   * process exit -- and only if this library was the one that brought it up.
   * These two flags are the only file-statics left, and they are the right
   * shape: they describe the PROCESS, not a tape.
   */
  static bool mpi_started_by_boltzmann = false;
  static bool mpi_shutdown_armed       = false;

  static void boltzmann_mpi_shutdown()
  {
    int done = 0;
    BZ_MPI( MPI_Finalized(&done) );

    if(!done && mpi_started_by_boltzmann){
      BZ_MPI( MPI_Finalize() );
    }
  }

  //brings MPI up if nobody has and arms the exit handler.  Idempotent.
  static void ensure_mpi()
  {
    int up = 0;
    BZ_MPI( MPI_Initialized(&up) );

    if(!up){
      BZ_MPI( MPI_Init(NULL,NULL) );
      mpi_started_by_boltzmann = true;
    }

    if(!mpi_shutdown_armed){
      std::atexit(boltzmann_mpi_shutdown);
      mpi_shutdown_armed = true;
    }
  }

  /*
   * Opening a tape while one is already open used to be silently tolerated:
   * get_proc_instance() handed back the SAME Process, so the second tape
   * inherited the first one's graph and counters, and the file-statics were
   * simply overwritten.  It is diagnosed and replaced now -- a fresh
   * initialize() means a fresh tape -- and replacing it is enough to release
   * it, because ~TapeState() owns everything a tape has.
   */
  static internals::TapeState * open_tape( largeint n , largeint m , largeint mem )
  {
    internals::TapeState * open = internals::current_tape();

    if(open){
      std::cerr << "boltzmann::initialize: a tape is already open -- closing it.  "
                   "Call finalize() before opening another.\n";
      internals::set_current_tape(NULL);
      delete open;
    }

    ensure_mpi();

    internals::TapeState * tp = new internals::TapeState(n,m,mem);
    internals::set_current_tape(tp);
    return tp;
  }
}//end of namespace

void boltzmann::initialize( largeint indep_size , largeint dep_size , largeint memory_available )
{
  open_tape(indep_size,dep_size,memory_available);
}

void boltzmann::initialize( largeint indep_size , largeint indep_rows , largeint indep_cols ,
 			largeint dep_size   , largeint dep_rows   , largeint dep_cols ,
			largeint memory_available )
{
  //the shadow copies are sized by indep_size / dep_size, so the 2-D shapes
  //handed over here must describe exactly those many elements
  assert(indep_rows*indep_cols==indep_size);
  assert(dep_rows*dep_cols==dep_size);

  TapeState * tp = open_tape(indep_size,dep_size,memory_available);

  tp->indep_rows = indep_rows;
  tp->indep_cols = indep_cols;
  tp->dep_rows   = dep_rows;
  tp->dep_cols   = dep_cols;
}

void boltzmann::finalize()
{
  TapeState * tp = current_tape();

  if(!tp){
    std::cerr << "boltzmann::finalize: no tape is open -- nothing to finalize.\n";
    return;
  }

  tp->proc.finalize();

  /*
   * THE LAST CHANCE TO SAY THE JACOBIAN IS WRONG.
   *
   * get_stale_reads() has been readable since the checkpoint contract was
   * enforced, but only a test ever looked.  A nonzero count means the section
   * read an active that did not survive a pass, so a derivative flowed
   * through a constant and the answer that just came out of harvest() is
   * wrong.  That is not something to leave to whoever remembers to ask, and
   * the tape is about to be deleted, so this is the last moment it can be
   * said at all.  The per-read message is printed once; this is printed once
   * per tape, with the total.
   */
  const largeint stale = tp->proc.get_stale_reads();

  if(stale){
    std::cerr << "boltzmann::finalize: this tape completed with " << stale
              << " read(s) of an active that did not survive checkpoint().\n"
              << "                     The Jacobian it produced is WRONG.\n";
  }

  /*
   * NOT MPI_Finalize().  See ensure_mpi() above: the runtime goes down at
   * process exit, and only if this library brought it up.  Ending the tape
   * and ending MPI are different events, and conflating them is what made a
   * second tape impossible.
   *
   * Everything a tape owns -- the graph, both arenas, the shadow copies and
   * the communicator -- goes with this delete.  There is nothing left to
   * reset by hand, which was the point of SVEGP-26: eleven file-statics had
   * to be put back one at a time and forgetting one poisoned the next tape.
   */
  set_current_tape(NULL);
  delete tp;
}

void boltzmann::independent( const active & x )
{
  TapeState * tp = current_tape();

  if(!tp){
    std::cerr << "boltzmann::independent: no tape is open -- call initialize() first.\n";
    return;
  }

  tp->proc.register_indep_vertex(x);

  //x.idx is 1-based; anything outside [1,tp->independent_size] would corrupt the heap
  assert(x.idx>0 && x.idx<=tp->independent_size);

  if(x.idx>0 && x.idx<=tp->independent_size){
    tp->indep_shadow_copy[x.idx-1] = x.val;
  }
}

void boltzmann::dependent( const active & x )
{
  TapeState * tp = current_tape();

  /*
   * Ported from Maxwell (SVEGP-18, SVEGP-26).  These two were the only entry
   * points that dereferenced the tape without checking it, so calling either
   * before initialize() -- or after finalize() -- was a null dereference
   * rather than a diagnostic.  Every other entry point in this file already
   * tested, so the inconsistency was the bug.
   */
  if(!tp){
    std::cerr << "boltzmann::dependent: no tape is open -- call initialize() first.\n";
    return;
  }

  tp->proc.register_dep_vertex(x);
}

void boltzmann::passive_op( const active & x )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.passive_op(x);
  }
}

void boltzmann::destructor( const active & x )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.destructor(x);
  }
  //after boltzmann::finalize() there is no current tape, so an active that
  //outlives the tape destructs harmlessly
}

void boltzmann::unary_op_ass( double dy_dx , const active & x )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.unary_op_ass(dy_dx,x); 
  }
}

void boltzmann::binary_op_ass( const active & x1 , double dy_dx1 , const active & x2, double dy_dx2 )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.binary_op_ass(x1,dy_dx1,x2,dy_dx2);
  }
}

void boltzmann::postfix_op( const active & x1 , const active & x2 )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.postfix_op(x1,x2);
  }
}

void boltzmann::unary_op( const active & x1 ,  double dy_dx1 , const active & x2 , bool overwrite )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.unary_op(x1,dy_dx1,x2,overwrite);
  }
}

void boltzmann::binary_op( const active & x1 , double dy_dx1 , const active & x2 , double dy_dx2  , const active & x3 )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.binary_op(x1,dy_dx1,x2,dy_dx2,x3);
  }
}

void boltzmann::harvest( largeint m , largeint n ,  double **& A )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.harvest(m,n,A);
  }
}

//Maxwell SVEGP-06
void boltzmann::harvest( largeint m , largeint n , double **& A , bool print_out )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.harvest(m,n,A,print_out);
  }
}

/*
 * Maxwell SVEGP-07.  Releases what harvest() handed out.  Tolerates A==NULL,
 * which is what every rank other than the assembling one has -- so a caller
 * can write the same three lines on every rank instead of guarding them.
 */
void boltzmann::free_jacobian( largeint m , double **& A )
{
  if(!A) return;

  for( largeint i=0 ; i<m ; i++ ){
    delete [] A[i];
  }

  delete [] A;
  A = NULL;
}

/*
 * Maxwell SVEGP-07/28 : the same Jacobian as a value.  Built on top of the
 * double** form deliberately -- one assembly path, so the two cannot
 * disagree.  Comes back empty on every rank that did not assemble it, which
 * is the question is_harvesting_rank() answers up front.
 */
Jacobian boltzmann::harvest( largeint m , largeint n , bool print_out )
{
  TapeState * tp = current_tape();

  Jacobian J;

  if(!tp) return J;

  double ** A = NULL;

  tp->proc.harvest(m,n,A,print_out);

  if(!A) return J;//not the assembling rank

  J = Jacobian(m,n);

  for( largeint i=0 ; i<m ; i++ ){
    for( largeint j=0 ; j<n ; j++ ){
      J(i,j) = A[i][j];
    }
  }

  boltzmann::free_jacobian(m,A);

  return J;
}

void boltzmann::set_input_array_dimension( largeint rows , largeint cols )
{
  TapeState * tp = current_tape();

  if(!tp){
    std::cerr << "boltzmann::set_input_array_dimension: no tape is open.\n";
    return;
  }

  tp->indep_rows = rows;
  tp->indep_cols = cols;
}

void boltzmann::set_output_array_dimension( largeint rows , largeint cols )
{
  TapeState * tp = current_tape();

  if(!tp){
    std::cerr << "boltzmann::set_output_array_dimension: no tape is open.\n";
    return;
  }

  tp->dep_rows = rows;
  tp->dep_cols = cols;
}

/* ----------------------------------------------------- the checkpoint loop --
 *
 * ONE IMPLEMENTATION, FOUR SIGNATURES.
 *
 * These four overloads were four copies of the same forty lines, differing
 * only in how the dependents are walked.  That is where SVEGP-17 lived: the
 * 2-D form indexed the shadow copy with the wrong dimension, and because no
 * other copy had a 2-D index there was nothing for it to disagree with.  The
 * shape-specific part is now three tiny helpers, each written once, and the
 * protocol is written once.
 *
 * The order inside is not free to change.  The dependents are saved before
 * finalize() because finalize() is what settles the partition count, and they
 * are restored after terminate() says the loop is over -- the caller gets its
 * dependent back with the idx it had when the tape closed, which is what the
 * Jacobian is harvested against.
 */
namespace
{
  using namespace boltzmann;
  using boltzmann::internals::TapeState;

  /* Save each dependent's value and index into the tape's shadow copy, and
   * tell the Process which vertex index the dependent sits on. */
  inline void save_dependents( TapeState * tp , active & y )
  {
    y.old_idx = y.idx;
    tp->dep_shadow_copy[0] = y.val;
    tp->proc.save_dependent_index(y.idx);
  }

  inline void save_dependents( TapeState * tp , active * y )
  {
    for( largeint k=0 ; k<tp->dependent_size ; k++ ){
      y[k].old_idx = y[k].idx;
      tp->dep_shadow_copy[k] = y[k].val;
      tp->proc.save_dependent_index(y[k].idx);
    }
  }

  inline void save_dependents( TapeState * tp , active ** y )
  {
    for( largeint i=0 ; i<tp->dep_rows ; i++ ){
      for( largeint j=0 ; j<tp->dep_cols ; j++ ){
        y[i][j].old_idx = y[i][j].idx;
        tp->dep_shadow_copy[tp->dep_cols*i+j] = y[i][j].val;
        tp->proc.save_dependent_index(y[i][j].idx);
      }
    }
  }

  /* Put them back as the loop ends.  The stride here and the stride above are
   * now adjacent and obviously the same expression; they were forty lines and
   * one function apart. */
  inline void restore_dependents( TapeState * tp , active & y )
  {
    y.idx = y.old_idx;
    y.val = tp->dep_shadow_copy[0];
  }

  inline void restore_dependents( TapeState * tp , active * y )
  {
    for( largeint k=0 ; k<tp->dependent_size ; k++ ){
      y[k].idx = y[k].old_idx;
      y[k].val = tp->dep_shadow_copy[k];
    }
  }

  inline void restore_dependents( TapeState * tp , active ** y )
  {
    for( largeint i=0 ; i<tp->dep_rows ; i++ ){
      for( largeint j=0 ; j<tp->dep_cols ; j++ ){
        y[i][j].idx = y[i][j].old_idx;
        y[i][j].val = tp->dep_shadow_copy[tp->dep_cols*i+j];
      }
    }
  }

  /*
   * Y BY REFERENCE, NOT BY VALUE, AND THE SUITE SAID SO IMMEDIATELY.
   *
   * Written as `Y y` first.  For the scalar overload Y deduces to `active`,
   * so the dependent was COPIED -- the protocol then saved and restored the
   * copy while the caller's y went untouched, and copying an active records
   * an operation on the tape as a side effect.  Every derivative came out
   * 0.0: fdcheck 36 failures, regress 45.  `Y & y` deduces active& for the
   * scalar, active*& and active**& for the array forms, which is what the
   * four hand-written bodies had.
   */
  template<class X, class Y>
  bool checkpoint_impl( X x , Y & y )
  {
    TapeState * tp = current_tape();

    /*
     * Maxwell SVEGP-26 : with no tape open this used to return true on the
     * first call -- tp->run_counter was 0 -- run the caller's section with
     * nothing recording, and then dereference a null tape on the second.  End
     * the loop before it starts instead.
     */
    if(!tp){
      std::cerr << "boltzmann::checkpoint: no tape is open -- call initialize() first.\n";
      return false;
    }

    if(!tp->run_counter){
      tp->run_counter++;
      return true;
    }

    if(tp->run_counter==1){
      save_dependents(tp,y);
      tp->proc.finalize();//top_owner_idx has been determined
      tp->proc.disable_profiling();
      /*
       * This used to print unconditionally, from inside checkpoint(), on
       * every rank -- so a p-rank job printed it p times, in the middle of
       * whatever the caller was doing, and an optimiser calling the library
       * in a loop printed it once per iteration per rank.  boltzmann::
       * get_partitions() is the way to ask (Maxwell SVEGP-22); a library
       * does not narrate.
       */
    }else if(tp->run_counter==2){
      tp->proc.max_rank_check_memory();
    }

    internals::restore_values(x,y);

    tp->proc.reinitialize();

    if(!(tp->proc.terminate())){
      tp->run_counter++;
      return true;
    }

    restore_dependents(tp,y);
    tp->proc.ignore_vertices();//make the destructor skip its vertex check
    return false;
  }
}

bool boltzmann::checkpoint( active * x , active & y )
{
  return checkpoint_impl(x,y);
}

bool boltzmann::checkpoint( active ** x , active & y )
{
  return checkpoint_impl(x,y);
}

bool boltzmann::checkpoint( active * x , active * y )
{
  return checkpoint_impl(x,y);
}

bool boltzmann::checkpoint( active ** x , active ** y )
{
  return checkpoint_impl(x,y);
}

int boltzmann::MPI_rank()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.MPI_rank();
  }

  return -1;
}

int boltzmann::MPI_size()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.MPI_size();
  }

  return -1;
}

bool boltzmann::is_harvesting_rank()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.is_final_rank_public();
  }

  return false;
}

largeint boltzmann::get_memory()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.get_memory();
  }

  return 0;
}

largeint boltzmann::get_heap()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.get_heap();
  }

  return 0;
}

//Maxwell SVEGP-22
largeint boltzmann::get_partitions()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.get_partitions();
  }

  return 0;
}

largeint boltzmann::get_stale_reads()
{
  TapeState * tp = current_tape();
  if(tp){ return tp->proc.get_stale_reads(); }
  return 0;
}

largeint boltzmann::get_total_partitions()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.get_total_partitions();
  }

  return 0;
}

largeint boltzmann::get_cost()
{
  TapeState * tp = current_tape();

  if(tp){
    return tp->proc.get_cost();
  }

  return 0;
}

//Maxwell SVEGP-32
void boltzmann::set_break_mode( break_t mode )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.set_break_mode(mode);
  }
}

//parallel-only; see pass_t in Typedefs.hpp
void boltzmann::set_pass_mode( pass_t mode )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.set_pass_mode(mode);
  }
}

double boltzmann::get_wall_time()
{
  struct timeval time;

  if(gettimeofday(&time,NULL)){
    return 0;
  }

  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double boltzmann::get_cpu_time()
{
  return (double)clock()/CLOCKS_PER_SEC;
}

void boltzmann::set_probe_frequency( int probe_freq )
{
  TapeState * tp = current_tape();

  if(tp){
    tp->proc.set_probe_frequency(probe_freq);
  }
}

void boltzmann::internals::restore_values( active * x , active & y )
{
  TapeState * tp = current_tape();

  /*
   * The pass boundary.  Everything restored below belongs to the new pass;
   * everything NOT restored below does not, and active::gen is how the library
   * tells the difference instead of dereferencing a freed vertex.
   */
  const largeint g = tp->proc.advance_pass();

  //y should be reinitialized before x because y could use the same program variable as one of x
  y.reachable = false;
  y.idx = 0;
  y.owner_idx = 0;
  y.val = 0.0;
  y.vtx = NULL;
  y.gen = g;

  for( largeint i=0 ; i<tp->independent_size ; i++ ){
    x[i].reachable = true;
    x[i].idx = i+1;
    x[i].owner_idx = 0;
    x[i].val = tp->indep_shadow_copy[i];
    x[i].vtx = NULL;
    x[i].gen = g;
  }
}

void boltzmann::internals::restore_values( active ** x , active & y )
{
  TapeState * tp = current_tape();

  /*
   * The pass boundary.  Everything restored below belongs to the new pass;
   * everything NOT restored below does not, and active::gen is how the library
   * tells the difference instead of dereferencing a freed vertex.
   */
  const largeint g = tp->proc.advance_pass();

  //y should be reinitialized before x because y could use the same program variable as one of x
  y.reachable = false;
  y.idx = 0;
  y.owner_idx = 0;
  y.val = 0.0;
  y.vtx = NULL;
  y.gen = g;

  for( largeint i=0 ; i<tp->independent_size ; i++ ){
    (*x[i]).reachable = true;
    (*x[i]).idx = i+1;
    (*x[i]).owner_idx = 0;
    (*x[i]).val = tp->indep_shadow_copy[i];
    (*x[i]).vtx = NULL;
    (*x[i]).gen = g;
  }
}

void boltzmann::internals::restore_values( active * x , active * y )
{
  TapeState * tp = current_tape();

  /*
   * The pass boundary.  Everything restored below belongs to the new pass;
   * everything NOT restored below does not, and active::gen is how the library
   * tells the difference instead of dereferencing a freed vertex.
   */
  const largeint g = tp->proc.advance_pass();

  for( largeint i=0 ; i<tp->dependent_size ; i++ ){
    y[i].reachable = false;
    y[i].idx = 0;
    y[i].owner_idx = 0;
    y[i].val = 0.0;
    y[i].vtx = NULL;
    y[i].gen = g;
  }

  for( largeint i=0 ; i<tp->independent_size ; i++ ){
    x[i].reachable = true;
    x[i].idx = i+1;
    x[i].owner_idx = 0;
    x[i].val = tp->indep_shadow_copy[i];
    x[i].vtx = NULL;
    x[i].gen = g;
  }
}

void boltzmann::internals::restore_values( active ** x , active ** y )
{
  TapeState * tp = current_tape();

  /*
   * The pass boundary.  Everything restored below belongs to the new pass;
   * everything NOT restored below does not, and active::gen is how the library
   * tells the difference instead of dereferencing a freed vertex.
   */
  const largeint g = tp->proc.advance_pass();

  for( largeint i=0 ; i<tp->dep_rows ; i++ ){
    for( largeint j=0 ; j<tp->dep_cols ; j++ ){
      y[i][j].reachable = false;
      y[i][j].idx = 0;
      y[i][j].owner_idx = 0;
      y[i][j].val = 0.0;
      y[i][j].vtx = NULL;
      y[i][j].gen = g;
    }
  }

  for( largeint i=0 ; i<tp->indep_rows ; i++ ){
    for( largeint j=0 ; j<tp->indep_cols ; j++ ){
      x[i][j].reachable = true;
      x[i][j].idx = tp->indep_cols*i+j+1;//row-major; indices start from 1
      x[i][j].owner_idx = 0;
      x[i][j].val = tp->indep_shadow_copy[tp->indep_cols*i+j];
      x[i][j].vtx = NULL;
      x[i][j].gen = g;
    }
  }
}


