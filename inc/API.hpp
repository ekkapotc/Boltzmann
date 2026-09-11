#ifndef INCLUDE_API_HPP
#define INCLUDE_API_HPP

/*
 * The public free-function surface.
 *
 * Ported from Maxwell's header-hygiene step, and it matters more here.  This
 * header used to pull in Vertex.hpp, Edge.hpp, Process.hpp, <mpi.h> and
 * <sys/time.h>, and then say `using namespace boltzmann;` at global scope.  So a
 * caller who wrote `#include "boltzmann.hpp"` got the whole graph representation,
 * the MPI headers and a POSIX header in their translation unit -- which meant
 * the graph could not be changed without recompiling, and often rewriting,
 * every program that used the library.  It has been changed twice in this
 * round alone (the edge arena, then sorted adjacency) and no caller noticed,
 * which is the point.
 *
 * <mpi.h> in particular is not the caller's business.  A program that only
 * differentiates does not need MPI's declarations, and a program that DOES
 * use MPI should include <mpi.h> itself rather than acquire it by accident
 * from an AD library.
 */
#include "Typedefs.hpp"
#include "Active.hpp"
#include "Jacobian.hpp"

namespace boltzmann
{ 
  void initialize( largeint n , largeint m , largeint mem_avail );

  void initialize( largeint n , largeint nrows , largeint ncols , largeint m , largeint mrows , largeint mcols, largeint mem_avail );

  void finalize();

  void independent( const active & x );

  void dependent( const active & x );

  void passive_op( const active & x );

  void destructor( const active & x );

  void unary_op_ass( double dy_dx , const active & x );

  void binary_op_ass( const active & x1 , double dy_dx1 , const active & x2 , double dy_dx2 );

  void postfix_op( const active & x1 , const active & x2 );

  void unary_op( const active & x1 ,  double dy_dx1 , const active & x2 , bool overwrite );

  void binary_op( const active & x1 , double dy_dx1 , const  active & x2 , double dy_dx2 , const active & x3 );

  void harvest( largeint m , largeint n , double **& A );

  /*
   * Ported from Maxwell (SVEGP-06) : the same call with the FJAC listing
   * suppressed.  harvest() used to print unconditionally, which is unusable
   * inside an optimiser loop -- and on p ranks it is p times as unusable.
   */
  void harvest( largeint m , largeint n , double **& A , bool print_out );

  /*
   * Ported from Maxwell (SVEGP-07) : the documented counterpart of harvest().
   * A is caller-owned and the library previously offered no way to release
   * it, so every harvested Jacobian leaked.
   */
  void free_jacobian( largeint m , double **& A );

  /*
   * Ported from Maxwell (SVEGP-07, SVEGP-28) : the Jacobian as a VALUE.
   * There is nothing to remember to free and copying one is a copy rather
   * than a second owner.  Empty on every rank that did not assemble it --
   * see Jacobian.hpp, which is where the parallel version of that rule is
   * written down.
   */
  Jacobian harvest( largeint m , largeint n , bool print_out=false );

  void set_input_array_dimension( largeint indep_x_dim , largeint indep_y_dim );
  
  void set_output_array_dimension( largeint dep_x_dim , largeint dep_y_dim );

  bool checkpoint( active * x , active & y );

  bool checkpoint( active ** x , active & y );

  bool checkpoint( active * x , active * y );

  bool checkpoint( active ** x , active ** y ); 

  int MPI_rank();

  int MPI_size();

  /*
   * True on the one rank that assembled the Jacobian.  Every harvest-shaped
   * call returns something meaningful only there, and until now a caller had
   * no way to ask except by testing the double** for NULL after the fact.
   */
  bool is_harvesting_rank();

  largeint get_memory();

  /*
   * Ported from Maxwell.  get_memory() is what the BUDGET is spent against --
   * the graph's payload.  get_heap() is what the arenas actually hold.  The
   * ratio is allocator granularity and nothing else; examples/invariants
   * watches it so it cannot drift unnoticed.
   */
  largeint get_heap();

  /*
   * Ported from Maxwell (SVEGP-22) : how many partitions THIS RANK recorded.
   * Nothing exposed this before, so an example could stop exercising the
   * replay path with no way to notice.  Summing it over the communicator
   * gives the tape's total partition count; the spread across ranks is the
   * pipeline's load balance, which nothing could measure before either.
   */
  largeint get_partitions();

  /*
   * How many partitions the whole tape was broken into, across every rank.
   *
   * NEEDS NO COMMUNICATION: the profiling pass is identical on every rank, so
   * every one of them has known this number since before the first productive
   * pass.  get_partitions() is this rank's share of it, and summing that over
   * the communicator gives the same answer the expensive way.
   *
   * Zero until the profiling pass has finished, i.e. until the checkpoint
   * loop has run at least twice.
   */
  largeint get_total_partitions();

  /*
   * How many times this rank's section read an active that did not survive
   * checkpoint().  Nonzero means the Jacobian is wrong and the library said so
   * on stderr; see active::gen.  Zero is the normal state, and the assertion a
   * test should make.
   */
  largeint get_stale_reads();

  //the multiply count of vertex elimination performed on this rank
  largeint get_cost();

  /*
   * Ported from Maxwell (SVEGP-32) : how a pass ends once its target
   * partition is recorded -- by throwing BreakException (BREAK_ON_TARGET, the
   * default and what the library has always done) or by letting the section
   * run to completion (RUN_TO_END).  Same partitions, same Jacobian, same
   * elimination cost.
   *
   * RUN_TO_END is a PREREQUISITE for a section that talks to other ranks --
   * under BREAK_ON_TARGET the rank that throws leaves the collective and the
   * others block in it forever -- but it is not on its own sufficient.  See
   * pass_t in Typedefs.hpp for what else is needed and what is still
   * outstanding.
   *
   * Must be called after initialize() and before the checkpoint loop.
   */
  void set_break_mode( break_t mode );

  /*
   * PARALLEL-ONLY, and not something Maxwell has to think about.  The
   * checkpoint loop is not collective by default -- a rank leaves it as soon
   * as it has recorded its own last partition, so with 8 partitions over 3
   * ranks one rank runs three passes of your section and two run four.  Any
   * collective inside the section is then either a deadlock or, worse, a
   * silent mispairing.
   *
   * PASSES_COLLECTIVE makes every rank run ceil(P/p) passes, so checkpoint()
   * returns the same value on every rank.  Needs no message and changes no
   * Jacobian; it costs at most one extra PASSIVE execution of the section on
   * at most p-1 ranks.
   *
   * It does NOT on its own make a collective inside the section safe, and
   * neither does RUN_TO_END, and neither do both: the library itself blocks
   * in MPI_Probe from inside your operators.  Read pass_t in Typedefs.hpp
   * before putting any MPI in the section.
   *
   * Must be called after initialize() and before the checkpoint loop.
   */
  void set_pass_mode( pass_t mode );

  double get_wall_time();

  double get_cpu_time();

  void set_probe_frequency( int probe_freq );
}

#endif

