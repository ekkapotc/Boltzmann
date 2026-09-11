/*
 * Frozen invariants.
 *
 * Ported from Maxwell, where this is the net under the graph-representation
 * work.  Two kinds of number are pinned, and the difference matters:
 *
 *   STRUCTURAL -- the elimination cost and a bit-exact hash of the harvested
 *       Jacobian.  Properties of the GRAPH and the order it is eliminated in,
 *       so they do not depend on how edges are stored.  Replacing std::map
 *       adjacency with an arena must not move either.  If it does, an edge
 *       was dropped, a fill-in was counted twice, or the accumulation order
 *       changed.
 *
 *   ACCOUNTING -- the partition count.  EXPECTED to move when get_memory()
 *       stops under-counting every edge by its two adjacency entries, because
 *       the same byte budget then buys a different amount of tape.  Pinned so
 *       the move is deliberate rather than noticed later.
 *
 * AND ONE INVARIANT MAXWELL DOES NOT HAVE, WHICH IS WHY THIS FILE IS SHAPED
 * DIFFERENTLY FROM ITS ANCESTOR.  The Jacobian must not depend on how many
 * ranks computed it -- but it CANNOT be bit-exact across rank counts, and
 * finding that out is what this file was for.
 *
 * Measured, np=1 against np=4 on the lattice kernel at a budget that gives 18
 * partitions: the two Jacobians differ by 8.8e-17 relative, one ulp, in some
 * entries and not at all in others.  The np=1 answer is bit-identical to the
 * unchunked answer; the np=4 answer is the one that moves, and every rank
 * count above one agrees with every other.  That is exactly what Maxwell's
 * DESIGN-NOTES 5 predicts: accumulation across eliminations is ordered by
 * ELIMINATION order, and splitting the graph over a pipeline changes which
 * vertex is eliminated where and therefore the order in which fill-in
 * accumulates into a surviving edge.  It is reassociation of a floating-point
 * sum, not a lost or double-counted edge.
 *
 * So there are two different assertions here and conflating them would make
 * one of them a lie:
 *
 *   AGREEMENT (checked at every rank count).  For each kernel, the Jacobian
 *       computed at a tight budget must agree with the same kernel computed
 *       at a budget so large it never chunks, to within a few ulp.  This is
 *       the real correctness statement about the checkpoint machinery, the
 *       ring serialisation and the cross-rank vertex merge, and it holds at
 *       any p.  `make ranks` runs it at 1, 2, 3 and 4.
 *
 *   PINNED (checked at np=1 only).  A bit-exact hash and an exact elimination
 *       cost, which are reproducible only when the partitioning is.  These
 *       are the refactoring change-detector: replacing std::map adjacency
 *       with an arena must not move either by one bit.  They are skipped
 *       above one rank rather than pinned per rank count, because pinning
 *       seven rows times four rank counts would be pinning the schedule, and
 *       the schedule is allowed to change.
 *
 * The hash is a change-detector, not a portable golden value: a different
 * compiler, libm or flag set will move it legitimately.  That is why this
 * compiles with -ffp-contract=off.  Re-pin with
 *
 *     mpirun -np 1 build/bin/invariants pin
 */

#include <mpi.h>

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "../../boltzmann.hpp"

using namespace boltzmann;

/* ------------------------------------------------------------------ hash -- */

static unsigned long long hash_jac( const Jacobian & J )
{
  unsigned long long h = 1469598103934665603ULL;         //FNV-1a, 64-bit

  for( largeint i=0 ; i<J.rows() ; i++ ){
    for( largeint j=0 ; j<J.cols() ; j++ ){
      double d = J(i,j);
      if(d==0.0) d = 0.0;                                //canonicalise -0.0
      unsigned char b[sizeof(double)];
      std::memcpy(b,&d,sizeof(double));
      for( size_t k=0 ; k<sizeof(double) ; k++ ){
        h ^= (unsigned long long)b[k];
        h *= 1099511628211ULL;
      }
    }
  }
  return h;
}

/* -------------------------------------------------------------- scratch --
 *
 * A destructor around `new active[n]`, and it is not decoration.
 *
 * The kernels below allocate working arrays INSIDE the checkpoint section,
 * and under the default BREAK_ON_TARGET that section is abandoned from
 * inside an operator, once per pass.  Written with a bare new[]/delete[]
 * this file leaked 49 432 bytes in 167 allocations at np=1 -- LeakSanitizer,
 * measured, and every byte of it the test's own, none the library's.  That
 * is exactly the failure Maxwell's SVEGP-32 note describes: everything the
 * section allocated before the throw is stranded, once per pass, and the
 * number of passes is the partition count, so tightening the budget silently
 * multiplies the leak.
 *
 * Two ways out, and this file wants the first: hold scratch in something with
 * a destructor, so unwinding cleans up.  The other is set_break_mode
 * (RUN_TO_END), which removes the unwind -- but this file has to exercise the
 * DEFAULT path, so it cannot use it.  examples/hybrdj1 is the one that does.
 */
struct scratch
{
  active * p;

  explicit scratch( largeint n ): p(new active[size_t(n)]) {}

  ~scratch(){ delete [] p; }

  active & operator[]( largeint i ){ return p[size_t(i)]; }

  operator active*(){ return p; }

private:
  scratch( const scratch & );
  scratch & operator=( const scratch & );
};

/* --------------------------------------------------------------- kernels -- */

enum kernel_t
{
  K_LATTICE,   //deep and narrow: a 3-point stencil marched NT steps.  Long
               //critical path, heavy fill-in, large predecessor out-degrees.
  K_WIDE,      //wide and shallow: many independents, one dependent.
  K_CHAIN,     //a single scalar through a long chain of elementals.
  K_NKERNEL
};

static const char * kernel_name( int k )
{
  switch(k){
    case K_LATTICE: return "lattice";
    case K_WIDE:    return "wide";
    case K_CHAIN:   return "chain";
  }
  return "?";
}

//the enumerator, so `invariants pin` emits code that compiles
static const char * kernel_enum( int k )
{
  switch(k){
    case K_LATTICE: return "K_LATTICE";
    case K_WIDE:    return "K_WIDE";
    case K_CHAIN:   return "K_CHAIN";
  }
  return "K_?";
}

struct result_t
{
  unsigned long long hash;
  largeint partitions;        //summed over the communicator
  largeint cost;              //summed over the communicator
  std::vector<double> values; //the Jacobian itself, on every rank
};

static result_t run_kernel( int kernel , largeint budget )
{
  result_t r;
  r.hash = 0; r.partitions = 0; r.cost = 0;

  const int NX = 6;
  const int NT = 8;
  const int NW = 24;
  const int NC = 40;

  largeint n = 0, m = 0;

  if(kernel==K_LATTICE){ n = NX; m = NX; }
  if(kernel==K_WIDE)   { n = NW; m = 1;  }
  if(kernel==K_CHAIN)  { n = 1;  m = 1;  }

  initialize(n,m,budget);

  /*
   * Plain arrays, not std::vector<active>.  Maxwell's DESIGN-NOTES 2 measured
   * a growing std::vector<active> costing a 5.5x tape, because every
   * reallocation copy records an identity vertex and edge.  A sized-once
   * vector does not reallocate and would be safe, but a file whose whole
   * purpose is to pin exact tape sizes should not be the one place that
   * depends on knowing that.
   */
  scratch x(n);
  scratch y(m);

  for( largeint i=0 ; i<n ; i++ ){
    x[size_t(i)] = 0.25 + 0.1*double(i);
    independent(x[size_t(i)]);
  }

  while( checkpoint(x.p,y.p) ){
    try{
      if(kernel==K_LATTICE){
        scratch u(n);
        for( largeint i=0 ; i<n ; i++ ) u[size_t(i)] = x[size_t(i)];
        for( int t=0 ; t<NT ; t++ ){
          scratch v(n);
          for( int i=0 ; i<NX ; i++ ){
            const int im = (i+NX-1)%NX;
            const int ip = (i+1)%NX;
            v[size_t(i)] = 0.5*u[size_t(i)] + 0.25*u[size_t(im)] + 0.25*u[size_t(ip)]
                         + 0.1*sin(u[size_t(i)]);
          }
          for( largeint i=0 ; i<n ; i++ ) u[size_t(i)] = v[size_t(i)];
        }
        for( largeint i=0 ; i<m ; i++ ) y[size_t(i)] = u[size_t(i)];
      }
      else if(kernel==K_WIDE){
        active s = 0.0;
        for( largeint i=0 ; i<n ; i++ ) s += x[size_t(i)]*x[size_t(i)] + sin(x[size_t(i)]);
        y[0] = s;
      }
      else{
        active a = x[0];
        for( int i=0 ; i<NC ; i++ ) a = 0.5*a + sin(a) + 0.25*exp(-a*a);
        y[0] = a;
      }
    }catch( BreakException const & ){}
  }

  for( largeint i=0 ; i<m ; i++ ) dependent(y[size_t(i)]);

  Jacobian J = harvest(m,n);

  unsigned long long h = J.empty() ? 0ULL : hash_jac(J);
  long parts = (long)get_partitions();
  long cost  = (long)get_cost();

  /* the values, broadcast off the assembling rank so that the agreement
   * check below can run on every rank rather than only where they landed */
  std::vector<double> vals(size_t(m*n),0.0);

  if(!J.empty()){
    for( largeint i=0 ; i<m ; i++ )
      for( largeint j=0 ; j<n ; j++ )
        vals[size_t(i*n+j)] = J(i,j);
  }

  finalize();

  /* the hash lives on one rank; the counts are per-rank and are summed, which
   * makes them an assertion about the whole pipeline rather than about
   * whichever rank happened to print */
  unsigned long long hsum = 0;
  long psum = 0, csum = 0;
  MPI_Allreduce(&h,&hsum,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&parts,&psum,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&cost,&csum,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  r.values.assign(size_t(m*n),0.0);
  MPI_Allreduce(&vals[0],&r.values[0],int(m*n),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  r.hash = hsum;//exactly one rank contributed a non-zero
  r.partitions = (largeint)psum;
  r.cost = (largeint)csum;
  return r;
}

/* largest relative difference between two Jacobians of the same shape */
static double max_rel( const std::vector<double> & a , const std::vector<double> & b )
{
  double worst = 0.0;

  if(a.size()!=b.size()) return 1.0e30;

  for( size_t i=0 ; i<a.size() ; i++ ){
    const double d = std::fabs(a[i]-b[i])/(1.0+std::fabs(b[i]));
    if(d>worst) worst = d;
  }

  return worst;
}

/* ------------------------------------------------------------- the table -- */

struct pin_t
{
  int kernel;
  largeint budget;
  unsigned long long hash;
  largeint partitions;
  largeint cost;
};

/* pinned 2026-09-11, g++ 13.3, -ffp-contract=off, after the arena +
 * honest-accounting port.  Re-pin with `invariants pin` at np=1. */
static pin_t PINNED[] = {
  { K_LATTICE ,   800000ULL ,  9419240632285471329ULL ,    1 ,     3078 },
  { K_LATTICE ,     4000ULL ,  9419240632285471329ULL ,   18 ,     3078 },
  { K_LATTICE ,     1000ULL ,  9419240632285471329ULL ,   68 ,     3078 },
  { K_WIDE    ,   800000ULL ,  6447923763444002627ULL ,    1 ,      143 },
  { K_WIDE    ,     2000ULL ,  6447923763444002627ULL ,    8 ,      143 },
  { K_CHAIN   ,   800000ULL , 14109583268233184124ULL ,    1 ,      481 },
  { K_CHAIN   ,     1500ULL , 14109583268233184124ULL ,   35 ,      481 }
};

static const int NPINNED = (int)(sizeof(PINNED)/sizeof(PINNED[0]));

int main( int argc , char ** argv )
{
  const bool pinning = (argc>1 && std::strcmp(argv[1],"pin")==0);

  std::vector<result_t> got((size_t)NPINNED);

  for( int i=0 ; i<NPINNED ; i++ ){
    got[(size_t)i] = run_kernel(PINNED[i].kernel,PINNED[i].budget);
  }

  /*
   * boltzmann::MPI_rank() answers only while a tape is open -- it goes through
   * the Process -- and every tape here is already closed.  Ask MPI, which is
   * still up because finalize() no longer tears it down.
   */
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  if(pinning){
    if(!rank){
      std::printf("static pin_t PINNED[] = {\n");
      for( int i=0 ; i<NPINNED ; i++ ){
        std::printf("  { %-9s , %8lluULL , %20lluULL , %4llu , %8llu }%s\n",
          kernel_enum(PINNED[i].kernel),
          (unsigned long long)PINNED[i].budget,
          got[(size_t)i].hash,
          (unsigned long long)got[(size_t)i].partitions,
          (unsigned long long)got[(size_t)i].cost,
          i+1<NPINNED?",":"");
      }
      std::printf("};\n");
    }
    return 0;
  }

  int bad = 0;

  /* ---- AGREEMENT: every rank count, every budget, against the unchunked
   * answer for the same kernel.  Eight ulp of the largest entry is the
   * tolerance; the measured worst across this table is one. */
  const double TOL   = 8.0*2.220446049250313e-16;
  double       worst = 0.0;

  for( int i=0 ; i<NPINNED ; i++ ){

    /* the unchunked row for this kernel is the first one listed for it */
    int ref = -1;
    for( int k=0 ; k<NPINNED ; k++ ){
      if(PINNED[k].kernel==PINNED[i].kernel){ ref = k; break; }
    }

    const double d = max_rel(got[(size_t)i].values,got[(size_t)ref].values);

    if(d>worst) worst = d;

    if(d>TOL){
      bad++;
      if(!rank) std::printf("  FAIL %-8s budget=%-8llu differs from unchunked by %.3e (> %.3e)\n",
        kernel_name(PINNED[i].kernel),(unsigned long long)PINNED[i].budget,d,TOL);
    }
  }

  /* ---- PINNED: only where the schedule is reproducible, i.e. one rank. */
  if(size==1){
    for( int i=0 ; i<NPINNED ; i++ ){
      const result_t & g = got[(size_t)i];
      const pin_t    & p = PINNED[i];

      /* a zero in the table means "not pinned yet" */
      if(p.hash && g.hash!=p.hash){
        bad++;
        if(!rank) std::printf("  FAIL %-8s budget=%-8llu hash %llu != pinned %llu\n",
          kernel_name(p.kernel),(unsigned long long)p.budget,g.hash,p.hash);
      }
      if(p.partitions && g.partitions!=p.partitions){
        bad++;
        if(!rank) std::printf("  FAIL %-8s budget=%-8llu partitions %llu != pinned %llu\n",
          kernel_name(p.kernel),(unsigned long long)p.budget,
          (unsigned long long)g.partitions,(unsigned long long)p.partitions);
      }
      if(p.cost && g.cost!=p.cost){
        bad++;
        if(!rank) std::printf("  FAIL %-8s budget=%-8llu cost %llu != pinned %llu\n",
          kernel_name(p.kernel),(unsigned long long)p.budget,
          (unsigned long long)g.cost,(unsigned long long)p.cost);
      }
    }
  }

  if(!rank){
    std::printf("invariants: np=%d, %d rows, %d failures, worst chunking error %.3e%s\n",
                size,NPINNED,bad,worst,
                size==1?"  (pinned hash+cost checked)":"  (pinned hash+cost skipped: np>1)");
    for( int i=0 ; i<NPINNED ; i++ ){
      std::printf("    %-8s budget=%-8llu hash=%020llu partitions=%-4llu cost=%llu\n",
        kernel_name(PINNED[i].kernel),(unsigned long long)PINNED[i].budget,
        got[(size_t)i].hash,
        (unsigned long long)got[(size_t)i].partitions,
        (unsigned long long)got[(size_t)i].cost);
    }
  }

  return bad ? 1 : 0;
}
