/*
 * One test per defect found, so that none of them can come back quietly.
 *
 * Ported from Maxwell's examples/regress, which exists because six defects
 * were fixed and nothing in the tree would have noticed if any of them were
 * reintroduced.  The list below is this library's own, from the 2026-09-11
 * review and from the Maxwell adaptation that followed it.
 *
 * The misuse section at the end is Maxwell's SVEGP-26 lesson: deleting a
 * global and letting the compiler find the uses has a blind spot -- a
 * function that used to read a file-static that was simply zero now
 * dereferences a pointer that is null, and that compiles.  So every public
 * entry point is called with no tape open and asserted not to crash.  That
 * is the test that would have caught it, and it did not exist until the bug
 * made the case for it.
 */

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../../boltzmann.hpp"

/* Internal, on purpose: M-10 pins a counter that is not part of the public
 * surface and should not become part of it just to be testable. */
#include "../../inc/TapeState.hpp"

using namespace boltzmann;

static int failures = 0;
static int checks   = 0;
static int rank     = 0;

static void ok( const char * what , bool cond )
{
  checks++;
  if(!cond){
    failures++;
    std::printf("  rank %d  FAIL  %s\n",rank,what);
  }
}

static void close_to( const char * what , double got , double want , double tol=1e-9 )
{
  ok(what, std::fabs(got-want) <= tol*(1.0+std::fabs(want)) );
}

/* ------------------------------------------------------------------------
 * R-01 .. R-06 : the six defects of the 2026-09-11 review.
 * ---------------------------------------------------------------------- */

/* R-01  six elementals had `x2.reachable = x2.reachable`, so asin, acos,
 *       atan, sinh, cosh and tanh silently dropped every derivative below
 *       them.  fdcheck covers all forty-eight; this pins the six that were
 *       actually wrong, and pins that the result is REACHABLE -- the bug was
 *       not a wrong number, it was a vertex that was never created. */
static double d_of( active (*f)( const active & ) , double x0 , largeint budget )
{
  initialize(1,1,budget);
  active x,y; x = x0; independent(x);
  while( checkpoint(&x,y) ){ try{ y = f(x); }catch( BreakException const & ){} }
  dependent(y);
  Jacobian J = harvest(1,1);
  double d = J.empty()?0.0:J(0,0);
  finalize();
  double s=0.0; MPI_Allreduce(&d,&s,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return s;
}

static active f_asin(const active&x){return asin(x);}
static active f_acos(const active&x){return acos(x);}
static active f_atan(const active&x){return atan(x);}
static active f_sinh(const active&x){return sinh(x);}
static active f_cosh(const active&x){return cosh(x);}
static active f_tanh(const active&x){return tanh(x);}
static active f_erf (const active&x){return erf(x);}
static active f_erfc(const active&x){return erfc(x);}

static void r01_reachability()
{
  const double x = 0.4;
  close_to("R-01 asin", d_of(f_asin,x,400),  1.0/std::sqrt(1.0-x*x));
  close_to("R-01 acos", d_of(f_acos,x,400), -1.0/std::sqrt(1.0-x*x));
  close_to("R-01 atan", d_of(f_atan,x,400),  1.0/(1.0+x*x));
  close_to("R-01 sinh", d_of(f_sinh,x,400),  std::cosh(x));
  close_to("R-01 cosh", d_of(f_cosh,x,400),  std::sinh(x));
  close_to("R-01 tanh", d_of(f_tanh,x,400),  1.0-std::tanh(x)*std::tanh(x));
}

/* R-02  erf' used pi where it wanted sqrt(pi), and divided by the exponential
 *       instead of multiplying. */
static void r02_erf()
{
  const double x  = 0.7;
  const double pi = 3.14159265358979323846264338327950288;
  const double d  = (2.0/std::sqrt(pi))*std::exp(-x*x);
  close_to("R-02 erf" , d_of(f_erf ,x,400),  d);
  close_to("R-02 erfc", d_of(f_erfc,x,400), -d);
}

/* R-03/R-04  the 2-D checkpoint overload and restore_values() both flattened
 *       a grid with the wrong stride, so any NON-SQUARE grid walked off the
 *       end of the heap.  The original segfaulted on this; it is two rows by
 *       three columns precisely because square is the case that worked. */
static void r03_nonsquare_grid()
{
  const largeint R=2, C=3;

  initialize(R*C,R,C,R*C,R,C,600);

  /*
   * RUN_TO_END here is not incidental.  The section below allocates its
   * scratch grid with new[], and under BREAK_ON_TARGET the throw abandons it
   * mid-operator once per pass: LeakSanitizer measured 2 880 bytes stranded
   * across nine passes, in the TEST, not the library.  This test is about
   * the 2-D stride, so it takes the mode that lets the section reach its own
   * cleanup and leaves the break-mode question to M-01.
   */
  set_break_mode(RUN_TO_END);

  active ** x = new active*[R];
  for( largeint i=0 ; i<R ; i++ ){
    x[i] = new active[C];
    for( largeint j=0 ; j<C ; j++ ) x[i][j].val = 1.0 + 0.1*double(C*i+j);
  }

  for( largeint i=0 ; i<R ; i++ ) for( largeint j=0 ; j<C ; j++ ) independent(x[i][j]);

  while( checkpoint(x,x) ){
    try{
      active ** r = new active*[R];
      for( largeint i=0 ; i<R ; i++ ) r[i] = new active[C];
      for( largeint i=0 ; i<R ; i++ )
        for( largeint j=0 ; j<C ; j++ )
          r[i][j] = x[i][j]*x[i][j] + 2.0*x[(i+1)%R][j] + sin(x[i][(j+1)%C]);
      for( largeint i=0 ; i<R ; i++ )
        for( largeint j=0 ; j<C ; j++ ) x[i][j] = r[i][j];
      for( largeint i=0 ; i<R ; i++ ) delete [] r[i];
      delete [] r;
    }catch( BreakException const & ){}
  }

  for( largeint i=0 ; i<R ; i++ ) for( largeint j=0 ; j<C ; j++ ) dependent(x[i][j]);

  Jacobian J = harvest(R*C,R*C);

  /* the analytic Jacobian of that one step */
  double v[2][3];
  for( largeint i=0 ; i<R ; i++ ) for( largeint j=0 ; j<C ; j++ ) v[i][j] = 1.0+0.1*double(C*i+j);

  double ref[6][6];
  for( int a=0 ; a<6 ; a++ ) for( int b=0 ; b<6 ; b++ ) ref[a][b] = 0.0;

  for( largeint i=0 ; i<R ; i++ ){
    for( largeint j=0 ; j<C ; j++ ){
      const size_t row = size_t(C*i+j);
      ref[row][C*i+j]           += 2.0*v[i][j];
      ref[row][C*((i+1)%R)+j]   += 2.0;
      ref[row][C*i+((j+1)%C)]   += std::cos(v[i][(j+1)%C]);
    }
  }

  double worst = 0.0;

  if(!J.empty()){
    for( int a=0 ; a<6 ; a++ )
      for( int b=0 ; b<6 ; b++ )
        worst = std::max(worst,std::fabs(J(largeint(a),largeint(b))-ref[a][b]));
  }

  double gworst=0.0;
  MPI_Allreduce(&worst,&gworst,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  finalize();

  for( largeint i=0 ; i<R ; i++ ) delete [] x[i];
  delete [] x;

  ok("R-03 non-square 2-D grid matches the analytic Jacobian", gworst<1e-12);
}

/* R-05  operator>= was declared and never defined: a link failure for any
 *       program that used it.  This function existing is the test. */
static void r05_ge()
{
  active a(2.0), b(1.0);
  ok("R-05 operator>= links and compares", (a>=b) && (a>=1.0) && (2.0>=b) && !(b>=a) );
}

/* R-06  a second tape in one process.  initialize() called MPI_Init and
 *       finalize() called MPI_Finalize, so the second initialize() aborted
 *       the job.  Every test above this line is already a second tape, but
 *       pin it explicitly -- it is the change that made an assertion suite
 *       possible at all. */
static void r06_many_tapes()
{
  double last = 0.0;
  for( int k=0 ; k<8 ; k++ ) last = d_of(f_tanh,0.3,400);
  close_to("R-06 eight tapes in one process", last, 1.0-std::tanh(0.3)*std::tanh(0.3));
}

/* ------------------------------------------------------------------------
 * Adapted from Maxwell: behaviour that must hold.
 * ---------------------------------------------------------------------- */

/* M-01  RUN_TO_END records the same partitions and the same Jacobian as
 *       BREAK_ON_TARGET (Maxwell SVEGP-32). */
static void m01_break_mode()
{
  double a=0.0,b=0.0; largeint pa=0,pb=0;

  for( int mode=0 ; mode<2 ; mode++ ){
    initialize(1,1,400);
    if(mode) set_break_mode(RUN_TO_END);
    active x,y; x = 0.6; independent(x);
    while( checkpoint(&x,y) ){
      try{ active u=x; for(int i=0;i<12;i++) u = 0.5*u + sin(u); y = u; }
      catch( BreakException const & ){}
    }
    dependent(y);
    Jacobian J = harvest(1,1);
    double d = J.empty()?0.0:J(0,0);
    largeint p = get_partitions();
    finalize();
    double s=0.0; MPI_Allreduce(&d,&s,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    long lp=(long)p,tp=0; MPI_Allreduce(&lp,&tp,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    if(mode){ b=s; pb=(largeint)tp; } else { a=s; pa=(largeint)tp; }
  }

  ok("M-01 RUN_TO_END gives the identical Jacobian", a==b);
  ok("M-01 RUN_TO_END records the identical partitions", pa==pb);
}

/* M-02  PASSES_COLLECTIVE equalises the number of passes across ranks
 *       without changing the number of partitions recorded.  This is the
 *       parallel-only defect the port surfaced: at 8 partitions over 3 ranks
 *       one rank ran three passes of the caller's section and two ran four. */
static void m02_pass_mode()
{
  for( int mode=0 ; mode<2 ; mode++ ){

    initialize(1,1,400);
    if(mode) set_pass_mode(PASSES_COLLECTIVE);

    active x,y; x = 0.6; independent(x);

    long passes = 0;
    while( checkpoint(&x,y) ){
      passes++;
      try{ active u=x; for(int i=0;i<12;i++) u = 0.5*u + sin(u); y = u; }
      catch( BreakException const & ){}
    }
    dependent(y);
    Jacobian J = harvest(1,1);
    long parts = (long)get_partitions();
    finalize();

    long lo=0,hi=0,tot=0;
    MPI_Allreduce(&passes,&lo,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&passes,&hi,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&parts ,&tot,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

    if(mode){
      ok("M-02 PASSES_COLLECTIVE: every rank runs the same number of passes", lo==hi);
    }
    ok("M-02 partitions are recorded regardless of pass mode", tot>0);
  }
}

/* M-03  the Jacobian value type is empty on exactly the ranks that did not
 *       assemble it, and is_harvesting_rank() says so BEFORE harvest(). */
static void m03_jacobian_value()
{
  initialize(2,1,400000);
  active x[2], y;
  x[0]=1.3; x[1]=0.7;
  independent(x[0]); independent(x[1]);
  while( checkpoint(x,y) ){ try{ y = x[0]*x[1] + sin(x[0]); }catch( BreakException const & ){} }
  dependent(y);

  const bool claims = is_harvesting_rank();
  Jacobian J = harvest(1,2);
  const bool has = !J.empty();

  int cl = claims?1:0, hs = has?1:0, ncl=0, nhs=0;
  MPI_Allreduce(&cl,&ncl,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&hs,&nhs,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  double d0 = has?J(0,0):0.0, s0=0.0;
  MPI_Allreduce(&d0,&s0,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  finalize();

  ok("M-03 exactly one rank assembles the Jacobian", nhs==1);
  ok("M-03 is_harvesting_rank() agrees with Jacobian::empty()", ncl==nhs);
  close_to("M-03 the assembled value is right", s0, 0.7+std::cos(1.3));
}

/* M-04  free_jacobian() releases the double** form and tolerates the NULL
 *       every non-assembling rank has.  Run under `make sanitize` this is
 *       the leak test. */
static void m04_free_jacobian()
{
  initialize(2,1,400000);
  active x[2], y;
  x[0]=1.3; x[1]=0.7;
  independent(x[0]); independent(x[1]);
  while( checkpoint(x,y) ){ try{ y = x[0]*x[1]; }catch( BreakException const & ){} }
  dependent(y);
  double ** A = NULL;
  harvest(1,2,A,false);
  free_jacobian(1,A);
  ok("M-04 free_jacobian nulls the pointer", A==NULL);
  free_jacobian(1,A);//twice, on purpose
  finalize();
  ok("M-04 free_jacobian tolerates a second call and a NULL", true);
}

/* M-05  run_tape() owns the checkpoint loop's try/catch.  Under the default
 *       BREAK_ON_TARGET a hand-written loop without a catch terminates the
 *       process, so the assertion here is simply that this function returns
 *       -- and that it agrees with the hand-written loop to the last bit. */
static void m05_run_tape()
{
  double hand = 0.0, owned = 0.0;
  largeint passes = 0, parts_hand = 0, parts_owned = 0;

  /* the hand-written loop, catch and all */
  {
    initialize(1,1,400);
    active x,y; x = 0.55; independent(x);
    while( checkpoint(&x,y) ){
      try{ active u=x; for(int i=0;i<14;i++) u = 0.5*u + sin(u)*cos(u); y = u; }
      catch( BreakException const & ){}
    }
    dependent(y);
    Jacobian J = harvest(1,1);
    if(!J.empty()) hand = J(0,0);
    parts_hand = get_partitions();
    finalize();
  }

  /* the same section, no catch anywhere in sight */
  {
    initialize(1,1,400);
    active x,y; x = 0.55; independent(x);
    passes = run_tape( &x , y , [&]{
      active u=x; for(int i=0;i<14;i++) u = 0.5*u + sin(u)*cos(u); y = u;
    });
    dependent(y);
    Jacobian J = harvest(1,1);
    if(!J.empty()) owned = J(0,0);
    parts_owned = get_partitions();
    finalize();
  }

  double h=0.0,o=0.0;
  MPI_Allreduce(&hand ,&h,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&owned,&o,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  ok("M-05 run_tape gives the identical Jacobian", h==o);
  ok("M-05 run_tape records the identical partitions", parts_hand==parts_owned);
  ok("M-05 run_tape reports 1 profiling pass + 1 per partition",
     passes==parts_owned+1);
}

/* M-06  many tapes in one process, each independent of the last.  The
 *       singleton made this impossible in two separate ways -- MPI_Finalize
 *       in finalize(), and a Process that the second initialize() handed
 *       back unchanged -- and SVEGP-26 is what makes the second one go.
 *       Differentiate three different functions in a row and check none of
 *       them inherited the previous tape's graph or counters. */
static void m06_tape_independence()
{
  const double a = 0.8;

  double d1 = d_of(f_tanh,a,400);
  double d2 = d_of(f_atan,a,400);
  double d3 = d_of(f_tanh,a,400);

  close_to("M-06 tape 1", d1, 1.0-std::tanh(a)*std::tanh(a));
  close_to("M-06 tape 2", d2, 1.0/(1.0+a*a));
  ok      ("M-06 tape 3 is bit-identical to tape 1, not polluted by tape 2", d1==d3);

  /* and the graph is empty between tapes: no tape, no accounting */
  ok("M-06 no tape means no memory reported", get_memory()==0 && get_heap()==0);
}

/* M-07  get_total_partitions() is the tape's partition count known LOCALLY,
 *       without a message: the profiling pass is identical on every rank, so
 *       every rank has the number before the first productive pass.  Assert
 *       it against the expensive answer -- the reduced sum of each rank's
 *       own share -- at whatever rank count this is running on. */
static void m07_total_partitions()
{
  initialize(1,1,400);

  active x,y; x = 0.45; independent(x);

  run_tape( &x , y , [&]{
    active u=x; for(int i=0;i<16;i++) u = 0.5*u + sin(u)*cos(u); y = u;
  });

  dependent(y);
  Jacobian J = harvest(1,1);

  const largeint lib   = get_total_partitions();
  const long     mine  = (long)get_partitions();
  const int      ranks = MPI_size();

  finalize();

  long summed = 0;
  MPI_Allreduce(&mine,&summed,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  int sz=0; MPI_Comm_size(MPI_COMM_WORLD,&sz);

  ok("M-07 get_total_partitions equals the reduced per-rank sum",
     (long)lib == summed);
  ok("M-07 the tape really did chunk, so the check means something", summed > 1);
  ok("M-07 MPI_size agrees with the communicator", ranks == sz);
}

/* ------------------------------------------------------------------------
 * Misuse: every public entry point with no tape open.  None may crash.
 * ---------------------------------------------------------------------- */
static void misuse_no_tape()
{
  active x(1.0), y(2.0);
  double ** A = NULL;

  independent(x);
  dependent(y);
  passive_op(x);
  destructor(x);
  unary_op_ass(1.0,x);
  binary_op_ass(x,1.0,y,1.0);
  postfix_op(x,y);
  unary_op(x,1.0,y,false);
  binary_op(x,1.0,y,1.0,y);
  harvest(1,1,A);
  harvest(1,1,A,false);
  Jacobian J = harvest(1,1);
  free_jacobian(1,A);
  set_input_array_dimension(1,1);
  set_output_array_dimension(1,1);
  set_break_mode(RUN_TO_END);
  set_pass_mode(PASSES_COLLECTIVE);
  set_probe_frequency(10);
  (void)get_total_partitions();

  ok("misuse: harvest with no tape is empty", J.empty());
  ok("misuse: MPI_rank with no tape is -1", MPI_rank()==-1);
  ok("misuse: MPI_size with no tape is -1", MPI_size()==-1);
  ok("misuse: get_memory with no tape is 0", get_memory()==0);
  ok("misuse: get_heap with no tape is 0", get_heap()==0);
  ok("misuse: get_partitions with no tape is 0", get_partitions()==0);
  ok("misuse: get_total_partitions with no tape is 0", get_total_partitions()==0);
  ok("misuse: get_cost with no tape is 0", get_cost()==0);
  ok("misuse: is_harvesting_rank with no tape is false", !is_harvesting_rank());

  /* the one that used to run the caller's section against nothing and then
   * null-deref on the second call */
  int trips = 0;
  while( checkpoint(&x,y) ){ trips++; if(trips>2) break; }
  ok("misuse: checkpoint with no tape ends the loop immediately", trips==0);

  /* and the loop the library owns must reach the same conclusion */
  int body_ran = 0;
  const largeint passes = run_tape( &x , y , [&]{ body_ran++; } );
  ok("misuse: run_tape with no tape runs the section zero times",
     passes==0 && body_ran==0);

  finalize();//with no tape open
  ok("misuse: finalize with no tape is survivable", true);
}


/* M-08  An active that did not survive checkpoint() is caught, not followed.
 *
 *       Found by giving examples/pde a budget that chunks: it built its
 *       Crank-Nicholson solver, and so its active work arrays, once outside
 *       the checkpoint loop.  Safe only while the example never chunked.
 *       AddressSanitizer caught the write path dereferencing a freed vertex
 *       in Vertex::kill(); the read path was worse, because it did not crash
 *       -- the derivative simply came out zero.  Maxwell had the same hole.
 *
 *       Two runs of one program differing only in where one active is
 *       declared.  The correct one must stay correct and stay quiet; the
 *       wrong one must be counted rather than dereferenced.
 */
static double stale_case( bool outside , largeint * reads )
{
  initialize(1,1,4000);

  active x; x = 1.0; independent(x);
  active y, carry;
  int pass = 0;

  run_tape( &x , y , [&]{
    active u = x, local;
    if(outside){ if(pass==0) carry = u*2.0; }
    else       { local = u*2.0; }
    for(int i=0;i<400;i++) u = u*1.0001 + 0.5;
    active z = (outside ? carry : local) + u;
    y = u;
    (void)z;
    pass++;
  });

  dependent(y);
  Jacobian J = harvest(1,1);

  const double   d = J.empty() ? 0.0 : J(0,0);
  const largeint r = get_stale_reads();

  finalize();

  double sum = 0.0;
  MPI_Allreduce(&d,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);//one contributor

  long lr = (long)r, tr = 0;
  MPI_Allreduce(&lr,&tr,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);//it can land on any rank
  *reads = (largeint)tr;

  return sum;
}

static void m08_stale_active()
{
  const double exact = std::pow(1.0001,400.0);

  largeint r_in = 1, r_out = 0;

  const double in = stale_case(false,&r_in);
  ok("M-08 declared inside the section: the derivative is right",
     std::fabs(in-exact) <= 1.0e-9*exact);
  ok("M-08 declared inside the section: nothing is reported",
     r_in==0);

  const double out = stale_case(true,&r_out);
  ok("M-08 declared outside: the stale read is counted, not dereferenced",
     r_out>0);
  ok("M-08 declared outside: the tape still finishes and harvests",
     out==out);//not NaN: the run completed
  ok("M-08 get_stale_reads() separates a clean tape from a dirty one",
     r_in==0 && r_out>0);
}


/* M-09  harvest(double**&) leaves the caller's pointer NULL on every rank
 *       that did not assemble.  free_jacobian()'s comment has always claimed
 *       this -- "tolerates A==NULL, which is what every rank other than the
 *       assembling one has" -- but A was only ever written inside the
 *       is_final_rank() branch, so elsewhere it kept whatever the caller had
 *       declared.  Every legacy example declares it `double ** A;`.
 */
static void m09_harvest_null_off_rank()
{
  initialize(1,1,100000);

  active x; x = 0.7; independent(x);
  active y;

  run_tape( &x , y , [&]{ active u=x; for(int i=0;i<6;i++) u = sin(u)+0.5*u; y = u; });

  dependent(y);

  double ** A = (double**)0xDEADBEEF;//poison: a NULL here has to be written, not inherited
  harvest(1,1,A,false);

  const bool assembling = is_harvesting_rank();

  ok("M-09 harvest gives the assembling rank a matrix",
     !assembling || A!=NULL);
  ok("M-09 harvest NULLs the pointer on every other rank",
     assembling || A==NULL);

  free_jacobian(1,A);
  ok("M-09 free_jacobian is then safe on every rank", A==NULL);

  finalize();
}


/* M-10  The partition counter advances exactly once per boundary crossed.
 *
 *       check_memory() incremented next_owner_idx inside the recording branch
 *       AND in the shared tail.  Under BREAK_ON_TARGET the throw skips the
 *       tail, so that was one increment; under RUN_TO_END nothing skips it and
 *       the rank advanced twice at its own target boundary.  The two modes are
 *       meant to differ only in whether the passive suffix runs.
 *
 *       Invisible at the time -- is_proc() is already permanently false by
 *       then -- but it made the counter a function of WHICH RANK OWNED WHAT
 *       rather than of the operation sequence, and 6.4 needs the second before
 *       a Revolve snapshot can restore it.  Measured on this tape, the value
 *       at the end of the section used to be 37 on some passes and 38 on
 *       others; it is now 37 on all of them, at every rank count.
 *
 *       Only RUN_TO_END can be checked: under BREAK_ON_TARGET the end of the
 *       section is never reached on a pass that records anything.
 */
static void m10_partition_counter()
{
  initialize(1,1,3000);
  set_break_mode(RUN_TO_END);

  active x; x = 0.7; independent(x);
  active y;

  largeint lo = (largeint)-1, hi = 0;
  int samples = 0;

  run_tape( &x , y , [&]{
    active u = x;
    for(int i=0;i<200;i++) u = 0.5*u + sin(u);
    y = u;

    const largeint c = internals::current_tape()->proc.owner_index();
    if(c<lo) lo = c;
    if(c>hi) hi = c;
    samples++;
  });

  dependent(y);
  Jacobian J = harvest(1,1);

  const largeint parts = get_total_partitions();

  finalize();

  long a=(long)lo, b=(long)hi, n=(long)samples, amin=0, bmax=0, nsum=0;
  MPI_Allreduce(&a,&amin,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&b,&bmax,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&n,&nsum,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  ok("M-10 the section ran to the end on every pass (RUN_TO_END)", nsum>0);
  ok("M-10 the partition counter is the same at the end of every pass",
     amin==bmax);
  ok("M-10 ...and equals the tape's partition count, on every rank",
     amin==(long)parts);
}

int main( int argc , char ** argv )
{
  (void)argc; (void)argv;

  misuse_no_tape();//before anything opens a tape

  r01_reachability();
  r02_erf();
  r03_nonsquare_grid();
  r05_ge();
  r06_many_tapes();

  m01_break_mode();
  m02_pass_mode();
  m03_jacobian_value();
  m04_free_jacobian();
  m05_run_tape();
  m06_tape_independence();
  m07_total_partitions();
  m08_stale_active();
  m09_harvest_null_off_rank();
  m10_partition_counter();

  misuse_no_tape();//and again after every tape has been closed

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int total = 0;
  MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(!rank) std::printf("regress: %d checks, %d failures -- %s\n",
                        checks,total,total?"FAIL":"PASS");

  return total ? 1 : 0;
}
