/*
 * 1-D heat equation: the gradient of a misfit functional.
 *
 * WHAT THIS IS FOR.  Two things, and the second is why it was rewritten.
 *
 *   1  It is the SHAPE of a data-assimilation problem: many inputs, one
 *      output.  The independents are the whole initial temperature profile,
 *      the dependent is a single scalar misfit against an observation, and
 *      what comes out is a gradient -- a 1 x (nx+1) Jacobian.  That is the
 *      opposite aspect ratio to lorenz (3 x 3) and to fdcheck, and vertex
 *      elimination behaves differently on a graph that funnels down to one
 *      sink than on one that fans out.
 *
 *   2  It CHUNKS, and it says so.  A time-stepping loop makes a tape that is
 *      deep and narrow, so a modest byte budget breaks it into dozens of
 *      partitions dealt round the ranks.  The old version of this example ran
 *      with a budget of it*215208 bytes -- hundreds of megabytes -- which
 *      bought exactly one partition, printed a column of numbers nobody
 *      checked, and returned 0 whatever happened.  It could not have failed.
 *
 * THE PROBLEM.  u_t = c u_xx on [0,1], Dirichlet at both ends, marched with
 * explicit Euler on nx intervals for nt steps of dt = 1/nt:
 *
 *     u[j] <- u[j] + k ( u[j+1] - 2 u[j] + u[j-1] ) ,   k = nx^2 c dt
 *
 * interior only; u[0] and u[nx] are held.  The cost is the squared misfit
 * against a fixed observed profile,
 *
 *     J(u0) = sum_{j=1}^{nx-1} ( u_nt[j] - obs[j] )^2 ,
 *
 * and the example computes dJ/du0[i] for every i, boundary nodes included --
 * they are not inert, because u[0] feeds u[1] at every step.
 *
 * STABILITY.  k <= 1/2 or explicit Euler diverges and the gradient is
 * enormous, correct, and meaningless.  main() checks k and refuses to run
 * rather than reporting a "pass" on a blown-up solution.
 *
 * TWO INDEPENDENT CHECKS, both exact in exact arithmetic.
 *
 *   1  THE DISCRETE ADJOINT.  The time stepping is LINEAR: u_nt = A^nt u_0
 *      with A the stencil matrix.  So the gradient is
 *
 *          dJ/du0 = (A^T)^nt w ,   w[j] = 2 ( u_nt[j] - obs[j] )
 *
 *      which this example builds by transposing the stencil BY HAND and
 *      sweeping it backwards in plain double arithmetic.  No AD anywhere in
 *      it, and it is a reverse sweep, whereas the library gets the same
 *      numbers by eliminating vertices of a graph in an order it chooses for
 *      itself.  Two genuinely different algorithms; agreement to machine
 *      precision is real evidence.
 *
 *   2  CENTRAL DIFFERENCES.  Usually a weak check, because the O(h^2)
 *      truncation error swamps everything.  Not here: J is a QUADRATIC
 *      function of u0, its third derivative is identically zero, and the
 *      central difference is therefore EXACT up to rounding.  That makes it
 *      a second independent reference good to ~1e-10 rather than the ~1e-6
 *      a finite difference is normally worth, and it needs no knowledge of
 *      the adjoint at all -- it would catch an error I had made in BOTH the
 *      library and the hand-written transpose.
 *
 * See the CHUNKING section at the bottom of main() for what is asserted about
 * the partitioning itself.
 */
/* The library does not leak <mpi.h> through boltzmann.hpp, on purpose.  This
 * example reduces its own counts across the ranks, so it includes it itself. */
#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "../../boltzmann.hpp"

using namespace boltzmann;

/* One explicit-Euler step, in whatever arithmetic.  The SAME source is what
 * the tape records and what the reference runs, so the two cannot drift
 * apart through an edit to one of them. */
template<class T>
static void heat_step( T * u , T * old , int nx , double k )
{
  for( int j=0 ; j<=nx ; j++ ) old[j] = u[j];

  for( int j=1 ; j<nx ; j++ )
    u[j] = old[j] + k*( old[j+1] - 2.0*old[j] + old[j-1] );
}

/* March, then form the misfit.  Again shared between active and passive. */
template<class T>
static void heat_run( T * u , int nx , int nt , double k , const double * obs , T & cost )
{
  std::vector<T> old(nx+1);

  for( int n=0 ; n<nt ; n++ ) heat_step<T>( u , &old[0] , nx , k );

  cost = 0.0;
  for( int j=1 ; j<nx ; j++ ){
    T d = u[j] - obs[j];
    cost += d*d;
  }
}

/* w <- A^T w for one step.  The transpose of
 *     new[j] = old[j] + k(old[j+1] - 2 old[j] + old[j-1])   (interior)
 *     new[j] = old[j]                                       (j = 0, nx)
 * scattered rather than gathered, which is what makes it an adjoint and not
 * just the same stencil written again. */
static void adjoint_step( double * w , double * t , int nx , double k )
{
  for( int j=0 ; j<=nx ; j++ ) t[j] = 0.0;

  t[0]  += w[0];
  t[nx] += w[nx];

  for( int j=1 ; j<nx ; j++ ){
    t[j]   += (1.0-2.0*k)*w[j];
    t[j+1] += k*w[j];
    t[j-1] += k*w[j];
  }

  for( int j=0 ; j<=nx ; j++ ) w[j] = t[j];
}

static double cost_at( const double * u0 , int nx , int nt , double k , const double * obs )
{
  std::vector<double> u( u0 , u0+nx+1 );
  double c = 0.0;
  heat_run<double>( &u[0] , nx , nt , k , obs , c );
  return c;
}

static void adjoint_gradient( const double * u0 , int nx , int nt , double k ,
                              const double * obs , double * g , double & cost_out )
{
  std::vector<double> u( u0 , u0+nx+1 ) , t(nx+1) , w(nx+1,0.0);

  for( int n=0 ; n<nt ; n++ ) heat_step<double>( &u[0] , &t[0] , nx , k );

  double c = 0.0;
  for( int j=1 ; j<nx ; j++ ){
    const double d = u[j]-obs[j];
    c   += d*d;
    w[j] = 2.0*d;
  }
  cost_out = c;

  for( int n=0 ; n<nt ; n++ ) adjoint_step( &w[0] , &t[0] , nx , k );

  for( int j=0 ; j<=nx ; j++ ) g[j] = w[j];
}

int main( int argc , char ** argv )
{
  /*
   * The defaults are the Makefile's RUNARGS, so a bare ./heat exercises the
   * pipeline rather than running one fat partition.  nt is what buys the
   * chunking: the tape is nt stencil sweeps long, and the budget decides how
   * many of them fit between two breaks.
   */
  const int      nx     = (argc>1) ? std::atoi(argv[1]) : 40;
  const int      nt     = (argc>2) ? std::atoi(argv[2]) : 400;
  const largeint budget = (argc>3) ? (largeint)std::atol(argv[3]) : 150000;

  if( nx<3 || nt<1 ){ std::printf("heat: need nx>=3 and nt>=1\n"); return 1; }

  const double c  = 0.001;
  const double dt = 1.0/double(nt);
  const double k  = double(nx)*double(nx)*c*dt;

  if( !(k<=0.5) ){
    std::printf("heat: k = nx^2 c dt = %g > 1/2 -- explicit Euler is unstable here.\n",k);
    std::printf("      raise nt or lower nx; the gradient would be correct and useless.\n");
    return 1;
  }

  std::vector<double> u0(nx+1) , obs(nx+1);
  u0[0] = 2.0;
  for( int i=1 ; i<=nx ; i++ ) u0[i] = 0.0;
  for( int i=0 ; i<=nx ; i++ ) obs[i] = 2.0 - double(i)/100.0;

  /* ---- the AD gradient ------------------------------------------------ */

  std::vector<double> gad(nx+1,0.0);

  initialize( (largeint)(nx+1) , 1 , budget );

  std::vector<active> x(nx+1);
  active y;

  for( int i=0 ; i<=nx ; i++ ){
    x[i] = u0[i];
    independent(x[i]);
  }

  /*
   * run_tape() owns the try/catch (inc/Run.hpp).  The body is abandoned from
   * inside an operator once per partition, and everything it allocates --
   * the two std::vector<active> here and the one inside heat_run -- is
   * released by unwinding, which is the whole point of the RAII surface.
   */
  const largeint passes = run_tape( &x[0] , y , [&]{
    std::vector<active> u(nx+1);
    for( int i=0 ; i<=nx ; i++ ) u[i] = x[i];
    heat_run<active>( &u[0] , nx , nt , k , &obs[0] , y );
  });

  dependent(y);

  Jacobian J = harvest( 1 , (largeint)(nx+1) );

  if(!J.empty()){
    for( int i=0 ; i<=nx ; i++ ) gad[i] = J(0,i);
  }

  const largeint my_parts  = get_partitions();
  const largeint my_cost   = get_cost();
  const largeint lib_total = get_total_partitions();//no message; see API.hpp
  const int      lib_ranks = MPI_size();

  finalize();

  /* ---- gather: the gradient lives on one rank, the counts on all ------ */

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  std::vector<double> g(nx+1,0.0);
  MPI_Allreduce(&gad[0],&g[0],nx+1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);//one contributor

  long lp=(long)my_parts, total_parts=0, min_parts=0, max_parts=0, ranks_working=0;
  long one = (my_parts>0) ? 1 : 0;
  MPI_Allreduce(&lp ,&total_parts  ,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&min_parts    ,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&max_parts    ,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&one,&ranks_working,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  long lc=(long)my_cost, total_cost=0, min_cost=0, max_cost=0;
  MPI_Allreduce(&lc,&total_cost,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&lc,&min_cost  ,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lc,&max_cost  ,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

  /* ---- reference 1: the discrete adjoint ------------------------------ */

  std::vector<double> gref(nx+1,0.0);
  double cost_ref = 0.0;
  adjoint_gradient( &u0[0] , nx , nt , k , &obs[0] , &gref[0] , cost_ref );

  double scale = 0.0;
  for( int i=0 ; i<=nx ; i++ )
    if( std::fabs(gref[i])>scale ) scale = std::fabs(gref[i]);

  double worst_adj = 0.0;
  int    where_adj = 0;
  for( int i=0 ; i<=nx ; i++ ){
    const double a = std::fabs(g[i]-gref[i]);
    if( a>worst_adj ){ worst_adj=a; where_adj=i; }
  }
  const double rel_adj = (scale>0.0) ? worst_adj/scale : worst_adj;

  /* ---- reference 2: central differences, exact for a quadratic -------- */

  double worst_fd = 0.0;
  int    where_fd = 0;
  const double h = 1.0e-4;

  for( int i=0 ; i<=nx ; i++ ){
    std::vector<double> up(u0), um(u0);
    up[i] += h;
    um[i] -= h;
    const double d = ( cost_at(&up[0],nx,nt,k,&obs[0])
                     - cost_at(&um[0],nx,nt,k,&obs[0]) ) / (2.0*h);
    const double a = std::fabs(g[i]-d);
    if( a>worst_fd ){ worst_fd=a; where_fd=i; }
  }
  const double rel_fd = (scale>0.0) ? worst_fd/scale : worst_fd;

  /* ---- verdicts ------------------------------------------------------- */

  /* Same arithmetic, different algorithm: machine precision, with headroom
   * over the ~1e-15 observed at every rank count. */
  const double ADJ_TOL = 1.0e-11;

  /* Exact in exact arithmetic, so the only error is cancellation in the
   * difference of two O(1) costs divided by 2h: eps*J/h ~ 1e-12 relative. */
  const double FD_TOL  = 1.0e-8;

  int bad = 0;

  if( !(rel_adj <= ADJ_TOL) ){
    bad++;
    if(!rank) std::printf("  FAIL adjoint: %.3e relative at node %d, tol %.1e\n",
                          rel_adj,where_adj,ADJ_TOL);
  }
  if( !(rel_fd <= FD_TOL) ){
    bad++;
    if(!rank) std::printf("  FAIL central differences: %.3e relative at node %d, tol %.1e\n",
                          rel_fd,where_fd,FD_TOL);
  }

  /* The library's own total against the reduced truth -- see lorenz. */
  if( lib_total != (largeint)total_parts ){
    bad++;
    if(!rank) std::printf("  FAIL library reports %lu partitions, the ranks recorded %ld\n",
                          (unsigned long)lib_total,total_parts);
  }
  if( lib_ranks != size ){
    bad++;
    if(!rank) std::printf("  FAIL library reports %d ranks, MPI reports %d\n",lib_ranks,size);
  }

  /* ---- CHUNKING: what this example exists to exercise ------------------ */

  if( total_parts < 8 ){
    bad++;
    if(!rank) std::printf("  FAIL only %ld partitions -- the tape is not being chunked.  "
                          "Lower the budget or raise nt.\n",total_parts);
  }
  if( size>1 && ranks_working<2 ){
    bad++;
    if(!rank) std::printf("  FAIL %ld partitions all landed on one rank of %d -- "
                          "the pipeline is not being exercised.\n",total_parts,size);
  }
  if( passes != my_parts+1 ){
    bad++;
    if(!rank) std::printf("  FAIL rank %d ran %lu passes for %lu partitions\n",
                          rank,(unsigned long)passes,(unsigned long)my_parts);
  }

  /* ---- report --------------------------------------------------------- */

  if(!rank){
    std::printf("heat: %lu partitions over %d MPI rank%s\n",
                (unsigned long)lib_total,lib_ranks,(lib_ranks==1)?"":"s");
    std::printf("  problem      nx=%d, nt=%d explicit-Euler steps, k=%g, budget=%lu bytes\n",
                nx,nt,k,(unsigned long)budget);
    std::printf("  gradient     1 x %d, misfit J = %.10e\n",nx+1,cost_ref);
    std::printf("  partitions   %ld recorded in total, %ld..%ld per rank, "
                "%ld of %d ranks working\n",
                total_parts,min_parts,max_parts,ranks_working,size);
    std::printf("  elim cost    %ld total, %ld..%ld per rank\n",
                total_cost,min_cost,max_cost);
    std::printf("  vs discrete adjoint        %.3e relative  (tol %.1e)\n",rel_adj,ADJ_TOL);
    std::printf("  vs central differences     %.3e relative  (tol %.1e)\n",rel_fd,FD_TOL);
    std::printf("heat: %s\n",(bad==0)?"PASS":"FAIL");
  }

  return (bad==0) ? 0 : 1;
}
