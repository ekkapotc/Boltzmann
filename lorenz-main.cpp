/*
 * Lorenz-63: the tangent-linear propagator of a chaotic flow.
 *
 * WHAT THIS IS FOR.  Every other example in this tree either does not chunk
 * at all or chunks a handful of times.  This one is built to drive the
 * machinery the library exists for: a long time integration makes a tape that
 * is DEEP and NARROW, so a modest byte budget breaks it into dozens or
 * hundreds of partitions, and those partitions are dealt round the ranks.  It
 * asserts that this happened rather than hoping -- see the CHUNKING section
 * at the bottom of main().
 *
 * THE PROBLEM.  Lorenz's 1963 convection model,
 *
 *     x' = s(y - x)
 *     y' = x(r - z) - y
 *     z' = xy - bz              s = 10, r = 28, b = 8/3
 *
 * marched with classical RK4.  The Jacobian of the final state with respect
 * to the initial state,
 *
 *     M(T) = d u(T) / d u(0)                                   (3x3)
 *
 * is the TANGENT-LINEAR PROPAGATOR.  It is the object at the centre of
 * variational data assimilation (4D-Var solves with M and its adjoint), of
 * Lyapunov exponent computation, and of any sensitivity study of a chaotic
 * system.  It is also a genuinely hard thing to compute by finite differences
 * -- the system is chaotic, so a perturbation grows like e^(0.9 t) and the
 * step size that is small enough to be linear is too small to be accurate.
 * That is the practical argument for AD here, and it is why this example
 * checks itself against something better than finite differences.
 *
 * TWO INDEPENDENT CHECKS, neither of them a finite difference.
 *
 *   1  THE VARIATIONAL EQUATION.  Integrate dM/dt = J(u(t)) M alongside the
 *      state, in plain double arithmetic, with no AD anywhere.  For a
 *      Runge-Kutta method the tangent-linear of the scheme IS the scheme
 *      applied to the variational equation using the same stage values, so
 *      this agrees with the AD Jacobian to MACHINE PRECISION, not to
 *      discretisation order.  Any disagreement above ~1e-13 is a bug in the
 *      library, not a modelling difference.
 *
 *   2  LIOUVILLE'S THEOREM.  The divergence of the Lorenz vector field is
 *
 *          div f = -s - 1 - b = -13.666...
 *
 *      a CONSTANT, so the exact flow contracts phase-space volume at a fixed
 *      rate and
 *
 *          det M(T) = exp( -(s + 1 + b) T )    exactly.
 *
 *      This is an analytic identity about the whole 3x3 Jacobian at once: it
 *      constrains all nine entries jointly, and no amount of agreeing with
 *      another numerical propagator would satisfy it by accident.
 *
 *      IT IS ALSO SEVERELY ILL-CONDITIONED, AND THE EXAMPLE SAYS SO RATHER
 *      THAN PRETENDING OTHERWISE.  The determinant is the product of the
 *      three singular values of M.  One of them grows like e^(0.9 T) while
 *      the product shrinks like e^(-13.67 T), so the smallest is doing all
 *      the work and the determinant is a tiny number formed by cancellation
 *      between entries that are not tiny.  Computing it in double precision
 *      costs roughly
 *
 *          relative error  ~  eps * ||M||^3 / |det M|  =  eps * e^(3L + 13.67)T
 *
 *      which is 5e-14 at T=0.25, 1e-11 at T=0.5, 4e-7 at T=1 and larger than
 *      ONE at T=2.  Measured on this example at T=2: 8e-2 at one rank and 2.0
 *      at four, while the tangent-linear check -- which is well conditioned --
 *      held at 3e-15 at every rank count.  The determinant disagreeing across
 *      rank counts there is not the library getting a different answer; it is
 *      the one-ulp reassociation of Section 5 amplified by 10^14.
 *
 *      So the example estimates that bound from the Jacobian it actually
 *      computed, adds the RK4 discretisation error (the discrete propagator
 *      satisfies Liouville only to the order of the scheme, O(dt^4) per unit
 *      time), and asserts the identity only where the arithmetic can resolve
 *      it.  Past that point it reports the check as unresolvable and says
 *      why.  A test that quietly widens its tolerance until it passes is
 *      worse than no test.
 *
 * The leading Lyapunov exponent is reported as well, from the growth of the
 * largest singular value of M.  It is an indicator, not an assertion: a
 * converged value (0.906 for these parameters) needs repeated
 * renormalisation over a much longer trajectory than one tape.
 *
 *   ./lorenz [steps] [budget_bytes] [dt]
 */

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "../../boltzmann.hpp"

using namespace boltzmann;

static const double SIGMA = 10.0;
static const double RHO   = 28.0;
static const double BETA  = 8.0/3.0;

/* ----------------------------------------------------------- the model --
 *
 * One template over the scalar type, so the ACTIVE and PASSIVE integrations
 * are demonstrably the same arithmetic rather than two transcriptions of it.
 * boltzmann::active supports everything used here.
 */
template<class T>
static void rhs( const T & x , const T & y , const T & z , T & fx , T & fy , T & fz )
{
  fx = SIGMA*(y - x);
  fy = x*(RHO - z) - y;
  fz = x*y - BETA*z;
}

template<class T>
static void rk4_step( T & x , T & y , T & z , double dt )
{
  T k1x,k1y,k1z, k2x,k2y,k2z, k3x,k3y,k3z, k4x,k4y,k4z;

  rhs<T>( x , y , z , k1x,k1y,k1z );
  rhs<T>( x+0.5*dt*k1x , y+0.5*dt*k1y , z+0.5*dt*k1z , k2x,k2y,k2z );
  rhs<T>( x+0.5*dt*k2x , y+0.5*dt*k2y , z+0.5*dt*k2z , k3x,k3y,k3z );
  rhs<T>( x+dt*k3x , y+dt*k3y , z+dt*k3z , k4x,k4y,k4z );

  x = x + (dt/6.0)*(k1x + 2.0*k2x + 2.0*k3x + k4x);
  y = y + (dt/6.0)*(k1y + 2.0*k2y + 2.0*k3y + k4y);
  z = z + (dt/6.0)*(k1z + 2.0*k2z + 2.0*k3z + k4z);
}

/* ------------------------------------------- the variational reference --
 *
 * The state and the 3x3 propagator marched together with the SAME RK4, in
 * plain doubles.  No active, no tape, no AD: this is the independent answer.
 */
static void jac( const double * u , double J[3][3] )
{
  J[0][0] = -SIGMA;        J[0][1] =  SIGMA;  J[0][2] =  0.0;
  J[1][0] =  RHO - u[2];   J[1][1] = -1.0;    J[1][2] = -u[0];
  J[2][0] =  u[1];         J[2][1] =  u[0];   J[2][2] = -BETA;
}

//dz/dt for the augmented system z = (u, M)
static void aug_rhs( const double * u , const double * M , double * du , double * dM )
{
  rhs<double>( u[0],u[1],u[2] , du[0],du[1],du[2] );

  double J[3][3];
  jac(u,J);

  for( int i=0 ; i<3 ; i++ ){
    for( int j=0 ; j<3 ; j++ ){
      double s = 0.0;
      for( int k=0 ; k<3 ; k++ ) s += J[i][k]*M[3*k+j];
      dM[3*i+j] = s;
    }
  }
}

static void variational( double * u , double * M , int steps , double dt )
{
  for( int i=0 ; i<3 ; i++ )
    for( int j=0 ; j<3 ; j++ )
      M[3*i+j] = (i==j) ? 1.0 : 0.0;

  double k1u[3],k1M[9], k2u[3],k2M[9], k3u[3],k3M[9], k4u[3],k4M[9];
  double tu[3], tM[9];

  for( int n=0 ; n<steps ; n++ ){

    aug_rhs(u,M,k1u,k1M);

    for(int i=0;i<3;i++) tu[i]=u[i]+0.5*dt*k1u[i];
    for(int i=0;i<9;i++) tM[i]=M[i]+0.5*dt*k1M[i];
    aug_rhs(tu,tM,k2u,k2M);

    for(int i=0;i<3;i++) tu[i]=u[i]+0.5*dt*k2u[i];
    for(int i=0;i<9;i++) tM[i]=M[i]+0.5*dt*k2M[i];
    aug_rhs(tu,tM,k3u,k3M);

    for(int i=0;i<3;i++) tu[i]=u[i]+dt*k3u[i];
    for(int i=0;i<9;i++) tM[i]=M[i]+dt*k3M[i];
    aug_rhs(tu,tM,k4u,k4M);

    for(int i=0;i<3;i++) u[i] += (dt/6.0)*(k1u[i]+2.0*k2u[i]+2.0*k3u[i]+k4u[i]);
    for(int i=0;i<9;i++) M[i] += (dt/6.0)*(k1M[i]+2.0*k2M[i]+2.0*k3M[i]+k4M[i]);
  }
}

/* --------------------------------------------------------------- misc -- */

static double det3( const double * M )
{
  return M[0]*(M[4]*M[8]-M[5]*M[7])
       - M[1]*(M[3]*M[8]-M[5]*M[6])
       + M[2]*(M[3]*M[7]-M[4]*M[6]);
}

//largest singular value, by power iteration on M^T M.  3x3: a few sweeps do.
static double sigma_max( const double * M )
{
  double A[9];

  for( int i=0 ; i<3 ; i++ )
    for( int j=0 ; j<3 ; j++ ){
      double s=0.0;
      for( int k=0 ; k<3 ; k++ ) s += M[3*k+i]*M[3*k+j];
      A[3*i+j]=s;
    }

  double v[3] = {1.0,1.0,1.0}, w[3];

  for( int it=0 ; it<200 ; it++ ){
    for( int i=0 ; i<3 ; i++ ){
      double s=0.0;
      for( int k=0 ; k<3 ; k++ ) s += A[3*i+k]*v[k];
      w[i]=s;
    }
    double n = std::sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
    if(n<=0.0) return 0.0;
    for( int i=0 ; i<3 ; i++ ) v[i]=w[i]/n;
  }

  double s=0.0;
  for( int i=0 ; i<3 ; i++ ){
    double t=0.0;
    for( int k=0 ; k<3 ; k++ ) t += A[3*i+k]*v[k];
    s += v[i]*t;
  }

  return std::sqrt(s>0.0?s:0.0);
}

int main( int argc , char ** argv )
{
  /*
   * The defaults are the Makefile's RUNARGS, so a bare ./lorenz shows the
   * example at its best rather than in a corner.  T = dt*steps = 0.5 is
   * chosen so that BOTH checks are live: much past T = 1 the Liouville
   * identity stops being resolvable in double precision (see above) and the
   * example would run with only one of its two checks doing any work.
   */
  const int      steps = (argc>1) ? std::atoi(argv[1]) : 2000;
  const largeint budget= (argc>2) ? (largeint)std::atol(argv[2]) : 180000;
  const double   dt    = (argc>3) ? std::atof(argv[3]) : 0.00025;

  const double T = dt*double(steps);

  const double u0[3] = { 1.0 , 1.0 , 20.0 };//on the attractor, near enough

  /* ---- the AD Jacobian ------------------------------------------------ */

  double Mad[9];
  for( int i=0 ; i<9 ; i++ ) Mad[i]=0.0;

  initialize(3,3,budget);

  active x[3], y[3];

  for( int i=0 ; i<3 ; i++ ){
    x[i] = u0[i];
    independent(x[i]);
  }

  /*
   * run_tape() owns the try/catch (inc/Run.hpp), so this section -- which is
   * abandoned from inside an operator dozens of times over -- needs no
   * exception handling of its own.  It allocates nothing, so there is
   * nothing for the throw to strand either.
   */
  const largeint passes = run_tape( x , y , [&]{
    active a=x[0], b=x[1], c=x[2];
    for( int n=0 ; n<steps ; n++ ) rk4_step<active>( a , b , c , dt );
    y[0]=a; y[1]=b; y[2]=c;
  });

  for( int i=0 ; i<3 ; i++ ) dependent(y[i]);

  Jacobian J = harvest(3,3);

  if(!J.empty()){
    for( int i=0 ; i<3 ; i++ )
      for( int j=0 ; j<3 ; j++ )
        Mad[3*i+j] = J(i,j);
  }

  const largeint my_parts = get_partitions();
  const largeint my_cost  = get_cost();

  finalize();

  /* ---- gather: the Jacobian lives on one rank, the counts on all ------ */

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  double M[9];
  MPI_Allreduce(Mad,M,9,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);//exactly one contributor

  long lp=(long)my_parts, total_parts=0, min_parts=0, max_parts=0, ranks_working=0;
  long one = (my_parts>0) ? 1 : 0;
  MPI_Allreduce(&lp ,&total_parts,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&min_parts  ,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&max_parts  ,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&one,&ranks_working,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  long lc=(long)my_cost, total_cost=0, min_cost=0, max_cost=0;
  MPI_Allreduce(&lc,&total_cost,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&lc,&min_cost  ,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lc,&max_cost  ,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);

  /* ---- the reference, and the two checks ------------------------------ */

  double uref[3] = { u0[0],u0[1],u0[2] }, Mref[9];
  variational(uref,Mref,steps,dt);

  double worst = 0.0, scale = 0.0;
  for( int i=0 ; i<9 ; i++ ){
    const double a=std::fabs(M[i]-Mref[i]);
    if(a>worst) worst=a;
    if(std::fabs(Mref[i])>scale) scale=std::fabs(Mref[i]);
  }
  const double rel = (scale>0.0) ? worst/scale : worst;

  const double det_ad   = det3(M);
  const double det_liou = std::exp(-(SIGMA+1.0+BETA)*T);
  const double det_rel  = std::fabs(det_ad-det_liou)/det_liou;

  const double smax = sigma_max(M);
  const double lyap = (T>0.0 && smax>0.0) ? std::log(smax)/T : 0.0;

  /* ---- verdicts -------------------------------------------------------
   *
   * TANGENT-LINEAR:  machine precision, because the tangent-linear of a
   * Runge-Kutta method is that method applied to the variational equation
   * with the same stage values.  1e-11 leaves four orders of headroom over
   * what is observed, which is ~3e-15 at every rank count.
   */
  const double TL_TOL = 1.0e-11;

  /*
   * LIOUVILLE:  two error sources, and the example works out both rather
   * than carrying a hand-tuned constant.
   *
   *   cancel   the conditioning of det for THIS Jacobian: eps*||M||^3/|det|.
   *            Grows like e^(3L+13.67)T and passes 1 somewhere near T=2.
   *   disc     RK4's discrete propagator satisfies the continuous identity
   *            only to the order of the scheme: O(dt^4) per unit time.
   *
   * If the cancellation floor alone is above 1e-4 the identity simply cannot
   * be resolved in double precision at this T, and the example says so
   * instead of asserting something arithmetic cannot decide.
   */
  double nrm2 = 0.0;
  for( int i=0 ; i<9 ; i++ ) nrm2 += M[i]*M[i];
  const double nrmF = std::sqrt(nrm2);

  const double EPS    = 2.220446049250313e-16;
  const double cancel = EPS*nrmF*nrmF*nrmF/det_liou;
  const double disc   = 200.0*T*dt*dt*dt*dt;

  const double LIOU_TOL       = 100.0*cancel + disc + 1.0e-14;
  const bool   liou_resolvable = (cancel < 1.0e-4);

  int bad = 0;

  if( !(rel <= TL_TOL) ){
    bad++;
    if(!rank) std::printf("  FAIL tangent-linear: %.3e relative, tol %.1e\n",rel,TL_TOL);
  }
  if( liou_resolvable && !(det_rel <= LIOU_TOL) ){
    bad++;
    if(!rank) std::printf("  FAIL Liouville: %.3e relative, tol %.1e\n",det_rel,LIOU_TOL);
  }

  /* ---- CHUNKING: what this example exists to exercise ------------------
   *
   * Maxwell's SVEGP-22 point, sharpened for a pipeline.  An example that
   * silently stops chunking stops testing anything, and one that chunks but
   * leaves every partition on one rank is not testing the pipeline at all.
   * Both are assertions here, not observations.
   */
  if( total_parts < 8 ){
    bad++;
    if(!rank) std::printf("  FAIL only %ld partitions -- the tape is not being chunked.  "
                          "Lower the budget or raise the step count.\n",total_parts);
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
    std::printf("lorenz: %d RK4 steps, dt=%g, T=%g, budget=%lu bytes, np=%d\n",
                steps,dt,T,(unsigned long)budget,size);
    std::printf("  partitions   %ld total, %ld..%ld per rank, %ld of %d ranks recorded\n",
                total_parts,min_parts,max_parts,ranks_working,size);
    std::printf("  elim cost    %ld total, %ld..%ld per rank\n",
                total_cost,min_cost,max_cost);
    std::printf("  tangent-linear vs variational   %.3e relative  (tol %.1e)\n",rel,TL_TOL);
    if(liou_resolvable){
      std::printf("  det M vs Liouville exp(-%.4f T) %.3e relative  (tol %.1e)\n",
                  SIGMA+1.0+BETA,det_rel,LIOU_TOL);
    }else{
      std::printf("  det M vs Liouville exp(-%.4f T) NOT RESOLVABLE at T=%g:\n",
                  SIGMA+1.0+BETA,T);
      std::printf("               |det| = %.2e against ||M|| = %.2e, so cancellation\n",
                  det_liou,nrmF);
      std::printf("               alone costs %.1e relative.  Observed %.1e.  Shorten T.\n",
                  cancel,det_rel);
    }
    std::printf("  finite-time Lyapunov, over T=%g   %.4f  (transient; the asymptotic\n",T,lyap);
    std::printf("               value is 0.906 and needs renormalisation over t ~ 100)\n");
    std::printf("lorenz: %s\n", bad?"FAIL":"PASS");
  }

  return bad ? 1 : 0;
}
