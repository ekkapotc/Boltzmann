/*
 * Every elementary function against central differences.
 *
 * Ported from Maxwell, where this net is what caught the two defects that
 * silently zeroed six derivatives and got erf's constant wrong.  This tree
 * had no equivalent: not one example asserted a numeric result, they all
 * printed a Jacobian and returned 0, so a wrong partial was something a human
 * had to notice in a column of numbers.
 *
 * Runs on every rank and asserts on every rank, because a wrong partial that
 * only appears when the operation lands in rank 2's partition is exactly the
 * kind of bug this library can have and Maxwell cannot.  It also runs each
 * check at a budget small enough to force several partitions, so the
 * checkpoint/replay path is under test rather than bypassed.
 *
 * Exit status is non-zero on any failure, so `make test` stops.
 */

#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

#include "../../boltzmann.hpp"

using namespace boltzmann;

static int failures = 0;
static int checks   = 0;
static int my_rank  = 0;

/* ------------------------------------------------------------------------ */

typedef active (*unary_fn)( const active & );
typedef double (*unary_ref)( double );

typedef active (*binary_fn)( const active & , const active & );
typedef double (*binary_ref)( double , double );

static const largeint TIGHT_BUDGET = 400;//forces the tape to chunk

static void report( const char * name , double got , double want , double tol )
{
  checks++;

  const double err = std::fabs(got-want);
  const double rel = err/(1.0+std::fabs(want));

  if( rel > tol || got!=got ){
    failures++;
    std::printf("  rank %d  FAIL  %-10s ad=% .12e  fd=% .12e  rel=%.3e\n",
                my_rank,name,got,want,rel);
  }
}

/* d/dx f(x) at x, taped, harvested on the assembling rank and broadcast so
 * that every rank can assert on it. */
static double ad_unary( unary_fn f , double x0 )
{
  initialize(1,1,TIGHT_BUDGET);

  active x, y;
  x = x0;
  independent(x);

  while( checkpoint(&x,y) ){
    try{ y = f(x); }catch( BreakException const & ){}
  }

  dependent(y);

  Jacobian J = harvest(1,1);

  double d = J.empty() ? 0.0 : J(0,0);
  int    who = J.empty() ? 0 : 1;

  finalize();

  /* exactly one rank assembled it; tell the others what it was */
  double sum = 0.0; int owners = 0;
  MPI_Allreduce(&d,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&who,&owners,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(owners!=1){
    failures++;
    std::printf("  rank %d  FAIL  %d ranks claimed the Jacobian, expected exactly 1\n",
                my_rank,owners);
  }

  return sum;
}

static void check_unary( const char * name , unary_fn f , unary_ref ref , double x0 , double tol=1e-6 )
{
  const double h  = 1e-6*(1.0+std::fabs(x0));
  const double fd = (ref(x0+h)-ref(x0-h))/(2.0*h);

  report(name,ad_unary(f,x0),fd,tol);
}

/* ------------------------------------------------------------------------ */

static double r_sin(double x){return std::sin(x);}     static active a_sin(const active&x){return sin(x);}
static double r_cos(double x){return std::cos(x);}     static active a_cos(const active&x){return cos(x);}
static double r_tan(double x){return std::tan(x);}     static active a_tan(const active&x){return tan(x);}
static double r_asin(double x){return std::asin(x);}   static active a_asin(const active&x){return asin(x);}
static double r_acos(double x){return std::acos(x);}   static active a_acos(const active&x){return acos(x);}
static double r_atan(double x){return std::atan(x);}   static active a_atan(const active&x){return atan(x);}
static double r_sinh(double x){return std::sinh(x);}   static active a_sinh(const active&x){return sinh(x);}
static double r_cosh(double x){return std::cosh(x);}   static active a_cosh(const active&x){return cosh(x);}
static double r_tanh(double x){return std::tanh(x);}   static active a_tanh(const active&x){return tanh(x);}
static double r_asinh(double x){return ::asinh(x);}    static active a_asinh(const active&x){return asinh(x);}
static double r_acosh(double x){return ::acosh(x);}    static active a_acosh(const active&x){return acosh(x);}
static double r_atanh(double x){return ::atanh(x);}    static active a_atanh(const active&x){return atanh(x);}
static double r_exp(double x){return std::exp(x);}     static active a_exp(const active&x){return exp(x);}
static double r_log(double x){return std::log(x);}     static active a_log(const active&x){return log(x);}
static double r_log10(double x){return std::log10(x);} static active a_log10(const active&x){return log10(x);}
static double r_sqrt(double x){return std::sqrt(x);}   static active a_sqrt(const active&x){return sqrt(x);}
static double r_cbrt(double x){return ::cbrt(x);}      static active a_cbrt(const active&x){return cbrt(x);}
static double r_fabs(double x){return std::fabs(x);}   static active a_fabs(const active&x){return fabs(x);}
static double r_erf(double x){return ::erf(x);}        static active a_erf(const active&x){return erf(x);}
static double r_erfc(double x){return ::erfc(x);}      static active a_erfc(const active&x){return erfc(x);}
static double r_neg(double x){return -x;}              static active a_neg(const active&x){return -x;}
static double r_pos(double x){return +x;}              static active a_pos(const active&x){return +x;}

/* composites, so the binary operators and the *_op_ass path are covered too */
static double r_mul(double x){return x*x*x;}           static active a_mul(const active&x){return x*x*x;}
static double r_div(double x){return (x+1.0)/(x*x+2.0);}
static active a_div(const active&x){return (x+1.0)/(x*x+2.0);}
static double r_pow(double x){return std::pow(x,2.7);} static active a_pow(const active&x){return pow(x,2.7);}
static double r_powa(double x){return std::pow(x,x);}  static active a_powa(const active&x){return pow(x,x);}
static double r_pow2(double x){return std::pow(2.3,x);}static active a_pow2(const active&x){return pow(2.3,x);}
static double r_at2(double x){return std::atan2(x,1.7);}
static active a_at2(const active&x){return atan2(x,1.7);}
static double r_at2b(double x){return std::atan2(1.7,x);}
static active a_at2b(const active&x){return atan2(1.7,x);}
static double r_at2c(double x){return std::atan2(x,x*x+1.0);}
static active a_at2c(const active&x){return atan2(x,x*x+1.0);}
static double r_hyp(double x){return ::hypot(x,2.1);}  static active a_hyp(const active&x){return hypot(x,2.1);}
static double r_hypb(double x){return ::hypot(2.1,x);} static active a_hypb(const active&x){return hypot(2.1,x);}
static double r_hypc(double x){return ::hypot(x,x+1.0);}
static active a_hypc(const active&x){return hypot(x,x+1.0);}
static double r_min(double x){return std::min(x,5.0);} static active a_min(const active&x){return min(x,5.0);}
static double r_minb(double x){return std::min(5.0,x);}static active a_minb(const active&x){return min(5.0,x);}
static double r_minc(double x){return std::min(x,x*x);}static active a_minc(const active&x){return min(x,x*x);}
static double r_max(double x){return std::max(x,0.1);} static active a_max(const active&x){return max(x,0.1);}
static double r_maxb(double x){return std::max(0.1,x);}static active a_maxb(const active&x){return max(0.1,x);}
static double r_maxc(double x){return std::max(x,x*x);}static active a_maxc(const active&x){return max(x,x*x);}
static double r_pluseq(double x){double a=1.5;a+=x;a+=x*x;return a;}
static active a_pluseq(const active&x){active a=1.5;a+=x;a+=x*x;return a;}
static double r_mineq(double x){double a=9.5;a-=x;a-=x*x;return a;}
static active a_mineq(const active&x){active a=9.5;a-=x;a-=x*x;return a;}
static double r_muleq(double x){double a=1.5;a*=x;a*=x+1.0;return a;}
static active a_muleq(const active&x){active a=1.5;a*=x;a*=x+1.0;return a;}
static double r_diveq(double x){double a=7.5;a/=x;a/=x+1.0;return a;}
static active a_diveq(const active&x){active a=7.5;a/=x;a/=x+1.0;return a;}
static double r_incr(double x){double a=x;++a;a++;return a*a;}
static active a_incr(const active&x){active a=x;++a;a++;return a*a;}
static double r_decr(double x){double a=x;--a;a--;return a*a;}
static active a_decr(const active&x){active a=x;--a;a--;return a*a;}
static double r_copy(double x){double a=x;double b=a;return b*b;}
static active a_copy(const active&x){active a=x;active b=a;return b*b;}
static double r_deep(double x){double a=x;for(int i=0;i<25;i++) a = a*0.5 + std::sin(a);return a;}
static active a_deep(const active&x){active a=x;for(int i=0;i<25;i++) a = a*0.5 + sin(a);return a;}

int main( int argc , char ** argv )
{
  (void)argc; (void)argv;

  check_unary("sin"   ,a_sin  ,r_sin  , 0.7);

  /*
   * boltzmann::MPI_rank() answers only while a tape is open -- it goes through
   * the Process, which finalize() destroys -- so it returns -1 here and
   * every rank would call itself rank 0 and print the summary.  Ask MPI,
   * which is still up because finalize() no longer tears it down.
   */
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

  check_unary("cos"   ,a_cos  ,r_cos  , 0.7);
  check_unary("tan"   ,a_tan  ,r_tan  , 0.7);
  check_unary("asin"  ,a_asin ,r_asin , 0.4);
  check_unary("acos"  ,a_acos ,r_acos , 0.4);
  check_unary("atan"  ,a_atan ,r_atan , 0.7);
  check_unary("sinh"  ,a_sinh ,r_sinh , 0.7);
  check_unary("cosh"  ,a_cosh ,r_cosh , 0.7);
  check_unary("tanh"  ,a_tanh ,r_tanh , 0.7);
  check_unary("asinh" ,a_asinh,r_asinh, 0.7);
  check_unary("acosh" ,a_acosh,r_acosh, 1.9);
  check_unary("atanh" ,a_atanh,r_atanh, 0.4);
  check_unary("exp"   ,a_exp  ,r_exp  , 0.7);
  check_unary("log"   ,a_log  ,r_log  , 1.7);
  check_unary("log10" ,a_log10,r_log10, 1.7);
  check_unary("sqrt"  ,a_sqrt ,r_sqrt , 1.7);
  check_unary("cbrt"  ,a_cbrt ,r_cbrt , 1.7);
  check_unary("fabs+" ,a_fabs ,r_fabs , 1.7);
  check_unary("fabs-" ,a_fabs ,r_fabs ,-1.7);
  check_unary("erf"   ,a_erf  ,r_erf  , 0.7);
  check_unary("erfc"  ,a_erfc ,r_erfc , 0.7);
  check_unary("neg"   ,a_neg  ,r_neg  , 0.7);
  check_unary("pos"   ,a_pos  ,r_pos  , 0.7);
  check_unary("mul"   ,a_mul  ,r_mul  , 0.7);
  check_unary("div"   ,a_div  ,r_div  , 0.7);
  check_unary("pow_d" ,a_pow  ,r_pow  , 1.3);
  check_unary("pow_aa",a_powa ,r_powa , 1.3);
  check_unary("pow_da",a_pow2 ,r_pow2 , 1.3);
  check_unary("atan2a",a_at2  ,r_at2  , 0.7);
  check_unary("atan2b",a_at2b ,r_at2b , 0.7);
  check_unary("atan2c",a_at2c ,r_at2c , 0.7);
  check_unary("hypota",a_hyp  ,r_hyp  , 0.7);
  check_unary("hypotb",a_hypb ,r_hypb , 0.7);
  check_unary("hypotc",a_hypc ,r_hypc , 0.7);
  check_unary("min_ad",a_min  ,r_min  , 0.7);
  check_unary("min_da",a_minb ,r_minb , 0.7);
  check_unary("min_aa",a_minc ,r_minc , 0.7);
  check_unary("max_ad",a_max  ,r_max  , 0.7);
  check_unary("max_da",a_maxb ,r_maxb , 0.7);
  check_unary("max_aa",a_maxc ,r_maxc , 0.7);
  check_unary("op+="  ,a_pluseq,r_pluseq,0.7);
  check_unary("op-="  ,a_mineq,r_mineq, 0.7);
  check_unary("op*="  ,a_muleq,r_muleq, 0.7);
  check_unary("op/="  ,a_diveq,r_diveq, 0.7);
  check_unary("prefix",a_incr ,r_incr , 0.7);
  check_unary("postfx",a_decr ,r_decr , 0.7);
  check_unary("copy"  ,a_copy ,r_copy , 0.7);
  check_unary("deep"  ,a_deep ,r_deep , 0.7, 1e-5);

  /*
   * MPI is still up: finalize() ends a tape, not the runtime (see the note in
   * src/API.cpp).  Before that change this reduction was impossible and so
   * was this whole file -- forty-eight tapes in one process could not exist.
   */
  int total_failures = 0;
  MPI_Allreduce(&failures,&total_failures,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(!my_rank){
    std::printf("fdcheck: %d checks, %d failures\n",checks,total_failures);
  }

  return total_failures ? 1 : 0;
}
