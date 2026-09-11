/* The library does not leak <mpi.h> through boltzmann.hpp, on purpose.  This
 * example reduces its own counts across the ranks, so it includes it itself. */
#include <mpi.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <memory>

#include "../../boltzmann.hpp"

using namespace boltzmann;
using namespace std;

template<typename ATYPE>
void dgttrf(int n, ATYPE *dl, ATYPE *d, ATYPE *du, ATYPE *du2, int *ipiv, int *info)
{
    /* System generated locals */
    int i__1;

    ATYPE d__1, d__2;

    /* Local variables */
    int i__;

    ATYPE fact, temp;


    /*  -- LAPACK routine (version 3.2) -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  DGTTRF computes an LU factorization of a real tridiagonal matrix A */
    /*  using elimination with partial pivoting and row interchanges. */

    /*  The factorization has the form */
    /*     A = L * U */
    /*  where L is a product of permutation and unit lower bidiagonal */
    /*  matrices and U is upper triangular with nonzeros in only the main */
    /*  diagonal and first two superdiagonals. */

    /*  Arguments */
    /*  ========= */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A. */

    /*  DL      (input/output) DOUBLE PRECISION array, dimension (N-1) */
    /*          On entry, DL must contain the (n-1) sub-diagonal elements of */
    /*          A. */

    /*          On exit, DL is overwritten by the (n-1) multipliers that */
    /*          define the matrix L from the LU factorization of A. */

    /*  D       (input/output) DOUBLE PRECISION array, dimension (N) */
    /*          On entry, D must contain the diagonal elements of A. */

    /*          On exit, D is overwritten by the n diagonal elements of the */
    /*          upper triangular matrix U from the LU factorization of A. */

    /*  DU      (input/output) DOUBLE PRECISION array, dimension (N-1) */
    /*          On entry, DU must contain the (n-1) super-diagonal elements */
    /*          of A. */

    /*          On exit, DU is overwritten by the (n-1) elements of the first */
    /*          super-diagonal of U. */

    /*  DU2     (output) DOUBLE PRECISION array, dimension (N-2) */
    /*          On exit, DU2 is overwritten by the (n-2) elements of the */
    /*          second super-diagonal of U. */

    /*  IPIV    (output) INTEGER array, dimension (N) */
    /*          The pivot indices; for 1 <= i <= n, row i of the matrix was */
    /*          interchanged with row IPIV(i).  IPIV(i) will always be either */
    /*          i or i+1; IPIV(i) = i indicates a row interchange was not */
    /*          required. */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -k, the k-th argument had an illegal value */
    /*          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization */
    /*                has been completed, but the factor U is exactly */
    /*                singular, and division by zero will occur if it is used */
    /*                to solve a system of equations. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /* Parameter adjustments */
    --ipiv;
    --du2;
    --du;
    --d;
    --dl;

    /* Function Body */
    *info = 0;
    if (n < 0) {
        *info = -1;
        return ;
    }

    /*     Quick return if possible */

    if (n == 0) {
        return ;
    }

    /*     Initialize IPIV(i) = i and DU2(I) = 0 */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ipiv[i__] = i__;
        /* L10: */
    }
    i__1 = n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
        du2[i__] = 0.;
        /* L20: */
    }

    i__1 = n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if ((d__1 = d[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {

            /*           No row interchange required, eliminate DL(I) */

            if (d[i__] != 0.) {
                fact = dl[i__] / d[i__];
                dl[i__] = fact;
                //d[i__ + 1] -= fact * du[i__];
		d[i__ + 1] = d[i__ + 1] -  fact * du[i__];
            }
        } else {

            /*           Interchange rows I and I+1, eliminate DL(I) */

            fact = d[i__] / dl[i__];
            d[i__] = dl[i__];
            dl[i__] = fact;
            temp = du[i__];
            du[i__] = d[i__ + 1];
            d[i__ + 1] = temp - fact * d[i__ + 1];
            du2[i__] = du[i__ + 1];
            du[i__ + 1] = -fact * du[i__ + 1];
            ipiv[i__] = i__ + 1;
        }
        /* L30: */
    }
    if (n > 1) {
        i__ = n - 1;
        if ((d__1 = d[i__], abs(d__1)) >= (d__2 = dl[i__], abs(d__2))) {
            if (d[i__] != 0.) {
                fact = dl[i__] / d[i__];
                dl[i__] = fact;
                //d[i__ + 1] -= fact * du[i__];
		d[i__ + 1] = d[i__ + 1] - fact * du[i__];
            }
        } else {
            fact = d[i__] / dl[i__];
            d[i__] = dl[i__];
            dl[i__] = fact;
            temp = du[i__];
            du[i__] = d[i__ + 1];
            d[i__ + 1] = temp - fact * d[i__ + 1];
            ipiv[i__] = i__ + 1;
        }
    }

    /*     Check for a zero on the diagonal of U. */

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (d[i__] == 0.) {
            *info = i__;
            goto L50;
        }
        /* L40: */
    }
L50:

    return ;

}

template<typename ATYPE>
void dgtts2(int itrans, int n, int nrhs,
            ATYPE *dl, ATYPE *d, ATYPE *du, ATYPE *du2,
            int *ipiv, ATYPE *b, int ldb)
{
    /* System generated locals */
    int b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    int i__, j, ip;

    ATYPE temp;


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  DGTTS2 solves one of the systems of equations */
    /*     A*X = B  or  A**T*X = B, */
    /*  with a tridiagonal matrix A using the LU factorization computed */
    /*  by DGTTRF. */

    /*  Arguments */
    /*  ========= */

    /*  ITRANS  (input) INTEGER */
    /*          Specifies the form of the system of equations. */
    /*          = 0:  A * X = B  (No transpose) */
    /*          = 1:  A**T* X = B  (Transpose) */
    /*          = 2:  A**T* X = B  (Conjugate transpose = Transpose) */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A. */

    /*  NRHS    (input) INTEGER */
    /*          The number of right hand sides, i.e., the number of columns */
    /*          of the matrix B.  NRHS >= 0. */

    /*  DL      (input) DOUBLE PRECISION array, dimension (N-1) */
    /*          The (n-1) multipliers that define the matrix L from the */
    /*          LU factorization of A. */

    /*  D       (input) DOUBLE PRECISION array, dimension (N) */
    /*          The n diagonal elements of the upper triangular matrix U from */
    /*          the LU factorization of A. */

    /*  DU      (input) DOUBLE PRECISION array, dimension (N-1) */
    /*          The (n-1) elements of the first super-diagonal of U. */

    /*  DU2     (input) DOUBLE PRECISION array, dimension (N-2) */
    /*          The (n-2) elements of the second super-diagonal of U. */

    /*  IPIV    (input) INTEGER array, dimension (N) */
    /*          The pivot indices; for 1 <= i <= n, row i of the matrix was */
    /*          interchanged with row IPIV(i).  IPIV(i) will always be either */
    /*          i or i+1; IPIV(i) = i indicates a row interchange was not */
    /*          required. */

    /*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
    /*          On entry, the matrix of right hand side vectors B. */
    /*          On exit, B is overwritten by the solution vectors X. */

    /*  LDB     (input) INTEGER */
    /*          The leading dimension of the array B.  LDB >= max(1,N). */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /*     Quick return if possible */

    /* Parameter adjustments */
    --dl;
    --d;
    --du;
    --du2;
    --ipiv;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    if (n == 0 || nrhs == 0) {
        return ;
    }

    if (itrans == 0) {

        /*        Solve A*X = B using the LU factorization of A, */
        /*        overwriting each right hand side vector with its solution. */

        if (nrhs <= 1) {
            j = 1;
L10:

            /*           Solve L*x = b. */

            i__1 = n - 1;
            for (i__ = 1; i__ <= i__1; ++i__) {
                ip = ipiv[i__];
                temp = b[i__ + 1 - ip + i__ + j * b_dim1] - dl[i__] * b[ip +
                        j * b_dim1];
                b[i__ + j * b_dim1] = b[ip + j * b_dim1];
                b[i__ + 1 + j * b_dim1] = temp;
                /* L20: */
            }

            /*           Solve U*x = b. */

            //b[n + j * b_dim1] /= d[n];
            b[n + j * b_dim1] = b[n + j * b_dim1]/d[n];
            if (n > 1) {
                b[n - 1 + j * b_dim1] = (b[n - 1 + j * b_dim1] - du[n - 1]
                                         * b[n + j * b_dim1]) / d[n - 1];
            }
            for (i__ = n - 2; i__ >= 1; --i__) {
                b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__
                                       + 1 + j * b_dim1] - du2[i__] * b[i__ + 2 + j * b_dim1]
                                      ) / d[i__];
                /* L30: */
            }
            if (j < nrhs) {
                ++j;
                goto L10;
            }
        } else {
            i__1 = nrhs;
            for (j = 1; j <= i__1; ++j) {

                /*              Solve L*x = b. */

                i__2 = n - 1;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    if (ipiv[i__] == i__) {
                        //b[i__ + 1 + j * b_dim1] -= dl[i__] * b[i__ + j *b_dim1];
                        b[i__ + 1 + j * b_dim1] =  b[i__ + 1 + j * b_dim1] - dl[i__] * b[i__ + j *b_dim1];                                 
                    } else {
                        temp = b[i__ + j * b_dim1];
                        b[i__ + j * b_dim1] = b[i__ + 1 + j * b_dim1];
                        b[i__ + 1 + j * b_dim1] = temp - dl[i__] * b[i__ + j *
                                                  b_dim1];
                    }
                    /* L40: */
                }

                /*              Solve U*x = b. */

                //b[n + j * b_dim1] /= d[n];
 	        b[n + j * b_dim1] = b[n + j * b_dim1]/d[n];
                if (n > 1) {
                    b[n - 1 + j * b_dim1] = (b[n - 1 + j * b_dim1] - du[n
                                             - 1] * b[n + j * b_dim1]) / d[n - 1];
                }
                for (i__ = n - 2; i__ >= 1; --i__) {
                    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[
                                               i__ + 1 + j * b_dim1] - du2[i__] * b[i__ + 2 + j *
                                                       b_dim1]) / d[i__];
                    /* L50: */
                }
                /* L60: */
            }
        }
    } else {

        /*        Solve A**T * X = B. */

        if (nrhs <= 1) {

            /*           Solve U**T*x = b. */

            j = 1;
L70:
            //b[j * b_dim1 + 1] /= d[1];
	    b[j * b_dim1 + 1]  =  b[j * b_dim1 + 1]/d[1];
            if (n > 1) {
                b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * b_dim1
                                     + 1]) / d[2];
            }
            i__1 = n;
            for (i__ = 3; i__ <= i__1; ++i__) {
                b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] * b[
                                           i__ - 1 + j * b_dim1] - du2[i__ - 2] * b[i__ - 2 + j *
                                                   b_dim1]) / d[i__];
                /* L80: */
            }

            /*           Solve L**T*x = b. */

            for (i__ = n - 1; i__ >= 1; --i__) {
                ip = ipiv[i__];
                temp = b[i__ + j * b_dim1] - dl[i__] * b[i__ + 1 + j * b_dim1]
                       ;
                b[i__ + j * b_dim1] = b[ip + j * b_dim1];
                b[ip + j * b_dim1] = temp;
                /* L90: */
            }
            if (j < nrhs) {
                ++j;
                goto L70;
            }

        } else {
            i__1 = nrhs;
            for (j = 1; j <= i__1; ++j) {

                /*              Solve U**T*x = b. */

                //b[j * b_dim1 + 1] /= d[1];
	        b[j * b_dim1 + 1] = b[j * b_dim1 + 1] / d[1];
                if (n > 1) {
                    b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j *
                                         b_dim1 + 1]) / d[2];
                }
                i__2 = n;
                for (i__ = 3; i__ <= i__2; ++i__) {
                    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] *
                                           b[i__ - 1 + j * b_dim1] - du2[i__ - 2] * b[i__ -
                                                   2 + j * b_dim1]) / d[i__];
                    /* L100: */
                }
                for (i__ = n - 1; i__ >= 1; --i__) {
                    if (ipiv[i__] == i__) {
                        //b[i__ + j * b_dim1] -= dl[i__] * b[i__ + 1 + j * b_dim1];
                        b[i__ + j * b_dim1] =  b[i__ + j * b_dim1] - dl[i__] * b[i__ + 1 + j * b_dim1];                                
                    } else {
                        temp = b[i__ + 1 + j * b_dim1];
                        b[i__ + 1 + j * b_dim1] = b[i__ + j * b_dim1] - dl[
                                                      i__] * temp;
                        b[i__ + j * b_dim1] = temp;
                    }
                    /* L110: */
                }
                /* L120: */
            }
        }
    }

}

template<typename ATYPE>
void dgttrs(const char *trans, int n, int nrhs,
            ATYPE *dl, ATYPE *d, ATYPE *du, ATYPE *du2,
            int *ipiv, ATYPE *b, int ldb, int *info)
{
    /* System generated locals */
    int b_dim1, b_offset;

    /* Local variables */
    int itrans;
    bool notran;


    /*  -- LAPACK routine (version 3.2) -- */
    /*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
    /*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  DGTTRS solves one of the systems of equations */
    /*     A*X = B  or  A**T*X = B, */
    /*  with a tridiagonal matrix A using the LU factorization computed */
    /*  by DGTTRF. */

    /*  Arguments */
    /*  ========= */

    /*  TRANS   (input) CHARACTER*1 */
    /*          Specifies the form of the system of equations. */
    /*          = 'N':  A * X = B  (No transpose) */
    /*          = 'T':  A**T* X = B  (Transpose) */
    /*          = 'C':  A**T* X = B  (Conjugate transpose = Transpose) */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A. */

    /*  NRHS    (input) INTEGER */
    /*          The number of right hand sides, i.e., the number of columns */
    /*          of the matrix B.  NRHS >= 0. */

    /*  DL      (input) DOUBLE PRECISION array, dimension (N-1) */
    /*          The (n-1) multipliers that define the matrix L from the */
    /*          LU factorization of A. */

    /*  D       (input) DOUBLE PRECISION array, dimension (N) */
    /*          The n diagonal elements of the upper triangular matrix U from */
    /*          the LU factorization of A. */

    /*  DU      (input) DOUBLE PRECISION array, dimension (N-1) */
    /*          The (n-1) elements of the first super-diagonal of U. */

    /*  DU2     (input) DOUBLE PRECISION array, dimension (N-2) */
    /*          The (n-2) elements of the second super-diagonal of U. */

    /*  IPIV    (input) INTEGER array, dimension (N) */
    /*          The pivot indices; for 1 <= i <= n, row i of the matrix was */
    /*          interchanged with row IPIV(i).  IPIV(i) will always be either */
    /*          i or i+1; IPIV(i) = i indicates a row interchange was not */
    /*          required. */

    /*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
    /*          On entry, the matrix of right hand side vectors B. */
    /*          On exit, B is overwritten by the solution vectors X. */

    /*  LDB     (input) INTEGER */
    /*          The leading dimension of the array B.  LDB >= max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -i, the i-th argument had an illegal value */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /* Parameter adjustments */
    --dl;
    --d;
    --du;
    --du2;
    --ipiv;
    b_dim1 =ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = *(unsigned char *)trans == 'N' || *(unsigned char *)trans == 'n';
    if (! notran && ! (*(unsigned char *)trans == 'T' || *(unsigned char *)
                       trans == 't') && ! (*(unsigned char *)trans == 'C' || *(unsigned
                                           char *)trans == 'c')) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (nrhs < 0) {
        *info = -3;
    } else if (( ldb < n) && (ldb<1)) {
        *info = -10;
    }
    if (*info != 0) {
        return ;
    }

    /*     Quick return if possible */

    if (n == 0 || nrhs == 0) {
        return;
    }

    /*     Decode TRANS */

    if (notran) {
        itrans = 0;
    } else {
        itrans = 1;
    }

    /*     Determine the number of right-hand sides to solve at a time. */

    dgtts2(itrans, n, nrhs, &dl[1], &d[1], &du[1], &du2[1], &ipiv[1],
           &b[b_offset], ldb);
}

template<typename ATYPE>
struct LocalVolSurface
{
    int m, n;

    vector<ATYPE> a, b;

    explicit LocalVolSurface() : m(0), n(0) { }

    LocalVolSurface(istream &in) 
   {
        in >> m >> n;
       
        a.resize(m+1);
        b.resize(n+1);
       
        for(int i=0; i<=m; i++){
	  in >> a[i];
        }

        for(int i=0; i<=n; i++){
	  in >> b[i];
        }
        
    }

    // squared volatilities; usually through interpolation
    // from measured market data
    ATYPE operator()(const ATYPE& x, const ATYPE& t) const
   { 
      // Horner
      ATYPE ax=a[m];
      
      for(int i=m-1; i>=0; i--) ax = a[i] + x*ax;
   
      ATYPE bx=b[n];
  
      for(int i=n-1;i>=0;i--) bx = b[i] + x*bx;
 
      return t*ax/bx;
    }
};

// collection of independent variables
template<typename ATYPE>
struct ACTIVE_INPUTS
{
  ATYPE S0,r,K,T; 
  
  LocalVolSurface<ATYPE> sigmaSq;
  
  ACTIVE_INPUTS(const double& S0, 
		const double& r, 
		const double& K,
		const double& T,
		istream& in) 
	  : S0(S0), r(r), K(K), T(T), sigmaSq(in) 
 { 

 }
	  
  ACTIVE_INPUTS(const ATYPE& S0, 
		const ATYPE& r, 
		const ATYPE& K,
		const ATYPE& T,
		const LocalVolSurface<ATYPE>& sigmaSq) 
	  : S0(S0), r(r), K(K), T(T), sigmaSq(sigmaSq) 
  { 

  }
	  
  ~ACTIVE_INPUTS()
  {
  
  }
};

// collection of passive parameters
struct PASSIVE_INPUTS
{
  int N;
  int M;
  mutable double rngseed[6];
  
  PASSIVE_INPUTS(
	 	const int& N, 
		const int& M 
		) : N(N), M(M) 
   {
       for(int i=0; i<6; i++) rngseed[i] = i+1; 
   }
   
  ~PASSIVE_INPUTS() { }
  
  void reseed() 
  {
    for(int i=0; i<6; i++) rngseed[i] = i+1; 
  }
};

// collection of dependent variables
template<typename ATYPE>
struct ACTIVE_OUTPUTS 
{
  ATYPE V;
};

// collection of dependent variables
struct PASSIVE_OUTPUTS
{
  double ci;
};

template<typename ATYPE>
class CrankNicholson 
{
  public:

    const ACTIVE_INPUTS<ATYPE>& X;
    
    const PASSIVE_INPUTS& XP;

    ATYPE *LHSj_dl,*LHSj_d,*LHSj_du,*RHSj;

    ATYPE xmin, xmax, tmax, *Vprev, *Vcurr;

    double tmin;

    CrankNicholson(const ACTIVE_INPUTS<ATYPE>& X,
		             const PASSIVE_INPUTS& XP) :
 	X(X), 
        XP(XP), 
        LHSj_dl(new ATYPE[XP.N-1]),
        LHSj_d(new ATYPE[XP.N]),
        LHSj_du(new ATYPE[XP.N-1]),
        Vprev( new ATYPE[XP.N+2]), 
        Vcurr( new ATYPE[XP.N+2]) 
       {
        /*    std::unique_ptr<ATYPE []> LHSj_dl(new ATYPE[XP.N-1]);
	    std::unique_ptr<ATYPE []> LHSj_d(new ATYPE[XP.N]);
	    std::unique_ptr<ATYPE []> LHSj_du(new ATYPE[XP.N-1]);
            std::unique_ptr<ATYPE []> Vprev(new ATYPE[XP.N+2]);
	    std::unique_ptr<ATYPE []> Vcurr(new ATYPE[XP.N+2]);*/
      }
    
    ~CrankNicholson()
    {
      delete [] LHSj_dl;
      delete [] LHSj_d;
      delete [] LHSj_du;
      delete [] Vcurr;
      delete [] Vprev;
    }

    void swapVs() {
      ATYPE * tmp = Vcurr;
      Vcurr = Vprev;
      Vprev = tmp;
    }

  private:
    CrankNicholson();
    CrankNicholson(const CrankNicholson&);
    CrankNicholson operator=(const CrankNicholson&);

  public:
    void setBoundaryConditions(int ts, ATYPE *V);
    void prepareLHS(int ts); 
    void prepareRHS(int); 
    void solveTridiagonalSystem();
};

template<typename ATYPE>
void CrankNicholson<ATYPE>::setBoundaryConditions(int ts, ATYPE * V) 
{
    const ATYPE dx = (xmax-xmin)/(XP.N+1);

    const ATYPE dt = X.T/(XP.M+1);

    if(ts==XP.M+1) 
    {
        for(int i=1; i<XP.N+1; i++) 
        {
            const ATYPE xi = xmin + i*dx;
            const ATYPE Si = exp(xi);
        
            if(Si-X.K > 0.0)
                V[i] = Si-X.K;
            else
                V[i] = 0.0;
        }
   }

   V[0] = 0.0;

   const ATYPE Smax = exp(xmax);
   const ATYPE t_j = ts*dt;

   V[XP.N+1] = Smax - exp(-X.r*(X.T-t_j))*X.K;
}

template<typename ATYPE>
void CrankNicholson<ATYPE>::prepareLHS(const int j_time)
{
    const int j = j_time;
    const ATYPE dt = (tmax-tmin)/(XP.M+1);
    const ATYPE t = tmin + j*dt;
    const ATYPE dx = (xmax-xmin)/(XP.N+1);
    const ATYPE alpha = dt/(dx*dx);
    const ATYPE z = 1.0 + X.r*dt;

    int row = 0;
    
    ATYPE x = xmin + (row+1)*dx;
    ATYPE u_ij, l_ij, c_ij;
    ATYPE vhat = X.sigmaSq(x,t);
    u_ij = 0.5*alpha*( vhat + dx*(X.r - 0.5*vhat) ); 
    c_ij = -alpha*vhat;                            
    LHSj_du[row] =   - 0.5*u_ij;
    LHSj_d[row]  = z - 0.5*c_ij;
    
    for(row=1; row<XP.N-1; row++) 
    {
        x = xmin + (row+1)*dx;
        vhat = X.sigmaSq(x,t);
        u_ij = 0.5*alpha*( vhat + dx*(X.r - 0.5*vhat) ); 
        l_ij = 0.5*alpha*( vhat - dx*(X.r - 0.5*vhat) ); 
        c_ij = -alpha*vhat;                            
        LHSj_dl[row-1] =   - 0.5*l_ij;
        LHSj_d[row]    = z - 0.5*c_ij;
        LHSj_du[row]   =   - 0.5*u_ij;
    }

    row = XP.N-1;
    x = xmin + (row+1)*dx;
    vhat = X.sigmaSq(x,t);
    l_ij = 0.5*alpha*( vhat - dx*(X.r - 0.5*vhat) ); 
    c_ij = -alpha*vhat;                            
    LHSj_dl[row-1] =   - 0.5*l_ij;
    LHSj_d[row]    = z - 0.5*c_ij;
}

template<typename ATYPE>
void CrankNicholson<ATYPE>::prepareRHS(const int j_time)
{
    const int j = j_time;
    const ATYPE dt = (tmax-tmin)/(XP.M+1);
    const ATYPE t = tmin + j*dt;
    const ATYPE dx = (xmax-xmin)/(XP.N+1);
    const ATYPE alpha = dt/(dx*dx);

    //allocate heap variables locally
    std::unique_ptr<ATYPE []> RHSj_d(new ATYPE[XP.N]);
    std::unique_ptr<ATYPE []> RHSj_du(new ATYPE[XP.N-1]);
    std::unique_ptr<ATYPE []> RHSj_dl(new ATYPE[XP.N-1]);

    //ATYPE* RHSj_d   = new ATYPE[XP.N];
    //ATYPE* RHSj_du  = new ATYPE[XP.N-1];
    //ATYPE* RHSj_dl  = new ATYPE[XP.N-1];

    int row = 0;

    ATYPE c_ij, u_ij, l_ij;
    ATYPE x = xmin + (row+1)*dx;
    ATYPE vhat = X.sigmaSq(x,t);
    u_ij = 0.5*alpha*( vhat + dx*(X.r - 0.5*vhat) ); 
    c_ij = -alpha*vhat;                            
    RHSj_du[row] = 0.5*u_ij;
    RHSj_d[row]  = 0.5*c_ij + 1.0;

    for(row=1; row<XP.N-1; row++) 
    {
        x = xmin + (row+1)*dx;
        vhat = X.sigmaSq(x,t);
        u_ij = 0.5*alpha*( vhat + dx*(X.r - 0.5*vhat) ); 
        l_ij = 0.5*alpha*( vhat - dx*(X.r - 0.5*vhat) ); 
        c_ij = -alpha*vhat;                            
        RHSj_dl[row-1] = 0.5*l_ij;
        RHSj_d[row]    = 0.5*c_ij + 1.0;
        RHSj_du[row]   = 0.5*u_ij;
    }

    row = XP.N-1;
    x = xmin + (row+1)*dx;
    vhat = X.sigmaSq(x,t);
    l_ij = 0.5*alpha*( vhat - dx*(X.r - 0.5*vhat) ); 
    c_ij = -alpha*vhat;                            
    RHSj_dl[row-1] = 0.5*l_ij;
    RHSj_d[row]    = 0.5*c_ij + 1.0;

    row = 0;
    RHSj[row] = RHSj_d[row]*Vprev[1] + RHSj_du[row]*Vprev[2];
    x = xmin + (row+1)*dx;
    vhat = X.sigmaSq(x,t);
    l_ij = 0.5*alpha*( vhat - dx*(X.r - 0.5*vhat) ); 
    RHSj[row] += l_ij*0.5*( Vcurr[0] + Vprev[0] );

    for(row=1; row<XP.N-1; row++) 
        RHSj[row] = RHSj_dl[row-1]*Vprev[row] + RHSj_d[row]*Vprev[row+1] + RHSj_du[row]*Vprev[row+2];

    row = XP.N-1;
    x = xmin + (row+1)*dx;
    RHSj[row] = RHSj_dl[row-1]*Vprev[row] + RHSj_d[row]*Vprev[row+1];
    vhat = X.sigmaSq(x,t);
    u_ij = 0.5*alpha*( vhat + dx*(X.r - 0.5*vhat) ); 
    RHSj[row] += u_ij*0.5*( Vcurr[XP.N+1] + Vprev[XP.N+1] );

    //deallocate heap variables locally
    //delete[] RHSj_d;
    //delete[] RHSj_dl;
    //delete[] RHSj_du;
}

template<typename ATYPE>
void
CrankNicholson<ATYPE>::solveTridiagonalSystem () 
{
    //int * LHSj_ipiv = new int[XP.N];
    //ATYPE * LHSj_du2 = new ATYPE[XP.N-2];

    std::unique_ptr<int []>  LHSj_ipiv(new int[XP.N]);
    std::unique_ptr<ATYPE []> LHSj_du2(new ATYPE[XP.N-2]);

    int info = 0;

    dgttrf(XP.N, LHSj_dl, LHSj_d, LHSj_du, LHSj_du2.get(), LHSj_ipiv.get(), &info);
    dgttrs("N", XP.N, 1, LHSj_dl, LHSj_d, LHSj_du, LHSj_du2.get(), LHSj_ipiv.get(), RHSj, XP.N, &info);

    //delete[] LHSj_ipiv;
    //delete[] LHSj_du2;
}

template<typename ATYPE>
void price(
                CrankNicholson<ATYPE> & CN,
                const ACTIVE_INPUTS<ATYPE>& X,
                const PASSIVE_INPUTS& XP,
                ACTIVE_OUTPUTS<ATYPE>& Y) 
{
    //CrankNicholson<ATYPE> CN(X,XP);
    const double tmin = 0.0; CN.tmin=tmin;
    const ATYPE tmax = X.T; CN.tmax=tmax;
    const ATYPE logS0 = log(X.S0);
    const ATYPE atmvol = sqrt( X.T*X.sigmaSq(logS0,X.T));

    const int C=10;
    const ATYPE dx = 2*C*atmvol/(XP.N+1);
    const ATYPE xmin = logS0 - XP.N/2*dx; CN.xmin=xmin;
    const ATYPE xmax = xmin + (XP.N+1)*dx; CN.xmax=xmax;

    CN.setBoundaryConditions(XP.M+1, CN.Vprev);
    
    for(int j=XP.M; j>=0; j--) 
   {
        CN.setBoundaryConditions(j, CN.Vcurr);
        CN.prepareLHS(j);
        CN.RHSj=&CN.Vcurr[1];
        CN.prepareRHS(j);
        CN.solveTridiagonalSystem();
        CN.swapVs();
    }

    Y.V=CN.Vprev[XP.N/2];
}

int N=40;
int M=20;
const double S0=70.0;
const double r=0.07;
const double K=40.0;
const double T=1.0;

/*
 * THE BUILT-IN SCENARIO.  Byte for byte scenario_3.in, which is what this
 * example ran with before: a local volatility surface as a ratio of two
 * degree-14 polynomials in log-spot, a(x)/b(x), scaled by t.  It is embedded
 * so that `./pde` with no arguments is a complete run -- the example used to
 * exit 1 with "Please specify input scenario file", which made it useless to
 * `make test` and invisible to `make ranks`.  A file named on the command
 * line still wins, and the six scenario_*.in files are still there.
 */
static const char DEFAULT_SCENARIO[] =
  "14 14\n"
  /* a[0..14] */
  "0.05\n"
  "0\n"
  "0.1375\n"
  "0\n"
  "0.135577\n"
  "0\n"
  "0.0639022\n"
  "0\n"
  "0.0154429\n"
  "0\n"
  "0.00184878\n"
  "0\n"
  "9.39685e-05\n"
  "0\n"
  "1.29428e-06\n"
  /* b[0..14] */
  "1\n"
  "0\n"
  "1.75\n"
  "0\n"
  "1.21154\n"
  "0\n"
  "0.420673\n"
  "0\n"
  "0.076486\n"
  "0\n"
  "0.00688374\n"
  "0\n"
  "0.000254953\n"
  "0\n"
  "2.27637e-06\n";

/*
 * One tape, one budget, one gradient.
 *
 * Called TWICE from main() with different byte budgets.  That is the check:
 * the Jacobian of a program is a property of the program, not of how much
 * memory the tape was allowed, so a run broken into a hundred partitions and
 * a run broken into a handful must agree to rounding.  It is the only
 * reference this example can have -- there is no closed form for a
 * Crank-Nicholson price under a rational local-volatility surface -- but it
 * is a real one, because the two runs eliminate DIFFERENT graphs: chunking
 * changes which vertices survive a partition boundary and therefore the whole
 * elimination order.
 */
static void run_one( const std::string & scen ,
                     int Mv , int Nv , largeint budget ,
                     std::vector<double> & g ,
                     largeint & parts , largeint & total ,
                     largeint & cost  , largeint & passes ,
                     int & ranks )
{
  istringstream in(scen);

  ACTIVE_INPUTS<active> X(S0,r,K,T,in);
  PASSIVE_INPUTS        XP(Nv,Mv);
  ACTIVE_OUTPUTS<active> Y;

  const unsigned int xmsz = 4+X.sigmaSq.a.size()+X.sigmaSq.b.size();

  initialize( (largeint)xmsz , 1 , budget );

  std::vector<active*> XM(xmsz);
  XM[0]=&X.S0; XM[1]=&X.r; XM[2]=&X.K; XM[3]=&X.T;

  for( unsigned int i=0 ; i<X.sigmaSq.a.size() ; i++ )
    XM[i+4]=&X.sigmaSq.a[i];

  for( unsigned int i=0 ; i<X.sigmaSq.b.size() ; i++ )
    XM[i+4+X.sigmaSq.a.size()]=&X.sigmaSq.b[i];

  for( unsigned int i=0 ; i<xmsz ; i++ ) independent(*XM[i]);

  active & YM = Y.V;

  /*
   * CN IS CONSTRUCTED INSIDE THE PASS, AND THAT IS NOT A STYLE CHOICE.
   *
   * It used to be built once, outside the checkpoint loop, and that was safe
   * only for as long as this example never chunked: with a 221 MB budget
   * there was one partition, one pass, and no boundary to survive.  Give it a
   * budget that chunks and it corrupts the heap -- AddressSanitizer reports a
   * use-after-free in Vertex::kill(), reached from active::operator= inside
   * price(), on memory released by the previous partition's elimination.
   *
   * The rule it breaks: an `active` that is not an independent must not live
   * from one pass to the next.  checkpoint() restores the independents and
   * starts a new tape; every other active is left holding a pointer into the
   * graph that was just freed.  CN's Vprev, Vcurr and LHSj_* arrays are
   * exactly that -- and worse, they are READ before they are written on a
   * later pass, so the stale graph node is spliced into the new graph.  The
   * same mistake in a program whose arrays happen to be fully overwritten
   * does not crash; it silently returns the wrong derivative.
   *
   * Everything a pass needs, a pass builds.
   */
  passes = run_tape( &XM[0] , YM , [&]{
    CrankNicholson<active> CN(X,XP);
    price(CN,X,XP,Y);
  });

  dependent(YM);

  /* The two-argument harvest, not the three-argument one: the old call went
   * through the printing path and dumped the whole gradient plus an
   * intmed_map size to the terminal, once per rank, whatever you wanted. */
  Jacobian J = harvest( 1 , (largeint)xmsz );

  g.assign(xmsz,0.0);
  if(!J.empty()){
    for( unsigned int i=0 ; i<xmsz ; i++ ) g[i] = J(0,i);
  }

  parts  = get_partitions();
  total  = get_total_partitions();
  cost   = get_cost();
  ranks  = MPI_size();//while the tape is still open: it is a tape query

  finalize();
}

int main( int argc , char* argv[] )
{
  /*
   * Everything is optional now.  `./pde` runs the built-in scenario; a file
   * name still overrides it, and "-" means "the built-in one" so that M, N
   * and the budget can be given without naming a file.
   *
   *     ./pde                                  built-in, M=50, N=100
   *     ./pde scenario_5.in 100 200            a file, bigger grid
   *     ./pde - 50 100 40000                   built-in, tighter budget
   */
  const char *   file   = (argc>1 && std::string(argv[1])!="-") ? argv[1] : 0;
  const int      Mv     = (argc>2) ? atoi(argv[2]) : 20;
  const int      Nv     = (argc>3) ? atoi(argv[3]) : 40;
  const largeint budget = (argc>4) ? (largeint)atol(argv[4]) : 300000;

  if( Mv<1 || Nv<4 ){ cerr << "pde: need M>=1 and N>=4\n"; return 1; }

  M = Mv;
  N = Nv;

  std::string scen;

  if(file){
    ifstream f(file);
    if(f.fail()){ cerr << "pde: cannot open scenario file '" << file << "'\n"; return 1; }
    scen.assign( istreambuf_iterator<char>(f) , istreambuf_iterator<char>() );
  }else{
    scen = DEFAULT_SCENARIO;
  }

  /* ---- the run under test: chunked ------------------------------------ */

  std::vector<double> gc;
  largeint parts=0, total=0, cost=0, passes=0;
  int      lib_ranks=0;
  run_one( scen , Mv , Nv , budget , gc , parts , total , cost , passes , lib_ranks );

  /* ---- the same program, told it has far more memory ------------------ */

  std::vector<double> gs;
  largeint sparts=0, stotal=0, scost=0, spasses=0;
  int      sranks=0;
  run_one( scen , Mv , Nv , budget*64 , gs , sparts , stotal , scost , spasses , sranks );

  /* ---- gather --------------------------------------------------------- */

  int rank=0, size=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  const int nin = (int)gc.size();

  std::vector<double> g(nin,0.0), gref(nin,0.0);
  MPI_Allreduce(&gc[0],&g[0]   ,nin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);//one contributor
  MPI_Allreduce(&gs[0],&gref[0],nin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  long lp=(long)parts, total_parts=0, min_parts=0, max_parts=0, ranks_working=0;
  long one = (parts>0) ? 1 : 0;
  MPI_Allreduce(&lp ,&total_parts  ,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&min_parts    ,1,MPI_LONG,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&lp ,&max_parts    ,1,MPI_LONG,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&one,&ranks_working,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  long lc=(long)cost, total_cost=0;
  MPI_Allreduce(&lc,&total_cost,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  long ls=(long)sparts, few_parts=0;
  MPI_Allreduce(&ls,&few_parts,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  /* ---- the check ------------------------------------------------------ */

  double scale=0.0, worst=0.0;
  int    where=0;
  for( int i=0 ; i<nin ; i++ ) if(fabs(gref[i])>scale) scale=fabs(gref[i]);
  for( int i=0 ; i<nin ; i++ ){
    const double d = fabs(g[i]-gref[i]);
    if(d>worst){ worst=d; where=i; }
  }
  const double rel = (scale>0.0) ? worst/scale : worst;

  /* Two eliminations of two different graphs reaching the same numbers: the
   * difference is accumulated rounding, not method error.  1e-10 is five
   * orders above what is observed. */
  const double TOL = 1.0e-10;

  int bad = 0;

  if( !(rel<=TOL) ){
    bad++;
    if(!rank) printf("  FAIL chunked vs unchunked: %.3e relative at input %d, tol %.1e\n",
                     rel,where,TOL);
  }
  if( total != (largeint)total_parts ){
    bad++;
    if(!rank) printf("  FAIL library reports %lu partitions, the ranks recorded %ld\n",
                     (unsigned long)total,total_parts);
  }
  if( lib_ranks != size ){
    bad++;
    if(!rank) printf("  FAIL library reports %d ranks, MPI reports %d\n",lib_ranks,size);
  }
  if( total_parts < 8 ){
    bad++;
    if(!rank) printf("  FAIL only %ld partitions -- the tape is not being chunked.  "
                     "Lower the budget or raise M.\n",total_parts);
  }
  if( total_parts <= few_parts ){
    bad++;
    if(!rank) printf("  FAIL the two budgets gave %ld and %ld partitions -- "
                     "the comparison is not comparing anything.\n",total_parts,few_parts);
  }
  if( size>1 && ranks_working<2 ){
    bad++;
    if(!rank) printf("  FAIL %ld partitions all landed on one rank of %d -- "
                     "the pipeline is not being exercised.\n",total_parts,size);
  }
  if( passes != parts+1 ){
    bad++;
    if(!rank) printf("  FAIL rank %d ran %lu passes for %lu partitions\n",
                     rank,(unsigned long)passes,(unsigned long)parts);
  }

  /* ---- report --------------------------------------------------------- */

  if(!rank){
    printf("pde: %lu partitions over %d MPI rank%s\n",
           (unsigned long)total,lib_ranks,(lib_ranks==1)?"":"s");
    printf("  problem      Crank-Nicholson, M=%d time steps, N=%d log-spot nodes, %s\n",
           Mv,Nv,file?file:"built-in scenario");
    printf("  gradient     1 x %d (S0, r, K, T and %d surface coefficients)\n",nin,nin-4);
    printf("  partitions   %ld recorded in total, %ld..%ld per rank, "
           "%ld of %d ranks working\n",
           total_parts,min_parts,max_parts,ranks_working,size);
    printf("  elim cost    %ld total\n",total_cost);
    printf("  budget %lu vs %lu bytes: %ld partitions vs %ld, gradients agree to %.3e  (tol %.1e)\n",
           (unsigned long)budget,(unsigned long)(budget*64),total_parts,few_parts,rel,TOL);
    printf("  dV/dS0 = %.12e   dV/dr = %.12e\n",g[0],g[1]);
    printf("  dV/dK  = %.12e   dV/dT = %.12e\n",g[2],g[3]);
    printf("pde: %s\n",(bad==0)?"PASS":"FAIL");
  }

  return (bad==0) ? 0 : 1;
}
