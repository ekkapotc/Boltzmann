#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
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

int N=2000; 
int M=50;
const double S0=70.0;
const double r=0.07;
const double K=40.0;
const double T=1.0;

int main( int argc , char* argv[] ) 
{
   if(argc!=4) 
   { 
     cerr << "Please specify input scenario file, e.g. scenario_1.in\n"; 
     return 1; 
   }

   cout.precision(15);

   ifstream vols(argv[1]); 
    
   if(vols.fail())
   { 
     cerr << "Cannot open scenario file '" << argv[1] << "'\n"; 
     return 1;
   }
 
   M = atoi(argv[2]);//M is varied
   N = atoi(argv[3]); //N=2000
  
   ACTIVE_INPUTS<active> X(S0,r,K,T,vols); 

   vols.close();

   PASSIVE_INPUTS XP(N,M);
 
   ACTIVE_OUTPUTS<active> Y;

   std::cout << "m = " << X.sigmaSq.m << std::endl;
   std::cout << "n = " << X.sigmaSq.n << std::endl;

   unsigned int xmsz = 4+X.sigmaSq.a.size()+X.sigmaSq.b.size();

   initialize(xmsz,1,221175912);//N=1000, M is varied
   //initialize(xmsz,1,442143912);//N=2000,M is varied

   active **XM=new active*[xmsz];

   XM[0]=&X.S0; XM[1]=&X.r; XM[2]=&X.K; XM[3]=&X.T;
   
   for (unsigned int i=0;i<X.sigmaSq.a.size();i++)
     XM[i+4]=&X.sigmaSq.a[i];
   
   for (unsigned int i=0;i<X.sigmaSq.b.size();i++)
     XM[i+4+X.sigmaSq.a.size()]=&X.sigmaSq.b[i];
   
  for (unsigned int i=0;i<xmsz;i++) independent(*XM[i]);   

  active & YM = Y.V;//dependent variable  

  CrankNicholson<active> CN(X,XP);

  while(checkpoint(XM,YM))
  {
    try
    {
        price(CN,X,XP,Y);
    }catch(BreakException const & e)
    {

    }
  }

  dependent(YM);

  double ** A;
  harvest(1,xmsz,A);
	
  finalize(); 
 
  delete [] XM;
   
  return 0;
}

