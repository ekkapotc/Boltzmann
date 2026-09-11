#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>

#include "../../boltzmann.hpp"

using namespace std;
using namespace boltzmann;

inline double nextU01(double *s)
{
    const double norm = 2.328306549295728e-10;
    const double m1 =  4294967087.0, m2=4294944443.0, a12=1403580.0, a13n=810728.0, a21=527612.0, a23n=1370589.0;
    const double m1inv = 1.0/m1, m2inv = 1.0/m2;
    long k;
    double p1, p2, u;;
    /* Component 1 */
    p1 = a12 * s[1] - a13n * s[0];
    k = long(p1 * m1inv);
    p1 -= k*m1;
    if (p1 < 0.0) p1 += m1;
    s[0] = s[1];
    s[1] = s[2];
    s[2] = p1;
    /* Component 2 */
    p2 = a21 * s[3+2] - a23n * s[3+0];
    k = long(p2 * m2inv);
    p2 -= k * m2;
    if (p2 < 0.0) p2 += m2;
    s[3+0] = s[3+1];
    s[3+1] = s[3+2];
    s[3+2] = p2;

    u = (  (p1 <= p2) ? (p1 - p2 + m1) * norm : (p1 - p2) * norm  );
    return u;
}


template<typename A>
void randNormal(int n, double *state, A * normal)
{
    const double pi = 3.14159265358979323846;
    double y1, y2, u;
    int i;
    for(i=0; i<=n-2; i+=2) {
        u = nextU01(state);
        double x2 = sqrt(-2.0*log(u));
        u = nextU01(state);
        double x = sin( 2.0*pi*u );

        y2 = cos( 2.0*pi*u );
        y2 *= x2;
        y1 = x*x2;
        normal[i] = y1;
        normal[i+1] = y2;
    }
    if(i<n) {
        u = nextU01(state);
        double x2 = sqrt(-2.0*log(u));
        u = nextU01(state);
        double x = sin( 2.0*pi*u );

        y2 = cos( 2.0*pi*u );
        y2 *= x2;
        y1 = x*x2;
        normal[i] = y1;
        i++;
        if(i<n) {
            normal[i] = y2;
        }
    }
}


template<typename A>
void randUniform(int n, double *state, A * P)
{
    for(int i=0; i<n; i++) {
        P[i] = nextU01(state);
    }
}

template<typename ATYPE>
struct LocalVolSurface {
    
    int m, n;
   
    vector<ATYPE> a, b;
    
   explicit LocalVolSurface() : m(0), n(0) { }
    
    LocalVolSurface(istream &in) {
        in >> m >> n;
        a.resize(m+1);
        b.resize(n+1);
        for(int i=0; i<=m; i++) in >> a[i];
        for(int i=0; i<=n; i++) in >> b[i];
    }

    // squared volatilities; usually through interpolation
    // from measured market data
    ATYPE operator()(const ATYPE& x, const ATYPE& t) const { 
      // Horner
      ATYPE ax=a[m];

      for(int i=m-1 ; i>=0 ; i--) ax = a[i] + x*ax;

      ATYPE bx=b[n];

      for(int i=n-1 ; i>=0 ; i--) bx = b[i] + x*bx;

      return t*ax/bx;
    }
};

// collection of independent variables
template<typename ATYPE>
struct ACTIVE_INPUTS {

  ATYPE S0,r,K,T; 

  LocalVolSurface<ATYPE> sigmaSq;

  ACTIVE_INPUTS(const double& S0, 
		const double& r, 
		const double& K,
		const double& T,
		istream& in) 
	  : S0(S0), r(r), K(K), T(T), sigmaSq(in) { }

  ACTIVE_INPUTS(const ATYPE& S0, 
		const ATYPE& r, 
		const ATYPE& K,
		const ATYPE& T,
		const LocalVolSurface<ATYPE>& sigmaSq) 
	  : S0(S0), r(r), K(K), T(T), sigmaSq(sigmaSq) { }

  ~ACTIVE_INPUTS() {}

};

// collection of passive parameters
struct PASSIVE_INPUTS {
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

  void reseed() {
    for(int i=0; i<6; i++) rngseed[i] = i+1; 
  }
};

// collection of dependent variables
template<typename ATYPE>
struct ACTIVE_OUTPUTS {
  ATYPE V;
};

// collection of dependent variables
struct PASSIVE_OUTPUTS {
  double ci;
};

template<typename ATYPE>
inline ATYPE mcpath(
    const ACTIVE_INPUTS<ATYPE> &X,
    const PASSIVE_INPUTS &XP) 
{
  const ATYPE mcdt=X.T/(XP.M-1);
  ATYPE logS=log(X.S0);
  ATYPE t=0;

  std::unique_ptr<double []> Z(new double[XP.M]);
  //static std::vector<double> Z(XP.M);

  randNormal(XP.M, XP.rngseed, Z.get());

  for(int i=0;i<XP.M;i++)
  {
    t+=mcdt;
    ATYPE volS=sqrt(X.sigmaSq(logS,t));
    logS+=(X.r-0.5*volS*volS)*mcdt+volS*sqrt(mcdt)*Z[i];
  }

  ATYPE ST=exp(logS);
  
  return ( ST<X.K ? ATYPE(0) : exp(-X.r*X.T)*(ST-X.K) );
}

template<typename ATYPE>
void price(
    const ACTIVE_INPUTS<ATYPE> &X,
    PASSIVE_INPUTS &XP,
    ACTIVE_OUTPUTS<ATYPE>& Y, 
    PASSIVE_OUTPUTS& YP
) {
    XP.reseed();

    Y.V=0;
    ATYPE sumsq=0;

    for(int p=0; p<XP.N; p++){
 
      //largeint mem_before = get_memory();

      ATYPE y=mcpath(X,XP);
      Y.V+=y;
      sumsq+=y*y;

      //largeint mem_after = get_memory();

      //std::cout << "mem_after-mem_before = " << mem_after-mem_before<< std::endl;
    }

    Y.V=Y.V/XP.N;

    YP.ci = passive_value( (ATYPE)(sqrt((sumsq/XP.N - Y.V*Y.V)/XP.N)) );
}

int N = 1000;     //number of MC paths
int M = 2000;     //number of Euler steps
const double S0=70.0; // asset
const double r=0.07;  // interest
const double K=40.0;  // strike
const double T=1.0;   // maturity

int main(int argc, char* argv[])
{
  if(argc!=5) { cerr << "Please specify input scenario file, e.g. scenario_1.in\n"; return 1; }
  cout.precision(15);
  ifstream vols(argv[1]); 
  if(vols.fail()) { cerr << "Cannot open scenario file '" << argv[1] << "'\n"; return 1; }
  
  N = atoi(argv[2]);
  M = atoi(argv[3]);
 
  int path_size = atoi(argv[4]);

  ACTIVE_INPUTS<active> X(S0,r,K,T,vols); vols.close();
  PASSIVE_INPUTS XP(N,M);
  ACTIVE_OUTPUTS<active> Y;
  PASSIVE_OUTPUTS YP;

  unsigned int xmsz = 4+X.sigmaSq.a.size()+X.sigmaSq.b.size();

  //initialize(xmsz,1,1515576);
  //initialize(xmsz,1,212113680);//scenario 6,M=2000
  //initialize(xmsz,1,2.5*212113680);//scenario 6, M=5000
   
  initialize(xmsz,1,path_size);

  active **XM=new active*[xmsz];
  
  XM[0]=&X.S0; XM[1]=&X.r; XM[2]=&X.K; XM[3]=&X.T;
  
  for (unsigned int i=0;i<X.sigmaSq.a.size();i++)
    XM[i+4]=&X.sigmaSq.a[i];
  for (unsigned int i=0;i<X.sigmaSq.b.size();i++)
    XM[i+4+X.sigmaSq.a.size()]=&X.sigmaSq.b[i];

  for (unsigned int i=0;i<xmsz;i++) independent(*XM[i]);   

  active & YM = Y.V;

  while(checkpoint(XM,YM))
  {
    try
    {
      price(X,XP,Y,YP);
    }catch(BreakException const & e)
    {

    }
  }

  dependent(YM);
  
  cout << "Y = " << YM.val << endl;

  double ** A;
  harvest(1,xmsz,A);
	
  finalize(); 
 
  delete [] XM;

  //cout << "Y=" << Y.V << " in (" << Y.V-3*YP.ci << ", " << Y.V+3*YP.ci << ")" << endl;

  return 0;
}
