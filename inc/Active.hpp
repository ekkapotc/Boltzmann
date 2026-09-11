#ifndef ACTIVE_INCLUDE
#define ACTIVE_INCLUDE

#include "Typedefs.hpp"

#include <iostream>

namespace boltzmann
{ 
  class active;//forward declaration of class active
}

namespace boltzmann
{
  namespace internals//delcare subnamespace internals
  {
    class Vertex;//forward declaration of class Vertex
  }
}

namespace boltzmann
{  
  bool operator!( const active & x );

  double passive_value( const active & x );

  std::ostream& operator<< ( std::ostream & out, const active & x );
  std::istream& operator>> ( std::istream & in, const active & x );

  const active operator+( const active & x );
  const active operator-( const active & x);

  const active operator+(const active & x1, const active & x2); 
  const active operator+(const active & x1, double x2);
  const active operator+(double x1,  const active & x2);
   
  const active operator-(const active & x1, const active & x2);      
  const active operator-(double x1,  const active & x2);
  const active operator-(const active & x1, double x2);

  const active operator*(const active & x1, const active & x2);
  const active operator*(const active & x1, double x2);
  const active operator*(double x1,  const active & x2);
 
  const active operator/( const active & x1, const active & x2);
  const active operator/( double x1, const active & x2);
  const active operator/( const active & x1, double x2);

  const active sin( const active & x );
  const active cos( const active & x );
  const active tan( const active & x );

  const active asin( const active & x );
  const active acos( const active & x );
  const active atan( const active & x );

  const active sinh( const active & x );
  const active cosh( const active & x );
  const active tanh( const active & x );

  const active asinh( const active & x );
  const active acosh( const active & x );
  const active atanh( const active & x );
 
  const active atan2( const active & x1 , const active & x2 );
  const active atan2( const active & x1 , double x2 );
  const active atan2( double x1 , const active & x2 );

  const active hypot( const active & x1 , const active & x2 );
  const active hypot( const active & x1 , double x2 );
  const active hypot( double x1 , const active & x2 );

  const active exp( const active & x);

  const active pow( const active & x1 , const active & x2 );
  const active pow( const active & x1 , double x2 );
  const active pow( double x1 , const active & x2 );

  const active sqrt( const active & x );
  const active cbrt( const active & x );

  const active log( const active & x );
  const active log10( const active & x );

  const active min( const active & x1 , const active & x2 );
  const active min( const active & x1 , double x2 );
  const active min( double x1 , const active & x2 );

  const active max( const active & x1 , const active & x2 );
  const active max( const active & x1 , double x2 );
  const active max( double x1 , const active & x2 );

  const active fabs( const active & x );

  const active erf( const active & x );
  const active erfc( const active & x );

  double abs( const active  & x );

  bool operator==(const active & x1,const active & x2);
  bool operator==(const active & x1,double  x2);
  bool operator==(double x1,const active & x2);

  bool operator!=(const active & x1,const active & x2);
  bool operator!=(const active & x1,double x2);
  bool operator!=(double x1,const active & x2);

  bool operator<(const active & x1, const active & x2);
  bool operator<(const active & x1, double x2);
  bool operator<(double x1, const active & x2);

  bool operator>(const active & x1,const active & x2);
  bool operator>(const active & x1,double x2);
  bool operator>(double x1,const active & x2);
		
  bool operator<=(const active & x1,const active & x2);
  bool operator<=(const active & x1,double x2);
  bool operator<=(double x1,const active & x2);

  bool operator>=(const active & x1,const active & x2);
  bool operator>=(const active & x1,double x2);
  bool operator>=(double x1,const active & x2);
}

namespace boltzmann
{

/*
 * Ported from Maxwell's header-hygiene step.  This header used to say
 *
 *     using namespace boltzmann::internals;
 *
 * at GLOBAL scope, so every translation unit that included boltzmann.hpp had the
 * library's internal namespace dumped into it -- a name-collision generator,
 * and it made the namespace decorative.  The one name it existed to supply is
 * qualified below instead.  Say `using namespace boltzmann;` in your own
 * translation unit, or qualify: boltzmann::active, boltzmann::sin.
 */
class active
{
public:

  mutable bool  reachable;
  mutable largeint idx;
  mutable largeint  owner_idx;	
  mutable largeint old_idx;
  mutable double val;
  mutable internals::Vertex * vtx;

  /*
   * WHICH PASS THIS active BELONGS TO.  Ported from Maxwell, where it was
   * added after chunking examples/pde here exposed the hole in both.
   *
   * checkpoint() restores the independents and the dependents and frees the
   * whole graph; every other active is left holding idx and vtx that name
   * vertices which no longer exist, and the new pass renumbers from the same
   * base, so a stale idx can even collide with a live one.  Reading such an
   * active spliced a freed vertex into the new graph -- AddressSanitizer
   * caught it as a use-after-free in Vertex::kill(), and in the cases that did
   * not crash the derivative silently came out zero.
   *
   * Process stamps this whenever it gives an active a vertex; restore_values()
   * re-stamps the independents and the dependents at the top of every pass.
   * An out-of-date stamp means the active did not survive the checkpoint.
   * 0 is "never recorded", which is a freshly constructed active.
   */
  mutable largeint gen;
		
public:		

  active();
  active( double  a);
  active( const active & x );		
  
  ~active();

  active & operator=(double x);  
  active & operator=( const active & x);

  active & operator+=( double x );   		
  active & operator+=( const active & x );
 
  active & operator-=( double x );   		
  active & operator-=( const active & x );

  active & operator*=( double x );      		
  active & operator*=( const active & x );

  active & operator/=( double x );    		
  active & operator/=( const active & x );

  active & operator++();
  active & operator--();

  active operator++(int);
  active operator--(int); 
};

}

#endif


