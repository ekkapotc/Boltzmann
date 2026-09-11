#include "../inc/Active.hpp"
#include "../inc/API.hpp"

#include <cmath>
#include <iostream>

using namespace boltzmann;

//namespace boltzmann
//{

active::active():
reachable(false),
idx(0),
owner_idx(0),
old_idx(0),
val(0),
vtx(NULL),
gen(0)
{ 
  //cout << " active() called...(val = " << val << ")" << endl;
}

active::active( double  a):
reachable(false),
idx(0),
owner_idx(0),
old_idx(0),
val(a),
vtx(NULL),
gen(0)
{
  //cout << "active( double ) called...(val = " << val << ")" << endl;	
}

active::active( const active & x ):
reachable(false),
idx(0),
owner_idx(0),
old_idx(0),
val(x.val),
vtx(NULL),
gen(0)
{
  //if(!internals::skip_mode){
    reachable = x.reachable;
    unary_op( x , 1.0 , *this , false );//overwrite=false
    //cout << "active( active & ) called...(val = " << val << ")" << endl;
  //}
}

active::~active()
{
  //if(!internals::skip_mode){
    destructor( *this );
    //cout << "~active() called..." << endl;
  //}
}

active & active::operator=( const active & x )
{
  if(this!=&x)//prevent self-assignment
  {
    //if(!internals::skip_mode){
      val = x.val; 
      reachable = x.reachable;
      unary_op( x , 1.0 , *this , true );//overwrite = true
      //cout << "operator=( active & ) called..." << endl;
    //}
  }

  return *this;
}

active & active::operator=( double x )
{
  //if(!internals::skip_mode){
    val = x; 
    reachable=false;
    passive_op(*this);
    //cout << "operator=( double ) called..." << endl;
  //}

  return *this;
}

active & active::operator+=( double x )
{
  //if(!internals::skip_mode){
    val+=x;//update val
    //reachable remains as it is
    unary_op_ass( 1.0 , *this );//A = A + x
    //cout << "operator+=(double) called..." << endl;
  //}

  return *this;
}

active & active::operator+=( const active & x )
{
  //if(!internals::skip_mode){
    val+=x.val;//update val
    reachable = reachable || x.reachable;
    binary_op_ass( x , 1.0 , *this , 1.0 );//A = A + B
    //cout << "operator+=(active&) called..." << endl;
  //}

  return *this;
}

active & active::operator-=( double x )
{
  //if(!internals::skip_mode){
    val-=x;//update val
    //reachable remains as it is
    unary_op_ass( 1.0 , *this );//A = A - x
    //cout << "operator-=(double) called..." << endl;
  //}
 
  return *this;
}

active & active::operator-=( const active & x )
{
  //if(!internals::skip_mode){
    val-=x.val;
    reachable = reachable || x.reachable;
    binary_op_ass( x , -1.0 , *this , 1.0 );
    //cout << "operator-=(active&) called..." << endl;
  //}

  return *this;
}

active & active::operator*=( double x )
{
  //if(!internals::skip_mode){
    val*=x;//update val
    //reachable remains as it is
    unary_op_ass( x , *this );//A = A*x
    //cout << "operator*=(double) called..." << endl;
  //}

  return *this;
}

active & active::operator*=( const active & x )
{
  //if(!internals::skip_mode){
    double old_val1 = val;
    double old_val2 = x.val;
    val*=x.val;
    reachable = reachable || x.reachable;
    binary_op_ass( x , old_val1 , *this , old_val2 );
    //cout << "operator*=(active&) called..." << endl;
  //}

  return *this;
}

active & active::operator/=( double x )
{
  //if(!internals::skip_mode){
    val/=x;//update val
    //reachable remains as it is
    unary_op_ass( (1.0/x) , *this );//A = A/x
    //cout << "operator/=(double) called..." << endl;
  //}

  return *this;
}

active & active::operator/=( const active & x )
{
  //if(!internals::skip_mode){
    double old_val1 = val;
    double old_val2 = x.val;
    val/=x.val;
    reachable = reachable || x.reachable;
    binary_op_ass( x , -(old_val1/(old_val2*old_val2)) , *this , 1.0/old_val2 );
    //cout << "operator/=(active&) called..." << endl;
  //}

  return *this;
}

active & active::operator++()
{
  //if(!internals::skip_mode){
    val += 1.0;
    //reachbale remains as it is
    unary_op_ass( 1.0 , *this );
    //cout << "operator++() called..." << endl;
    return *this;
  //}
}

active & active::operator--()
{
  //if(!internals::skip_mode){
    val -= 1.0;
    //reachbale remains as it is
    unary_op_ass( 1.0 , *this );
    //cout << "operator--() called..." << endl;
  //} 
 
  return *this;
}

active active::operator++(int)
{
  active x;

  //if(!internals::skip_mode){
    x.val = val;
    x.reachable = reachable;
    postfix_op( *this , x );
    ++(*this);//prefix operator
  //}

  return x;//return saved state
}

active active::operator--(int)
{
  active x;

  //if(!internals::skip_mode){
    x.val = val;
    x.reachable = reachable;
    postfix_op( *this , x );
    --(*this);//prefix operator
  //}

  return x;//return saved state
}

bool boltzmann::operator!( const active & x )
{
  return x.val==0;
}

double boltzmann::passive_value( const active & x )
{
  return x.val;
}

const active boltzmann::operator+( const active & x1 )
{
  active x2;
  
  //if(!internals::skip_mode){
    x2.val = x1.val;
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0 , x2 , false );
    //cout << "operator+(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::operator-( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = -x1.val;
    x2.reachable = x1.reachable;
    unary_op( x1 , -1.0 , x2 , false );
    //cout << "operator-(active&) called..." << endl;
  //} 
 
  return x2;
}

const active boltzmann::operator+( const active & x1 , const active & x2 )
{
  active x3;
 
  //if(!internals::skip_mode){  
    x3.val = x1.val + x2.val;
    x3.reachable = x1.reachable || x2.reachable;
    binary_op( x1 , 1.0 , x2 , 1.0 , x3 );
    //cout << "operator+(active&,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator+( const active & x1 ,  double  x2 )
{
  active x3;

  //if(!internals::skip_mode){ 
    x3.val = x1.val+x2;
    x3.reachable=x1.reachable;
    unary_op( x1 , 1.0 , x3 , false );
    //cout << "operator+(active&,double) called..." << endl;	
  //}

  return x3;
}

const active boltzmann::operator+( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = x1 + x2.val;
    x3.reachable = x2.reachable;
    unary_op( x2 , 1.0 , x3 , false );
    //cout << "operator+(double,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator-( const active & x1 ,  const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){ 
    x3.val = x1.val-x2.val;
    x3.reachable = x1.reachable || x2.reachable;
    binary_op( x1 , 1.0 , x2 , -1.0 , x3 );
    //cout << "operator-(active&,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator-( const active & x1 , double x2 )
{
  active x3;
 
  //if(!internals::skip_mode){
    x3.val = x1.val -x2;
    x3.reachable=x1.reachable;
    unary_op( x1 , 1.0 , x3 , false );
    //cout << "operator-(active&,double) called..." << endl;	
  //}

  return x3;
}

const active boltzmann::operator-( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = x1-x2.val;
    x3.reachable = x2.reachable;
    unary_op( x2 , -1.0 , x3 , false );
    //cout << "operator-(double,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator*( const active & x1 , const active & x2 )
{
  active x3; 
 
  //if(!internals::skip_mode){
    x3.val = x1.val*x2.val;
    x3.reachable = x1.reachable || x2.reachable;
    binary_op( x1 , x2.val , x2 , x1.val , x3 );	
    //cout << "operator*(active&,active&) called..." << endl;
  //}
  
  return x3;
}

const active boltzmann::operator*( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){ 
    x3.val = x1.val*x2;
    x3.reachable = x1.reachable;
    unary_op( x1 , x2 , x3 , false );
    //cout << "operator*(active&,double) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator*( double x1 , const active & x2 )
{
  active x3;
  
  //if(!internals::skip_mode){
    x3.val = x1*x2.val;
    x3.reachable = x2.reachable;
    unary_op( x2 , x1 , x3 , false );
    //cout << "operator*(double,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::operator/( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = x1.val / x2.val;
    x3.reachable = x1.reachable || x2.reachable;
    binary_op( x1 , 1.0/x2.val , x2 , -x1.val/std::pow(x2.val,2) , x3 );
    //cout << "operator/(active&,active&) called..." << endl;
  //}
  
  return x3;
}

const active boltzmann::operator/( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = x1.val/x2;
    x3.reachable = x1.reachable;
    unary_op( x1 , 1.0/x2 , x3 , false );
    //cout << "operator/(active&,double) called..." << endl;
  //} 

  return x3;
}

const active boltzmann::operator/( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = x1/x2.val;
    x3.reachable = x2.reachable;
    unary_op( x2 , -x1/std::pow(x2.val,2) , x3 , false );
    //cout << "operator/(double,active&) called..." << endl;	
  //}

  return x3;
}

const active boltzmann::sin( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::sin(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , std::cos(x1.val) , x2 , false );
    //cout << "sin(active&) called..." << endl;	
  //}

  return x2; 
}

const active boltzmann::cos( const active & x1 )
{
  active x2;
  
  //if(!internals::skip_mode){
    x2.val = std::cos(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , -std::sin(x1.val) , x2 , false );
    //cout << "cos(active&) called..." << endl;	
  //}

  return x2;
}

const active boltzmann::tan( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::tan(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0+(x2.val*x2.val) , x2 , false );
    //cout << "tan(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::asin( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::asin(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(std::sqrt(1.0-(x1.val*x1.val))) , x2 , false );
    //cout << "asin(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::acos( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::acos(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , -1.0/(std::sqrt(1.0-(x1.val*x1.val))) , x2 , false );
    //cout << "acos(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::atan( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::atan(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(1.0+(x1.val*x1.val)) , x2 , false );
    //cout << "atan(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::sinh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::sinh(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , std::cosh(x1.val) , x2 , false );
    //cout << "sinh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::cosh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::cosh(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , std::sinh(x1.val) , x2 , false );
    //cout << "cosh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::tanh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::tanh(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0-(x2.val*x2.val) , x2 , false );
    //cout << "tanh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::asinh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = ::asinh(x1.val);//c++11 std::asinh
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(std::sqrt(1.0+(x1.val*x1.val))), x2 , false );
    //cout << "asinh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::acosh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = ::acosh(x1.val);//c++11 std::asinh
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(std::sqrt((x1.val*x1.val)-1.0)), x2 , false );
    //cout << "acosh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::atanh( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = ::atanh(x1.val);//c++11 std::asinh
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(1.0-(x1.val*x1.val)) , x2 , false );
    //cout << "atanh(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::atan2( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::atan2(x1.val,x2.val);
    x3.reachable = x1.reachable || x2.reachable;
    binary_op( x1 , x2.val/((x1.val*x1.val)+(x2.val*x2.val)) , x2 , -x1.val/((x1.val*x1.val)+(x2.val*x2.val)) , x3 );
    //cout << "atan2(active&,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::atan2( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::atan2(x1.val,x2);
    x3.reachable = x1.reachable;
    unary_op( x1 , x2/((x1.val*x1.val)+(x2*x2)) , x3 , false );
    //cout << "atan2(active&,double) called..." << endl;
  //}

  return x3;
}

const active boltzmann::atan2( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::atan2(x1,x2.val);
    x3.reachable = x2.reachable;
    unary_op( x2 , (-x1)/((x1*x1)+(x2.val*x2.val)) , x3 , false );
    //cout << "atan2(double,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::hypot( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = ::hypot(x1.val,x2.val);//c++11 std::hypot
    x3.reachable = x1.reachable || x2.reachable;
    double hypot_val = x3.val;
    if(x1.val==0.0 && x2.val==0.0){
      hypot_val += 1.0e-12;//stablize
    }
    binary_op( x1 , x1.val/hypot_val , x2 , x2.val/hypot_val , x3 );
    //cout << "hypot(active&,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::hypot( const active & x1, double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = ::hypot(x1.val,x2);//c++11 std::hypot
    x3.reachable = x1.reachable;
    double hypot_val = x3.val;
    if(x1.val==0.0 && x2==0.0){
      hypot_val += 1.0e-12;//stablize
    }
    unary_op( x1 , x1.val/hypot_val , x3 , false );
    //cout << "hypot(active&,double) called..." << endl;
  //}

  return x3;
}

const active boltzmann::hypot( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = ::hypot(x1,x2.val);//c++11 std::hypot
    x3.reachable = x2.reachable;

    double hypot_val = x3.val;

    if(x1==0.0 && x2.val==0.0){
      hypot_val += 1.0e-12;//stablize
    }

    unary_op( x2 , x2.val/hypot_val , x3 , false );
    //cout << "hypot(double,active&) called..." << endl;
  //}
    
  return x3;
}

const active boltzmann::exp( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::exp(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , x2.val , x2 , false );
    //cout << "exp(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::pow( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::pow(x1.val,x2.val);
    x3.reachable = x1.reachable || x2.reachable;

    double dy_dx1 = (x2.val)*std::pow(x1.val,(x2.val)-1.0);//a*x^(x-1)
    double dy_dx2 = 0.0;

    if(x1.val>0.0){
      dy_dx2 = std::log(x1.val)*std::pow(x1.val,x2.val);//(a^x)*ln(a)
    }

    binary_op( x1 , dy_dx1 , x2 , dy_dx2 , x3 );
    //cout << "pow(active&,active&) called..." << endl;
 //}

  return x3;
}

const active boltzmann::pow( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::pow(x1.val,x2);
    x3.reachable = x1.reachable;
    unary_op( x1 , x2*std::pow(x1.val,x2-1.0) , x3 , false );
    //cout << "pow(active&,double) called..." << endl;
  //}

  return x3;
}

const active boltzmann::pow( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){  
    x3.val = std::pow(x1,x2.val);
    x3.reachable = x2.reachable;
    
    double dy_dx2 = 0.0;

    if(x1>0.0){
      dy_dx2 = std::log(x1)*std::pow(x1,x2.val);//(a^x)*ln(a)
    }

    unary_op( x2 , dy_dx2 , x3 , false );
    //cout << "pow(double,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::sqrt( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::sqrt(x1.val);
    x2.reachable = x1.reachable;
  
    if(x1.val!=0){
      unary_op( x1 , 1.0/(2.0*std::sqrt(x1.val)) , x2 , false );
    }else{
      unary_op( x1 , 1.0/(2.0*(std::sqrt(x1.val)+1.0e-12)) , x2 , false );//stablize
    }
    //cout << "sqrt(active&) called..." << endl;	
  //}  

  return x2;
}

const active boltzmann::cbrt( const active & x1 )//c++11 std::cbrt
{
  active x2;
  
  //if(!internals::skip_mode){
    x2.val = ::cbrt(x1.val);
    x2.reachable = x1.reachable;

    double denom_val = 3.0*(x2.val*x2.val);

    if(x2.val==0.0){
      denom_val += 1.0e-12;//stablize
    }

    unary_op( x1 , 1.0/denom_val , x2 , false );
    //cout << "cbrt(active&) called..." << endl;	
  //}
 
  return x2;
}

const active boltzmann::log( const active & x1 )
{
  active x2;
  
  //if(!internals::skip_mode){
    x2.val = std::log(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/x1.val , x2 , false );
    //cout << "log(active&) called..." << endl;	
  //}

  return x2;
}

const active boltzmann::log10( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::log10(x1.val);
    x2.reachable = x1.reachable;
    unary_op( x1 , 1.0/(x1.val*std::log(10.0)) , x2 , false );
    //cout << "log10(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::min( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::min(x1.val,x2.val);

    if(x1.val==x2.val){
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }else if(x1.val>x2.val){  
      x3.reachable = x2.reachable;
      unary_op( x2 , 1.0 , x3 , false );
    }else{
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }
    //cout << "min(active&,active&) called..." << endl;
  //}

  return x3;
}

const active boltzmann::min( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::min(x1.val,x2);

    if(x1.val<=x2){
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }else{
      x3.reachable = false;
      passive_op( x3 );
    }
    //cout << "min(active&,double) called..." << endl;
  //}  

  return x3;
}

const active boltzmann::min( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::min(x1,x2.val);

    if(x2.val<=x1){
      x3.reachable = x2.reachable;
      unary_op( x2 , 1.0 , x3 , false );
    }else{
      x3.reachable = false;
      passive_op( x3 );
    }
    //cout << "min(double,active&) called..." << endl;
  //}
  
  return x3;
}

const active boltzmann::max( const active & x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::max(x1.val,x2.val);

    if(x1.val==x2.val){
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }else if(x1.val>x2.val){  
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }else{
      x3.reachable = x2.reachable;
      unary_op( x2 , 1.0 , x3 , false );
    }
    //cout << "max(active&,active&) called..." << endl;
  //}  

  return x3;
}

const active boltzmann::max( const active & x1 , double x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::max(x1.val,x2);

    if(x1.val>=x2){
      x3.reachable = x1.reachable;
      unary_op( x1 , 1.0 , x3 , false );
    }else{
      x3.reachable = false;
      passive_op( x3 );
    }
  
    //cout << "max(active&,double) called..." << endl;
  //}

  return x3;
}

const active boltzmann::max( double x1 , const active & x2 )
{
  active x3;

  //if(!internals::skip_mode){
    x3.val = std::max(x1,x2.val);

    if(x2.val>=x1){
      x3.reachable = x2.reachable;
      unary_op( x2 , 1.0 , x3 , false );
    }else{
      x3.reachable = false;
      passive_op( x3 );
    }
    //cout << "max(double,active&) called..." << endl;
  //}
 
  return x3;
}

const active boltzmann::fabs( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = std::fabs(x1.val);
    x2.reachable = x1.reachable;

    if(x1.val<0){
      unary_op( x1 , -1.0 , x2 , false );
    }else{
      unary_op( x1 , 1.0 , x2 , false );
    }
    //cout << "fabs(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::erf( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = ::erf(x1.val);//c++11 std::erf
    x2.reachable = x1.reachable;
    unary_op( x1 , (2.0/std::sqrt(3.14159265358979323846264338327950288))*std::exp(-x1.val*x1.val) , x2 , false );
    //cout << "erf(active&) called..." << endl;
  //}

  return x2;
}

const active boltzmann::erfc( const active & x1 )
{
  active x2;

  //if(!internals::skip_mode){
    x2.val = ::erfc(x1.val);//c++11 std::erfc
    x2.reachable = x1.reachable;
    unary_op( x1 , -(2.0/std::sqrt(3.14159265358979323846264338327950288))*std::exp(-x1.val*x1.val) , x2 , false );
  //}
   
  return x2;
}

std::ostream & boltzmann::operator<<(std::ostream & os, const active & x )
{
  os << "value = " << x.val;
  os << " , ";
  os << "idx = " << x.idx;
  //cout << "operator<< called..." << endl;
  return os;
}

std::istream & boltzmann::operator>>(std::istream & in, const active & x )
{
  in >> x.val;
  x.reachable = false;
  passive_op(x);
  //cout << "operator>> called..." << endl;
  return in;
}

double boltzmann::abs( const active & x )
{
  return std::abs(x.val);
}

bool boltzmann::operator==(const active & x1,const active & x2) { return x1.val==x2.val; }
bool boltzmann::operator==(const active & x1,double  x2)  { return x1.val==x2; }
bool boltzmann::operator==(double x1,const active & x2)   { return x1==x2.val; }

bool boltzmann::operator!=(const active & x1,const active & x2) { return x1.val!=x2.val; }
bool boltzmann::operator!=(const active & x1,double x2)   { return x1.val!=x2; }
bool boltzmann::operator!=(double x1,const active & x2)   { return x1!=x2.val; }

bool boltzmann::operator<(const active & x1, const active & x2) { return x1.val < x2.val; }
bool boltzmann::operator<(const active & x1, double x2)   { return x1.val < x2; }
bool boltzmann::operator<(double x1, const active & x2)   { return x1<x2.val; }

bool boltzmann::operator>(const active & x1,const active & x2)  { return x1.val > x2.val; }
bool boltzmann::operator>(const active & x1,double x2)    { return x1.val > x2; }
bool boltzmann::operator>(double x1,const active & x2)    { return x1 > x2.val; }
		 
bool boltzmann::operator<=(const active & x1,const active & x2) { return x1.val <= x2.val; }
bool boltzmann::operator<=(const active & x1,double x2)   { return x1.val <= x2; }
bool boltzmann::operator<=(double x1,const active & x2)   { return x1 <= x2.val; }

bool boltzmann::operator>=(const active & x1,const active & x2) { return x1.val >= x2.val; }
bool boltzmann::operator>=(const active & x1,double x2)   { return x1.val >= x2; }
bool boltzmann::operator>=(double x1,const active & x2)   { return x1 >= x2.val; }

//}//end of namespace boltzmann

