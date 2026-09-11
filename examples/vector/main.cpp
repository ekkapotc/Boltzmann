#include <iostream>
#include <cstdlib>

#include "../../boltzmann.hpp"

using namespace boltzmann;

template<class T>
class vector3
{
public:
  vector3():
   _x(0),
   _y(0),
   _z(0)
  {

  }

  vector3( T x , T y , T z ):
  _x(x),
  _y(y),
  _z(z)
  { 
    
  }

  vector3( const vector3<T> & v ):
  _x(v._x),
  _y(v._y),
  _z(v._z)
  {

  }

  vector3<T> & operator=( const vector3<T> & v )
  {
    if(this!=&v){
      _x = v._x;
      _y = v._y;
      _z = v._z;
    }
    return *this;
  }

  friend vector3<T> operator+( const vector3<T> & v1 , const vector3<T> & v2 )//vector sum
  {
    vector3<T> v3;
    v3._x = v1._x + v2._x;
    v3._y = v1._y + v2._y;
    v3._z = v1._z + v2._z;
    return v3;
  }

  friend vector3<T> operator-( const vector3<T> & v1 , const vector3<T> & v2 )//vector sum
  {
    vector3<T> v3;
    v3._x = v1._x - v2._x;
    v3._y = v1._y - v2._y;
    v3._z = v1._z - v2._z;
    return v3;
  }

  friend vector3<T> operator*( const vector3<T> & v1 , double c )//scalar mult
  {
    vector3<T> v3;
    v3._x = c*v1._x;
    v3._y = c*v1._y;
    v3._z = c*v1._z;
    return v3;
  }

  friend vector3<T> operator*( double c , const vector3<T> & v1 )//scalar mult
  {
    vector3<T> v3;
    v3._x = c*v1._x;
    v3._y = c*v1._y;
    v3._z = c*v1._z;
    return v3;
  }

  friend T operator&( const vector3<T> & v1 , const vector3<T> & v2 )//inner product
  {
    return v1._x*v2._x+v1._y*v2._y+v1._z*v2._z;
  }

  friend vector3<T> sin( const vector3<T> & v1 )//sin
  {
    vector3<T> v3;
    v3._x = sin(v1._x);
    v3._y = sin(v1._y);
    v3._z = sin(v1._z);
    return v3;
  }

  friend vector3<T> cos( const vector3<T> & v1 )//cos
  {
    vector3<T> v3;
    v3._x = cos(v1._x);
    v3._y = cos(v1._y);
    v3._z = cos(v1._z);
    return v3;
  }

  friend T min( const vector3<T> & v1 )//min
  {
    return min(v1._x,min(v1._y,v1._z));
  }

  friend T max( const vector3<T> & v1 )//max
  {
    return max(v1._x,max(v1._y,v1._z));
  }

private:
  T _x;
  T _y;
  T _z;
};

void f(
       int nx, 
       int nt, 
       active * temp,
       active & cost
       ) 
{
  cost = 0.0;

  for( int i=0 ; i<nt ; i++ ){

    //double mem_before = get_memory();

    for( int j=0 ; j<nx ; j++ )
    {
      vector3<active> v( temp[3*j] , temp[3*j+1] , temp[3*j+2] );
      cost += (v&v);
    }

    //double mem_after = get_memory();
    //std::cout << "mem_after-mem_before = " << mem_after-mem_before << std::endl;
  }
}

int main( int argc , char ** argv )
{
  int nx = atoi(argv[1]);
  int nt = atoi(argv[2]);
  double it = atoi(argv[3]);
  
  //int k=0;

  //while(k<5){

  active * temp=new active[3*nx];
  active cost;
 	
  for( int i=0 ; i<3*nx ; i++ ){
    temp[i] = 1.0+i*0.5;
  }

  initialize(3*nx,1,it*33600);//nx=20 //0->1

  for( int i=0 ; i<3*nx ; i++ ){
    independent(temp[i]);//1->2
  }

  while(checkpoint(temp,cost)){
    try{
      f(nx,nt,temp,cost);//2->3
    }catch(BreakException const & e){

    }
  }

  dependent(cost);//3->4

  //std::cout << "cost.idx = " << cost.idx << std::endl;

  double ** A;
  harvest(1,3*nx,A);//A is an 1 x 3*nx matrix //4->5
	
  finalize();//5->0
	
  delete [] temp;

  //k++;

  //}

  return 0;
}

