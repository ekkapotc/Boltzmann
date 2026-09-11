#ifndef VERTEX_INCLUDE
#define VERTEX_INCLUDE

#include "Active.hpp"
#include "Typedefs.hpp"
#include "Edge.hpp"

#include <vector>

namespace boltzmann
{
namespace internals
{

/*
 * Ported from Maxwell (SVEGP-30).
 *
 * Adjacency used to be std::map<largeint,Edge*>, twice per vertex.  That cost
 * 96 of the vertex's bytes before a single edge existed, and every edge bought
 * two red-black tree nodes -- so an edge was one allocation for itself plus
 * two more for its tree nodes, each with malloc's own header, scattered across
 * the heap.
 *
 * It is a sorted array now, binary-searched: one contiguous block per vertex
 * per direction, carved from an arena.  The members are called first/second so
 * that the elimination code, which is full of it->second, reads exactly as it
 * did against std::map, and insert() keeps std::map's semantics deliberately:
 * an existing key is left alone, not overwritten.
 *
 * WHAT IS DIFFERENT HERE.  Maxwell's graph is built once and eliminated once.
 * This one is also SERIALISED: Process::send() walks every surviving vertex's
 * in-list and writes a cji_t per edge onto the wire.  A sorted array is a
 * contiguous scan for that walk where a std::map was a pointer-chasing tree
 * traversal, so the gain shows up twice -- once in elimination and once in the
 * pack loop that feeds MPI_Send.  It also makes the serialised edge ORDER a
 * deterministic function of the graph rather than of tree-node addresses,
 * which is what lets the receiving rank's accumulation be reproducible.
 */
struct AdjEntry
{
  largeint first;   //the neighbouring vertex's index
  Edge *   second;  //the edge itself, owned by the EdgeArena
};

/*
 * Size-classed storage for adjacency blocks.  Capacities are 1, 2, 4, 8, ...
 * one class per power of two -- most vertices on a tape have in-degree one or
 * two and out-degree one, so a larger minimum overshoots the commonest case on
 * the two blocks every single vertex carries.
 *
 * A block is carved from a slab of its class and, when released, goes on that
 * class's free list to be handed straight back out.  Slabs are never returned
 * to the allocator until the whole graph is torn down.
 */
class AdjArena
{

private:

  static const largeint NCLASS   = 22;//caps 1 .. 2 097 152
  static const largeint SLAB_MIN = 32;
  static const largeint SLAB_CAP = 4096;

  largeint slab_max;//per-class ceiling, set from the budget

  std::vector<AdjEntry*> slabs;               //every allocation, freed at clear()
  std::vector<AdjEntry*> recycled[NCLASS];    //blocks handed back, per class

  AdjEntry * cursor[NCLASS];                  //current slab for the class
  largeint   left[NCLASS];                    //entries still uncarved in it
  largeint   next_slab[NCLASS];               //entries the class's next slab holds

  largeint entries;                           //total entries across every slab

  AdjArena( const AdjArena & );//an arena owns raw storage: not copyable
  AdjArena & operator=( const AdjArena & );

  static largeint class_of( largeint want );

public:

  static largeint cap_of( largeint cls ){ return largeint(1)<<cls; }

  AdjArena();
  ~AdjArena();

  void set_budget( largeint bytes );

  //hands back a block of at least want entries; cap is what it actually holds
  AdjEntry * acquire( largeint want , largeint & cap );

  void release( AdjEntry * p , largeint cap );

  //what the arena actually holds, slabs and free lists included
  largeint bytes() const;

  void clear();

};//end of class

class Adjacency
{

private:

  AdjEntry * a;     //arena block, or NULL while empty
  largeint   n;     //entries in use
  largeint   cap;   //entries the block holds

public:

  typedef AdjEntry * iterator;

  Adjacency(): a(NULL), n(0), cap(0) {}

  iterator begin(){ return a; }

  iterator end(){ return a+n; }

  largeint size() const { return n; }

  //first entry whose key is >= k
  iterator lower( largeint k )
  {
    largeint lo = 0;
    largeint hi = n;

    while( lo < hi ){
      const largeint mid = lo + (hi-lo)/2;
      if( a[mid].first < k ) lo = mid+1;
      else                   hi = mid;
    }

    return a + lo;
  }

  iterator find( largeint k )
  {
    iterator it = lower(k);
    if( it!=end() && it->first==k ) return it;
    return end();
  }

  //std::map::insert: a key already present is left as it is
  void insert( largeint k , Edge * e , AdjArena & arena )
  {
    iterator it = lower(k);
    if( it!=end() && it->first==k ) return;

    const largeint at = largeint(it - a);

    if( n==cap ){
      largeint new_cap = 0;
      AdjEntry * grown = arena.acquire( cap ? cap*2 : 1 , new_cap );

      for( largeint i=0 ; i<n ; i++ ) grown[i] = a[i];

      if(a) arena.release(a,cap);

      a   = grown;
      cap = new_cap;
    }

    for( largeint i=n ; i>at ; i-- ) a[i] = a[i-1];

    a[at].first  = k;
    a[at].second = e;

    n++;
  }

  void erase( largeint k )
  {
    iterator it = lower(k);
    if( it==end() || it->first!=k ) return;

    const largeint at = largeint(it - a);

    for( largeint i=at ; i+1<n ; i++ ) a[i] = a[i+1];

    n--;
  }

  /*
   * Gives the block back.  An Adjacency cannot clean up after itself -- it
   * does not know its arena, deliberately, a back-pointer per vertex being
   * most of what the arena saves.  Process owns the arena and retires
   * vertices through it.
   */
  void clear( AdjArena & arena )
  {
    if(a) arena.release(a,cap);
    a   = NULL;
    n   = 0;
    cap = 0;
  }

  largeint bytes() const { return cap*sizeof(AdjEntry); }

};//end of class

class Vertex 
{
public :
  
  bool alive;
  bool dep_vertex;
  largeint idx;
  largeint owner_idx;

  Adjacency in_edges;
  Adjacency out_edges;

public:

  Vertex( );

  Vertex ( largeint idx , largeint owner_idx );

  /*
   * Deliberately empty.  The adjacency blocks belong to the arena and are
   * given back by Process::retire(), which is the only place a Vertex is
   * destroyed.
   */
  ~Vertex();
	
  bool is_dead();

  void kill();

  bool is_dep_vertex();

  largeint get_idx();	

  largeint get_owner_idx();

  largeint in_degree();

  largeint out_degree();

  Edge * to( Vertex * tgt );

  Edge * from( Vertex * src );

  Edge * add_in_edge( Vertex * src , double partial , EdgeArena & arena , AdjArena & adj );

  largeint eliminate( EdgeArena & arena , AdjArena & adj );

};//end of class

}//end of namespace internals	
}//end of namespace boltzmann

#endif
