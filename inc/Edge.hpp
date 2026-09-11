#ifndef EDGE_INCLUDE
#define EDGE_INCLUDE

#include "Typedefs.hpp"

#include <vector>

namespace boltzmann
{
namespace internals
{

class Vertex;//only ever pointed to from here

class Edge
{

public :

  Vertex * src;
  Vertex * tgt;
  double   eval;

public:

  Edge();

  Edge( Vertex * src , Vertex * tgt , double eval );

  Vertex * get_src();

  Vertex * get_tgt();

  double get_partial();

};//end of class

/*
 * Ported from Maxwell (SVEGP-29).
 *
 * Edges used to be one heap allocation each -- new in add_in_edge(), delete in
 * eliminate() -- on top of the two red-black tree nodes the adjacency maps
 * allocate for every edge.  Elimination churns edges hard, so that was several
 * million malloc/free round trips for a graph that never held more than a
 * fraction of them at once.
 *
 * They come from blocks now, and freed slots go on a free list to be handed
 * straight back out.  The behaviour is deliberately identical to what malloc
 * was doing, recycling included.
 *
 * WHY THIS MATTERS MORE HERE THAN IN MAXWELL.  Every rank in the pipeline
 * builds a graph, eliminates it, serialises what survives and then TEARS THE
 * WHOLE THING DOWN, once per partition, N/p times per rank.  Maxwell does that
 * churn once; a p-rank job does it p times concurrently against one kernel
 * allocator, and malloc's arena lock is per-process, not per-rank -- but the
 * ranks are separate processes, so what they contend for is the page cache and
 * the kernel's mmap path.  Removing the allocations removes the contention
 * along with the time.
 *
 * The block size follows the memory budget (Maxwell SVEGP-31): a fixed block
 * put a floor under the whole checkpointing scheme, since the arena claimed
 * the same capacity whatever the budget said, and tightening the budget then
 * bought nothing.
 */
class EdgeArena
{

private:

  std::vector<Edge*> blocks;
  std::vector<Edge*> free_list;

  largeint block_used;//slots taken from the newest block
  largeint slots;     //total slots across every block
  largeint live;      //slots currently handed out

  static const largeint BLOCK_MAX = 8192;
  static const largeint BLOCK_MIN = 64;

  largeint block;//slots per block; set from the budget before anything is carved

  EdgeArena( const EdgeArena & );//an arena owns raw storage: not copyable
  EdgeArena & operator=( const EdgeArena & );

public:

  EdgeArena();
  ~EdgeArena();

  /*
   * Sized from the per-partition budget by Process::initialize(), before a
   * single edge exists.  Ignored once storage has been carved.
   */
  void set_budget( largeint bytes );

  largeint block_slots() const;

  Edge * acquire( Vertex * src , Vertex * tgt , double eval );

  void release( Edge * e );

  largeint live_edges() const;

  //what the arena actually holds, blocks and free list included
  largeint bytes() const;

  void clear();

};//end of class

}//end of namespace internals
}//end of namespace boltzmann

#endif
