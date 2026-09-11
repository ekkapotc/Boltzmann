#include "../inc/Edge.hpp"
#include "../inc/Vertex.hpp"

using namespace boltzmann;
using namespace boltzmann::internals;

Edge::Edge():
src(NULL),
tgt(NULL),
eval(0.0)
{

}

Edge::Edge( Vertex * s , Vertex * t , double e ):
src( s ),
tgt( t ),
eval( e )
{

}

Vertex * Edge::get_src()
{
  return src;
}

Vertex * Edge::get_tgt()
{
  return tgt;
}

double Edge::get_partial()
{
  return eval;
}

/* ------------------------------------------------------------ EdgeArena -- */

EdgeArena::EdgeArena():
block_used(0),
slots(0),
live(0),
block(BLOCK_MAX)
{

}

EdgeArena::~EdgeArena()
{
  clear();
}

/*
 * A chunk of `bytes` holds at most bytes/EDGE_BYTES edges, so that is what one
 * block needs to cover.  BLOCK_MAX keeps the old fixed value as a ceiling, so
 * a large or absent budget behaves exactly as before; BLOCK_MIN stops a
 * pathologically small budget turning the arena back into one malloc per edge,
 * which is what the arena exists to remove.
 */
void EdgeArena::set_budget( largeint bytes )
{
  if(!blocks.empty()) return;

  //sizeof(Edge) plus the two adjacency entries that reference it; see
  //Process.cpp's EDGE_BYTES, which this deliberately mirrors
  const largeint per_edge = largeint(sizeof(Edge)) + 2*largeint(2*sizeof(void*));

  largeint want = (bytes>0) ? (bytes/per_edge) : BLOCK_MAX;

  if( want < BLOCK_MIN ) want = BLOCK_MIN;
  if( want > BLOCK_MAX ) want = BLOCK_MAX;

  block = want;
}

largeint EdgeArena::block_slots() const
{
  return block;
}

Edge * EdgeArena::acquire( Vertex * src , Vertex * tgt , double eval )
{
  Edge * e = NULL;

  if(!free_list.empty()){
    e = free_list.back();
    free_list.pop_back();
  }else{

    if( blocks.empty() || block_used==block ){
      Edge * b = new Edge[size_t(block)];
      blocks.push_back(b);
      slots     += block;
      block_used = 0;
    }

    e = blocks.back() + block_used;
    block_used++;
  }

  e->src  = src;
  e->tgt  = tgt;
  e->eval = eval;

  live++;
  return e;
}

void EdgeArena::release( Edge * e )
{
  if(!e) return;

  free_list.push_back(e);

  if(live) live--;
}

largeint EdgeArena::live_edges() const
{
  return live;
}

largeint EdgeArena::bytes() const
{
  return slots*sizeof(Edge)
       + blocks.capacity()*sizeof(Edge*)
       + free_list.capacity()*sizeof(Edge*);
}

void EdgeArena::clear()
{
  for( std::vector<Edge*>::iterator it=blocks.begin() ; it!=blocks.end() ; it++ ){
    delete [] (*it);
  }

  blocks.clear();
  free_list.clear();

  block_used = 0;
  slots      = 0;
  live       = 0;
}
