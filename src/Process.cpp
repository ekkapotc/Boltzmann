#include "../inc/Process.hpp"
#include "../inc/API.hpp"
#include "../inc/BreakException.hpp"

#include <cassert>
#include <cstdio>
#include <set>//destroy_graph(), so no vertex can be deleted twice

#define DEFAULT_PROBE_FREQUENCY 50

/*
 * Ported from Maxwell (SVEGP-24) : what a recorded vertex and a recorded edge
 * actually cost the graph, in bytes.  These are the numbers check_memory()
 * spends the user's budget against, so they decide how many partitions the
 * tape is broken into -- and, in this library, how much work each rank of the
 * pipeline is given.
 *
 * The edge figure is the one that was wrong, and had been since the
 * beginning.  It counted sizeof(Edge) alone -- but an edge is not just its
 * Edge object: it is also the two ADJACENCY ENTRIES that reference it, one in
 * the target's in-list and one in the source's out-list.  Before the arena
 * those two were a pair of red-black tree nodes and went uncounted; now they
 * are a pair of AdjEntry and are counted.  An edge costs
 * sizeof(Edge) + 2*sizeof(AdjEntry) of payload, not sizeof(Edge).
 *
 * WHY AN HONEST BUDGET MATTERS MORE HERE.  In Maxwell an under-counted budget
 * means more partitions than the user asked for, which costs replay time.
 * Here the partition count also decides how the work is DISTRIBUTED: rank r
 * records partitions r, r+p, r+2p, ... so a budget that is wrong by 2.4x
 * changes the number of partitions and therefore changes which rank owns
 * which piece of the tape, the length of the pipeline, and the size of every
 * message on the ring.  Load balance across ranks is downstream of this
 * number being right.
 *
 * What these deliberately do NOT include is allocator overhead, std::vector
 * capacity slack and the arenas' block granularity -- properties of the
 * allocator rather than of the graph, and not knowable during the profiling
 * pass, when the budget decision is made and nothing has been allocated yet.
 * get_heap() reports what is actually held so the ratio cannot drift
 * unnoticed.
 */
static const boltzmann::largeint VERTEX_BYTES = sizeof(boltzmann::internals::Vertex);
static const boltzmann::largeint EDGE_BYTES   = sizeof(boltzmann::internals::Edge) + 2*sizeof(boltzmann::internals::AdjEntry);

using namespace boltzmann;
using namespace boltzmann::internals;

Process::Process():
comm(MPI_COMM_NULL),
break_mode(BREAK_ON_TARGET),//Maxwell SVEGP-32: what the library has always done
pass_mode(PASSES_PER_RANK),//what the library has always done
profiling(true),
throwable(false),
dtor_ignore_vertex(false),
has_vector(false),
mpi_rank(0),
mpi_size(1),
probe_freq(DEFAULT_PROBE_FREQUENCY),
MPI_EDGE(MPI_DATATYPE_NULL),
mem_size(0),
prev_mem_usage(0),
next_vertex_idx(1),
next_owner_idx(0),
top_owner_idx(0),
topmost_owner_idx(0),
indep_count(0),
intmed_count(0),
dep_count(0),
edge_count(0),
elim_cost(0),
partition_count(0),
pass_gen(1),
stale_reads(0),
productive_pass(0),
dependent_index(0)
{
 
}

/*
 * Ported from Maxwell (SVEGP-27) : the teardown, in one place, callable twice.
 *
 * intmed_vec and intmed_map never hold the same vertex at a moment a caller
 * can observe, but a teardown is the wrong place to lean on that -- a Process
 * abandoned mid-section (the checkpoint loop left by an exception, a rank that
 * threw out of a collective) is precisely the case where the two containers
 * are least likely to be in the state the happy path leaves them in.  The
 * pointers go through a set and each is retired exactly once.
 *
 * dep_vec is not retired from: register_dep_vertex() only ever pushes
 * vertices it has also put in intmed_map.  Clearing it drops the now-dangling
 * pointers instead of leaving them to be followed.
 *
 * Edges are not walked at all -- the arena owns every one of them.
 */
void Process::destroy_graph()
{
  std::set<Vertex*> owned;

  for( std::map<largeint,Vertex*>::iterator it=intmed_map.begin() ; it!=intmed_map.end() ; it++ ){
    if(it->second) owned.insert(it->second);
  }

  for( std::vector<Vertex*>::iterator it=intmed_vec.begin() ; it!=intmed_vec.end() ; it++ ){
    if(*it) owned.insert(*it);
  }

  for( std::vector<Vertex*>::iterator it=alive_intmed_vec.begin() ; it!=alive_intmed_vec.end() ; it++ ){
    if(*it) owned.insert(*it);
  }

  for( std::set<Vertex*>::iterator it=owned.begin() ; it!=owned.end() ; it++ ){
    retire(*it);
  }

  intmed_map.clear();
  intmed_vec.clear();
  alive_intmed_vec.clear();
  dep_vec.clear();
  cji_vec.clear();

  edge_arena.clear();
  adj_arena.clear();
}

Process::~Process()
{
  destroy_graph();

  if(MPI_EDGE!=MPI_DATATYPE_NULL){
    MPI_Type_free(&MPI_EDGE);//released before MPI_Finalize()
    MPI_EDGE = MPI_DATATYPE_NULL;
  }
}


//private member functions

/*
 * THE CHECKPOINT CONTRACT, ENFORCED.  Ported from Maxwell; found here, by
 * giving examples/pde a budget that actually chunks.
 *
 * checkpoint() restores the independents and the dependents and frees the
 * whole graph.  Every other active keeps an idx and a vtx naming vertices that
 * were destroyed, and the next pass renumbers from the same base, so a stale
 * idx can collide with a live one.  Before this, reading such an active
 * spliced a freed vertex into the new graph: AddressSanitizer caught it as a
 * use-after-free in Vertex::kill() from set_vertex_dead(), and where the
 * allocator happened not to reuse the block the derivative simply came out
 * zero with no complaint at all.
 */
bool Process::stale( const active & x ) const
{
  return x.gen != pass_gen;
}

void Process::adopt( const active & x ) const
{
  x.gen     = pass_gen;
  x.idx     = 0;
  x.old_idx = 0;
  x.vtx     = NULL;
}

largeint Process::advance_pass()
{
  return ++pass_gen;
}

largeint Process::generation() const
{
  return pass_gen;
}

largeint Process::get_stale_reads() const
{
  return stale_reads;
}

/*
 * Once, not once per operator, and prefixed with the rank: a section that does
 * this does it thousands of times on every rank at once, and the first message
 * is the one that tells you where to look.  The count is readable with
 * get_stale_reads() so a test can assert it.
 */
void Process::report_stale_read()
{
  if(!stale_reads){
    std::cerr <<
      "boltzmann: an active that did not survive checkpoint() has been read.\n"
      "           Only the independents and the dependents handed to\n"
      "           checkpoint() are restored between passes; everything else is\n"
      "           left pointing at a graph that has been freed.  It is being\n"
      "           treated as a constant, so any derivative that flows through\n"
      "           it will be WRONG.  Build the section's own variables inside\n"
      "           the section.\n"
      "           (reported once per rank; see get_stale_reads() for the count)\n";
  }
  stale_reads++;
}

Vertex * Process::vertex_on_rhs( const active & x )
{
  /*
   * A READ.  This is the one that was silently wrong, so it is the one that
   * reports.  NULL means "no vertex": add_edge() already drops a null operand,
   * so the stale value is used as a constant and contributes no derivative --
   * which is what it was doing anyway, now defined and audible.
   */
  if( stale(x) ){
    report_stale_read();
    adopt(x);
    return NULL;
  }

  if(is_proc())
  {
    if(x.owner_idx==next_owner_idx)
    {
      return x.vtx;
    }else
    {
      if(x.vtx)
      {
        return x.vtx;
      }else
      {
        Vertex * vPtr = new Vertex( x.idx , x.owner_idx );
        
	vertex_counter++;//debug

        if(x.idx==x.old_idx)
	{ 
          vPtr->dep_vertex = true;
        }

        x.vtx = vPtr;
        intmed_vec.push_back(vPtr);
        return vPtr;
      }
    }
  }
  
  return NULL;
}

Vertex * Process::vertex_on_lhs( const active & x )
{
  x.gen = pass_gen;//a write: this active now belongs to this pass

  intmed_count++;
 
  if(is_proc())
  {
    Vertex * vPtr = new Vertex( x.idx , x.owner_idx ); 
		
    vertex_counter++;//debug

    if(x.idx==x.old_idx)
    { 
      vPtr->dep_vertex = true;
    }

    x.vtx = vPtr;
    intmed_vec.push_back(vPtr);
    return vPtr;
  }

  return NULL;
}

void Process::set_vertex_dead( const active & x )
{
  /*
   * A WRITE.  Overwriting an active that did not survive the checkpoint is not
   * an error -- it is how a scratch variable declared outside the section gets
   * reused -- but the vertex it used to name is gone.  Drop it.  This is the
   * line that dereferenced freed memory in pde.
   */
  if( stale(x) ){
    adopt(x);
    return;
  }

  if(is_proc())
  {
    if(x.idx)
    { 
      if(x.idx<=next_vertex_idx)//bug fix on 08.07.2016
      {
        if(x.owner_idx==next_owner_idx)
        {
          Vertex * vPtr = x.vtx;
          vPtr->kill();
        }
       }
    }//end of x.idx

  }//end of is_proc()

  x.idx = 0;//reset idx to 0
  x.vtx = NULL;//reset vtx to NULL
}

Vertex * Process::set_vertex_dead_for_unary_op_ass( const active & x )
{
  /*
   * x op= c is a READ as well as a write -- the new value depends on the old
   * one -- so a stale x here loses a derivative path and is reported, unlike
   * the plain overwrite in set_vertex_dead().
   */
  if( stale(x) ){
    report_stale_read();
    adopt(x);
    return NULL;
  }

  if(is_proc())
  {
    if(x.idx)
    { 
      if(x.idx<=next_vertex_idx)//bug fix on 08.07.2016
      {
        if(x.owner_idx==next_owner_idx)
        {
          Vertex * vPtr = x.vtx;
          vPtr->kill();
          x.idx = 0;
          x.vtx = NULL;
          return vPtr;
        }else
        {
          if(x.vtx)
          {
            Vertex * vPtr = x.vtx;
            x.idx = 0;
            x.vtx = NULL;
            return vPtr;
          }else
          {
            Vertex * vPtr = new Vertex( x.idx , x.owner_idx );
  
 	    vertex_counter++;//debug

            if(x.idx==x.old_idx)
            { 
              vPtr->dep_vertex = true; 
            }

            intmed_vec.push_back(vPtr);
            x.idx = 0;
            x.vtx = NULL;
            return vPtr;
          }
        }
      }else{//bug fix on 08.07.2016
	x.idx = 0;
	x.vtx = NULL;
	return NULL;
      }
    }
  }

  x.idx = 0;//reset idx to 0
  x.vtx = NULL;//reset vtx to NULL
  return NULL;
}

Vertex * Process::set_vertex_dead_for_binary_op_ass( const active & x , bool & meaningful )
{
  if( stale(x) ){//x op= y reads x; see set_vertex_dead_for_unary_op_ass()
    report_stale_read();
    adopt(x);
    meaningful = false;
    return NULL;
  }

  if(is_proc())
  {
    if(x.idx)
    { 
      if(x.idx<=next_vertex_idx)
      {
        if(x.owner_idx==next_owner_idx)
        {
          Vertex * vPtr = x.vtx;
          vPtr->kill();
          x.vtx = NULL;
          return vPtr;
        }else
        {
          if(x.vtx)
	  {
            Vertex * vPtr = x.vtx;
            x.vtx = NULL;
            return vPtr;
          }else
  	  {
            Vertex * vPtr = new Vertex( x.idx , x.owner_idx );

	    vertex_counter++;//debug
      
            if(x.idx==x.old_idx){
              vPtr->dep_vertex = true;
            }

            intmed_vec.push_back(vPtr);
            x.vtx = NULL;
            return vPtr;
          }
        }
      }else{
        meaningful = false;
        x.vtx = NULL;
	return NULL;
      }
    }
  }

  if(x.idx>next_vertex_idx) 
    meaningful = false;

  x.vtx = NULL;//reset vtx to NULL
  return NULL;
}

bool Process::is_proc()
{
  if(next_owner_idx==mpi_rank+(top_owner_idx-mpi_size+1))
  {
    return true;
  }

  return false;
}

bool Process::is_max_rank()
{
 return mpi_rank+1==mpi_size;
}

bool Process::is_final_rank()
{
   return (largeint)mpi_rank==(mpi_size-1-(topmost_owner_idx%mpi_size)+1)%mpi_size;
}

void Process::move_edges( Vertex * from_vertex , Vertex * to_vertex )
{
  for( Adjacency::iterator it=from_vertex->out_edges.begin() ; it!=from_vertex->out_edges.end() ; it++ ){

    Edge * et = it->second;
    Adjacency::iterator dup = to_vertex->out_edges.find(it->first);

    if(dup!=to_vertex->out_edges.end()){
      //an edge to the same target already exists: partials accumulate,
      //std::map::insert() would have thrown this contribution away
      dup->second->eval += et->eval;
      et->tgt->in_edges.erase(to_vertex->idx);
      et->tgt->in_edges.insert(to_vertex->idx,dup->second,adj_arena);
      edge_arena.release(et);
      edge_counter--;//debug
    }else{
      et->src = to_vertex;//update src pointer
      to_vertex->out_edges.insert( it->first , et , adj_arena );
    }
  }

  for( Adjacency::iterator it=from_vertex->in_edges.begin() ; it!=from_vertex->in_edges.end() ; it++ ){

    Edge * et = it->second;
    Adjacency::iterator dup = to_vertex->in_edges.find(it->first);

    if(dup!=to_vertex->in_edges.end()){
      dup->second->eval += et->eval;
      et->src->out_edges.erase(to_vertex->idx);
      et->src->out_edges.insert(to_vertex->idx,dup->second,adj_arena);
      edge_arena.release(et);
      edge_counter--;//debug
    }else{
      et->tgt = to_vertex;//update tgt pointer
      to_vertex->in_edges.insert( it->first , et , adj_arena );
    }
  }

  from_vertex->out_edges.clear(adj_arena);
  from_vertex->in_edges.clear(adj_arena);
}

void Process::local_vertex_elimination_vector( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to )
{
  for( std::vector<Vertex*>::reverse_iterator it=from.rbegin() ; it!=from.rend() ; it++ ){

    Vertex * vPtr = (*it);
    
    if(vPtr->in_degree() && vPtr->out_degree()){
      if(!(vPtr->dep_vertex)){
        //std::cout << "*Rank " << MPI_rank() << " eliminated " << vPtr->idx << std::endl;
        elim_cost += vPtr->eliminate(edge_arena,adj_arena);
        retire(vPtr);
      }else{
        to.insert(std::pair<largeint,Vertex*>(vPtr->idx,vPtr));
      }
    }else if(vPtr->in_degree() && !(vPtr->out_degree())){
      if(!(vPtr->dep_vertex)){
        //std::cout << "*Rank " << MPI_rank() << " eliminated " << vPtr->idx << std::endl;
        elim_cost += vPtr->eliminate(edge_arena,adj_arena);
        retire(vPtr);
      }else{
        to.insert(std::pair<largeint,Vertex*>(vPtr->idx,vPtr));
      }
    }else{
      to.insert(std::pair<largeint,Vertex*>(vPtr->idx,vPtr));
    }
  }
  
  from.clear();
}

void Process::local_vertex_elimination_normal( std::vector<Vertex*> & from , std::vector<cji_t> & mid , std::map<largeint,Vertex*> & to )
{
  int ready = 0;
  int mpi_src_rank = (mpi_rank+1)%mpi_size;
  int cji_count;
  largeint pos = 0;

  for( std::vector<Vertex*>::reverse_iterator it=from.rbegin() ; it!=from.rend() ; it++ ){

    Vertex * vPtr = (*it);

    if(vPtr->is_dead()){
      //std::cout << "Rank " << MPI_rank() << " eliminated " << vPtr->idx << std::endl;
      elim_cost += vPtr->eliminate(edge_arena,adj_arena);
      retire(vPtr);
    }else{
      to.insert( std::pair<largeint,Vertex*>( vPtr->idx , vPtr ) );
    }

    pos++;
 
    if((pos%probe_freq)==0)
    {
      MPI_Iprobe( mpi_src_rank ,MPI_COMM_DATA_MESSAGE , comm , &ready , &status );
    
      if(ready){
        break;
      }
    }
  }

  if(!ready){
     MPI_Probe( mpi_src_rank , MPI_COMM_DATA_MESSAGE , comm , &status );
  }

  MPI_Get_count( &status , MPI_EDGE , &cji_count );

  if(cji_count<0){
    cji_count = 0;
  }

  mid.resize(cji_count);//resize mid to accomodate incoming message

  //&mid.front() is undefined on an empty vector
  MPI_Recv( cji_count?&(mid.front()):NULL , cji_count , MPI_EDGE , mpi_src_rank , MPI_COMM_DATA_MESSAGE , comm , &status );

  //std::cout << "Rank " << mpi_rank << " received " << cji_count-1 << " edges from Rank " << mpi_src_rank << std::endl;

  //now take care of front: pos vertices were consumed from the back
  largeint remaining = (pos<from.size()) ? (from.size()-pos) : 0;//unsigned: never wrap

  for( largeint i=0 ; i<remaining ; i++ ){
    to.insert( std::pair<largeint,Vertex*>( from[i]->idx , from[i] ) );
  }

  from.clear();

  //build from mid
  build( mid , to );
}

void Process::move_vertices( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to )
{
  for( std::vector<Vertex*>::iterator it=from.begin() ; it!=from.end() ; it++ ){

    Vertex * new_vertex = (*it);

    std::map<largeint,Vertex*>::iterator map_it = to.find(new_vertex->idx);

    if(map_it!=to.end()){
      Vertex * old_vertex = map_it->second;
      move_edges(new_vertex,old_vertex);
      retire(new_vertex);
    }else{
      to.insert( std::pair<largeint,Vertex*>( new_vertex->idx,new_vertex) );
    }
  }

  from.clear();
}

void Process::move_vertices_back( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to )
{
  for( std::vector<Vertex*>::iterator it=from.begin() ; it!=from.end() ; it++ ){
    to.insert( std::pair<largeint,Vertex*>( (*it)->idx , (*it) ) );
  }

  from.clear();
}

void Process::global_vertex_elimination( std::map<largeint,Vertex*> & from , std::vector<Vertex*> & to )
{
  for( std::map<largeint,Vertex*>::reverse_iterator it=from.rbegin() ; it!=from.rend() ; it++ ){

    Vertex * vPtr = it->second;
    
    if(vPtr->in_degree() && vPtr->out_degree()){
      if(!(vPtr->dep_vertex)){
        elim_cost += vPtr->eliminate(edge_arena,adj_arena);
        retire(vPtr);
      }else{
        to.push_back(vPtr);
      }
    }else if(vPtr->in_degree() && !(vPtr->out_degree())){
      if(!(vPtr->dep_vertex)){
        elim_cost += vPtr->eliminate(edge_arena,adj_arena);
        retire(vPtr);
      }else{
        to.push_back(vPtr);
      }
    }else{
      to.push_back(vPtr);
    }
  }

  from.clear();
}

void Process::send( std::vector<Vertex*> & from , std::vector<cji_t> & to )
{
  to.push_back( cji_t( 0 , 0 , 0 ) );//dummy

  for( std::vector<Vertex*>::iterator it=from.begin() ; it!=from.end() ; it++ )
  {
    Vertex * vPtr = (*it);
 
    for( Adjacency::iterator eit=vPtr->in_edges.begin() ; eit!=vPtr->in_edges.end() ; eit++ )
    {      
      to.push_back( cji_t(
                             eit->second->src->idx,
                             eit->second->tgt->idx,
                             eit->second->eval ) );

      edge_arena.release(eit->second);
      edge_counter--;//debug
    }
  }

  for( std::vector<Vertex*>::iterator it=from.begin() ; it!=from.end() ; it++ )
  {
    retire(*it);
  }

  from.clear();

  int mpi_dst_rank = (mpi_rank+mpi_size-1)%mpi_size;
  MPI_Send( to.empty()?NULL:&(to.front()) , (int)to.size() , MPI_EDGE , mpi_dst_rank , MPI_COMM_DATA_MESSAGE , comm );

  //std::cout << "Rank " << mpi_rank << " sent " << to.size()-1 << " edges to Rank " << mpi_dst_rank << std::endl;

  to.clear();
}

void Process::build( std::vector<cji_t> & from , std::map<largeint,Vertex*> & to )
{
  Vertex * srcVptr;
  Vertex * tgtVptr;
  bool dummy_found = false;

  for( std::vector<cji_t>::iterator from_it=from.begin() ; from_it!=from.end() ; from_it++ ){
    if(!dummy_found)
    {
      dummy_found = true;
    }else
    {
      largeint src_idx = (*from_it).src;
      largeint tgt_idx = (*from_it).tgt;
      double cji_val  = (*from_it).cji;
    
      std::map<largeint,Vertex*>::iterator to_it = to.find(tgt_idx);//look up for target vertex

      if(to_it!=to.end()){
        tgtVptr = to_it->second;
      }else{
        to.insert( std::pair<largeint,Vertex*>( tgt_idx , tgtVptr = new Vertex(tgt_idx,0) ) );

        vertex_counter++;//debug

        //check if tgtVptr is dep vertex
        if(dep_count==1){
          if(tgt_idx==dependent_index){
            tgtVptr->dep_vertex = true;
          }
        }else{
          if(dep_index_map.find(tgt_idx)!=dep_index_map.end()){
            tgtVptr->dep_vertex = true;
          }
        }
      }
   
      to_it = to.find(src_idx);//look up for source vertex
     
      if(to_it!=to.end()){
        srcVptr = to_it->second;
      }else{
        to.insert( std::pair<largeint,Vertex*>( src_idx , srcVptr = new Vertex(src_idx,0) ) );

        vertex_counter++;//debug

        //check if srcVptr is dep vertex
        if(dep_count==1){
          if(src_idx==dependent_index){
            srcVptr->dep_vertex = true;
          }
        }else{
          if(dep_index_map.find(src_idx)!=dep_index_map.end()){
            srcVptr->dep_vertex = true;
          }
        }
      }

      //an edge src->tgt may already be present: accumulate instead of
      //dropping the contribution on the floor (std::map::insert() keeps the old one)
      Edge * existing = tgtVptr->from( srcVptr );

      if(existing){
        existing->eval += cji_val;
      }else{
        tgtVptr->add_in_edge( srcVptr , cji_val , edge_arena , adj_arena );//add edge
      }
    }
  }
  
  from.clear();
}

void Process::serial_run()
{
  move_vertices( intmed_vec , intmed_map );

  global_vertex_elimination( intmed_map , alive_intmed_vec );

  move_vertices_back( alive_intmed_vec , intmed_map );
}

void Process::vector_run()
{
  local_vertex_elimination_vector( intmed_vec , intmed_map );
    
  global_vertex_elimination( intmed_map , alive_intmed_vec );

  if(next_owner_idx!=1){
    send( alive_intmed_vec , cji_vec );
  }else{
    move_vertices_back( alive_intmed_vec , intmed_map );
  }
  
  has_vector = false;
}

void Process::normal_run()
{
  local_vertex_elimination_normal( intmed_vec , cji_vec , intmed_map );

  global_vertex_elimination( intmed_map , alive_intmed_vec );

  if(next_owner_idx!=1){
    send( alive_intmed_vec , cji_vec );
  }else{
    move_vertices_back( alive_intmed_vec , intmed_map );
  }
}

void Process::check_memory()
{
  largeint part_size = (intmed_count*VERTEX_BYTES + edge_count*EDGE_BYTES);

  if(part_size>=mem_size){

    if(!profiling){

      if(is_proc()){

        if(mpi_size>1){
          if(has_vector){
	    vector_run();
          }else{
	    normal_run();
          }
        }else{
          serial_run();
        }

        top_owner_idx = top_owner_idx-mpi_size;//update top

	intmed_count = 0;
        edge_count = 0;
        //prev_mem_usage = cur_mem_usage;//update memory here
        next_owner_idx++;
        //the same value of prev_mem_usage here could be again used in max_rank_check_memory()
        //again if it is the last partition so we have to update it here so that when
        //max_rank_check_mememory() is called part_size will equal zero.
       
        dtor_ignore_vertex = true;
        partition_count++;

        /*
         * Ported from Maxwell (SVEGP-32).  The target partition is recorded
         * and shipped, so this pass has nothing left to do.  BREAK_ON_TARGET
         * abandons the rest of it; RUN_TO_END returns and lets the section
         * finish on its own.
         *
         * Falling through is safe and is not a new code path.  top_owner_idx
         * has just been decremented by mpi_size and next_owner_idx
         * incremented, so this rank's target is now strictly BELOW
         * next_owner_idx, which only ever increases -- is_proc() cannot
         * become true again before reinitialize() resets them.  The rest of
         * the pass therefore executes exactly as the partitions BEFORE the
         * target already do on every pass: the counters advance so the
         * boundaries stay put, vertex_on_lhs() returns NULL, and nothing is
         * allocated or sent.
         *
         * THE PARALLEL ARGUMENT FOR PREFERRING IT.  A throw leaves the
         * caller's section at a point that depends on the budget and appears
         * nowhere in the source, so the ranks that throw and the ranks that
         * do not make different sequences of calls within one pass.  That is
         * one of the two things standing between this library and a section
         * that may contain MPI; pass_t in Typedefs.hpp has the other, and the
         * measurement showing that neither is sufficient yet.
         */
        if(break_mode==BREAK_ON_TARGET){
          throw BreakException();
        }
      }
    }

    next_owner_idx++;
    intmed_count = 0;
    edge_count = 0;
  }
}

void Process::max_rank_check_memory()
{
  //largeint cur_mem_usage = (intmed_count*sizeof(Vertex) + edge_count*sizeof(Edge));
  //largeint part_size = cur_mem_usage - prev_mem_usage;

  largeint part_size = (intmed_count*VERTEX_BYTES + edge_count*EDGE_BYTES);

  if(part_size){

    if(!profiling){

      if(is_max_rank()){

        if(mpi_size>1){
          vector_run();
        }else{
          serial_run();
        }

        top_owner_idx = top_owner_idx-mpi_size;//update top

        dtor_ignore_vertex = true;
        partition_count++;
        //cannot throw here because max_rank_check_memory() is called inside checkpoint()
        //and there is no need to throw here because it is the last partition anyway.
      }
    }

    //prev_mem_usage = cur_mem_usage;
    next_owner_idx++;
    intmed_count = 0;
    edge_count = 0;
  }
}

void Process::add_edge( Vertex * src , Vertex * tgt , double eval )
{
  if(src && tgt){
    tgt->add_in_edge(src,eval,edge_arena,adj_arena);
  }
}

//Maxwell SVEGP-30: the only way a Vertex is destroyed
void Process::retire( Vertex * v )
{
  if(!v) return;

  v->in_edges.clear(adj_arena);
  v->out_edges.clear(adj_arena);

  delete v;
  vertex_counter--;//debug
}

void Process::initialize_mpi_env()
{
  MPI_Comm_rank( comm , &mpi_rank );
  MPI_Comm_size( comm , &mpi_size );

  struct cji_t e(0,0,0.0);//dummy element
  
  int lengths[3] = { 1 , 1 , 1 };

  MPI_Aint offsets[3];
  MPI_Aint base_addr, addr;
  MPI_Get_address(&e,&base_addr);	
  MPI_Get_address(&e.src,&addr);
  offsets[0] = addr - base_addr;
  MPI_Get_address(&e.tgt,&addr);
  offsets[1] = addr - base_addr;
  MPI_Get_address(&e.cji,&addr);
  offsets[2] = addr - base_addr; 

  /*
   * MPI_UINT64_T, not MPI_UNSIGNED_LONG: idx_t is std::uint64_t now, and on
   * an LLP64 platform (Windows) unsigned long is four bytes.  A heterogeneous
   * job would then have disagreed about the width of a vertex index and read
   * every edge after the first at the wrong offset.
   */
  MPI_Datatype types[3] = { MPI_UINT64_T , MPI_UINT64_T , MPI_DOUBLE };

  /*
   * ...and the committed type is RESIZED to sizeof(cji_t).  Without this the
   * type's extent is the offset of its last member plus that member's size,
   * which equals sizeof(cji_t) only while the struct has no tail padding.
   * Sending an ARRAY of them -- which is the only way this type is ever used
   * -- strides by the extent, so the day a field is added or reordered and
   * padding appears, every element after the first is read from the wrong
   * address.  Pinning the extent to the real object size makes that
   * impossible rather than unlikely.
   */
  MPI_Datatype packed;
  MPI_Type_create_struct(3,lengths,offsets,types,&packed);
  MPI_Type_create_resized(packed,0,(MPI_Aint)sizeof(cji_t),&MPI_EDGE);
  MPI_Type_free(&packed);
  MPI_Type_commit(&MPI_EDGE);
}

void Process::initialize( largeint indep_count , largeint dep_count , largeint mem_size , MPI_Comm comm )
{
  this->comm = comm;

  initialize_mpi_env();

  this->indep_count = indep_count;
  this->dep_count = dep_count;
  this->mem_size = mem_size;

  /*
   * Ported from Maxwell (SVEGP-31) : the arenas are told the budget before
   * anything is recorded.  They are the only two places where the library
   * claims storage in units unrelated to the budget, and until they were
   * told, a tight budget could not actually make a rank small -- the graph
   * obeyed the budget and the arena capacity under it did not.  On a cluster
   * that is the difference between a per-rank budget that fits in cache and
   * one that does not.
   */
  this->edge_arena.set_budget(mem_size);
  this->adj_arena.set_budget(mem_size);

  if(is_max_rank()){//initially, max rank owns the vector
    this->has_vector = true;
  }
}

void Process::reinitialize()
{
  if(!profiling)
  {
    this->throwable = true;
    this->dtor_ignore_vertex = false;
    this->next_vertex_idx = this->indep_count+1;//starts from indep_count+1
    this->next_owner_idx = 1;//starts from next_owner_idx=1
    this->intmed_count = 0;
    this->edge_count = 0;
    this->prev_mem_usage = 0;
    this->productive_pass++;//one execution of the caller's section
  }
}

void Process::finalize()
{
  if(profiling)
  {
    //largeint cur_mem_usage = (intmed_count*sizeof(Vertex) + edge_count*sizeof(Edge));
    //largeint part_size = cur_mem_usage - prev_mem_usage;
   
    largeint part_size = (intmed_count*VERTEX_BYTES + edge_count*EDGE_BYTES);
 
    if(part_size){
      //prev_mem_usage = cur_mem_usage;
      next_owner_idx++;
      intmed_count = 0;
      edge_count = 0;
    }

    top_owner_idx = next_owner_idx-1;
    topmost_owner_idx = top_owner_idx;
  }else
  {
    /*
     * Maxwell SVEGP-27 : this used to BE the teardown, and only on the final
     * rank, which is why every OTHER rank leaked its whole graph on every
     * run -- and why a rank that died any other way leaked it too.  It
     * delegates now, ~Process() calls the same function, and running both is
     * harmless: the second finds the containers empty.
     *
     * Doing it on every rank is the fix, not a widening.  A non-final rank
     * has already sent its surviving vertices away and had send() retire
     * them, so what destroy_graph() finds there is whatever the pipeline left
     * behind -- usually nothing, and on an aborted run, everything.
     */
    destroy_graph();

    dtor_ignore_vertex = true;

    //std::cout << "Rank " << mpi_rank << " , vertex_counter = " << vertex_counter << std::endl;
    //std::cout << "Rank " << mpi_rank << " , edge_counter   = " << edge_counter << std::endl;
  }
}

void Process::register_indep_vertex( const active & x )
{
  if(profiling)
  {
    x.gen = pass_gen;//an independent belongs to every pass; restore_values re-stamps it
    x.reachable = true;
    x.idx = next_vertex_idx;
    //x.owner_idx = next_owner_idx;//redundant ? 

    if(next_vertex_idx==indep_count)
    {
      next_owner_idx++;
    }

    next_vertex_idx++;
  }
}

void Process::register_dep_vertex( const active & x )
{
  if(!profiling)
  {
    if(is_final_rank())
    {
      std::map<largeint,Vertex*>::iterator map_it = intmed_map.find(x.idx);

      if(map_it!=intmed_map.end()){
        dep_vec.push_back(map_it->second);
      }else{
        if(x.idx>0 && x.idx<=indep_count)
        {
          Vertex * vPtr = new Vertex(x.idx,0);
	  
	  vertex_counter++;//debug	

          dep_vec.push_back(vPtr);
          intmed_map.insert( std::pair<largeint,Vertex*>(vPtr->idx,vPtr) );
        }else
        {
	  dep_vec.push_back(NULL);
        }
      }  
    }
  } 
}

void Process::unary_op( const active & x1 ,  double dy_dx1 , const active & x2 , bool overwrite )
{
  if(profiling){
    if(x2.reachable){
      x2.idx = next_vertex_idx++;
      intmed_count++;
      edge_count++;
    }else{
      x2.idx = 0;
    }
  }else{
    Vertex * v1 = NULL , * v2 = NULL;

    if(overwrite){//assignment operator
      set_vertex_dead(x2);
    }

    x2.owner_idx = next_owner_idx;

    if(x2.reachable){    
      x2.idx = next_vertex_idx++;
      v1 = vertex_on_rhs(x1);
      v2 = vertex_on_lhs(x2);
      add_edge(v1,v2,dy_dx1);
      edge_count++;
    }
  }

  check_memory();  
}

void Process::binary_op(  const active & x1 , double dy_dx1 , const active & x2 , double dy_dx2  , const active & x3 )
{
  if(profiling){
    if(x3.reachable){

      x3.idx = next_vertex_idx++;
      intmed_count++;

      if(&x1==&x2){
        edge_count++;
      }else{

        if(x1.reachable){
          edge_count++;
        }
     
        if(x2.reachable){
          edge_count++;
        }
      }
    }else{
      x3.idx = 0;//as in unary_op(): an unreachable result carries no vertex
    }
  }else{
    Vertex * v1 = NULL , * v2 = NULL , * v3 = NULL;

    if(x3.reachable){
      if(x1.reachable){
        v1 = vertex_on_rhs(x1);
      }

      if(x2.reachable){
        v2 = vertex_on_rhs(x2);
      } 
    }

    x3.owner_idx = next_owner_idx;

    if(x3.reachable){
   
      x3.idx = next_vertex_idx++;
      v3 = vertex_on_lhs(x3);

      if(&x1==&x2){
        add_edge(v1,v3,dy_dx1+dy_dx2);
        edge_count++;	
      }else{

        if(x1.reachable){		
          add_edge(v1,v3,dy_dx1);
	  edge_count++;
        }

        if(x2.reachable){ 	
          add_edge(v2,v3,dy_dx2);
	  edge_count++;	
        }		
      }
    }else{
      x3.idx = 0;
    }
  }

  check_memory();
}

void Process::unary_op_ass( double dy_dx , const active & x )
{
  if(profiling)
  {
    if(x.reachable){
      x.idx = next_vertex_idx++;
      intmed_count++;
      edge_count++;
    }else{
      x.idx = 0;
    }
  }else
  {
    Vertex * v1 = NULL , * v2 = NULL;
  
    v1 = set_vertex_dead_for_unary_op_ass(x);//x = x op a , extract x.vtx and save it in v1
    x.owner_idx = next_owner_idx;

    if(x.reachable){
      x.idx = next_vertex_idx++;
      v2 = vertex_on_lhs(x);
      add_edge(v1,v2,dy_dx);
      edge_count++;
    }
  }

  check_memory();
}

void Process::binary_op_ass( const active & x1 , double dy_dx1 , const active & x2 , double dy_dx2 )
{
  if(profiling){
    if(x2.reachable){

      intmed_count++;

      if(&x1==&x2){
        edge_count++;
      }else{

        if(x1.reachable){
          edge_count++;
        }

        if(x2.idx){
          edge_count++;
        }
      }

      x2.idx = next_vertex_idx++;//update x2.idx
    }else{
      x2.idx = 0;
    }
  }else
  {
    Vertex *v1 = NULL , *v2 = NULL , *v3 = NULL;
  
    largeint old_x2_idx = 0;
    bool meaningful = true;

    v2 = set_vertex_dead_for_binary_op_ass(x2,meaningful);//v2 could be NULL

    if(meaningful)    
      old_x2_idx = x2.idx;//save x2.idx before it gets overwritten
    else
      old_x2_idx = 0;

    x2.owner_idx = next_owner_idx;

    if(x2.reachable){

      x2.idx = next_vertex_idx++;//x2.idx gets overwritten here
      v3 = vertex_on_lhs(x2);

      if(&x1==&x2){
         add_edge(v2,v3,dy_dx1+dy_dx2);
         edge_count++;
      }else{

        if(x1.reachable){
          v1 = vertex_on_rhs(x1);
          add_edge(v1,v3,dy_dx1);
          edge_count++;
        }
       
        if(old_x2_idx)//check if x2 was previously reachable before the overwrite
        {
          add_edge(v2,v3,dy_dx2);
          edge_count++;
        }
      }
    }else{
      x2.idx = 0;
    }
  }

  check_memory();
}

void Process::postfix_op( const active & x1 , const active & x2 )
{
  if(profiling){
    if(x2.reachable){
      x2.idx = next_vertex_idx++;
      intmed_count++;
      edge_count++;
    }else{
      x2.idx = 0;
    }
  }else
  {
    Vertex * v1 = NULL , * v2 = NULL;

    if(x1.reachable){
      v1 = vertex_on_rhs(x1);
    }

    x2.owner_idx = next_owner_idx;

    if(x2.reachable){
      x2.idx = next_vertex_idx++;
      v2 = vertex_on_lhs(x2);
      add_edge(v1,v2,1.0);
      edge_count++;
    }else{
      x2.idx = 0;
    }
  }

  check_memory();
}

void Process::passive_op( const active & x )
{
  if(profiling){
    x.idx = 0;
  }else{
    set_vertex_dead(x);
    x.owner_idx = next_owner_idx;
  }
}

void Process::destructor( const active & x )//to be made inline
{
  if(!profiling){
    if(is_proc()){
      if(!dtor_ignore_vertex){
        if(x.idx){
          Vertex * vPtr = x.vtx;
          if(vPtr){
	    if(x.owner_idx==next_owner_idx){
	      vPtr->kill();
            }
          }
        }
      }
    }
  }
}

void Process::harvest( largeint nrows , largeint ncols ,  double **& A , bool print_out )//default : print_out = true
{
  if(!profiling)
  {
    if(is_final_rank())
    {
       A = new double*[nrows];
  
       for( largeint i=0 ; i<nrows ; i++ ){
         A[i] = new double[ncols]; 
           for( largeint j=0 ; j<ncols ; j++ ){
             A[i][j] = 0.0;
           }
       }

       const char *FMT = "%4lu %4lu %15.8E\n";
       largeint dep_pos = 0;

       /*
        * Maxwell's audit left this open as a known defect: harvest() printed
        * intmed_map.size() to stderr and the FJAC banner to stdout
        * UNCONDITIONALLY, straight through the quiet path that print_out
        * exists to provide.  Here it was p times as loud, and it made the
        * assertion suites unreadable before they had asserted anything.
        */
       if(print_out){
         std::cerr << "intmed_map.size() = " << this->intmed_map.size() << std::endl;
         std::cout << "FJAC:" << std::endl;
       }

       for( std::vector<Vertex*>::iterator dep_it=this->dep_vec.begin() ; dep_it!=this->dep_vec.end() ; dep_it++ )
       {
         Vertex * dep_vPtr = (*dep_it);

         if(dep_vPtr)
         {
           for( largeint i=0 ; i<indep_count ; i++ )
           {
             std::map<largeint,Vertex*>::iterator indep_it = this->intmed_map.find(i+1);
      
             if(indep_it!=intmed_map.end())
             {
               Vertex * indep_vPtr = indep_it->second;
               Edge * e = dep_vPtr->from(indep_vPtr); 
 
               largeint col = indep_vPtr->get_idx()-1;

               if(e)
               {
                 if(print_out)
                   printf( FMT , (unsigned long)e->get_src()->get_idx() , (unsigned long)(dep_pos+1) , e->get_partial() );
	
	         if(dep_pos<nrows && col<ncols){
	           A[dep_pos][col] = e->get_partial();
	         }
               }
               else if(indep_vPtr==dep_vPtr)
               {
                 if(print_out)
	           printf( FMT , (unsigned long)indep_vPtr->get_idx() , (unsigned long)(dep_pos+1) , 1.0 );

                 if(dep_pos<nrows && col<ncols){
                   A[dep_pos][col] = 1.0;
                 }
               }
             }
          }
        }

        dep_pos++;//move on

        if(print_out)
          std::printf("\n");
      }
    }
  }
}

largeint Process::get_memory()
{
  return intmed_count*VERTEX_BYTES + edge_count*EDGE_BYTES;
}

/*
 * What the graph actually holds, as opposed to what the budget is spent
 * against.  The difference is allocator granularity -- arena blocks, slabs
 * and free lists -- and examples/invariants watches the ratio so that it
 * cannot drift unnoticed.
 */
largeint Process::get_heap()
{
  largeint verts = 0;

  for( std::map<largeint,Vertex*>::iterator it=intmed_map.begin() ; it!=intmed_map.end() ; it++ ) verts++;

  verts += intmed_vec.size() + alive_intmed_vec.size();

  return verts*sizeof(Vertex) + edge_arena.bytes() + adj_arena.bytes();
}

largeint Process::get_cost()
{
  return elim_cost;
}

//Maxwell SVEGP-22: how many partitions THIS rank recorded
largeint Process::get_partitions()
{
  return partition_count;
}

/*
 * How many partitions the whole tape was broken into, across every rank --
 * AND IT COSTS NO MESSAGE.
 *
 * topmost_owner_idx is set by finalize() at the end of the profiling pass,
 * and the profiling pass is identical on every rank: it builds no graph, only
 * counts, and the counting is a deterministic function of the operation
 * sequence.  So every rank already knows the total, and has since before the
 * first productive pass.  terminate() and the PASSES_COLLECTIVE branch both
 * rely on exactly that.
 *
 * Which means a caller asking "how big did this tape get?" should not have to
 * MPI_Allreduce get_partitions() -- every example that wanted the number was
 * reducing something the library could have handed over.  Zero before the
 * profiling pass has finished, like get_partitions().
 */
largeint Process::get_total_partitions()
{
  return topmost_owner_idx;
}

void Process::set_break_mode( break_t mode )
{
  break_mode = mode;
}

void Process::set_pass_mode( pass_t mode )
{
  pass_mode = mode;
}

void Process::disable_profiling()
{
  if(profiling){
    profiling = false;
  }
}

void Process::save_dependent_index( largeint idx )
{
  if(dep_count==1){
    dependent_index = idx;
  }else{
    dep_index_map.insert( std::pair<largeint,largeint>( idx , idx ) );
  }
}

void Process::ignore_vertices()
{
  dtor_ignore_vertex = true;
}

bool Process::terminate()
{
  /*
   * PASSES_COLLECTIVE : every rank runs the same number of passes, which is
   * the maximum any rank needs -- ceil(P/p), P being the tape's partition
   * count and p the communicator size.  Both are known to every rank without
   * communicating: topmost_owner_idx is set from the profiling pass, which
   * every rank runs identically.  So this costs no message, only up to one
   * extra passive execution of the caller's section on at most p-1 ranks.
   *
   * A rank with nothing left to record has already driven top_owner_idx past
   * its own target, so is_proc() is false for the whole of the extra pass:
   * the counters advance, vertex_on_lhs() returns NULL, and neither the
   * arena nor the ring is touched.
   */
  if(pass_mode==PASSES_COLLECTIVE)
  {
    const largeint needed = (topmost_owner_idx + (largeint)mpi_size - 1)/(largeint)mpi_size;

    return productive_pass > needed;
  }

  largeint target_owner_idx = mpi_rank+(top_owner_idx-mpi_size+1);
 
  if(!target_owner_idx)
  {
    return true;
  }
  else if(target_owner_idx>topmost_owner_idx)//underflow
  {
    return true;
  }

  return false;
}

bool Process::is_final_rank_public()
{
  return is_final_rank();
}

int Process::MPI_rank()
{
  return mpi_rank;
}

int Process::MPI_size()
{
  return mpi_size;
}

void Process::set_probe_frequency( int freq )
{
  probe_freq = freq;
}

