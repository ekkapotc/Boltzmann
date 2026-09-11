#ifndef INCLUDE_PROCESS_HPP
#define INCLUDE_PROCESS_HPP

#include "Active.hpp"
#include "Vertex.hpp"
#include "Edge.hpp"

#include <mpi.h>
#include <map>
#include <vector>

namespace boltzmann
{

namespace internals
{

extern largeint vertex_counter;//debug
extern largeint edge_counter;  //debug

class Process
{
public:
  
  /*
   * SVEGP-26 : this used to be a SINGLETON -- a static instance with
   * get_proc_instance()/del_proc_instance() -- sitting beside a dozen
   * file-statics in API.cpp that held everything else about the recording.
   * Two global states that had to agree with each other, and the 2-D stride
   * bug fixed on 2026-09-11 lived in the second one, in code that had no
   * business existing apart from the tape it described.
   *
   * A Process is now an ordinary object owned by one TapeState, and the
   * constructor is public because there is no singleton left to enforce.
   */

  /*
   * THE COMMUNICATOR IS A MEMBER, which is the parallel half of SVEGP-26 and
   * has no counterpart in Maxwell.
   *
   * It was a file-global, duplicated from MPI_COMM_WORLD by initialize().
   * With one tape per process that is indistinguishable from a member; with
   * two, the second initialize() overwrote the first tape's communicator and
   * the first tape's ring traffic then went out on a communicator that
   * finalize() would later free underneath it.  The graph state and the
   * communicator have to be owned by the same object for the same reason the
   * shadow copies do.
   */
  MPI_Comm comm;

  /*
   * Ported from Maxwell (SVEGP-29, SVEGP-30).  The graph's storage: edges come
   * from one arena and adjacency blocks from another, so an edge costs no
   * allocation of its own and a vertex costs one rather than three.
   */
  EdgeArena edge_arena;
  AdjArena  adj_arena;

  break_t break_mode;//Maxwell SVEGP-32
  pass_t  pass_mode;//parallel-only; see pass_t in Typedefs.hpp

  bool profiling;
  bool throwable;
  bool dtor_ignore_vertex;
  bool has_vector;

  //mpi-related fields
  int mpi_rank;
  int mpi_size;
  int probe_freq;
  MPI_Datatype MPI_EDGE;
  MPI_Status status;

  //memory-relevant fields
  largeint mem_size;
  largeint prev_mem_usage;
  
  //index fields
  largeint next_vertex_idx;
  largeint next_owner_idx;
  largeint top_owner_idx;
  largeint topmost_owner_idx;
  
  //counter fields
  largeint indep_count;
  largeint intmed_count;
  largeint dep_count;
  largeint edge_count;
  largeint elim_cost;//Maxwell: the multiply count of vertex elimination
  largeint partition_count;//Maxwell SVEGP-22: partitions this rank recorded
  largeint productive_pass;//passes of the caller's section run since profiling ended
 
  std::vector<Vertex*> intmed_vec;
  std::vector<Vertex*> alive_intmed_vec;
  std::vector<Vertex*> dep_vec;
  std::vector<cji_t>   cji_vec;

  largeint dependent_index;//for a problem with one output

  std::map<largeint,Vertex*> intmed_map;
  std::map<largeint,largeint> dep_index_map;

  Process();

private:

  //give back every Vertex and both arenas' storage (Maxwell SVEGP-27)
  void destroy_graph();

  inline Vertex * vertex_on_rhs( const active & x );

  inline Vertex * vertex_on_lhs( const active & x );

  inline void set_vertex_dead( const active & x );

  inline Vertex * set_vertex_dead_for_unary_op_ass( const active & x );

  inline Vertex * set_vertex_dead_for_binary_op_ass( const active & x , bool & meaningful );

  inline bool is_proc();

  inline bool is_max_rank();

  inline bool is_final_rank();

  inline void check_memory();

  inline void add_edge( Vertex * src , Vertex * tgt , double eval );

  /*
   * Ported from Maxwell (SVEGP-30) : the only way a Vertex is destroyed.  Its
   * two adjacency blocks belong to adj_arena and have to go back before the
   * vertex does; a Vertex cannot do it itself, having no arena, and giving it
   * one would have cost a back-pointer per vertex -- most of what the arena
   * saves.
   */
  inline void retire( Vertex * v );

  inline void move_edges( Vertex * from_vertex , Vertex * to_vertex );

  inline void serial_run();

  inline void vector_run();

  inline void normal_run();

  inline void initialize_mpi_env();

  inline void move_vertices( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to );

  inline void move_vertices_back( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to );

  inline void send( std::vector<Vertex*> & from , std::vector<cji_t> & to );

  inline void build( std::vector<cji_t> & from , std::map<largeint,Vertex*> & to );

  inline void local_vertex_elimination_vector( std::vector<Vertex*> & from , std::map<largeint,Vertex*> & to );

  inline void local_vertex_elimination_normal( std::vector<Vertex*> & from , std::vector<cji_t> & mid , std::map<largeint,Vertex*> & to );

  inline void global_vertex_elimination( std::map<largeint,Vertex*> & from , std::vector<Vertex*> & to );

public:

  /*
   * Ported from Maxwell (SVEGP-27) : a Process owns every Vertex in
   * intmed_map/intmed_vec and the storage behind both arenas, so destroying
   * one has to give all of it back.  finalize() delegates here rather than
   * being the only place the teardown exists.
   */
  ~Process();

  void initialize( largeint indep_count , largeint dep_count , largeint mem_size , MPI_Comm comm );

  void reinitialize();

  void finalize();

  void register_indep_vertex( const active & x );

  void register_dep_vertex( const active & x );

  void unary_op( const active & x1 ,  double dy_dx1 , const active & x2 , bool overwrite );

  void binary_op( const active & x1 , double dy_dx1 ,  const active & x2 , double dy_dx2  , const active & x3 );
  
  void unary_op_ass( double dy_dx , const active & x );

  void binary_op_ass( const active & x1 , double dy_dx1 , const active & x2 , double dy_dx2 );
 
  void postfix_op( const active & x1 , const active & x2 );
 
  void passive_op( const active & x );

  void destructor( const active & x );
  
  void harvest( largeint m , largeint n ,  double **& A , bool print_out=true );

  largeint get_memory();

  //Maxwell: what the graph actually stores, arenas included
  largeint get_heap();

  largeint get_cost();

  largeint get_partitions();

  //Maxwell SVEGP-32: end a pass by throwing, or by running the section to its end
  void set_break_mode( break_t mode );

  //parallel-only: make the checkpoint loop collective
  void set_pass_mode( pass_t mode );

  void max_rank_check_memory();  

  void disable_profiling();

  void save_dependent_index( largeint dep_idx );

  void ignore_vertices();

  bool terminate();

  int MPI_rank();

  int MPI_size();

  //the rank that assembles the Jacobian; is_final_rank() is private
  bool is_final_rank_public();

  void set_probe_frequency( int probe_freq );
};

}//end of namespace internals
}//end of namespace boltzmann

#endif

