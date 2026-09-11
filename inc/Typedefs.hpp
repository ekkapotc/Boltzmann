#ifndef INCLUDE_TYPEDEF_HPP
#define INCLUDE_TYPEDEF_HPP

#include <cstddef>
#include <cstdint>

namespace boltzmann
{

/*
 * Ported from Maxwell (SVEGP-08).  "unsigned long" is 64-bit on Linux/macOS
 * but only 32-bit on Windows/MSVC (LLP64).  largeint is used both as a vertex
 * index and as the MEMORY BUDGET type, so on Windows the budget silently
 * capped at 4 GB and the vertex counter wrapped after ~4e9 recorded
 * operations.
 *
 * PARALLEL-SPECIFIC CONSEQUENCE, which Maxwell does not have.  idx_t crosses
 * the wire: cji_t is serialised into an MPI_Datatype and shipped between
 * ranks.  With "unsigned long" the MPI type had to be MPI_UNSIGNED_LONG, and
 * that is a DIFFERENT WIDTH on a heterogeneous job -- a Windows rank and a
 * Linux rank would disagree about how many bytes a vertex index is and the
 * receiver would read the edge list misaligned.  A fixed-width type has a
 * fixed-width MPI counterpart (MPI_UINT64_T), so the wire format is now
 * pinned by the language rather than by the platform.
 */
typedef unsigned char flag_t;
typedef std::uint64_t idx_t;
typedef idx_t         largeint;

/*
 * Message tags.  These were #defines in a header, which is a macro leaking
 * into every translation unit that includes it; they are ordinary constants
 * now and carry the namespace.
 */
enum mpi_tag_t
{
  MPI_COMM_METADATA_MESSAGE = 1,
  MPI_COMM_DATA_MESSAGE     = 2
};

/*
 * Ported from Maxwell (SVEGP-32) : how a pass ENDS once its target partition
 * has been recorded.
 *
 * BREAK_ON_TARGET  throw BreakException and abandon the rest of the pass.
 *                  What the library has always done.
 *
 * RUN_TO_END       return normally and let the section run to completion.  The
 *                  suffix executes passively for nothing -- roughly twice the
 *                  passive work -- and in exchange the caller's section is
 *                  never abandoned mid-operator.
 *
 * Maxwell's argument for RUN_TO_END is that a throw out of the middle of an
 * operator strands everything the section allocated before it, once per pass.
 * That argument applies here unchanged.  THIS LIBRARY HAS A SHARPER ONE.
 *
 * A throw leaves the caller's section at a point that depends on the memory
 * budget and appears nowhere in the source.  If that section contains ANY
 * collective -- an MPI_Allreduce over a residual, an MPI_Bcast of a parameter,
 * a barrier -- then the rank that throws leaves the collective and the ranks
 * that do not throw block in it forever.  The job does not crash, it HANGS,
 * and it hangs only for some budgets and some rank counts.  There is no
 * exception-safety discipline the caller can adopt that fixes this, because
 * the problem is not that the unwind is unsafe locally; it is that the unwind
 * is not collective.  RUN_TO_END is the only mode under which a section
 * containing a collective is well defined at all.
 */
enum break_t
{
  BREAK_ON_TARGET,
  RUN_TO_END
};

/*
 * FOUND WHILE PORTING MAXWELL'S RUN_TO_END, AND NOT PRESENT IN MAXWELL AT ALL.
 *
 * The checkpoint loop is not collective.  A rank leaves it as soon as IT has
 * recorded its last partition, and with P partitions over p ranks the ranks
 * get ceil/floor shares -- so at P=8, p=3 rank 0 runs three passes of the
 * caller's section while ranks 1 and 2 run four.  Measured, not argued.
 *
 * PASSES_PER_RANK    each rank stops as soon as its own work is done.  The
 *                    default, and what the library has always done.
 *
 * PASSES_COLLECTIVE  every rank runs ceil(P/p) passes, the maximum over
 *                    ranks.  A rank with nothing left to record runs the
 *                    extra passes passively -- is_proc() is false throughout,
 *                    so nothing is allocated, nothing is eliminated and
 *                    nothing goes on the wire; it costs one extra execution
 *                    of the caller's section on at most p-1 ranks, needs no
 *                    message, and does not move the Jacobian by one bit.  In
 *                    exchange checkpoint() returns the same value on every
 *                    rank, so the LOOP is collective.
 *
 * READ THIS BEFORE ASSUMING IT MAKES YOUR SECTION SAFE.  It does not, and
 * neither does RUN_TO_END, and neither do the two together.  Equal pass
 * counts and no mid-operator throw are both NECESSARY for a section that
 * talks to other ranks, and they are still not SUFFICIENT, because of
 * something further down:
 *
 *     THE LIBRARY DOES BLOCKING POINT-TO-POINT MPI FROM INSIDE THE CALLER'S
 *     OPERATORS.  check_memory() runs at the end of every recorded operation,
 *     and when the budget is spent it calls normal_run(), which blocks in
 *     MPI_Probe waiting for the next rank's edge list.  So a rank can be
 *     sitting in MPI_Probe inside active::operator= while another rank is
 *     sitting in MPI_Allreduce inside the caller's own code.  Neither moves.
 *     Confirmed under gdb, two ranks, both stacks captured:
 *
 *       rank 0  PMPI_Probe <- local_vertex_elimination_normal <- normal_run
 *               <- check_memory <- Process::unary_op <- boltzmann::unary_op
 *               <- active::operator= <- main
 *       rank 1  PMPI_Allreduce <- main
 *
 *     Duplicating the communicator (which initialize() does) prevents message
 *     CROSSTALK; it cannot prevent this, because the problem is not which
 *     communicator the messages are on, it is that one rank is blocked in the
 *     library while the other is blocked in the user.
 *
 * So the contract, until the ring exchange is made non-blocking and drained
 * at checkpoint() boundaries instead of mid-operator, is: THE SECTION BETWEEN
 * checkpoint() CALLS MUST NOT CONTAIN ANY MPI.  That was always true and was
 * written down nowhere.  PASSES_COLLECTIVE and RUN_TO_END are the two
 * prerequisites for lifting it, and lifting it is a roadmap item, not a
 * setting.
 */
enum pass_t
{
  PASSES_PER_RANK,
  PASSES_COLLECTIVE
};

/*
 * One serialised edge of the partial Jacobian, as it crosses the ring from
 * rank r to rank r-1.  Kept trivially copyable and free of padding surprises:
 * two fixed-width indices and a double, which is exactly what the committed
 * MPI_Datatype describes.
 */
struct cji_t
{
  cji_t():
  src(0),
  tgt(0),
  cji(0.0)
  {

  }

  cji_t( largeint src , largeint tgt , double cji ):
  src(src),
  tgt(tgt),
  cji(cji)
  {

  }

  idx_t  src;
  idx_t  tgt;
  double cji;
};
  
}

#endif
