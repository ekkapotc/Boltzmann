#ifndef RUN_INCLUDE
#define RUN_INCLUDE

#include "Typedefs.hpp"
#include "Active.hpp"
#include "API.hpp"
#include "BreakException.hpp"

/*
 * Ported from Maxwell (SVEGP-28) : the library owns the checkpoint loop.
 *
 * The protocol is
 *
 *     while( checkpoint(x,y) ){
 *       try{ body(); }catch( BreakException const & ){}
 *     }
 *
 * and the try/catch is the caller's to write, every time, correctly.  That
 * imposes three requirements on the section which were written down nowhere:
 * it must be RE-RUNNABLE, it must be SIDE-EFFECT FREE, and -- under the
 * default BREAK_ON_TARGET -- it must be EXCEPTION-SAFE, because the section
 * is abandoned from inside an operator at a point that depends on the memory
 * budget and appears nowhere in the source.
 *
 * run_tape() takes the third off the caller entirely: the try/catch is inside
 * the library, so it cannot be forgotten or written to catch the wrong thing.
 * The first two are still the caller's, and are stated here because a
 * contract this sharp deserves to be written down.
 *
 * THIS IS NOT A STYLE PREFERENCE.  examples/hybrdj1 had no try/catch -- its
 * section was a bare call to fcn() -- and died with `terminate called after
 * throwing an instance of boltzmann::BreakException`, on every rank, the
 * first time anything in this tree actually ran it.  Two of the three
 * assertion suites written for the Maxwell port leaked (49 432 bytes in 167
 * allocations, and 2 880 in 9) through scratch arrays stranded by that same
 * throw.  Three out of three authors got it wrong, and all three were me.
 *
 * WHAT run_tape() DOES NOT FIX, because it cannot:
 *
 *   - It does not make the loop COLLECTIVE.  A rank leaves the loop as soon
 *     as it has recorded its own last partition, so at 8 partitions over 3
 *     ranks one rank runs three passes and two run four.  set_pass_mode
 *     (PASSES_COLLECTIVE) is what equalises that.
 *
 *   - It does not make MPI in the section safe.  The library blocks in
 *     MPI_Probe from inside the caller's operators, so a rank can sit in the
 *     library's probe while another sits in the user's collective.  See
 *     pass_t in Typedefs.hpp.
 *
 * Under set_break_mode(RUN_TO_END) nothing is thrown and the catch below
 * simply never fires, so RUN_TO_END and run_tape() are complements rather
 * than alternatives: the first removes the throw, the second makes it
 * harmless if you have not removed it.
 *
 * It returns the number of passes it ran -- one profiling pass plus one per
 * partition this rank recorded -- because that number was previously
 * unobservable and is the first thing you want when a budget behaves oddly.
 * boltzmann::get_partitions() is the other half.
 *
 * These are free functions rather than members of a Tape object because a
 * parallel RAII surface has to make COLLECTIVE misuse unrepresentable, not
 * just local misuse, and that wants designing before it is written.  When
 * that lands, this becomes Tape::run and these stay as the layer underneath.
 */

namespace boltzmann
{

/* one independent vector, one dependent scalar */
template<class Body>
largeint run_tape( active * x , active & y , Body body )
{
  largeint passes = 0;

  while( checkpoint(x,y) ){
    passes++;
    try{ body(); }catch( BreakException const & ){}
  }

  return passes;
}

/* a 2-D independent grid, one dependent scalar */
template<class Body>
largeint run_tape( active ** x , active & y , Body body )
{
  largeint passes = 0;

  while( checkpoint(x,y) ){
    passes++;
    try{ body(); }catch( BreakException const & ){}
  }

  return passes;
}

/* one independent vector, one dependent vector */
template<class Body>
largeint run_tape( active * x , active * y , Body body )
{
  largeint passes = 0;

  while( checkpoint(x,y) ){
    passes++;
    try{ body(); }catch( BreakException const & ){}
  }

  return passes;
}

/* a 2-D independent grid and a 2-D dependent grid.  set_output_array_dimension
 * must have been called, or given to the seven-argument initialize(), because
 * checkpoint() reads the dependent shape on its first pass. */
template<class Body>
largeint run_tape( active ** x , active ** y , Body body )
{
  largeint passes = 0;

  while( checkpoint(x,y) ){
    passes++;
    try{ body(); }catch( BreakException const & ){}
  }

  return passes;
}

}//end of namespace boltzmann

#endif
