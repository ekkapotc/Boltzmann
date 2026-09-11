#ifndef BOLTZMANN_INCLUDE
#define BOLTZMANN_INCLUDE

/*
 * The public surface: the active scalar, the free-function API, the Jacobian
 * value type, run_tape() -- which owns the checkpoint loop's try/catch so a
 * caller cannot forget it -- the typedefs those need, and the exception the
 * loop throws.
 *
 * Vertex, Edge, Process and the two arenas are internals and are deliberately
 * NOT included here -- nor is <mpi.h>.  Nothing outside src/ has ever used
 * them, and keeping them out means the graph representation can change
 * without breaking a single caller.  It already has, twice in one round: the
 * edge arena and then sorted adjacency.
 *
 * This header no longer opens the namespace for you.  Say
 *
 *     using namespace boltzmann;
 *
 * in your own translation unit, or qualify: boltzmann::active, boltzmann::sin.
 *
 * A program that also calls MPI itself should include <mpi.h> on its own
 * account.  The library duplicates MPI_COMM_WORLD at initialize() and works
 * on its own communicator, so it does not disturb yours.
 */

#include "inc/Typedefs.hpp"
#include "inc/Active.hpp"
#include "inc/Jacobian.hpp"
#include "inc/API.hpp"
#include "inc/BreakException.hpp"
#include "inc/Run.hpp"

#endif
