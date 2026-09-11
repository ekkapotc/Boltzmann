#ifndef INCLUDE_BREAK_EXCEPTION_HPP
#define INCLUDE_BREAK_EXCEPTION_HPP

#include <exception>

namespace boltzmann
{

/*
 * Ported from Maxwell.  This lived in boltzmann::internals, and the caller is
 * required to catch it -- so it was public in everything but name, reachable
 * only because a public header said `using namespace boltzmann::internals;` at
 * global scope.  Removing that line is what made the misplacement visible.
 * It is part of the public surface, so it lives in the public namespace.
 *
 * Under RUN_TO_END (Typedefs.hpp) nothing throws this at all and the caller's
 * catch becomes dead code rather than wrong, so an existing caller need not
 * remove it.
 */
class BreakException : public std::exception
{

public:

  virtual const char* what() const throw();

};//end of class

}//end of namespace boltzmann

#endif
