#ifndef JACOBIAN_INCLUDE
#define JACOBIAN_INCLUDE

#include "Typedefs.hpp"

#include <vector>

namespace boltzmann
{

/*
 * Ported from Maxwell (SVEGP-07, SVEGP-28).
 *
 * harvest() hands back a caller-owned double**, and until free_jacobian() was
 * added the library offered no way to release it, so every harvested Jacobian
 * leaked.  A value type retires the whole class of bug -- there is nothing to
 * remember to free, and copying one is a copy rather than a second owner.
 * Storage is one flat row-major block, so a Jacobian is as cheap to copy as
 * its data and has no interior pointers to get wrong.
 *
 * WHAT IS DIFFERENT IN A PARALLEL LIBRARY, and why this type is worth more
 * here than it is in Maxwell.
 *
 * Exactly ONE rank assembles the Jacobian -- the one at the end of the
 * pipeline, where the last partition's edges arrive.  Every other rank has
 * nothing to report.  With the double** interface that fact was invisible:
 * harvest() simply did not touch A on the other ranks, so a caller who wrote
 *
 *     double ** A;            // uninitialised
 *     harvest(m,n,A);
 *     use(A[0][0]);           // fine on one rank, garbage pointer on p-1
 *
 * got undefined behaviour on every rank but one, and got it silently, because
 * the single-rank run they debugged on was the one case that worked.  Every
 * example in this tree is written that way.
 *
 * A Jacobian cannot be in that state.  It is default-constructed empty, and
 * empty() is true on precisely the ranks that did not assemble it -- so the
 * question "am I the rank that has the answer?" has an answer the caller can
 * ask for, rather than one they have to know.  Broadcasting it, if they want
 * it everywhere, is then a well-defined operation over rows()*cols() doubles
 * at data(), which is one contiguous block precisely so that it can be.
 */
class Jacobian
{

public:

  Jacobian(): m_(0), n_(0) {}

  Jacobian( largeint rows , largeint cols ):
  m_(rows),
  n_(cols),
  a_(size_t(rows*cols),0.0)
  {

  }

  largeint rows() const { return m_; }

  largeint cols() const { return n_; }

  /*
   * True when this rank has no Jacobian: it is not the assembling rank, or
   * harvest() refused -- a shape mismatch, or a tape that never ran.
   */
  bool empty() const { return a_.empty(); }

  double operator()( largeint i , largeint j ) const { return a_[size_t(i*n_+j)]; }

  double & operator()( largeint i , largeint j ) { return a_[size_t(i*n_+j)]; }

  //one contiguous row-major block, so MPI_Bcast over it is well defined
  const double * data() const { return a_.empty() ? 0 : &a_[0]; }

  double * data() { return a_.empty() ? 0 : &a_[0]; }

  largeint size() const { return m_*n_; }

private:

  largeint m_;
  largeint n_;
  std::vector<double> a_;

};//end of class

}//end of namespace boltzmann

#endif
