/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Groebner bases with cofactors.
 *
 * Computes a Groebner basis of an ideal of Q[x1,...,xn] and, for every basis
 * element, the polynomials ("cofactors") that express it as a combination of
 * the generators the ideal was given by.
 *
 * Hooks into CoCoA, which reports the S-polynomials it builds and the
 * reduction steps it takes. The cofactors are reconstructed from these
 * reports, and every combination this class hands out is checked by
 * multiplying it out.
 */

#include "cvc5_private.h"

#ifdef CVC5_USE_COCOA

#ifndef CVC5__THEORY__ARITH__NL__IDEAL__COFACTOR_BASIS_H
#define CVC5__THEORY__ARITH__NL__IDEAL__COFACTOR_BASIS_H

#include <CoCoA/PPMonoid.H>
#include <CoCoA/ring.H>

#include <string>
#include <unordered_map>
#include <vector>

namespace cvc5::internal {
namespace theory {
namespace arith {
namespace nl {

class CofactorBasis
{
 public:
  /**
   * @param polyRing A polynomial ring over the rationals.
   * @param gens The generators of the ideal, which are non-zero elements of
   * polyRing.
   */
  CofactorBasis(const CoCoA::ring& polyRing,
                const std::vector<CoCoA::RingElem>& gens);
  ~CofactorBasis();

  /**
   * Compute the basis.
   *
   * @param timeLimit The time limit in seconds, where 0 means no limit.
   * @return false if the basis could not be computed within the limit, or if
   * the cofactors of one of its elements could not be reconstructed. No other
   * method may be called in this case.
   */
  bool compute(double timeLimit);

  /** Is the ideal the whole ring? */
  bool isTrivial() const;

  /**
   * Reduce p modulo the ideal.
   *
   * @param p The polynomial to reduce.
   * @param nf Set to the normal form of p with respect to the basis. Two
   * polynomials have the same normal form iff their difference is in the
   * ideal.
   * @param cofs Set to polynomials, one for each generator, such that
   * p = nf + sum_i cofs[i] * gens[i].
   * @return false if this identity does not hold for the result computed, in
   * which case nf and cofs must not be used.
   */
  bool reduce(const CoCoA::RingElem& p,
              CoCoA::RingElem& nf,
              std::vector<CoCoA::RingElem>& cofs) const;

 private:
  /** A monic polynomial of the ideal as a combination of the generators */
  struct Entry
  {
    /** The (monic) polynomial */
    CoCoA::RingElem d_poly;
    /** d_poly = sum_i d_cofs[i] * d_gens[i] */
    std::vector<CoCoA::RingElem> d_cofs;
  };

  /** Hook: CoCoA built the S-polynomial s of p and q. */
  void sPoly(CoCoA::ConstRefRingElem p,
             CoCoA::ConstRefRingElem q,
             CoCoA::ConstRefRingElem s);
  /** Hook: CoCoA starts reducing p. */
  void reductionStart(CoCoA::ConstRefRingElem p);
  /**
   * Hook: CoCoA cancels the term with power product t of the polynomial under
   * reduction, using a multiple of q.
   */
  void reductionStep(CoCoA::ConstRefRingElem q, CoCoA::ConstRefPPMonoidElem t);
  /** Hook: the reduction that is in progress resulted in r. */
  void reductionEnd(CoCoA::ConstRefRingElem r);

  /** Install the hooks above into CoCoA. */
  void setHandlers();
  /** Remove them. */
  void unsetHandlers();

  /**
   * CoCoA computes in a ring of its own choice (e.g. with integer
   * coefficients). Returns the same polynomial as an element of d_ring.
   */
  CoCoA::RingElem import(CoCoA::ConstRefRingElem p);
  /**
   * Returns the entry for the monic multiple of p, where p is in d_ring, or
   * nullptr if p is zero or not known.
   */
  const Entry* find(const CoCoA::RingElem& p);
  /**
   * Record that value = sum_i cofs[i] * d_gens[i], where CoCoA says that value
   * is a scalar multiple of expected.
   */
  void record(CoCoA::RingElem value,
              std::vector<CoCoA::RingElem> cofs,
              const CoCoA::RingElem& expected);
  /** Does poly = sum_i cofs[i] * d_gens[i] hold? */
  bool verify(const CoCoA::RingElem& poly,
              const std::vector<CoCoA::RingElem>& cofs) const;
  /** Note that the trace CoCoA gave us could not be followed. */
  void traceFailure(const char* reason);

  /** The polynomial ring, over the rationals */
  CoCoA::ring d_ring;
  /** The generators */
  std::vector<CoCoA::RingElem> d_gens;
  /** The polynomials we know cofactors for, by their string representation */
  std::unordered_map<std::string, Entry> d_known;
  /** The basis, after a successful compute() */
  std::vector<Entry> d_basis;
  /** The polynomial under reduction, as a scalar multiple of CoCoA's */
  CoCoA::RingElem d_cur;
  /** d_cur = sum_i d_curCofs[i] * d_gens[i] */
  std::vector<CoCoA::RingElem> d_curCofs;
  /** Is a reduction in progress that we could follow so far? */
  bool d_curValid;
  /** The number of times we could not follow the trace */
  size_t d_numFailures;
  /** Are our hooks installed? */
  bool d_handlersSet;
};

}  // namespace nl
}  // namespace arith
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5__THEORY__ARITH__NL__IDEAL__COFACTOR_BASIS_H */

#endif /* CVC5_USE_COCOA */
