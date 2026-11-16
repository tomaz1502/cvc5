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
 */

#ifdef CVC5_USE_COCOA

#include "theory/arith/nl/ideal/cofactor_basis.h"

#include <CoCoA/BigRat.H>
#include <CoCoA/CpuTimeLimit.H>
#include <CoCoA/SparsePolyIter.H>
#include <CoCoA/SparsePolyOps-RingElem.H>
#include <CoCoA/SparsePolyOps-ideal.H>
#include <CoCoA/SparsePolyRing.H>
#include <CoCoA/TmpGPoly.H>
#include <CoCoA/error.H>
#include <CoCoA/ideal.H>

#include <functional>
#include <sstream>

#include "base/check.h"
#include "base/output.h"

namespace cvc5::internal {
namespace theory {
namespace arith {
namespace nl {

namespace {

std::string toString(const CoCoA::RingElem& p)
{
  std::ostringstream o;
  o << p;
  return o.str();
}

/** Returns 1/c for the leading coefficient c of the non-zero polynomial p. */
CoCoA::BigRat inverseOfLeadingCoeff(const CoCoA::RingElem& p)
{
  CoCoA::BigRat c;
  bool isRat = CoCoA::IsRational(c, CoCoA::LC(p));
  AlwaysAssert(isRat) << "Expected rational coefficients";
  return CoCoA::BigRat(CoCoA::den(c), CoCoA::num(c));
}

}  // namespace

CofactorBasis::CofactorBasis(const CoCoA::ring& polyRing,
                             const std::vector<CoCoA::RingElem>& gens)
    : d_ring(polyRing),
      d_gens(gens),
      d_cur(CoCoA::zero(polyRing)),
      d_curValid(false),
      d_numFailures(0),
      d_handlersSet(false)
{
}

CofactorBasis::~CofactorBasis()
{
  // if an exception unwound the stack while the hooks were installed
  if (d_handlersSet)
  {
    unsetHandlers();
  }
}

void CofactorBasis::setHandlers()
{
  Assert(!CoCoA::handlersEnabled);
  CoCoA::sPolyHandler = [this](CoCoA::ConstRefRingElem p,
                               CoCoA::ConstRefRingElem q,
                               CoCoA::ConstRefRingElem s) { sPoly(p, q, s); };
  CoCoA::reductionStartHandler = [this](CoCoA::ConstRefRingElem p) {
    reductionStart(p);
  };
  CoCoA::reductionStepHandler = [this](CoCoA::ConstRefRingElem q,
                                       CoCoA::ConstRefPPMonoidElem t) {
    reductionStep(q, t);
  };
  CoCoA::reductionEndHandler = [this](CoCoA::ConstRefRingElem r) {
    reductionEnd(r);
  };
  CoCoA::handlersEnabled = true;
  d_handlersSet = true;
}

void CofactorBasis::unsetHandlers()
{
  CoCoA::handlersEnabled = false;
  CoCoA::sPolyHandler = {};
  CoCoA::reductionStartHandler = {};
  CoCoA::reductionStepHandler = {};
  CoCoA::reductionEndHandler = {};
  d_handlersSet = false;
}

bool CofactorBasis::compute(double timeLimit)
{
  const CoCoA::SparsePolyRing ring = d_ring;
  d_known.clear();
  d_basis.clear();
  d_numFailures = 0;
  d_curValid = false;
  // the generators are combinations of themselves
  for (size_t i = 0, n = d_gens.size(); i < n; ++i)
  {
    Assert(!CoCoA::IsZero(d_gens[i]));
    std::vector<CoCoA::RingElem> cofs(n, CoCoA::zero(d_ring));
    cofs[i] = CoCoA::one(d_ring);
    record(d_gens[i], cofs, d_gens[i]);
  }
  std::vector<CoCoA::RingElem> basis;
  try
  {
    setHandlers();
    const CoCoA::ideal ideal(d_gens);
    basis = timeLimit > 0
                ? CoCoA::GBasis(ideal, CoCoA::CpuTimeLimit(timeLimit))
                : CoCoA::GBasis(ideal);
    unsetHandlers();
  }
  catch (const CoCoA::ErrorInfo& e)
  {
    // includes CoCoA::TimeoutException
    unsetHandlers();
    Trace("nl-ideal-gb") << "CofactorBasis: no basis, CoCoA says: "
                         << CoCoA::message(e) << std::endl;
    return false;
  }
  Trace("nl-ideal-gb") << "CofactorBasis: " << d_gens.size()
                       << " generators, " << basis.size()
                       << " basis elements, " << d_numFailures
                       << " trace failures" << std::endl;
  for (const CoCoA::RingElem& b : basis)
  {
    const Entry* e = find(b);
    if (e == nullptr)
    {
      Trace("nl-ideal-gb") << "...no cofactors for " << b << std::endl;
      d_basis.clear();
      return false;
    }
    if (!verify(e->d_poly, e->d_cofs))
    {
      Trace("nl-ideal-gb") << "...wrong cofactors for " << b << std::endl;
      Assert(false) << "Wrong cofactors for " << b;
      d_basis.clear();
      return false;
    }
    d_basis.push_back(*e);
  }
  // only the basis is needed from now on
  d_known.clear();
  return true;
}

bool CofactorBasis::isTrivial() const
{
  return d_basis.size() == 1 && CoCoA::IsConstant(d_basis[0].d_poly);
}

bool CofactorBasis::reduce(const CoCoA::RingElem& p,
                           CoCoA::RingElem& nf,
                           std::vector<CoCoA::RingElem>& cofs) const
{
  const CoCoA::SparsePolyRing ring = d_ring;
  nf = CoCoA::zero(d_ring);
  cofs.assign(d_gens.size(), CoCoA::zero(d_ring));
  // The division algorithm: cancel the leading term of the remainder with a
  // basis element whenever possible, and move it to the normal form otherwise.
  CoCoA::RingElem rem = p;
  while (!CoCoA::IsZero(rem))
  {
    CoCoA::ConstRefPPMonoidElem lpp = CoCoA::LPP(rem);
    const Entry* reducer = nullptr;
    for (const Entry& e : d_basis)
    {
      if (CoCoA::IsDivisible(lpp, CoCoA::LPP(e.d_poly)))
      {
        reducer = &e;
        break;
      }
    }
    if (reducer == nullptr)
    {
      CoCoA::RingElem lt = CoCoA::monomial(ring, CoCoA::LC(rem), lpp);
      nf += lt;
      rem -= lt;
      continue;
    }
    // basis elements are monic
    CoCoA::RingElem mult =
        CoCoA::monomial(ring, CoCoA::LC(rem), lpp / CoCoA::LPP(reducer->d_poly));
    rem -= mult * reducer->d_poly;
    for (size_t i = 0, n = cofs.size(); i < n; ++i)
    {
      if (!CoCoA::IsZero(reducer->d_cofs[i]))
      {
        cofs[i] += mult * reducer->d_cofs[i];
      }
    }
  }
  return verify(p - nf, cofs);
}

bool CofactorBasis::verify(const CoCoA::RingElem& poly,
                           const std::vector<CoCoA::RingElem>& cofs) const
{
  Assert(cofs.size() == d_gens.size());
  CoCoA::RingElem sum = CoCoA::zero(d_ring);
  for (size_t i = 0, n = cofs.size(); i < n; ++i)
  {
    if (!CoCoA::IsZero(cofs[i]))
    {
      sum += cofs[i] * d_gens[i];
    }
  }
  return sum == poly;
}

CoCoA::RingElem CofactorBasis::import(CoCoA::ConstRefRingElem p)
{
  if (CoCoA::owner(p) == d_ring)
  {
    return CoCoA::RingElem(p);
  }
  const CoCoA::SparsePolyRing ring = d_ring;
  CoCoA::RingElem res = CoCoA::zero(d_ring);
  for (CoCoA::SparsePolyIter it = CoCoA::BeginIter(p); !CoCoA::IsEnded(it);
       ++it)
  {
    CoCoA::BigRat c;
    bool isRat = CoCoA::IsRational(c, CoCoA::coeff(it));
    AlwaysAssert(isRat) << "Expected rational coefficients";
    res += CoCoA::monomial(ring, c, CoCoA::exponents(CoCoA::PP(it)));
  }
  return res;
}

const CofactorBasis::Entry* CofactorBasis::find(const CoCoA::RingElem& p)
{
  if (CoCoA::IsZero(p))
  {
    return nullptr;
  }
  auto it = d_known.find(toString(p * inverseOfLeadingCoeff(p)));
  return it == d_known.end() ? nullptr : &it->second;
}

void CofactorBasis::record(CoCoA::RingElem value,
                           std::vector<CoCoA::RingElem> cofs,
                           const CoCoA::RingElem& expected)
{
  if (CoCoA::IsZero(value))
  {
    if (!CoCoA::IsZero(expected))
    {
      traceFailure("expected a non-zero polynomial");
    }
    return;
  }
  if (CoCoA::IsZero(expected))
  {
    traceFailure("expected zero");
    return;
  }
  // make it monic
  CoCoA::BigRat inv = inverseOfLeadingCoeff(value);
  value = value * inv;
  if (value != expected * inverseOfLeadingCoeff(expected))
  {
    traceFailure("unexpected polynomial");
    return;
  }
  std::string key = toString(value);
  if (d_known.find(key) != d_known.end())
  {
    return;
  }
  for (CoCoA::RingElem& c : cofs)
  {
    if (!CoCoA::IsZero(c))
    {
      c = c * inv;
    }
  }
  Entry& e = d_known[key];
  e.d_poly = value;
  e.d_cofs = std::move(cofs);
}

void CofactorBasis::traceFailure(const char* reason)
{
  ++d_numFailures;
  Trace("nl-ideal-gb") << "CofactorBasis: cannot follow CoCoA: " << reason
                       << std::endl;
}

void CofactorBasis::sPoly(CoCoA::ConstRefRingElem p,
                          CoCoA::ConstRefRingElem q,
                          CoCoA::ConstRefRingElem s)
{
  const Entry* ep = find(import(p));
  const Entry* eq = find(import(q));
  if (ep == nullptr || eq == nullptr)
  {
    traceFailure("S-polynomial of an unknown polynomial");
    return;
  }
  const CoCoA::SparsePolyRing ring = d_ring;
  // both are monic, hence the leading terms cancel
  CoCoA::PPMonoidElem l = CoCoA::lcm(CoCoA::LPP(ep->d_poly), CoCoA::LPP(eq->d_poly));
  CoCoA::RingElem mp = CoCoA::monomial(ring, 1, l / CoCoA::LPP(ep->d_poly));
  CoCoA::RingElem mq = CoCoA::monomial(ring, 1, l / CoCoA::LPP(eq->d_poly));
  std::vector<CoCoA::RingElem> cofs;
  for (size_t i = 0, n = d_gens.size(); i < n; ++i)
  {
    cofs.emplace_back(mp * ep->d_cofs[i] - mq * eq->d_cofs[i]);
  }
  record(mp * ep->d_poly - mq * eq->d_poly, std::move(cofs), import(s));
}

void CofactorBasis::reductionStart(CoCoA::ConstRefRingElem p)
{
  d_curValid = false;
  const Entry* e = find(import(p));
  if (e == nullptr)
  {
    traceFailure("reduction of an unknown polynomial");
    return;
  }
  d_cur = e->d_poly;
  d_curCofs = e->d_cofs;
  d_curValid = true;
}

void CofactorBasis::reductionStep(CoCoA::ConstRefRingElem q,
                                  CoCoA::ConstRefPPMonoidElem t)
{
  if (!d_curValid)
  {
    return;
  }
  const Entry* e = find(import(q));
  if (e == nullptr)
  {
    traceFailure("reduction by an unknown polynomial");
    d_curValid = false;
    return;
  }
  const CoCoA::SparsePolyRing ring = d_ring;
  // t as a power product of our ring
  CoCoA::PPMonoidElem tt(CoCoA::PPM(ring), CoCoA::exponents(t));
  // the term of d_cur that is cancelled
  CoCoA::RingElem c = CoCoA::zero(CoCoA::CoeffRing(ring));
  for (CoCoA::SparsePolyIter it = CoCoA::BeginIter(d_cur); !CoCoA::IsEnded(it);
       ++it)
  {
    if (CoCoA::PP(it) == tt)
    {
      c = CoCoA::coeff(it);
      break;
    }
  }
  if (CoCoA::IsZero(c) || !CoCoA::IsDivisible(tt, CoCoA::LPP(e->d_poly)))
  {
    traceFailure("reduction step does not apply");
    d_curValid = false;
    return;
  }
  // e->d_poly is monic
  CoCoA::RingElem mult = CoCoA::monomial(ring, c, tt / CoCoA::LPP(e->d_poly));
  d_cur -= mult * e->d_poly;
  for (size_t i = 0, n = d_curCofs.size(); i < n; ++i)
  {
    if (!CoCoA::IsZero(e->d_cofs[i]))
    {
      d_curCofs[i] -= mult * e->d_cofs[i];
    }
  }
}

void CofactorBasis::reductionEnd(CoCoA::ConstRefRingElem r)
{
  if (!d_curValid)
  {
    return;
  }
  d_curValid = false;
  record(d_cur, std::move(d_curCofs), import(r));
  d_curCofs.clear();
}

}  // namespace nl
}  // namespace arith
}  // namespace theory
}  // namespace cvc5::internal

#endif /* CVC5_USE_COCOA */
