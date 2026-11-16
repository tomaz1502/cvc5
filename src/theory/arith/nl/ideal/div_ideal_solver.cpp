/******************************************************************************
 * This file is part of the cvc5 project.
 *
 * Copyright (c) 2009-2026 by the authors listed in the file AUTHORS
 * in the top-level source directory and their institutional affiliations.
 * All rights reserved.  See the file COPYING in the top-level source
 * directory for licensing information.
 * ****************************************************************************
 *
 * Ideal membership reasoning for integer division and modulus.
 */

#include "theory/arith/nl/ideal/div_ideal_solver.h"

#include <algorithm>
#include <sstream>

#include "expr/skolem_manager.h"
#include "options/arith_options.h"
#include "theory/arith/inference_manager.h"
#include "theory/arith/nl/nl_model.h"
#include "theory/rewriter.h"
#include "theory/theory_state.h"
#include "theory/valuation.h"
#include "util/rational.h"

#ifdef CVC5_USE_COCOA
#include <CoCoA/BigInt.H>
#include <CoCoA/BigRat.H>
#include <CoCoA/RingQQ.H>
#include <CoCoA/SparsePolyIter.H>
#include <CoCoA/SparsePolyOps-RingElem.H>
#include <CoCoA/SparsePolyRing.H>
#include <CoCoA/error.H>
#include <CoCoA/symbol.H>

#include "theory/arith/nl/ideal/cofactor_basis.h"
#endif /* CVC5_USE_COCOA */

using namespace cvc5::internal::kind;

namespace cvc5::internal {
namespace theory {
namespace arith {
namespace nl {

namespace {

/** The most equalities we give to a Groebner basis computation */
const size_t s_maxEqualities = 25;
/** The time limit for a Groebner basis computation, in seconds */
const double s_basisTimeLimit = 2.0;

/** The variable that stands for the integer division t = (div s n) */
Node getQuotientVar(TNode t)
{
  Assert(t.getKind() == Kind::INTS_DIVISION_TOTAL);
  return SkolemManager::mkPurifySkolem(t);
}

/** The variable that stands for the quotient of t = (mod s n) */
Node getQuotientVarOfMod(NodeManager* nm, TNode t)
{
  Assert(t.getKind() == Kind::INTS_MODULUS_TOTAL);
  return getQuotientVar(nm->mkNode(Kind::INTS_DIVISION_TOTAL, t[0], t[1]));
}

#ifdef CVC5_USE_COCOA

/**
 * Converts between terms and polynomials of Q[x_0,...,x_k], where x_i stands
 * for the i-th of the given leaves.
 */
class PolyEncoder
{
 public:
  PolyEncoder(NodeManager* nm, const std::vector<Node>& leaves)
      : d_nm(nm), d_leaves(leaves), d_ring(mkRing(leaves.size()))
  {
    for (size_t i = 0, n = leaves.size(); i < n; ++i)
    {
      d_index[leaves[i]] = i;
    }
  }
  const CoCoA::ring& getRing() const { return d_ring; }
  /** The polynomial for a term all of whose leaves were given */
  CoCoA::RingElem encode(TNode n)
  {
    auto itc = d_cache.find(n);
    if (itc != d_cache.end())
    {
      return itc->second;
    }
    const CoCoA::SparsePolyRing ring = d_ring;
    CoCoA::RingElem res = CoCoA::zero(d_ring);
    switch (n.getKind())
    {
      case Kind::CONST_INTEGER:
        res = CoCoA::RingElem(
            d_ring,
            CoCoA::BigIntFromString(
                n.getConst<Rational>().getNumerator().toString()));
        break;
      case Kind::ADD:
        for (const Node& nc : n)
        {
          res += encode(nc);
        }
        break;
      case Kind::SUB: res = encode(n[0]) - encode(n[1]); break;
      case Kind::NEG: res = -encode(n[0]); break;
      case Kind::MULT:
      case Kind::NONLINEAR_MULT:
        res = CoCoA::one(d_ring);
        for (const Node& nc : n)
        {
          res *= encode(nc);
        }
        break;
      case Kind::INTS_MODULUS_TOTAL:
        res = encode(n[0])
              - encode(n[1]) * encodeLeaf(getQuotientVarOfMod(d_nm, n));
        break;
      case Kind::INTS_DIVISION_TOTAL: res = encodeLeaf(getQuotientVar(n)); break;
      default: res = encodeLeaf(n); break;
    }
    d_cache.emplace(n, res);
    return res;
  }
  /**
   * The term for a polynomial, or the null node if it has a coefficient that is
   * not an integer.
   */
  Node decode(const CoCoA::RingElem& p) const
  {
    std::vector<Node> sum;
    for (CoCoA::SparsePolyIter it = CoCoA::BeginIter(p); !CoCoA::IsEnded(it);
         ++it)
    {
      CoCoA::BigRat c;
      if (!CoCoA::IsRational(c, CoCoA::coeff(it)) || !CoCoA::IsOneDen(c))
      {
        return Node::null();
      }
      std::ostringstream cs;
      cs << CoCoA::num(c);
      std::vector<Node> factors;
      factors.push_back(d_nm->mkConstInt(Rational(Integer(cs.str(), 10))));
      std::vector<long> exps = CoCoA::exponents(CoCoA::PP(it));
      for (size_t i = 0, n = exps.size(); i < n; ++i)
      {
        for (long j = 0; j < exps[i]; ++j)
        {
          factors.push_back(d_leaves[i]);
        }
      }
      sum.push_back(factors.size() == 1 ? factors[0]
                                         : d_nm->mkNode(Kind::MULT, factors));
    }
    if (sum.empty())
    {
      return d_nm->mkConstInt(Rational(0));
    }
    return sum.size() == 1 ? sum[0] : d_nm->mkNode(Kind::ADD, sum);
  }

 private:
  static CoCoA::ring mkRing(size_t numLeaves)
  {
    std::vector<CoCoA::symbol> syms;
    for (size_t i = 0; i < numLeaves; ++i)
    {
      syms.push_back(CoCoA::symbol("x", static_cast<long>(i)));
    }
    return CoCoA::NewPolyRing(CoCoA::RingQQ(), syms);
  }
  CoCoA::RingElem encodeLeaf(const Node& n) const
  {
    auto it = d_index.find(n);
    Assert(it != d_index.end()) << "Unknown leaf " << n;
    return CoCoA::indet(d_ring, static_cast<long>(it->second));
  }
  NodeManager* d_nm;
  std::vector<Node> d_leaves;
  std::unordered_map<Node, size_t> d_index;
  CoCoA::ring d_ring;
  std::unordered_map<Node, CoCoA::RingElem> d_cache;
};

/** Are all coefficients of p integers? */
bool hasIntegerCoefficients(const CoCoA::RingElem& p)
{
  for (CoCoA::SparsePolyIter it = CoCoA::BeginIter(p); !CoCoA::IsEnded(it);
       ++it)
  {
    CoCoA::BigRat c;
    if (!CoCoA::IsRational(c, CoCoA::coeff(it)) || !CoCoA::IsOneDen(c))
    {
      return false;
    }
  }
  return true;
}

std::string toString(const CoCoA::RingElem& p)
{
  std::ostringstream o;
  o << p;
  return o.str();
}

#endif /* CVC5_USE_COCOA */

/** Add to dest the elements of src that it does not contain. */
void addNew(std::vector<Node>& dest, const std::vector<Node>& src)
{
  for (const Node& n : src)
  {
    if (std::find(dest.begin(), dest.end(), n) == dest.end())
    {
      dest.push_back(n);
    }
  }
}

}  // namespace

DivIdealSolver::DivIdealSolver(Env& env,
                               TheoryState& state,
                               InferenceManager& im,
                               NlModel& model)
    : EnvObj(env),
      d_astate(state),
      d_im(im),
      d_model(model),
      d_derivedSent(false)
{
  d_true = nodeManager()->mkConst(true);
  d_zero = nodeManager()->mkConstInt(Rational(0));
#ifndef CVC5_USE_COCOA
  if (options().arith.nlDivIdeal)
  {
    warning() << "--nl-div-ideal has no effect, since this build does not "
                 "include CoCoA. Configure with --cocoa --gpl."
              << std::endl;
  }
#endif
}

DivIdealSolver::~DivIdealSolver() {}

void DivIdealSolver::collectDivTerms(TNode n,
                                     std::unordered_set<Node>& visited,
                                     std::vector<Node>& quots)
{
  SkolemManager* sm = nodeManager()->getSkolemManager();
  std::vector<Node> visit{n};
  while (!visit.empty())
  {
    Node cur = visit.back();
    visit.pop_back();
    if (!visited.insert(cur).second)
    {
      continue;
    }
    if (cur.getKind() == Kind::SKOLEM)
    {
      auto it = d_divTerms.find(cur);
      if (it == d_divTerms.end())
      {
        DivTerm dt;
        if (sm->getId(cur) == SkolemId::PURIFY)
        {
          Node u = SkolemManager::getUnpurifiedForm(cur);
          // divisions by constants are linear
          if (u.getKind() == Kind::INTS_DIVISION_TOTAL && !u[1].isConst()
              && u[0].getType().isInteger())
          {
            dt.d_quot = cur;
            dt.d_num = u[0];
            dt.d_den = u[1];
          }
        }
        it = d_divTerms.emplace(cur, dt).first;
      }
      if (!it->second.d_quot.isNull())
      {
        quots.push_back(cur);
      }
      continue;
    }
    for (const Node& nc : cur)
    {
      visit.push_back(nc);
    }
  }
}

bool DivIdealSolver::collectLeaves(TNode n, std::vector<Node>& leaves)
{
  auto it = d_leaves.find(n);
  if (it == d_leaves.end())
  {
    std::pair<bool, std::vector<Node>> res(true, {});
    switch (n.getKind())
    {
      case Kind::CONST_INTEGER: break;
      case Kind::ADD:
      case Kind::SUB:
      case Kind::NEG:
      case Kind::MULT:
      case Kind::NONLINEAR_MULT:
      case Kind::INTS_MODULUS_TOTAL:
        for (const Node& nc : n)
        {
          res.first = res.first && collectLeaves(nc, res.second);
        }
        if (n.getKind() == Kind::INTS_MODULUS_TOTAL)
        {
          addNew(res.second, {getQuotientVarOfMod(nodeManager(), n)});
        }
        break;
      case Kind::INTS_DIVISION_TOTAL:
        res.second.push_back(getQuotientVar(n));
        break;
      default:
        // an integer term that we treat as a variable
        res.first = n.getType().isInteger();
        res.second.push_back(n);
        break;
    }
    it = d_leaves.emplace(n, res).first;
  }
  addNew(leaves, it->second.second);
  return it->second.first;
}

bool DivIdealSolver::isAsserted(const Node& lit)
{
  bool pol = lit.getKind() != Kind::NOT;
  TNode atom = pol ? lit : lit[0];
  Valuation& val = d_astate.getValuation();
  bool value;
  return val.isSatLiteral(atom) && val.hasSatValue(atom, value) && value == pol;
}

void DivIdealSolver::check(const std::vector<Node>& assertions)
{
  d_derived.clear();
  d_derivedSent = false;
#ifdef CVC5_USE_COCOA
  NodeManager* nm = nodeManager();
  // the divisions that occur in the assertions, by their divisor
  std::unordered_set<Node> visited;
  std::vector<Node> quots;
  for (const Node& a : assertions)
  {
    collectDivTerms(a, visited, quots);
  }
  if (quots.empty())
  {
    return;
  }
  Trace("nl-div-ideal") << "divs = \n";
  for (auto q: quots) {
    Trace("nl-div-ideal") << q << "\n";
  }
  Trace("nl-div-ideal") << "-------\n";
  std::sort(quots.begin(), quots.end());
  std::map<Node, std::vector<DivTerm>> groups;
  for (const Node& q : quots)
  {
    const DivTerm& dt = d_divTerms[q];
    groups[dt.d_den].push_back(dt);
  }
  // the asserted equalities between integer terms
  std::vector<Node> eqs;
  for (const Node& a : assertions)
  {
    if (a.getKind() == Kind::EQUAL && a[0].getType().isInteger())
    {
      eqs.push_back(a);
    }
  }
  Trace("nl-div-ideal") << "DivIdealSolver::check: " << quots.size()
                        << " divisions, " << groups.size() << " divisors, "
                        << eqs.size() << " equalities" << std::endl;
  for (const std::pair<const Node, std::vector<DivTerm>>& g : groups)
  {
    const Node& den = g.first;
    // The sign of the divisor in the model decides which of the premises
    // den > 0 and den < 0 we use.
    Node denv = d_model.computeAbstractModelValue(den);
    if (!denv.isConst() || denv.getConst<Rational>().sgn() == 0)
    {
      continue;
    }
    bool isPos = denv.getConst<Rational>().sgn() > 0;
    Node sign =
        rewrite(nm->mkNode(isPos ? Kind::GT : Kind::LT, den, d_zero));
    if (sign.isConst() && !sign.getConst<bool>())
    {
      continue;
    }
    // The terms we work with, as polynomials. We take the equalities that are
    // connected to the divisions via the variables they contain.
    std::vector<Node> leaves;
    std::vector<DivTerm> divs;
    if (!collectLeaves(den, leaves))
    {
      continue;
    }
    for (const DivTerm& dt : g.second)
    {
      std::vector<Node> dleaves;
      if (collectLeaves(dt.d_num, dleaves))
      {
        divs.push_back(dt);
        addNew(leaves, dleaves);
      }
    }
    if (divs.empty())
    {
      continue;
    }
    std::vector<Node> geqs;
    std::unordered_set<Node> vars(leaves.begin(), leaves.end());
    bool changed = true;
    while (changed && geqs.size() < s_maxEqualities)
    {
      changed = false;
      for (const Node& eq : eqs)
      {
        if (std::find(geqs.begin(), geqs.end(), eq) != geqs.end())
        {
          continue;
        }
        std::vector<Node> eleaves;
        if (!collectLeaves(eq[0], eleaves) || !collectLeaves(eq[1], eleaves))
        {
          continue;
        }
        if (std::none_of(eleaves.begin(), eleaves.end(), [&vars](const Node& l) {
              return vars.find(l) != vars.end();
            }))
        {
          continue;
        }
        geqs.push_back(eq);
        vars.insert(eleaves.begin(), eleaves.end());
        changed = true;
        if (geqs.size() >= s_maxEqualities)
        {
          break;
        }
      }
    }
    std::sort(geqs.begin(), geqs.end());
    // compute the inferences, unless we did so already
    std::vector<Node> keyc;
    keyc.push_back(den);
    for (const DivTerm& dt : divs)
    {
      keyc.push_back(dt.d_quot);
    }
    keyc.insert(keyc.end(), geqs.begin(), geqs.end());
    Node key = nm->mkNode(Kind::SEXPR, keyc);
    auto itc = d_cache.find(key);
    if (itc == d_cache.end())
    {
      itc = d_cache.emplace(key, computeInferences(den, divs, geqs)).first;
    }
    // send those that the model violates
    for (const Inference& inf : itc->second)
    {
      std::vector<Node> premises = inf.d_eqs;
      bool allAsserted = true;
      if (!sign.isConst())
      {
        premises.push_back(sign);
        allAsserted = isAsserted(sign);
      }
      if (!inf.d_rangeTerm.isNull())
      {
        // 0 <= t < |den| must hold in the model
        Node absDen = isPos ? den : nm->mkNode(Kind::NEG, den);
        Node lower = rewrite(nm->mkNode(Kind::GEQ, inf.d_rangeTerm, d_zero));
        Node upper = rewrite(nm->mkNode(Kind::LT, inf.d_rangeTerm, absDen));
        bool inRange = true;
        for (const Node& r : {lower, upper})
        {
          if (d_model.computeAbstractModelValue(r) != d_true)
          {
            inRange = false;
          }
          else if (!r.isConst())
          {
            premises.push_back(r);
            allAsserted = allAsserted && isAsserted(r);
          }
        }
        if (!inRange)
        {
          continue;
        }
      }
      Node lem = premises.empty()
                     ? inf.d_conclusion
                     : nm->mkNode(
                           Kind::IMPLIES, nm->mkAnd(premises), inf.d_conclusion);
      if (d_model.computeAbstractModelValue(inf.d_conclusion) != d_true)
      {
        Trace("nl-div-ideal") << "...lemma: " << lem << std::endl;
        d_im.addPendingLemma(lem, InferenceId::ARITH_NL_DIV_IDEAL);
      }
      else if (allAsserted)
      {
        Trace("nl-div-ideal")
            << "...derived fact: " << inf.d_conclusion << std::endl;
        d_derived.emplace_back(inf.d_conclusion, lem);
      }
    }
  }
#endif /* CVC5_USE_COCOA */
}

std::vector<DivIdealSolver::Inference> DivIdealSolver::computeInferences(
    const Node& den,
    const std::vector<DivTerm>& divs,
    const std::vector<Node>& eqs)
{
  std::vector<Inference> infs;
#ifdef CVC5_USE_COCOA
  NodeManager* nm = nodeManager();
  Trace("nl-div-ideal") << "DivIdealSolver: divisor " << den << " with "
                        << divs.size() << " divisions and " << eqs.size()
                        << " equalities" << std::endl;
  std::vector<Node> leaves;
  collectLeaves(den, leaves);
  for (const DivTerm& dt : divs)
  {
    collectLeaves(dt.d_num, leaves);
  }
  for (const Node& eq : eqs)
  {
    collectLeaves(eq[0], leaves);
    collectLeaves(eq[1], leaves);
  }
  try
  {
    PolyEncoder enc(nm, leaves);
    // The generators of the ideal: the divisor comes first, followed by the
    // polynomials that are zero. We remember the equality each one is from.
    std::vector<CoCoA::RingElem> gens;
    std::vector<Node> genEqs;
    gens.push_back(enc.encode(den));
    genEqs.push_back(Node::null());
    if (CoCoA::IsZero(gens[0]) || CoCoA::IsConstant(gens[0]))
    {
      return infs;
    }
    for (const Node& eq : eqs)
    {
      CoCoA::RingElem p = enc.encode(eq[0]) - enc.encode(eq[1]);
      if (!CoCoA::IsZero(p))
      {
        gens.push_back(p);
        genEqs.push_back(eq);
      }
    }
    CofactorBasis basis(enc.getRing(), gens);
    if (!basis.compute(s_basisTimeLimit) || basis.isTrivial())
    {
      Trace("nl-div-ideal") << "...no basis" << std::endl;
      return infs;
    }
    // reduce the numerators: nums[i] = nfs[i] + sum_j cofs[i][j] * gens[j]
    size_t ndivs = divs.size();
    std::vector<CoCoA::RingElem> nums;
    std::vector<CoCoA::RingElem> nfs(ndivs, CoCoA::zero(enc.getRing()));
    std::vector<std::vector<CoCoA::RingElem>> cofs(ndivs);
    std::vector<Node> rems;
    // the first division with a given normal form
    std::map<std::string, size_t> reps;
    for (size_t i = 0; i < ndivs; ++i)
    {
      const DivTerm& dt = divs[i];
      nums.push_back(enc.encode(dt.d_num));
      rems.push_back(nm->mkNode(
          Kind::SUB, dt.d_num, nm->mkNode(Kind::MULT, den, dt.d_quot)));
      if (!basis.reduce(nums[i], nfs[i], cofs[i]))
      {
        Assert(false) << "Failed to reduce " << dt.d_num;
        continue;
      }
      Trace("nl-div-ideal") << "  " << dt.d_num << " is congruent to "
                            << enc.decode(nfs[i]) << std::endl;
      if (TraceIsOn("nl-div-ideal"))
      {
        // the certificate: numerator = normal form + sum_j cofactor_j * gen_j
        for (size_t j = 0, ngens = gens.size(); j < ngens; ++j)
        {
          if (CoCoA::IsZero(cofs[i][j]))
          {
            continue;
          }
          Node c = enc.decode(cofs[i][j]);
          Trace("nl-div-ideal") << "    cofactor ";
          if (c.isNull())
          {
            Trace("nl-div-ideal") << toString(cofs[i][j]);
          }
          else
          {
            Trace("nl-div-ideal") << c;
          }
          Trace("nl-div-ideal") << " of ";
          if (j == 0)
          {
            Trace("nl-div-ideal") << "the divisor " << den;
          }
          else
          {
            Trace("nl-div-ideal") << genEqs[j];
          }
          Trace("nl-div-ideal") << std::endl;
        }
      }
      size_t rep = reps.emplace(toString(nfs[i]), i).first->second;
      Inference inf;
      if (rep != i)
      {
        // The numerators differ by k * den, up to polynomials that are zero.
        // If k is integral, the remainders are the same.
        if (!hasIntegerCoefficients(cofs[i][0] - cofs[rep][0]))
        {
          continue;
        }
        for (size_t j = 1, ngens = gens.size(); j < ngens; ++j)
        {
          if (cofs[i][j] != cofs[rep][j])
          {
            inf.d_eqs.push_back(genEqs[j]);
          }
        }
        inf.d_conclusion = rewrite(rems[i].eqNode(rems[rep]));
      }
      else
      {
        // The numerator differs from its normal form t by k * den, up to
        // polynomials that are zero. If k is integral and t is in the range of
        // the remainder, then t is the remainder.
        if (nfs[i] == nums[i] || !hasIntegerCoefficients(cofs[i][0]))
        {
          continue;
        }
        Node t = enc.decode(nfs[i]);
        if (t.isNull())
        {
          continue;
        }
        t = rewrite(t);
        for (size_t j = 1, ngens = gens.size(); j < ngens; ++j)
        {
          if (!CoCoA::IsZero(cofs[i][j]))
          {
            inf.d_eqs.push_back(genEqs[j]);
          }
        }
        if (t != d_zero)
        {
          inf.d_rangeTerm = t;
        }
        inf.d_conclusion = rewrite(rems[i].eqNode(t));
      }
      if (inf.d_conclusion == d_true)
      {
        continue;
      }
      Trace("nl-div-ideal") << "  infer " << inf.d_conclusion << " from "
                            << inf.d_eqs << ", range term " << inf.d_rangeTerm
                            << std::endl;
      infs.push_back(inf);
    }
  }
  catch (const CoCoA::ErrorInfo& e)
  {
    Trace("nl-div-ideal") << "...CoCoA says: " << CoCoA::message(e)
                          << std::endl;
    infs.clear();
  }
#endif /* CVC5_USE_COCOA */
  return infs;
}

void DivIdealSolver::getDerivedFacts(std::vector<Node>& assertions) const
{
  std::unordered_set<Node> has(assertions.begin(), assertions.end());
  for (const std::pair<Node, Node>& d : d_derived)
  {
    if (has.insert(d.first).second)
    {
      assertions.push_back(d.first);
    }
  }
}

void DivIdealSolver::sendDerivedLemmas()
{
  if (d_derivedSent)
  {
    return;
  }
  d_derivedSent = true;
  for (const std::pair<Node, Node>& d : d_derived)
  {
    d_im.addPendingLemma(d.second, InferenceId::ARITH_NL_DIV_IDEAL);
  }
}

}  // namespace nl
}  // namespace arith
}  // namespace theory
}  // namespace cvc5::internal
