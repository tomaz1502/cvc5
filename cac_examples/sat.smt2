(set-info :smt-lib-version 2.6)
(set-logic QF_NRA)
(set-info :source |
These benchmarks used in the paper:

  Dejan Jovanovic and Leonardo de Moura.  Solving Non-Linear Arithmetic.
  In IJCAR 2012, published as LNCS volume 7364, pp. 339--354.

The meti-tarski benchmarks are proof obligations extracted from the
Meti-Tarski project, see:

  B. Akbarpour and L. C. Paulson. MetiTarski: An automatic theorem prover
  for real-valued special functions. Journal of Automated Reasoning,
  44(3):175-205, 2010.

Submitted by Dejan Jovanovic for SMT-LIB.


|)
(set-info :category "industrial")
(set-info :status sat)
(declare-fun X () Real)
(declare-fun Y () Real)
(declare-fun Z () Real)
(declare-fun W () Real)
(assert (and (not (= (+ (- 1) (* X X)) Y)) (and (<= 0 Y) (and (<= 0 W) (and (<= 0 X) (and (<= 0 Z) (not (<= 1 Y))))))))
(check-sat)
(exit)
