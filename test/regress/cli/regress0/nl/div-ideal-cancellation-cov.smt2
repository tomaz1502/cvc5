; REQUIRES: cocoa
; REQUIRES: poly
; COMMAND-LINE: --nl-div-ideal --nl-cov --nl-ext=none
; EXPECT: unsat
; As div-ideal-cancellation.smt2, where the coverings solver is the only other
; non-linear solver.
(set-logic QF_NIA)
(set-info :status unsat)
(declare-const a Int)
(declare-const n Int)
(declare-const x Int)
(declare-const y Int)
(declare-const d Int)
(declare-const u Int)
(declare-const v Int)
(assert (> n 0))
(assert (= (- (* a y) (* a x)) (* n d)))
(assert (= (+ (* a u) (* n v)) 1))
(assert (not (= (mod (- y x) n) 0)))
(check-sat)
