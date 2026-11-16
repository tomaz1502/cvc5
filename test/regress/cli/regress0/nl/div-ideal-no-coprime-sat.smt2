; REQUIRES: cocoa
; COMMAND-LINE: --nl-div-ideal
; EXPECT: sat
; The cancellation law for congruences does not hold when a and n are not
; coprime, e.g. for a = n = 2, x = 0, y = 1.
(set-logic QF_NIA)
(set-info :status sat)
(declare-const a Int)
(declare-const n Int)
(declare-const x Int)
(declare-const y Int)
(declare-const d Int)
(assert (> n 0))
(assert (= (- (* a y) (* a x)) (* n d)))
(assert (not (= (mod (- y x) n) 0)))
(check-sat)
