; REQUIRES: cocoa
; COMMAND-LINE: --nl-div-ideal
; EXPECT: unsat
; Multiplication distributes over addition modulo p. The numerators of the two
; outermost applications of mod differ by a multiple of p.
(set-logic QF_NIA)
(set-info :status unsat)
(declare-const p Int)
(declare-const x Int)
(declare-const y Int)
(declare-const z Int)
(assert (> p 0))
(assert (not (= (mod (* x (mod (+ y z) p)) p)
                (mod (+ (mod (* x y) p) (mod (* x z) p)) p))))
(check-sat)
