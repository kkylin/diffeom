;;; This file loads the appropriate definitions for the numerical
;;; differentiation of real functions.

(declare (usual-integrations))
(load "manifolds/linear")
(load "manifolds/lshared")
(load "manifolds/richardson")


;;; Stolen from manifolds/misc.scm:

(define (pdiff i)
  (lambda (f)
    (let ((df (diff f)))
      (lambda (x)
	(let ((v (make-vector (vector-length x) 0)))
	  (vector-set! v i 1)
	  ((df x) v))))))
