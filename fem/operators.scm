;;; This file contains some ad-hoc definitions of simple differential
;;; operators, such as the two-dimensional Laplacian or the 1+1-dimensional
;;; d'Alembertian.  It also defines some important constructors and operations
;;; on differential operators that compute the adjoint and encapsulate
;;; integration-by-parts.  (See ELEMENT-MAKER in fem.scm.)

;;; opalg.scm contains the beginnings of a much more abstract (and complete)
;;; approach.

;;; This file uses various procedures from basis.scm.

(declare (usual-integrations))


;;; Simple constructor and methods for an operator structure:

(define (make-operator left-op right-op combine)
  (vector left-op right-op combine))

(define (operator:get-left-op operator)
  (vector-ref operator 0))

(define (operator:get-right-op operator)
  (vector-ref operator 1))

(define (operator:get-combine operator)
  (vector-ref operator 2))

(define (operator:get-local-form op)
  (let ((combine (operator:get-combine op))
	(left-op (operator:get-left-op op))
	(right-op (operator:get-right-op op)))
    (lambda (f g)
      (combine (left-op f) (right-op g)))))
