;;; This file defines some auxiliary data structures for the PDE code.

(declare (usual-integrations))


;;; Linear equations:

(define (make-equation node constant terms)
  (vector constant terms node))

(define (equation:get-constant equation)
  (vector-ref equation 0))

(define (equation:get-terms equation)
  (vector-ref equation 1))

(define (equation:get-node equation)
  (vector-ref equation 2))

(define (equation:get-id equation)
  (node:get-id (equation:get-node equation)))

(define (null-equation? equation)
  (let loop ((terms (equation:get-terms equation)))
    (if (null? terms)
	(let ((node (equation:get-node equation)))
	  (list (node:get-point node) (node:boundary? node)))
	(if (almost-zero? (term:get-coeff (car terms)))
	    (loop (cdr terms))
	    #f))))


;;; Terms in linear equations:

(define (make-term node value)
  (vector node value))

(define (term:get-node term)
  (vector-ref term 0))

(define (term:get-id term)
  (node:get-id (term:get-node term)))

(define (term:get-coeff term)
  (vector-ref term 1))
