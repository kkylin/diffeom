;;; This file defines nodes, which localize assembly computations.  Nodes use
;;; elements to compute source terms and integrals.  Linearity of the PDE is
;;; probabliy the only assumption here.  This also appears to be the only
;;; definition that restricts this code to two dimensions.

(declare (usual-integrations))


;;; Constructor:

(define (make-node x y . aux)
  (vector x y (if (null? aux) #f (car aux)) '() '() 37 0.))


;;; Access methods for nodes.  These are necessary for the Delaunay
;;; triangulation program, and are also useful for testing the algorithm, as
;;; the canonical coordinate functions are harmonic.

(define (node:get-x node) (vector-ref node 0))
(define (node:get-y node) (vector-ref node 1))
(define (node:get-coords node) (vector (node:get-x node) (node:get-y node)))
(define (node:get-value node) (vector-ref node 6))
(define (node:set-value! node val) (vector-set! node 6 val))
(define (node:get-id node) (vector-ref node 5))
(define (node:set-id! node id) (vector-set! node 5 id))
(define (node:boundary? node) (vector-ref node 2))
(define (node:set-boundary! node flag) (vector-set! node 2 flag))
(define (node:get-elements node) (vector-ref node 3))
(define (node:get-local-ids node) (vector-ref node 4))

(define (node:add-element node element index)
  (vector-set! node 3 (cons element (vector-ref node 3)))
  (vector-set! node 4 (cons index (vector-ref node 4))))
