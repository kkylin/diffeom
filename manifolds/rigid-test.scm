(load "load-ode")

(define p bad-init)
(define m (tangent:get-anchor p))
(define p0 (make-tangent c0 m (chart:push-forward p c0)))
(define p1 (make-tangent c1 m (chart:push-forward p c1)))
(define p2 (make-tangent c2 m (chart:push-forward p c2)))
(define p3 (make-tangent c3 m (chart:push-forward p c3)))

(define atlas (manifold:get-finite-atlas so3))

(define c0 (list-ref atlas 0))
(define c1 (list-ref atlas 1))
(define c2 (list-ref atlas 2))
(define c3 (list-ref atlas 3))

(define tc0 (make-tangent-chart c0))
(define tc1 (make-tangent-chart c1))
(define tc2 (make-tangent-chart c2))
(define tc3 (make-tangent-chart c3))

(define ttc0 (make-tangent-chart tc0))
(define ttc1 (make-tangent-chart tc1))
(define ttc2 (make-tangent-chart tc2))
(define ttc3 (make-tangent-chart tc3))

(define (evaluate-rigid-field chart p)
  (let ((x (chart:point->coords p chart)))
    (vector-append x ((make-rigid-body-field chart) x))))

(define v0 (evaluate-rigid-field tc0 p0))
(define v1 (evaluate-rigid-field tc1 p1))
(define v2 (evaluate-rigid-field tc2 p2))
(define v3 (evaluate-rigid-field tc3 p3))

(for-each
 (lambda (v chart)
   (write-line
    (vector:distance (chart:point->coords (chart:coords->point v chart) ttc0)
		     v0)))
 (list v0 v1 v2 v3)
 (list ttc0 ttc1 ttc2 ttc3))
