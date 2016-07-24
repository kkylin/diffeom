#|

This file is part of DIFFEOM, a system for solving
differential equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<kkylin@alum.mit.edu>

This program is free software; you can redistribute
it and/or modify it under the terms of the GNU
General Public License as published by the Free
Software Foundation; either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General
Public License along with this program; if not, write
to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.

|#

;;; This file defines nodes for solving PDEs on manifolds.  These are more
;;; complicated than nodes in the old FEM program because they need to keep
;;; track of corresponding structures on manifolds.

(declare (usual-integrations))


;;; The constructor now takes a chart.  (X should be *coordinates* in the
;;; chart, not the point on the manifold.)

(define (make-node x chart)
  (let ((p (chart:coords->point x chart))
	(b? (and (boundary-chart? chart)
		 (chart:range-boundary? x chart))))
    (vector x       ;; 0. Coordinates.
	    0.      ;; 1. Value.
	    37      ;; 2. ID.
	    b?      ;; 3. Boundary?
	    '()     ;; 4. Elements.
	    '()     ;; 5. Local IDs.
	    p       ;; 6. Point on manifold.
	    chart   ;; 7. Chart.
	    '()     ;; 8. Basis functions.
	    '()     ;; 9. Extra structures.
	    #f      ;; A. Is this node constrained to another chart?
	    #f      ;; B. Is this node on the boundary of the chart?
	    #t)))   ;; C. Is this guy still alive?


;;; Same old methods:

(define (node:get-x node) (vector-ref (node:get-coords node) 0))
(define (node:get-y node) (vector-ref (node:get-coords node) 1))
(define (node:get-z node) (vector-ref (node:get-coords node) 2))
(define (node:get-coords node) (vector-ref node 0))
(define (node:get-value node) (vector-ref node 1))
(define (node:set-value! node val) (vector-set! node 1 val))
(define (node:get-id node) (vector-ref node 2))
(define (node:set-id! node id) (vector-set! node 2 id))
(define (node:boundary? node) (vector-ref node 3))
(define (node:get-elements node) (vector-ref node 4))
(define (node:get-local-ids node) (vector-ref node 5))

(define (node:add-element node element index)
  (vector-set! node 4 (cons element (vector-ref node 4)))
  (vector-set! node 5 (cons index (vector-ref node 5))))


;;; Some new additions:

(define (node:get-real-x node) (vector-ref (node:get-point node) 0))
(define (node:get-real-y node) (vector-ref (node:get-point node) 1))
(define (node:get-real-z node) (vector-ref (node:get-point node) 2))

(define (node:get-point node)
  (vector-ref node 6))

(define (node:get-chart node)
  (vector-ref node 7))

(define (node:get-basis-functions node)
  (vector-ref node 8))

(define (node:add-basis-function node basis-function)
  (vector-set! node 8 (cons basis-function (vector-ref node 8))))

(define (node:install-extra node tag datum)
  (let ((result (assq (vector-ref node 9) tag)))
    (if result
	(set-cdr! result datum)
	(vector-set! node 9 (cons (cons tag datum) (vector-ref node 9))))))

(define (node:get-extra node tag)
  (let ((result (assq (vector-ref node 9) tag)))
    (if result
	(cdr result)
	#f)))

(define (node:constrained? node)
  (if (node:get-constraint node)
      #t
      #f))

(define (node:get-constraint node)
  (vector-ref node 10))

(define (node:set-constraint! node chart)
  (vector-set! node 10 chart))

;;; Copy a node to a different chart:

(define (node:copy node chart)
  (let ((new-node (make-node (chart:point->coords (node:get-point node) chart)
			     chart)))
    (node:set-constraint! new-node node)
    new-node))

;;; Recursively find the ID of the node to which a given node is constrained:

(define (node:get-real-id node)
  (let ((constraint (node:get-constraint node)))
    (if constraint
	(node:get-real-id constraint)
	(node:get-id node))))

;;; (Some of these properties may be obsolete.)

(define (node:local-boundary? node)
  (and (vector-ref node 11)
       (not (node:boundary? node))))

(define (node:set-local-boundary! node flag)
  (vector-set! node 11 flag))

(define (node:active? node)
  (vector-ref node 12))

(define (node:kill! node)
  (vector-set! node 12 #f))

(define (node:resurrect! node)
  (vector-set! node 12 #t))
