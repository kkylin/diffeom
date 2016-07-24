#|

This file is part of DIFFEOM ("DIFFerential Equations
On Manifolds"), a system for solving differential
equations on manifolds.

Copyright (C) 2016 by Kevin K Lin
<klin@math.arizona.edu>

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
