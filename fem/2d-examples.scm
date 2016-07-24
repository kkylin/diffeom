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

;;; This file defines some examples of PDEs over planar regions, particularly
;;; the unit square.  This file goes with 2d-domains.scm and 2d-basis.scm.


;;; Some methods for constructing elements:

;;; For Laplace's equation:

(define make-laplacian-element
  (element-maker laplacian
		 make-triangular-integrator
		 make-polynomial-basis-function))


;;; For the linear wave equation:

(define (make-wave-element-with-coeff c)
  (element-maker (make-wave-operator c)
		 make-triangular-integrator
		 make-polynomial-basis-function))

(define *wave-constant* 1.00001)

(define make-wave-element (make-wave-element-with-coeff *wave-constant*))


;;; For characteristic bending:

(define make-bent-element
  (element-maker (make-bent-operator *wave-constant* .5 1.)
		 make-triangular-integrator
		 make-polynomial-basis-function))

;;; For testing real functions:

(define make-real-laplacian-element
  (element-maker real-laplacian
		 (trapezoidal-integrator-maker 16)
		 make-real-basis-function))


;;; Constructors for various test cases:

;;; Estimate the solution in the center of a disk, using N nodes on the
;;; boundary.

(define make-circular-domain
  (domain-maker make-circular-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		delaunay-triangulation
		make-laplacian-element
		(predicate->make-boundary circular-boundary?)))


;;; Make a square consisting of MxN interior nodes.  M is the number of nodes
;;; along the x-axis, N is the number of nodes along the y-axis.

(define make-square-domain
  (domain-maker make-square-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		square-domain-triangulation
		make-laplacian-element
		(predicate->make-boundary dirichlet-boundary?)))

(define make-quadratic-domain
  (domain-maker make-square-domain-vertices
		make-midpoint-node
		make-no-interior-nodes
		square-domain-triangulation
		make-laplacian-element
		(predicate->make-boundary dirichlet-boundary?)))

(define make-wave-domain
  (domain-maker make-square-domain-vertices
		make-midpoint-node
		make-no-interior-nodes
		square-domain-triangulation
		make-wave-element
		(predicate->make-boundary cauchy-boundary?)))

(define make-dirichlet-wave-domain
  (domain-maker make-square-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		square-domain-triangulation
		make-wave-element
		(predicate->make-boundary dirichlet-boundary?)))

(define make-bent-domain
  (domain-maker make-square-domain-vertices
		make-midpoint-node
		make-no-interior-nodes
		square-domain-triangulation
		make-bent-element
		(predicate->make-boundary cauchy-boundary?)))

(define make-hat-domain
  (domain-maker make-square-domain-vertices
		make-midpoint-node
		make-no-interior-nodes
		square-domain-triangulation
		(make-wave-element-with-coeff 4.)
		(predicate->make-boundary hat-boundary?)))

(define make-true-hat-domain
  (domain-maker make-hat-vertices
		make-midpoint-node
		make-no-interior-nodes
		delaunay-triangulation
		make-wave-element
		(predicate->make-boundary true-hat-boundary?)))


;;; For testing real functions:

(define make-real-square-domain
  (domain-maker make-square-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		square-domain-triangulation
		make-real-laplacian-element
		(predicate->make-boundary dirichlet-boundary?)))


;;; A triangular domain:

(define make-triangular-domain
  (domain-maker make-triangular-domain-vertices
		make-midpoint-node
		make-no-interior-nodes
		delaunay-triangulation
		(make-wave-element-with-coeff 1.)
		(predicate->make-boundary triangular-boundary?)))


;; Now let's try a randomized distribution, using the Delaunay triangulation:

(define make-random-square-domain
  (domain-maker make-random-square-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		delaunay-triangulation
		make-laplacian-element
		(predicate->make-boundary dirichlet-boundary?)))


;;; The boundary, instead of the unit square, is the convex hull:

(define make-random-domain
  (domain-maker make-random-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		random-domain-triangulation
		make-laplacian-element
		do-nothing-to-nodes))


;;; Same thing, but the nodes shouldn't get too close:

(define make-not-so-random-domain
  (domain-maker make-not-so-random-domain-vertices
		make-no-edge-nodes
		make-no-interior-nodes
		random-domain-triangulation
		make-laplacian-element
		do-nothing-to-nodes))


;;; Some useful boundary/initial conditions:

(define potential
  (let* ((pi 3.141592653589793)
         (sinh (lambda (x) (/ (- (exp x) (exp (- x))) 2)))
         (A (/ (sinh pi))))
    (lambda (node)
      (* A (sinh (* pi (node:get-y node))) (sin (* pi (node:get-x node)))))))

(define (wave node)
  (let ((x (node:get-x node))
	(t (node:get-y node)))
    (cos (* 2 pi (- x (* t *wave-constant*))))))

(define (standing-wave node)
  (* (sin (* 2 pi (node:get-x node)))
     (sin (* 2 pi *wave-constant* (node:get-y node)))))


;;; A harmonic function for testing the programs:

(define (test-function node)
  (let ((x (node:get-x node))
	(y (node:get-y node)))
    (* 4 (- (square (- x .5)) (square (- y .5))))))
