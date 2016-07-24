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

;;; This file defines the structures necessary to support Lagrangian vector
;;; fields on configuration spaces of classical mechanical systems.

(declare (usual-integrations))


;;; The Lagrangian should be a smooth map from the tangent bundle of some
;;; manifold into the real line.

;;; This is very slow, as every evaluation of the field involves a matrix
;;; inversion.  Which is why Hamiltonians are *better*, even for comuptational
;;; purposes!

(define (lagrangian->v.field L)
  (let ((TM (smooth-map:get-domain L))
	(R (smooth-map:get-range L)))
    (lambda (p)
      (let ((U
	     (if (tangent? p)
		 (make-tangent-chart (tangent:get-chart p))
		 (manifold:find-best-chart TM p))))
	(let ((f (smooth-map:make-transition
		  L U (car (manifold:get-finite-atlas R))))
	      (x (chart:point->coords p U)))
	  (let ((v (vector-tail x (/ (vector-length x) 2))))
	    (let ((E-L (euler-lagrange-in-coords f x)))
	      (let ((A (car E-L))
		    (B (cadr E-L))
		    (c (caddr E-L)))
		(let ((accel (matrix:solve-linear-system
			      A
			      (vector:+ (apply-linear-transformation B v) c))))
		  (make-tangent U p (vector-append v accel)))))))))))


;;; Derive the Euler-Lagrange equations for f at x (in coordinates) in the form
;;; A*xdotdot = B*xdot + c.

(define (euler-lagrange-in-coords f x)
  (let* ((n (/ (vector-length x) 2))
	 (A (make-matrix n n))
	 (B (make-matrix n n))
	 (c (make-vector n 0)))

    (do ((i n (+ i 1))
	 (p 0 (+ p 1)))
	((>= p n))

      ;; First, compute the hessian of f with respect to the velocity part of
      ;; the independent variable:

      (matrix-set! A p p (vector-first (((pdiff i) ((pdiff i) f)) x)))

      (do ((j (+ i 1) (+ j 1))
	   (q (+ p 1) (+ q 1)))
	  ((>= q n))
	(let ((val (vector-first (((pdiff j) ((pdiff i) f)) x))))
	  (matrix-set! A p q val)
	  (matrix-set! A q p val)))

      ;; Next, compute the rest of the terms involving the partials of the
      ;; Lagrangian with respect to the positions (note the minus sign):

      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(let ((val (- (vector-first (((pdiff j) ((pdiff i) f)) x)))))
	  (matrix-set! B p j val)))

      ;; And then there's the term due to the derivative of the Lagrangian with
      ;; respect to the position variables:

      (vector-set! c p (vector-first (((pdiff p) f) x))))

    (list A B c)))


;;; In many mechanics problems, it's natural to check conservation laws:

(define (check-vector-conservation-law quantity ref-point)
  (let ((ref (quantity ref-point)))
    (lambda (chart tangent)
      (vector:distance (quantity (tangent:get-anchor tangent)) ref))))
