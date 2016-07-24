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

;;; This file defines some procedures that help extend the polynomial basis
;;; functions used in fem.scm.  Note that these functions limit the accuracy of
;;; computation because a polynomial may not stay a polynomial under coordinate
;;; transformations, and yet that is how we transform these guys between
;;; coordinate systems.

(declare (usual-integrations))


;;; A little wrapper that let's us keep track of basis functions:

(define (pde:make-poly-basis-function nodes i)
  (let ((basis (make-polynomial-basis-function nodes i)))
    (node:add-basis-function (list-ref nodes i) basis)
    basis))


;;; This procedure turns an operator on M, where M is represented as an
;;; imbedded submanifold of R^n, into an operator on functions on M.  It
;;; depends heavily on the fact that it's working with polynomial interpolants.

(define (operator:pull-back-poly-op left-op right-op combine)
  (lambda (chart nodes)

    ;; Take a polynomial basis function, pull back to the canonical coordinates
    ;; of the ambient space, interpolate by polynomial, then apply operator to
    ;; form a polynomial approximating the image of the original polynomial in
    ;; the original chart under the differential operator.

    (let* ((pl (map node:get-point nodes))
	   (pv (list->vector pl))
	   (cl (map node:get-coords nodes))
	   (cv (list->vector cl)))
      (lambda (f g)
	(let* ((f1 (basis-function->function f))
	       (f2 (vector->poly (poly:point-value->coeff
				  (list->vector (map f1 cl)) pv)))
	       (g1 (basis-function->function g))
	       (g2 (vector->poly (poly:point-value->coeff
				  (list->vector (map g1 cl)) pv)))

	       (h (basis-function->function
		   (combine (left-op f2) (right-op g2)))))

	  (vector->poly (poly:point-value->coeff
			 (list->vector (map h pl)) cv)))))))

;;; Here's a problem: If we integrate purely in local coordinates, then the
;;; integral is in fact using the *wrong* measure.  In order to perform the
;;; correct integration, we need to put in the Jacobian of the coordinate
;;; function.  Since integrators are given nodes (not coordinates), and nodes
;;; have charts attached to them, this could be done (very approximately):

(define (make-triangular-chart-integrator nodes)
  (let ((triangular-integrate (make-triangular-integrator nodes))
	(jacobian (abs (/ (apply triangle-area (map node:get-point nodes))
			  (apply triangle-area (map node:get-coords nodes))))))

    (lambda (f . rest)
      (* (apply triangular-integrate (cons f rest)) jacobian))))

(define (triangle-area a b c)
  (let ((x1 (- (vector-ref b 0) (vector-ref a 0)))
	(y1 (- (vector-ref b 1) (vector-ref a 1)))
	(x2 (- (vector-ref c 0) (vector-ref a 0)))
	(y2 (- (vector-ref c 1) (vector-ref a 1))))
    (abs (* 1/2 (- (* x1 y2) (* y1 x2))))))
