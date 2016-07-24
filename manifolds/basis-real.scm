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

;;; This file defines a class of basis functions based on real functions (not
;;; just polynomials).  At bottom, we still use polynomial basis functions, but
;;; these guys don't get truncated under changes of coordinates.  The drawback
;;; is that we actually need to use numerical integration, which is less
;;; accurate and a lot slower.

(declare (usual-integrations))


;;; Constructor:

(define (pde:make-real-basis-function nodes i)
  (let ((f (make-real-basis-function nodes i)))
    (node:add-basis-function (list-ref nodes i) f)
    f))


;;; Differential operators:

(define (operator:pull-back-real-op left-op right-op combine)
  (lambda (chart nodes)

    ;; Take a basis function, pull back, apply operator, and then send it back
    ;; to the chart.

    (let ((coord-map (chart:get-coord-map chart))
	  (inverse-map (chart:get-inverse-map chart)))

      (lambda (f g)
	(let ((f (proc->real (compose (basis:get-rep f) coord-map)))
	      (g (proc->real (compose (basis:get-rep g) coord-map))))
	  (let ((h (combine (left-op f) (right-op g))))
	    (proc->real (compose (basis:get-rep h) inverse-map))))))))


;;; This integrates with the wrong *measure*, though.  What is required is to
;;; take into account the Jacobian of the coordinate charts.  (See
;;; basis-poly.scm, where this is done very approximately.)  Of course, this
;;; particular approach assumes that the manifold is imbedded in some Euclidean
;;; space, which can be restrictive for some applications.  To fix this, we
;;; probably need some computational representation of Riemann metrics or
;;; differential forms on manifolds.

(define (trapezoidal-integrator-maker-on-charts count)
  (let ((make-integrator (trapezoidal-integrator-maker count)))
    (lambda (nodes)
      (let* ((integrate (make-integrator nodes))
	     (g (chart:get-inverse-map (node:get-chart (car nodes))))
	     (dg (proc->real (function->jacobian g))))
	(lambda (f . rest)
	  (apply integrate `(,dg ,f ,@rest)))))))


;;; Given that F is an imbedding of a subset of the plane in a
;;; higher-dimensional Euclidean space, how do we (efficiently) compute its
;;; Jacobian?

;;; This guy will currently work only if F goes from R^2 to R^2.  Needs fixing.

(define (function->jacobian f)
  (lambda (x)
    (let ((M (jacobian-matrix f x)))
      (abs (det M)))))

(define (jacobian-matrix f x)
  (let ((n (vector-length x))
	(m (vector-length (f x))))
    (let ((mat (make-matrix m n)))
      (do ((j (- n 1) (- j 1)))
	  ((< j 0) mat)
	(let ((v (((pdiff j) f) x)))
	  (do ((i 0 (+ i 1)))
	      ((>= i m))
	    (matrix-set! mat i j (vector-ref v i))))))))
