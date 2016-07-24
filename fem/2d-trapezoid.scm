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

;;; This file defines a simple numerical integrator over triangular subregions
;;; of the plane.  It uses the trapezoidal rule because that's the easiest
;;; thing to implement, and I'd just like to see if it improves FEM on
;;; manifolds.

(declare (usual-integrations))
;;; Heh heh...  ^^^^^^^^^^^

;;; We can follow the same idea as in 2d-basis.scm: Map the triangular region
;;; to a standard isoceles triangle by an affine transformation and apply the
;;; trapezoidal rule.  Note that its assumptions about basis functions are
;;; different.

(define (trapezoidal-integrator-maker count)

  ;; This parameter determines how many thingamajigs to use for integration.
  ;; It should really scale depending on the element, but for simplicity let's
  ;; keep it a constant (for now).

  (let* ((count-1 (- count 1))
	 (h (/ 1. count))
	 (area (/ (square h) 2)))

    (lambda (vertex-nodes)

      ;; We assume that there are three vertex nodes, and that the triangle
      ;; they form is the boundary of the element:

      (if (not (= (length vertex-nodes) 3))
	  (error (string-append "Error: Elements must have three vertex nodes."
				" -- MAKE-TRIANGULAR-INTEGRATOR")))

      (let ((p1 (car vertex-nodes))
	    (p2 (cadr vertex-nodes))
	    (p3 (caddr vertex-nodes)))

	;; Find the absolute value of the Jacobian of the affine transformation
	;; mapping the reference triangle {(0,0),(1,0),(0,1)} to this triangle.

	(let* ((A (list->matrix
		   2 2
		   (list
		    (- (node:get-x p2) (node:get-x p1))
		    (- (node:get-x p3) (node:get-x p1))
		    (- (node:get-y p2) (node:get-y p1))
		    (- (node:get-y p3) (node:get-y p1)))))
	       (b (node:get-coords p1))
	       (jacobian (abs (det A)))
	       (aff (lambda (x) (apply-affine-transformation A b x))))

	  (lambda (f . rest)
	    (let ((flist (map basis-function->function (cons f rest)))
		  (sum 0))
	      (do ((i 0 (+ i 1)))
		  ((>= i count))
		(let ((x1 (* i h))
		      (x2 (* (+ i 1) h)))
		  (do ((j (- count i 1) (- j 1)))
		      ((< j 0))
		    (let* ((y1 (* j h))
			   (y2 (* (+ j 1) h))

			   (ll (aff (vector x1 y1)))
			   (lr (aff (vector x2 y1)))
			   (ul (aff (vector x1 y2))))

		      (set! sum (+ sum (trapezoidal-average
					flist (list ll lr ul))))

		      (if (< (+ i j) count-1)
			  (let ((ur (aff (vector x2 y2))))
			    (set! sum (+ sum (trapezoidal-average
					      flist (list ul ur lr))))))))))
	      (* sum area jacobian))))))))

(define (trapezoidal-average flist plist)
  (let next-point ((plist plist) (sum 0) (count 0))
    (if (null? plist)
	(/ sum count)
	(let ((p (car plist)))
	  (let next-function ((flist flist) (prod 1))
	    (if (null? flist)
		(next-point (cdr plist) (+ sum prod) (+ count 1))
		(next-function (cdr flist) (* prod ((car flist) p)))))))))
