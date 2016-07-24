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

;;; This file defines basis function constructors and integrator codes.  What
;;; it provides is a way to handle functions over elements.  Note that the code
;;; in fem.scm operates independent of the representation we use here.

;;; There is an abuse of terms here.  By "basis function" we mean basis
;;; functions and their linear combinations.  Thus, functions over elements
;;; represented by sums of basis functions are also considered basis functions.

;;; This code is specific to polynomial basis functions of two variables.  The
;;; polynomial code we have still needs lots of work, so we won't use it here.
;;; Most of the procedures operate directly on vector representations of basis
;;; functions.

(declare (usual-integrations))


;;; Basic constructor:

(define (vector->poly v)
  (package-basis-function-methods
   '2d-poly-basis-function
   v
   (poly->function v)
   (make-2d-poly-adder v)
   (make-2d-poly-subtractor v)
   (make-2d-poly-multiplier v)
   (make-2d-poly-scalar-multiplier v)))

(define (make-polynomial-basis-function nodes center)
  (let* ((n (length nodes))
	 (vals (make-vector n))
	 (points (make-vector n)))

    (let loop ((nodes nodes) (i 0))
      (if (null? nodes)
	  (vector->poly (poly:point-value->coeff vals points))
	  (let ((node (car nodes)))
	    (if (= i center)
		(vector-set! vals i 1)
		(vector-set! vals i 0))
	    (vector-set! points i (vector (node:get-x node) (node:get-y node)))
	    (loop (cdr nodes) (+ i 1)))))))

(define (make-2d-poly-adder v)
  (lambda (w)
    (vector->poly (poly:+ v (basis:get-rep w)))))

(define (make-2d-poly-subtractor v)
  (lambda (w)
    (vector->poly (poly:- v (basis:get-rep w)))))

(define (make-2d-poly-multiplier v)
  (lambda (w)
    (vector->poly (poly:* v (basis:get-rep w)))))

(define (make-2d-poly-scalar-multiplier v)
  (lambda (a)
    (vector->poly (poly:scalar* a v))))


;;; A slightly different kind of constructor:

(define (function->poly f nodes)
  (let* ((n (length nodes))
	 (vals (make-vector n))
	 (points (make-vector n)))

    (let loop ((i 0) (nodes nodes))
      (if (null? nodes)
	  (vector->poly (poly:point-value->coeff vals points))
	  (let ((node (car nodes)))
	    (vector-set! points i (node:get-coords node))
	    (vector-set! vals i (f node))
	    (loop (+ i 1) (cdr nodes)))))))


;;; And its inverse:

(define (poly->function f)
  (lambda (x)
    (poly:evaluate f x)))

(define (poly:evaluate f x)
  (vector-first (poly:coeff->point-value f (vector x))))


;;; Operations on basis functions:

(define (poly:+ v w)
  (let ((m (vector-length v))
	(n (vector-length w)))

    (let ((m (max m n))
	  (n (min m n))
	  (v (if (>= m n) v w))
	  (w (if (< m n) v w)))
      (let ((result (make-vector m)))
	(do ((i 0 (+ i 1)))
	    ((>= i n))
	  (vector-set! result i (+ (vector-ref v i) (vector-ref w i))))
	(do ((i n (+ i 1)))
	    ((>= i m) result)
	  (vector-set! result i (vector-ref v i)))))))

(define (poly:scalar* a v)
  (let* ((n (vector-length v))
	 (w (make-vector n)))
    (do ((i 0 (+ i 1)))
	((>= i n) w)
      (vector-set! w i (* a (vector-ref v i))))))

(define (poly:* p1 p2)
  (let* ((n1 (vector-length p1))
	 (n2 (vector-length p2))
	 (degree (+ (poly:degree p1) (poly:degree p2)))
	 (n (choose (+ degree 2) 2))
	 (p (make-vector n 0)))

    (do ((i 0 (+ i 1)))
	((>= i n1) p)
      (let ((powers (zig-zag i))
	    (coeff (vector-ref p1 i)))
	(do ((j 0 (+ j 1)))
	    ((>= j n2))
	  (let ((k (apply inverse-zig-zag (map + powers (zig-zag j)))))
	    (vector-set! p k (+ (vector-ref p k)
				(* coeff (vector-ref p2 j))))))))))

(define (poly:- v w)
  (poly:+ v (poly:scalar* -1 w)))

(define (poly-basis:partial k)
  (if (not (or (= k 0) (= k 1)))
      (error "Only (partial 0) and (partial 1) exist (for now)! -- PARTIAL")
      (let ((select (if (= k 0) car cadr)))
	(lambda (v)
	  (let* ((v (basis:get-rep v))
		 (n (vector-length v))
		 (w (make-vector n 0)))
	    (do ((i 0 (+ i 1)))
		((>= i n) (vector->poly w))
	      (let ((powers (zig-zag i)))
		(if (> (select powers) 0)
		    (vector-set!
		     w (- i (apply + powers) k)
		     (* (select powers) (vector-ref v i)))))))))))


;;; Some useful definitions in two dimensions:

(define d/dx (poly-basis:partial 0))
(define d/dy (poly-basis:partial 1))
(define d/dt d/dy)


;;; Converting between point-value and coefficient representations; is there a
;;; higher-dimensional analog of the FFT trick?  Point-value representation is
;;; great for everything *except* differentiation...

(define (poly:coeff->point-value v sample-points)
  (let* ((n (vector-length sample-points))
	 (w (make-vector n))
	 (m (vector-length v)))

    (do ((i 0 (+ i 1)))
	((>= i n) w)

      (let* ((coords (vector-ref sample-points i))
	     (x (vector-ref coords 0))
	     (y (vector-ref coords 1)))

	(let loop ((j 0) (sum 0.))
	  (if (< j m)
	      (let ((powers (zig-zag j)))
		(loop (+ j 1) (+ sum (* (vector-ref v j)
					(expt x (car powers))
					(expt y (cadr powers))))))
	      (vector-set! w i sum)))))))

(define (poly:point-value->coeff w sample-points)
  (let* ((n (vector-length sample-points))
	 (A (make-matrix n (+ n 1))))

    (let next-row ((i 0))
      (if (< i n)
	  (let* ((coords (vector-ref sample-points i))
		 (x (vector-ref coords 0))
		 (y (vector-ref coords 1)))
	    (let next-column ((j 0))
	      (if (< j n)
		  (let* ((powers (zig-zag j))
			 (p (car powers))
			 (q (cadr powers)))
		    (matrix-set! A i j (* (expt x p) (expt y q)))
		    (next-column (+ j 1)))
		  (begin
		    (matrix-set! A i n (vector-ref w i))
		    (next-row (+ i 1))))))
	  (lu-solve A 'no-copy)))))

(define (poly:slow-make-sample-points n)
  (if (> n 0)
      (let* ((delta (exact->inexact (/ n)))
	     (n (choose (+ n 2) 2))
	     (v (make-vector n)))

	(do ((i 0 (+ i 1)))
	    ((>= i n) v)
	  (let ((powers (zig-zag i)))
	    (vector-set! v i (vector (* (car powers) delta)
				     (* (cadr powers) delta))))))
      (vector (vector 0 0))))

(define poly:make-sample-points
  (simple-memoize poly:slow-make-sample-points 10))


;;; Some useful operations on basis functions:

(define (poly:degree p)
  (apply + (zig-zag (- (vector-length p) 1))))

(define (poly:coeff->expr v)
  (let* ((v (basis:get-rep v))
	 (n (vector-length v)))
    (let loop ((expr '()) (i (- n 1)))
      (if (>= i 0)
	  (let ((powers (zig-zag i)))
	    (loop (cons `((x ,(car powers) y ,(cadr powers)) ,(vector-ref v i))
			expr)
		  (- i 1)))
	  expr))))


;;; The truly messy stuff: Integrals!  This needs to run a lot faster.  What
;;; about doing away with the coordinate transformations?

(define (make-triangular-integrator vertex-nodes)

  ;; We assume that there are three vertex nodes, and that the triangle they
  ;; form is the boundary of the element:

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
	   (jacobian (abs (det A))))

      (define (integrate f . rest)
	(let* ((f (apply basis:* (cons f rest)))
	       (degree (poly:degree f))
	       (reference (poly:make-sample-points degree))
	       (n (choose (+ degree 2) 2))
	       (real (make-vector n)))

	  (do ((i 0 (+ i 1)))
	      ((>= i n))
	    (vector-set! real i
			 (apply-affine-transformation
			  A b (vector-ref reference i))))

	  (* jacobian
	     (inner-product
	      (poly:point-value->coeff
	       (poly:coeff->point-value (basis:get-rep f) real) reference)
	      (make-reference-integrals degree)))))

      integrate)))

(define (slow-make-reference-integrals degree)
  (let* ((n (choose (+ degree 2) 2))
	 (integrals (make-vector n)))
    (do ((i 0 (+ i 1)))
	((>= i n) integrals)
      (vector-set! integrals i (apply reference-integral (zig-zag i))))))

(define make-reference-integrals
  (simple-memoize slow-make-reference-integrals 10))

(define (reference-integral m n)
  (let ((n+1 (+ n 1)))
    (let loop ((i 0) (sum 0.) (-1^i 1))
      (if (<= i n+1)
	  (loop (+ i 1)
		(+ sum (* (choose n+1 i) (/ (exact->inexact (+ i m 1))) -1^i))
		(* -1^i -1))
	  (/ sum n+1)))))


;;; Zig-zag across the two-dimensional square lattice:

(define (slow-zig-zag n)
  (let loop ((n n) (p 0) (q 0))
    (if (> n 0)
	(if (zero? p)
	    (loop (- n 1) (+ q 1) 0)
	    (loop (- n 1) (- p 1) (+ q 1)))
	(list p q))))

(define zig-zag (simple-memoize slow-zig-zag 20))

(define (inverse-zig-zag m n)
  (if (and (zero? m) (zero? n))
      0
      (+ (choose (+ m n 1) 2) n)))


;;; We need to define the gradient to help define the laplacian:

(define (poly-gradient f)
  (vector (d/dx f) (d/dy f)))
