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

(declare (usual-integrations sqrt))


;;; Some useful procedures for manipulating stereographic projections and such:

(define (vector:drop-coord v i)
  ;; Project v onto the orthogonal complement of the ith standard basis vector:
  (let* ((n (vector-length v))
	 (w (make-vector (- n 1))))
    (let loop ((j 0) (k 0))
      (if (< j n)
	  (if (= j i)
	      (loop (+ j 1) k)
	      (begin
		(vector-set! w k (vector-ref v j))
		(loop (+ j 1) (+ k 1))))
	  w))))

(define (vector:add-coord v i)
  ;; Do the opposite:
  (let* ((n (+ (vector-length v) 1))
	 (w (make-vector n)))
    (let loop ((j 0) (k 0))
      (if (< j n)
	  (if (= j i)
	      (begin
		(vector-set! w j 0.)
		(loop (+ j 1) k))
	      (begin
		(vector-set! w j (vector-ref v k))
		(loop (+ j 1) (+ k 1))))
	  w))))

(define (vector:basis dim i val)
  (let ((v (make-vector dim 0.)))
    (vector-set! v i val)
    v))


;;; Useful to return the last COUNT elements of a vector:

(define (vector-end v count)
  (vector-tail v (- (vector-length v) count)))


;;; Based on vector:dot:

(define (vector:magnitude^2 v)
  (vector:dot v v))

(define (vector:magnitude v)
  (sqrt (vector:magnitude^2 v)))

(define (vector:distance^2 v w)
  (vector:magnitude^2 (vector:- v w)))

(define (vector:distance v w)
  (sqrt (vector:distance^2 v w)))


;;; Least-squares approximation:

(define (project-onto-basis vlist vector)

  ;; Just do least-squares...

  (let* ((trans (list->vector vlist))
	 (basis (transpose trans))

	 ;; Taking this transpose, of course, is where we implicitly use the
	 ;; metric structure of the ambient Euclidean space.

	 (result (matrix:solve-linear-system (matrix:* trans basis)
					     (apply-linear-transformation
					      trans vector))))
    (apply-linear-transformation basis result)))


;;; Make an identity matrix:

(define (make-identity-matrix n)
  (let ((mat (make-matrix n n)))
    (do ((i 0 (+ i 1)))
	((>= i n) mat)
      (matrix-set! mat i i 1))))


;;; Compute the square root of the trace of a matrix multiplied by its own
;;; transpose:

(define (matrix:magnitude A)
  (let ((m (matrix-row-count A))
	(n (matrix-column-count A)))

    (let next-row ((i 0) (sum 0))
      (if (< i m)
	  (let next-col ((j 0) (sum sum))
	    (if (< j n)
		(next-col (+ j 1) (+ sum (square (matrix-ref A i j))))
		(next-row (+ i 1) sum)))
	  (sqrt sum)))))


;;; Find the maximum element of a matrix:

(define (matrix:max A)
  (let ((m (matrix-row-count A))
	(n (matrix-column-count A)))

    (let next-row ((i 0) (max 0))
      (if (< i m)
	  (let next-col ((j 0) (max max))
	    (if (< j n)
		(let ((mag (magnitude (matrix-ref A i j))))
		  (if (> mag max)
		      (next-col (+ j 1) mag)
		      (next-col (+ j 1) max)))
		(next-row (+ i 1) max)))
	  max))))


;;; Do a least-squares approximation:

(define (least-squares mat)
  (let* ((m (matrix-row-count mat))
	 (n+1 (matrix-column-count mat))
	 (n (- n+1 1))
	 (out (make-matrix n n+1)))

    (write-line '(preparing normal equations...))

    (do ((i 0 (+ i 1)))
	((>= i n))
      (do ((j 0 (+ j 1)))
	  ((> j n))
	(let loop ((k 0) (sum 0))
	  (if (< k m)
	      (loop (+ k 1) (+ sum (* (matrix-ref mat k i)
				      (matrix-ref mat k j))))
	      (matrix-set! out i j sum)))))

    (write-line '(solving normal equations using lu-decomposition...))
    (lu-solve mat)))


;;; Here's a matrix deconstructor:

(define (matrix:flatten A)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (size (* m n))
	 (v (make-vector size 0)))
    (do ((i 0 (+ i 1)))
	((>= i m))
      (do (( j 0 (+ j 1)))
	  ((>= j n))
	(vector-set! v (+ (* i n) j) (matrix-ref A i j))))
    v))

(define (vector->matrix m n v)
  (let ((A (make-matrix m n)))
    (do ((i 0 (+ i 1)))
	((>= i m))
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! A i j (vector-ref v (+ (* i n) j)))))
    A))
