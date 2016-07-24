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

;;; Using Scmutils' differentiation facilities to replace the numerical
;;; differentiation stuff.

(declare (usual-integrations))


;;; Ugly kludge, but sort of works.  For everything we care about, at any rate.

(define (diff f)
  (let ((df (derivative f)))
    (lambda (p)
      (let ((J (df p)))
	(if (number? J)
	    (set! J (vector (vector J)))
	    (if (not (vector? (vector-ref J 0)))
		(set! J (vector->row-matrix J))))
	(lambda (v)
	  (if (number? v)
	      (set! v (vector v)))
	  (apply-linear-transformation J v))))))


;;; Scmutils can't handle this, though:

;(define (f x)
;(if (> x 0)
;(exp (- (/ (square x))))
;0))

;(((derivative f) x) 1)


;;; Need to replace some linear algebra stuff, too:

(define (list->matrix m n l)
  (let ((v (list->vector l)))
    (generate-matrix
     m n
     (lambda (i j)
       (vector-ref v (+ (* i n) j))))))

(define vector:* vector:scalar*vector)
(define vector:dot vector:dot-product)
(define inner-product vector:dot)
(define apply-linear-transformation matrix:matrix*vector)


;;; Printing matrices:

(define matrix-row-count matrix:num-rows)
(define matrix-column-count matrix:num-cols)

(define (matrix-size M)
  (list (matrix-row-count M) (matrix-column-count M)))

(define (print-matrix matrix)
  (newline)
  (let ((m (matrix-row-count matrix))
	(n (matrix-column-count matrix)))
    (do ((i 0 (+ i 1)))
	((>= i m))
      (display (matrix-ref matrix i 0))

      (do ((j 1 (+ j 1)))
	  ((>= j n))
	(display #\tab)
	(display (matrix-ref matrix i j)))

      (newline))))

(define matrix-get-column matrix:nth-col)


;;; This can come in handy sometimes:

(define (make-matrix m n)
  (generate-matrix m n (lambda (i j) 0)))


;;; As can this:

(define det matrix:determinant)


;;; And to solve linear equations:

(define (lu-solve eqs . whatever)
  (let ((m (matrix-row-count eqs))
	(n (matrix-column-count eqs)))
    (if (= n (+ m 1))
	(let ((v (make-vector m))
	      (A (make-matrix m m)))
	  (do ((i 0 (+ i 1)))
	      ((>= i m))
	    (do ((j 0 (+ j 1)))
		((>= j m))
	      (matrix-set! A i j (matrix-ref eqs i j)))
	    (vector-set! v i (matrix-ref eqs i m)))
	  (matrix:solve-linear-system A v))
	(error "Input has incorrect dimensions. -- LU-SOLVE"))))


;;; And something else as well... (This is getting out of hand!)

(define (apply-affine-transformation A b v)
  (vector:+ (matrix:matrix*vector A v) b))
