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

;;; This file describes a data structure useful for describing sparse matrices.
;;; It is geared towards saving space, and is rather handy for performing SOR
;;; on large matrices or assembling finite-element equations.

(declare (usual-integrations))


;;; Basic data structure and associated methods:

(define (make-sparse-matrix nrows ncols)
  (list (make-vector nrows '()) nrows ncols))

(define sparse-matrix-row-count cadr)
(define sparse-matrix-column-count caddr)
(define sparse-matrix-size cdr)

(define (sparse-matrix-ref sm i j)
  (let ((m (sparse-matrix-row-count sm))
	(n (sparse-matrix-column-count sm)))

    (if (or (< i 0) (>= i m) (< j 0) (>= j n))
	(error "Matrix access out of bound. -- SPARSE-MATRIX-REF"))

    (let ((result (assq j (vector-ref (car sm) i))))
      (if result
	  (cadr result)
	  0))))

(define (sparse-matrix-set! sparse i j val)
  (let ((m (sparse-matrix-row-count sparse))
	(n (sparse-matrix-column-count sparse))
	(sm (car sparse)))

    (if (or (< i 0) (>= i m) (< j 0) (>= j n))
	(error "Matrix access out of bound. -- SPARSE-MATRIX-SET!"))

    (let ((result (assq j (vector-ref sm i))))
      (if (zero? val)
	  (if result
	      (vector-set! sm i (all-but (vector-ref sm i) result)))
	  (if result
	      (set-cdr! result (list val))
	      (vector-set! sm i (cons (list j val) (vector-ref sm i))))))))
	  

(define (sparse-matrix-get-row sm i)
  (let ((m (sparse-matrix-row-count sm)))

    (if (or (< i 0) (>= i m))
	(error "Matrix access out of bound. -- SPARSE-MATRIX-GET-ROW"))

    (vector-ref (car sm) i)))

(define sparse-matrix-get-rows car)

(define (sparse-matrix-get-column sm j)
  (let ((m (sparse-matrix-row-count sm)))
    (let next-row ((i 0) (result '()))
      (if (< i m)
	  (let next-term ((row (sparse-matrix-get-row sm i)))
	    (if (null? row)
		(next-row (+ i 1) result)
		(if (= (caar row) j)
		    (next-row (+ i 1) (cons (list i (cadar row)) result))
		    (next-term (cdr row)))))
	  result))))

(define (sparse-matrix-get-columns sm)
  (let ((m (sparse-matrix-row-count sm))
	(v (make-vector (sparse-matrix-column-count sm) '())))
    (let next-row ((i 0))
      (if (< i m)
	  (let next-term ((row (sparse-matrix-get-row sm i)))
	    (if (null? row)
		(next-row (+ i 1))
		(let ((j (caar row)))
		  (vector-set! v j (cons (list i (cadar row))
					 (vector-ref v j)))
		  (next-term (cdr row)))))
	  v))))

;;; A predicate that can come in handy:

(define (sparse-matrix? sm)
  (call-with-current-continuation
   (lambda (exit)
     (if (and (list? sm) (= (length sm) 3))
	 (let ((m (cadr sm))
	       (n (caddr sm))
	       (sm (car sm)))
	   (if (and (integer? m) (integer? n) (> m 0) (> n 0) (vector? sm))
	       (do ((i 0 (+ i 1)))
		   ((>= i m))
		 (let ((l (vector-ref sm i)))
		   (if (or (not (list? l))
			   (memq #f (map list? l)))
		       (exit #f))))
	       #f)
	   #t)
	 #f))))


;;; Converters:

(define (sparse->matrix sm)
  (let* ((m (sparse-matrix-row-count sm))
	 (n (sparse-matrix-column-count sm))
	 (matrix (make-matrix m n)))
    (do ((i 0 (+ i 1)))
	((>= i m) matrix)
      (for-each
       (lambda (pair)
	 (matrix-set! matrix i (car pair) (cadr pair)))
       (sparse-matrix-get-row sm i)))))

(define (matrix->sparse matrix)
  (let* ((m (matrix-row-count matrix))
	 (n (matrix-column-count matrix))
	 (sm (make-sparse-matrix m n)))
    (do ((i 0 (+ i 1)))
	((>= i m) sm)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(sparse-matrix-set! sm i j (matrix-ref matrix i j))))))


;;; Very useful routine:

(define (print-sparse-matrix matrix . argl)
  (if (null? argl)
      (set! argl (list (current-output-port))))

  (if (not (and (null? (cdr argl)) (output-port? (car argl))))
      (error "Invalid argument(s) -- PRINT-MATRIX"))

  (let ((port (car argl)))
    (newline port)
    (let ((m (sparse-matrix-row-count matrix))
          (n (sparse-matrix-column-count matrix)))
      (do ((i 0 (+ i 1)))
	  ((>= i m))
	(display (sparse-matrix-ref matrix i 0) port)

	(do ((j 1 (+ j 1)))
	    ((>= j n))
	  (display #\tab port)
	  (display (sparse-matrix-ref matrix i j) port))

	(newline port)))))


;;; Prepare for least squares on sparse matrices:

(define (sparse-normal-equations mat)
  (let* ((m (sparse-matrix-row-count mat))
	 (n+1 (sparse-matrix-column-count mat))
	 (n (- n+1 1))
	 (out (make-sparse-matrix n n+1))
	 (columns (sparse-matrix-get-columns mat)))

    ;; Compute the normal equations:

    (do ((i 0 (+ i 1)))
	((>= i n) out)
      (let ((ith-column (vector-ref columns i)))

	;; First, the diagonal:

	(sparse-matrix-set! out i i (sparse-dot ith-column ith-column))

	;; Next, the off-diagonal terms:

	(do ((j (+ i 1) (+ j 1)))
	    ((>= j n))
	  (let ((val (sparse-dot ith-column (vector-ref columns j))))
	    (sparse-matrix-set! out i j val)
	    (sparse-matrix-set! out j i val)))

	;; Finally, the RHS:

	(sparse-matrix-set! out i n (sparse-dot ith-column
						(vector-ref columns n)))))))

(define (sparse-dot a b)
  (let a-loop ((a a) (sum 0))
    (if (null? a)
	sum
	(let ((id (caar a)))
	  (let b-loop ((b b))
	    (if (null? b)
		(a-loop (cdr a) sum)
		(if (= id (caar b))
		    (a-loop (cdr a) (+ sum (* (cadar a) (cadar b))))
		    (b-loop (cdr b)))))))))
