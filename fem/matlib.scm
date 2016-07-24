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

;;; Some useful ideas from concrete linear algebra.  It is pretty poorly
;;; organized and some implementations need improving.

(declare (usual-integrations))


;;; Vector operations (more to come as needed):

(define (inner-product v1 v2)
  (if (and (vector? v1)
           (vector? v2)
           (= (vector-length v1) (vector-length v2)))
      (let ((n (vector-length v1)))
        (let loop ((sum 0) (i 0))
          (if (< i n)
              (loop (+ sum (* (vector-ref v1 i) (vector-ref v2 i))) (+ i 1))
              sum)))
      #f))

(define (apply-affine-transformation A b v)

  ;; Left-multiplication by a matrix:

  (let* ((n (vector-length v))
	 (w (make-vector n)))
    (do ((i 0 (+ i 1)))
	((>= i n) w)
      (let loop ((j 0) (sum 0.))
	(if (< j n)
	    (loop (+ j 1) (+ sum (* (matrix-ref A i j) (vector-ref v j))))
	    (vector-set! w i (+ sum (vector-ref b i))))))))

(define (apply-linear-transformation A v)
  (let ((m (matrix-row-count A))
	(n (matrix-column-count A))
	(p (vector-length v)))

    (if (not (= n p))
	(error "Wah!  A mistake! -- APPLY-LINEAR-TRANSFORMATION"))

    (let ((w (make-vector m)))
      (do ((i 0 (+ i 1)))
	  ((>= i m) w)
	(let loop ((j 0) (sum 0))
	  (if (< j n)
	      (loop (+ j 1) (+ sum (* (matrix-ref A i j) (vector-ref v j))))
	      (vector-set! w i sum)))))))


;;; Miscellaneous matrix operations:

(define (transpose A)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (At (make-matrix n m)))
    (do ((i 0 (+ i 1)))
	((>= i m) At)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! At j i (matrix-ref A i j))))))

(define (matrix:binary* A B)
  (let ((p (matrix-column-count A)))

    (if (not (= p (matrix-row-count B)))
	(error "Incompatible matrix sizes!"))

    (let* ((m (matrix-row-count A))
	   (n (matrix-column-count B))
	   (result (make-matrix m n)))
      (do ((i 0 (+ i 1)))
	  ((>= i m))
	(do ((j 0 (+ j 1)))
	    ((>= j n))
	  (let loop ((k 0) (sum 0))
	    (if (< k p)
		(loop (+ k 1)
		      (+ sum (* (matrix-ref A i k) (matrix-ref B k j))))
		(matrix-set! result i j sum)))))
      result)))

(define (matrix:* A . rest)
  (let loop ((A A) (l rest))
    (if (null? l)
	A
	(loop (matrix:binary* A (car l)) (cdr l)))))

(define (matrix:+ A . rest)
  (let loop ((A A) (l rest))
    (if (null? l)
	A
	(loop (matrix:binary+ A (car l)) (cdr l)))))

(define (matrix:binary+ A B)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (C (make-matrix m n)))

    (if (not (and (= m (matrix-row-count B))
		  (= n (matrix-column-count B))))
	(error "Cannot add matrices of different dimensions!"))

    (do ((i 0 (+ i 1)))
	((>= i m) C)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! C i j (+ (matrix-ref A i j) (matrix-ref B i j)))))))

(define (matrix:- A . rest)
  (if (null? rest)
      A
      (matrix:binary- A (apply matrix:+ rest))))

(define (matrix:binary- A B)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (C (make-matrix m n)))

    (if (not (and (= m (matrix-row-count B))
		  (= n (matrix-column-count B))))
	(error "Cannot subtract matrices of different dimensions!"))

    (do ((i 0 (+ i 1)))
	((>= i m) C)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! C i j (- (matrix-ref A i j) (matrix-ref B i j)))))))

(define (matrix:trace A)
  (let ((m (matrix-row-count A)))

    (if (not (= m (matrix-column-count A)))
	(error "Cannot compute the trace of a non-square matrix!"))

    (let loop ((i 0) (sum 0))
      (if (< i m)
	  (loop (+ i 1) (+ sum (matrix-ref A i i)))
	  sum))))


;;; Matrix constructors and methods:

(define (make-matrix m n)
  (let ((vector-of-rows (make-vector m)))
    (do ((i 0 (+ i 1)))
        ((>= i m))
      (vector-set! vector-of-rows i (make-vector n 0)))
    vector-of-rows))

(define (list->matrix m n l)
  (if (not (= (length l) (* m n)))
      (error "Incorrect dimensions -- LIST->MATRIX")
      (let ((A (make-matrix m n)))
        (do ((i 0 (+ i 1)))
	    ((>= i m) A)
	  (do ((j 0 (+ j 1)))
	      ((>= j n))
	    (matrix-set! A i j (car l))
	    (set! l (cdr l)))))))

(define (matrix-set! matrix i j newval)
  (vector-set! (vector-ref matrix i) j newval))

(define (matrix-ref matrix i j)
  (vector-ref (vector-ref matrix i) j))

(define matrix-row-count vector-length)

(define (matrix-column-count matrix)
  (vector-length (vector-ref matrix 0)))

(define (matrix-size matrix)
  (list (matrix-row-count matrix) (matrix-column-count matrix)))

(define matrix-dimensions matrix-size)

(define (matrix-copy A)
  (let* ((m (matrix-row-count A))
	 (n (matrix-column-count A))
	 (B (make-matrix m n)))
    (do ((i 0 (+ i 1)))
	((>= i m) B)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(matrix-set! B i j (matrix-ref A i j))))))

(define (matrix-get-row A i)
  (let* ((n (matrix-column-count A))
	 (result (make-vector n)))
    (do ((j 0 (+ j 1)))
	((>= j n) result)
      (vector-set! result j (matrix-ref A i j)))))

(define (matrix-get-column A j)
  (let* ((m (matrix-row-count A))
	 (result (make-vector m)))
    (do ((i 0 (+ i 1)))
	((>= i m) result)
      (vector-set! result i (matrix-ref A i j)))))


;;; Some predicates that might be useful:

(define (matrix? matrix)
  (and (vector? matrix)
       (> (vector-length matrix) 0)
       (vector? (vector-ref matrix 0))
       (let ((m (vector-length matrix))
             (n (vector-length (vector-ref matrix 0))))
         (let loop ((i 0))
           (if (< i m)
               (if (and (vector? (vector-ref matrix i))
                        (= (vector-length (vector-ref matrix i)) n))
                   (loop (+ i 1))
                   #f)
               #t)))))

(define (diag-dom? matrix)
  (let ((m (matrix-row-count matrix))
	(n (matrix-column-count matrix)))

    (if (> m n) (error "Matrix has more rows than columns - DIAG-DOM?"))

    (call-with-current-continuation
     (lambda (return)
       (let ((sum 0.))
         (do ((i 0 (+ i 1)))
	     ((>= i m))
	   (set! sum 0.)

	   (do ((j 0 (+ j 1)))
	       ((>= j m))
	     (set! sum (+ sum (abs (matrix-ref matrix i j)))))

	   (let ((diag (abs (matrix-ref matrix i i))))
	     (if (not (or (> diag (- sum diag))
			  (almost-zero? (- (* 2 diag) sum))))
		 (begin
		   (write-line `(- ,diag ,(- sum diag)))
		   (return #f))))))
       #t))))

(define (symmetric? matrix)
  (let ((m (matrix-row-count matrix))
	(n (matrix-column-count matrix)))
    (if (> m n)
	#f
	(call-with-current-continuation
	 (lambda (return)
	   (do ((i 1 (+ i 1)))
	       ((>= i m))
	     (do ((j 0 (+ j 1)))
		 ((>= j i))
	       (if (not (almost-zero? (- (matrix-ref matrix i j)
					 (matrix-ref matrix j i))))
		   (return #f))))
	   #t)))))


;;; LU-decomposition, done in a pretty primitive way.  It is almost directly
;;; lifted out of _Numerical Recipes_.  It can also be used to compute
;;; determinants, but it mutates its argument.  DET does not.
;;;
;;; Note that the matrix may be left in a partially modified state, since the
;;; procedure aborts on singular matrices.

(define (LU-decomp A . aux)
  (call-with-current-continuation
   (lambda (exit)
     (let ((m (matrix-row-count A))
	   (n (matrix-column-count A))
	   (det (if (and (not (null? aux))
			 (eq? (car aux) 'no-det))
		    #f
		    1)))

       ;; Compute the upper part:

       (do ((j 0 (+ j 1)))
	   ((>= j m))
        (do ((i 0 (+ i 1)))
	    ((>= i j))

	  (let ((sum (matrix-ref A i j)))
	    (do ((k 0 (+ k 1)))
		((>= k i))
	      (set! sum (- sum (* (matrix-ref A i k) (matrix-ref A k j)))))
	    (matrix-set! A i j sum)))

	;; Compute the lower portion with partial pivoting:

        (let ((pivot 0) (new-j -1))
          (do ((i j (+ i 1)))
	      ((>= i m))

	    (let ((sum (matrix-ref A i j)))
	      (do ((k 0 (+ k 1)))
		  ((>= k j))
		(set! sum (- sum (* (matrix-ref A i k) (matrix-ref A k j)))))
	      (matrix-set! A i j sum)

	      (if (>= (abs sum) (abs pivot))
		  (begin
		    (set! pivot sum)
		    (set! new-j i)))))

	  ;; Swap rows, if necessary.

          (if (> new-j j)
              (begin
		(if det (set! det (* det -1)))
		(do ((k 0 (+ k 1)))
		    ((>= k n))
		  (let ((val (matrix-ref A j k)))
		    (matrix-set! A j k (matrix-ref A new-j k))
		    (matrix-set! A new-j k val)))))

	  ;; If the matrix is singular, return 0 and leave matrix as is.

          (if (almost-zero? pivot) (exit 0))

          (do ((i (+ j 1) (+ i 1)))
              ((>= i m))
	    (matrix-set! A i j (/ (matrix-ref A i j) pivot)))))

       ;; Compute the determinant, if necessary.  Note that we kept track of
       ;; its sign during row swaps.

       (if det
	   (do ((k 0 (+ k 1)))
	       ((>= k m) det)

	     ;; Look out for underflows:

	     (if (or (almost-zero? det) (almost-zero? (matrix-ref A k k)))
		 (exit 0)
		 (set! det (* det (matrix-ref A k k)))))
	   (begin
	     (do ((k 0 (+ k 1)))
		 ((>= k m))
	       (if (almost-zero? (matrix-ref A k k))
		   (exit 0)))
	     (exit 1)))))))

;;; Use LU-decomposition to solve a linear system of equations (signals error
;;; if the system is singular).

(define (LU-solve A . aux)

  (if (sparse-matrix? A)
      (set! A (sparse->matrix A)))

  (let ((m (matrix-row-count A))
	(n (matrix-column-count A)))

    (if (or (null? aux)
	    (not (eq? (car aux) 'no-copy)))
	(set! A (matrix-copy A)))

    (rref A)

    ;; Form the result:

    (if (> n (+ m 1))
	(let ((result (make-matrix m (- m n))))
	  (do ((i 0 (+ i 1)))
	      ((>= i m) result)
	    (do ((j m (+ j 1)))
		((>= j n))
	      (matrix-set! result i (- j m) (matrix-ref A i j)))))
	(let ((result (make-vector m)))
	  (do ((i 0 (+ i 1)))
	      ((>= i m) result)
	    (vector-set! result i (matrix-ref A i m)))))))

(define (rref A)
  ;; Solve Ax = b.

  (let ((m (matrix-row-count A))
	(n (matrix-column-count A)))

    (if (>= m n) (error "Incorrect dimensions -- LU-SOLVE"))

    ;; Get A = LU.

    (let ((det (lu-decomp A 'no-det)))
      (if (almost-zero? det)
	  (error (string-append "Singular matrix! (Determinant = "
				(number->string det)
				") -- LU-SOLVE"))))

    ;; Forward substitution to solve Ly = b.

    (do ((j 0 (1+ j)))
	((>= j m))
      (do ((i (1+ j) (1+ i)))
	  ((>= i m))
	(do ((k m (1+ k)))
	    ((>= k n))
	  (matrix-set! A i k
		       (- (matrix-ref A i k)
			  (* (matrix-ref A j k) (matrix-ref A i j))))
	  (matrix-set! A i j 0))))

    ;; Backward substitution to solve Ux = y.

    (do ((i (-1+ m) (-1+ i)))
	((< i 0))

      (let ((diag (matrix-ref A i i)))
	(do ((k i (1+ k)))
	    ((>= k n))
	  (matrix-set! A i k (/ (matrix-ref A i k) diag))))

      (do ((j (-1+ i) (-1+ j)))
	  ((< j 0))
	(let ((factor (matrix-ref A j i)))
	  (do ((k i (1+ k)))
	      ((>= k n))
	    (matrix-set! A j k
			 (- (matrix-ref A j k)
			    (* factor (matrix-ref A i k))))))))))


;;; Compute determinants without mutating the argument:

(define (det A)
  (lu-decomp (matrix-copy A)))


;;; A very useful procedure to have around:

(define (print-matrix matrix . argl)
  (if (null? argl)
      (set! argl (list (current-output-port))))

  (if (not (and (null? (cdr argl)) (output-port? (car argl))))
      (error "Invalid argument(s) -- PRINT-MATRIX"))

  (let ((port (car argl)))
    (newline port)
    (let ((m (matrix-row-count matrix))
          (n (matrix-column-count matrix)))
      (do ((i 0 (+ i 1)))
	  ((>= i m))
	(display (matrix-ref matrix i 0) port)

	(do ((j 1 (+ j 1)))
	    ((>= j n))
	  (display #\tab port)
	  (display (matrix-ref matrix i j) port))

	(newline port)))))


;;; Also very useful, especially for printing large matrices in a Emacs buffer:

(define (round-matrix mat precision)
  (let* ((m (matrix-row-count mat))
	 (n (matrix-column-count mat))
	 (out (make-matrix m n)))
    (do ((i 0 (+ i 1)))
	((>= i m) out)
      (do ((j 0 (+ j 1)))
	  ((>= j n))
	(let ((val (matrix-ref mat i j)))
	  (matrix-set! out i j (* (round (/ val precision)) precision)))))))
