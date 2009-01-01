;;; This file plays around with Richardson extrapolation.  ORDER sould be the
;;; order of the error.

(declare (usual-integrations))


;;; Here is the slow way:

(define (richardson f order)
  (if (> order 0)
      (let ((f (richardson f (- order 1)))
	    (k (expt 2 order)))
	(lambda (h)
	  (/ (- (* k (f (/ h 2))) (f h)) (- k 1))))
      f))


;;; Here is a quicker way:

(define (richardson-coeffs order)
  (let ((v (make-vector order 0))
	(w (make-vector order 0)))

    (vector-set! v 0 1)

    (let loop ((i 1) (2^i 1) (from v) (to w))
      (if (< i order)
	  (let ((2^i+1 (* 2 2^i))
		(i-1 (- i 1)))
	    (vector-set! to 0 (/ (vector-ref from 0) (- 1 2^i+1)))

	    (do ((j 1 (+ j 1)))
		((> j i-1))
	      (vector-set! to j (/ (- (* 2^i+1 (vector-ref from (- j 1)))
				      (vector-ref from j))
				   (- 2^i+1 1))))

	    (vector-set! to i (* (/ 2^i+1 (- 2^i+1 1))
				 (vector-ref from i-1)))

	    (loop (+ i 1) 2^i+1 to from))
	  from))))

(define (quick-r f order)
  (let ((v (richardson-coeffs order)))
    (lambda (h)
      (let loop ((i 0) (h/2^i h) (sum 0.))
	(if (< i order)
	    (loop (+ i 1) (/ h/2^i 2) (+ sum (* (vector-ref v i) (f h/2^i))))
	    sum)))))


;;; Try some numerical differentiation:

;;; An observation that may make this better: When the derivative is computed
;;; using the central difference (instead of forward or backward difference),
;;; only the even-degree terms in the difference quotient survive.

;;; Also, try using roots of unity (complex arithmetic) to get rid of
;;; higher-order terms.  Better yet, use contour integrals!

(define (make-differentiator v+ v- v* h order)
  (let ((v (richardson-coeffs order)))
    (lambda (f)
      (let ((diff-quo
	     (lambda (x)
	       (lambda (h)
		 (v* (/ (* 2 h)) (v- (f (+ x h)) (f (- x h))))))))
	(lambda (x)
	  (let ((f (diff-quo x)))
	    (let loop ((i 1)
		       (h/2^i (/ h 2))
		       (sum (v* (vector-ref v 0) (f h))))
	      (if (< i order)
		  (loop (+ i 1) (/ h/2^i 2)
			(v+ sum (v* (vector-ref v i) (f h/2^i))))
		  sum))))))))

;(define derivative (make-differentiator + - * 1e-5 3))

;;; For now, let's use numerical differentiation with Richardson extrapolation,
;;; and use charts to represent tangent vectors.  Later, we should provide ways
;;; of automagically switching between different representations (e.g. among
;;; different charts and/or between chart-vector representations and imbedding
;;; representations (imbeddings *are* very useful)).

(define (make-numerical-differential-operator h n scale)
  (let ((d (make-differentiator vector:+ vector:- vector:* h n)))
    (lambda (f)
      (lambda (x)
	(lambda (v)
	  (let* ((change-factor (/ (vector:magnitude v) scale))
		 (w (if (zero? change-factor)
			v
			(vector:* (/ change-factor) v)))
		 (on-path (lambda (s) (f (vector:+ x (vector:* s w))))))
	    (vector:* change-factor ((d on-path) 0))))))))


;;; Is this basically equivalent to Ridder's method, described in _Numerical
;;; Recipes in Fortran, Second Edition_?


;;; Minimize a convex function:

(define (minimize-convex-function f a b n)
  (let loop ((n n) (a a) (b b))
    (if (> n 0)
	(let ((x1 (+ (* 2/3 a) (/ b 3)))
	      (x2 (+ (/ a 3) (* 2/3 b))))
	  (let ((y1 (f x1))
		(y2 (f x2)))
	    (cond ((< y1 y2) (loop (- n 1) a x2))
		  ((> y1 y2) (loop (- n 1) x1 b))
		  (else (loop (- n 1) x1 x2)))))
	(let ((result (/ (+ a b) 2)))
	  (list result (f result))))))


;;; Let's define the differential operator:

(define diff
  (let ((order 5)
	(min-step-size 1e-5)
	(max-step-size 1e-1)
	(search-depth 11)
	(scale 1))
    (make-numerical-differential-operator
     (let ((result
	    (minimize-convex-function
	     (lambda (h)
	       (let ((d (make-numerical-differential-operator h order scale)))
		 (abs (- (vector-ref (((d (lambda (x)
					    (vector (sqrt (vector-ref x 0)))))
				       (vector 1))
				      (vector 1))
				     0)
			 1/2))))
	     min-step-size max-step-size search-depth)))
       (write-line `(richardson order = ,order))
       (write-line `(h = ,(car result)))
       (write-line `(error = ,(cadr result)))
       (car result))
     order
     scale)))

;(define diff-old (make-numerical-differential-operator 1e-5 3))
