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

(declare (usual-integrations))


;;; Sometimes useful (stolen from nscmutils):

(define (make-comparator tol)
  (lambda (a b)
    (< (magnitude (- a b))
       (* .5 tol
	  (+ (magnitude a) (magnitude b) 2.0)))))

(define almost-equal? (make-comparator 1e-10))

(define (almost-zero? x)
  (almost-equal? x 0))


;;; Useful in making product manifolds:

(define (all-pairs l1 l2)
  (let loop1 ((l1 l1) (result '()))
    (if (null? l1)
	result
	(let ((obj (car l1)))
	  (let loop2 ((l2 l2) (result result))
	    (if (null? l2)
		(loop1 (cdr l1) result)
		(loop2 (cdr l2) (cons (list obj (car l2)) result))))))))


;;; How do you tell if an object is a vector in R^n?

(define (make-euclidean-test dim)
  (lambda (v)
    (and (vector? v)
	 (= (vector-length v) dim)
	 (let loop ((i 0))
	   (if (< i dim)
	       (and (real? (vector-ref v i))
		    (loop (+ i 1)))
	       #t)))))


;;; Always useful to memoize things:

(define (simple-memoize proc size)
  ;; Memoize a function whose argument is a non-negative integer:
  (let ((cache (make-vector size 'undefined)))
    (lambda (n)
      (if (>= n size)
	  (proc n)
	  (let ((val (vector-ref cache n)))
	    (if (eq? val 'undefined)
		(let ((val (proc n)))
		  (vector-set! cache n val)
		  val)
		val))))))


;;; Useful in our implementation of spherical coordinates:

(define (list-integers i)
  (let loop ((result '()) (i i))
    (if (< i 0)
	result
	(loop (cons i result) (- i 1)))))


;;; Surprisingly enough, this is not in MIT Scheme:

(define (cot z)
  (/ (tan z)))


;;; PARTIAL in ScmUtils doesn't do the right thing for functions of vector
;;; arguments, so we can't just use that and stub out the numerical equivalent:

(define (pdiff i)
  (lambda (f)
    (let ((df (diff f)))
      (lambda (x)
	(let ((v (make-vector (vector-length x) 0)))
	  (vector-set! v i 1)
	  ((df x) v))))))
