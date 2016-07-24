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

;;; Memoization is frequently useful:

(declare (usual-integrations))

;;; ALWAYS useful:

(define (square z)
  (* z z))

(define (cube z)
  (* z z z))


;;; Even more useful:

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

(define almost-zero?
  (let ((*tolerance* 1e-14))
    (lambda (z)
      (< (magnitude z) *tolerance*))))


;;; Why is this not defined elsewhere?

(define (compose f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (binary-compose f (car l)) (cdr l)))))

(define (binary-compose f g)
  (lambda (first . rest)
    (f (apply g (cons first rest)))))
