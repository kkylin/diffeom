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

;;; This file defines abstract operator algebra.  It turns out to be more
;;; general than what we need.  If we're going to be this abstract anyway, why
;;; don't we just use the polynomial code?  Anyway, this is the wrong thing, so
;;; we won't continue along this line of work.  Let's keep it around, though,
;;; just in case we need it someday.

(declare (usual-integrations))


;;; We use the multi-index notation:

(define (make-multi-index n)
  (make-vector n 0))

(define multi-index vector)

(define multi-index:length vector-length)

(define multi-index:get vector-ref)

(define multi-index:set! vector-set!)

(define (multi-index:sum multind)
  (let ((len (multi-inde:length multind)))
    (let loop ((sum 0) (i 0))
      (if (< i len)
	  (loop (+ sum (multi-index:get multind i)) (+ i 1))
	  sum))))

(define (mult-index:fact multind)
  (let ((len (multi-index:length multind)))
    (let loop ((prod 1) (i 0))
      (if (< i len)
	  (loop (* prod (factorial (mult-index:get multind i))) (+ i 1))
	  prod))))

(define (multi-index:+ mi1 mi2)
  (let* ((len (multi-index:length mi1))
	 (result (make-multi-index len)))
    (do ((i 0 (+ i 1)))
	((>= i len) result)
      (multi-index:set! result i
			(+ (multi-index:get mi1 i) (multi-index:get mi2 i))))))


;;; Make a partial differential operator with one summand:

(define (multi-index->diffop multind)
  `((1 ,multind)))


;;; Elementary operations:

(define (diffop:+ op1 op2)
  (append op1 op2))

(define (diffop:- op1 op2)
  (append op1 (diffop:negate op2)))

(define diffop:first-term car)

(define diffop:remaining-terms cdr)

(define (diffop:zero? op)
  (null? (diffop:simplify op)))


;;; Composing differential operators:

(define (diffop:compose op1 op2)
  (if (diffop:zero? op1)
      op2
      (diffop:compose (diffop:remaining-terms op1)
		      (let ((op1 (diffop:first-term op1)))
			(append-map
			 (lambda (op2)
			   (diffop:compose-terms op1 op2))
			 op2)))))
