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

;;; Let's define some differential operators so we have something to test.

(declare (usual-integrations))


;;; Differential operators (ours only act on scalar functions, not sections of
;;; vector bundles):

(define (make-operator M make-local-form)
  (vector M make-local-form #f '()))

(define (operator:get-manifold L)
  (vector-ref L 0))

;;; Some extra structures for working with FEM.  "Context" is a kluge that lets
;;; operators learn about the particular element they are working in, etc.

(define (operator:get-local-form operator)
  (let ((context (operator:get-context operator)))
    (if context
	(apply (vector-ref operator 1) context)
	#f)))

(define (operator:set-context! operator . contextual-data)
  (vector-set! operator 2 contextual-data))

(define (operator:get-context operator)
  (vector-ref operator 2))

;;; This might come in handy:

(define (operator:install-extra L tag datum)
  (let ((result (assq tag (vector-ref L 4))))
    (if result
	(set-cdr! result datum)
	(vector-set! L 4 (cons (cons tag datum) (vector-ref L 4))))))

(define (operator:get-extra L tag)
  (let ((result (assq tag (vector-ref L 4))))
    (if result
	(cdr result)
	#f)))
