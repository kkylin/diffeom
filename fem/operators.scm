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

;;; This file contains some ad-hoc definitions of simple differential
;;; operators, such as the two-dimensional Laplacian or the 1+1-dimensional
;;; d'Alembertian.  It also defines some important constructors and operations
;;; on differential operators that compute the adjoint and encapsulate
;;; integration-by-parts.  (See ELEMENT-MAKER in fem.scm.)

;;; opalg.scm contains the beginnings of a much more abstract (and complete)
;;; approach.

;;; This file uses various procedures from basis.scm.

(declare (usual-integrations))


;;; Simple constructor and methods for an operator structure:

(define (make-operator left-op right-op combine)
  (vector left-op right-op combine))

(define (operator:get-left-op operator)
  (vector-ref operator 0))

(define (operator:get-right-op operator)
  (vector-ref operator 1))

(define (operator:get-combine operator)
  (vector-ref operator 2))

(define (operator:get-local-form op)
  (let ((combine (operator:get-combine op))
	(left-op (operator:get-left-op op))
	(right-op (operator:get-right-op op)))
    (lambda (f g)
      (combine (left-op f) (right-op g)))))
