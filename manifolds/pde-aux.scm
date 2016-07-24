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

;;; This file defines some auxiliary data structures for the PDE code.

(declare (usual-integrations))


;;; Linear equations:

(define (make-equation node constant terms)
  (vector constant terms node))

(define (equation:get-constant equation)
  (vector-ref equation 0))

(define (equation:get-terms equation)
  (vector-ref equation 1))

(define (equation:get-node equation)
  (vector-ref equation 2))

(define (equation:get-id equation)
  (node:get-id (equation:get-node equation)))

(define (null-equation? equation)
  (let loop ((terms (equation:get-terms equation)))
    (if (null? terms)
	(let ((node (equation:get-node equation)))
	  (list (node:get-point node) (node:boundary? node)))
	(if (almost-zero? (term:get-coeff (car terms)))
	    (loop (cdr terms))
	    #f))))


;;; Terms in linear equations:

(define (make-term node value)
  (vector node value))

(define (term:get-node term)
  (vector-ref term 0))

(define (term:get-id term)
  (node:get-id (term:get-node term)))

(define (term:get-coeff term)
  (vector-ref term 1))
