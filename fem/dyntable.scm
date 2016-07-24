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

(define (make-dynamic-table . argl)
  (if (or (null? argl)
	  (> (length argl) 0))
      (vector 0 (make-vector 8))
      (vector 0 (make-vector (car argl)))))

(define (dynamic-table-size table)
  (vector-ref table 0))

(define (dynamic-table-add table new-element)
  (let ((size (dynamic-table-size table))
	(real-size (vector-length (vector-ref table 1))))
    (if (= size real-size)
	(let ((old-table (vector-ref table 1))
	      (new-table (make-vector (* 2 real-size))))
	  (do ((i 0 (+ i 1)))
	      ((>= i real-size))
	    (vector-set! new-table i (vector-ref old-table i)))
	  (vector-set! table 1 new-table)))

    (vector-set! (vector-ref table 1) size new-element)
    (vector-set! table 0 (+ size 1))))

(define (dynamic-table-fetch table i)
  (if (>= i (dynamic-table-size table))
      (error "Access out of bound -- DYNAMIC-TABLE-FETCH")
      (vector-ref (vector-ref table 1) i)))

(define (dynamic-table->list table)
  (let loop ((l '()) (n (- (dynamic-table-size table) 1)))
    (if (< n 0)
	l
	(loop (cons (dynamic-table-fetch table n) l) (- n 1)))))

(define (dtableq table element)
  (let loop ((n (- (dynamic-table-size table) 1)))
    (if (< n 0)
	#f
	(if (eq? (dynamic-table-fetch table n) element)
	    n
	    (loop (- n 1))))))

(define (dynamic-table-set! table i val)
  (if (>= i (dynamic-table-size table))
      (error "Access out of bound -- DYNAMIC-TABLE-SET!")
      (vector-set! (vector-ref table 1) i val)))
