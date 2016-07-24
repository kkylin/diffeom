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

;;; Some special structures that lets us tell simple *range* shapes.  This
;;; makes mesh generation for the PDE solver much easier.

(declare (usual-integrations))


;;; Spherical ranges are very useful (e.g. stereographic projection):

(define (make-spherical-range chart center radius)
  (chart:install-extra chart 'spherical-range (vector center radius)))

(define (spherical-range:get-structs chart)
  (chart:get-extra chart 'spherical-range))

(define (spherical-range:get-center chart)
  (let ((result (spherical-range:get-structs chart)))
    (if result
	(vector-ref result 0)
	#f)))

(define (spherical-range:get-radius chart)
  (let ((result (spherical-range:get-structs chart)))
    (if result
	(vector-ref result 1)
	#f)))

(define (chart:spherical-range? chart)
  (if (spherical-range:get-structs chart)
      #t
      #f))


;;; Ranges that are n-cells are also very useful (e.g. spherical coordinates).

(define (make-cell-range chart intervals)
  (chart:install-extra
   chart 'cell-range (if (list? intervals)
			 (list->vector intervals)
			 intervals)))

(define (cell-range:get-structs chart)
  (chart:get-extra chart 'cell-range))

(define (cell-range:get-interval chart i)
  (let ((result (cell-range:get-structs chart)))
    (if result
	(vector result i)
	#f)))

(define cell-range:get-intervals cell-range:get-structs)

(define (cell-range:get-interval-list chart)
  (let ((result (cell-range:get-intervals chart)))
    (if result
	(vector->list result)
	#f)))

(define (chart:cell-range? chart)
  (if (cell-range:get-structs chart)
      #t
      #f))


;;; The interval:

(define (make-interval a b)
  (vector a b))

(define (interval:inf interval)
  (vector-ref interval 0))

(define (interval:sup interval)
  (vector-ref interval 1))

(define (interval:member? x interval)
  (and (real? x)
       (< (interval:inf interval) x)
       (< x (interval:sup interval))))
