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


;;; Some definitions that are always useful:

(define (square z)
  (* z z))

(define (compose f . rest)
  (let loop ((f f) (l rest))
    (if (null? l)
	f
	(loop (binary-compose f (car l)) (cdr l)))))

(define (binary-compose f g)
  (lambda (first . rest)
    (f (apply g (cons first rest)))))

(define pi (* 4 (atan 1)))
(define -pi (- pi))

(define (identity x)
  x)

(define (relative-error val ref)
  (if (zero? ref)
      val
      (/ (- val ref) ref)))
