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

;;; This file loads the appropriate definitions for the numerical
;;; differentiation of real functions.

(declare (usual-integrations))
(load "manifolds/linear")
(load "manifolds/lshared")
(load "manifolds/richardson")


;;; Stolen from manifolds/misc.scm:

(define (pdiff i)
  (lambda (f)
    (let ((df (diff f)))
      (lambda (x)
	(let ((v (make-vector (vector-length x) 0)))
	  (vector-set! v i 1)
	  ((df x) v))))))
