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

(load "load-ode")

(define result
  (show-time
   (lambda ()
     (rigid-body-path singular-init 1.))))

(define e-list
  (show-time
   (lambda ()
     (map (compose vector-first (make-rigid-body-energy 1. (sqrt 2) 2.) cadr)
	  result))))

(define e-errors
  (let ((ref (car e-list)))
    (show-time
     (lambda ()
       (map (lambda (val) (relative-error val ref)) e-list)))))
