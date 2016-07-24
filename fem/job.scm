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

(load "load")

(define nodes (make-real-square-domain 3 3))
(define m (sparse->matrix (fem 0-function nodes node:get-x)))
(define v (lu-solve m))

(let ((port (open-output-file "test"))
      (n (vector-length v)))
  (print-matrix m port)
  (newline port)
  (newline port)
  (do ((i 0 (+ i 1)))
      ((>= i n))
    (write-line (vector-ref v i) port))
  (close-output-port port))
