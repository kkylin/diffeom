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

(load "rigid")


;;; Compute some expressions:

(for-each
 (lambda (make-sysder filename)
   (let ((port (open-output-file filename)))
     (pp (traditional->correct-order
	  (vector-tail
	   (show-time
	    (lambda ()
	      (*sysder-simplify*
	       ((make-sysder 'a 'b 'c) rigid-qqdot))))
	   1))
	 port)
     (close-output-port port)))
 (list rigid-sysder-0 rigid-sysder-1 rigid-sysder-2 rigid-sysder-3)
 (list "rigid-field-0" "rigid-field-1" "rigid-field-2" "rigid-field-3"))
