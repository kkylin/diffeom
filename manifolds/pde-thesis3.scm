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

;;; Turns out one of the bug fixes introduced an error.  The new code with the
;;; error got reloaded into the system, and...

(load "pde-test")


(test-5 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test5j")

(test-1 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test1k")
(test-2 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test2k")
(test-3 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test3k")
(test-4 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test4k")
(test-5 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test5k")
