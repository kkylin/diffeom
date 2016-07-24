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

;;; In fixing a previous bug, I introduced another bug.

(load "pde-test")


;;; Here goes again:

(test-1 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test1a")
(test-2 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test2a")

(test-1 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test1b")
(test-2 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test2b")

(test-1 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test1c")
(test-2 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test2c")

(test-1 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test1d")
(test-2 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test2d")

(test-1 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test1e")
(test-2 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test2e")

(test-1 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test1f")
(test-2 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test2f")

(test-1 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test1g")
(test-2 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test2g")

(test-1 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test1h")
(test-2 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test2h")

(test-1 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test1i")
(test-2 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test2i")

(test-1 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test1j")
(test-2 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test2j")

(test-1 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test1k")
(test-2 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test2k")
