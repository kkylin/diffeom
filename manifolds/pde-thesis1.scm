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

;;; Use this file to collect data for theis work.  Based on pde-collect.scm.

(load "pde-test")


;;; Run the experiments (sorted by size, not test):

(test-6 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test6a")
(test-7 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test7a")
(test-8 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test8a")

(test-6 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test6b")
(test-7 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test7b")
(test-8 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test8b")

(test-6 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test6c")
(test-7 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test7c")
(test-8 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test8c")

(test-6 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test6d")
(test-7 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test7d")
(test-8 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test8d")

(test-6 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test6e")
(test-7 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test7e")
(test-8 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test8e")

(test-6 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test6f")
(test-7 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test7f")
(test-8 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test8f")

(test-6 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test6g")
(test-7 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test7g")
(test-8 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test8g")

(test-6 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test6h")
(test-7 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test7h")
(test-8 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test8h")

(test-6 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test6i")
(test-7 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test7i")
(test-8 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test8i")

(test-6 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test6j")
(test-7 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test7j")
(test-8 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test8j")

(test-6 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test6k")
(test-7 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test7k")
(test-8 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test8k")
