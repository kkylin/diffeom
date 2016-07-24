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

;;; Use this file to collect data for theis work.  Based on pde-collect.scm.

(load "pde-test")


;;; Run the experiments (sorted by size, not test):

(test-1 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test1a")
(test-2 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test2a")
(test-3 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test3a")
(test-4 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test4a")
(test-5 '(rectangular 10 5) '(spherical 5 10) 10000 1.9 "Data/thesis/test5a")

(test-1 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test1b")
(test-2 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test2b")
(test-3 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test3b")
(test-4 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test4b")
(test-5 '(rectangular 14 7) '(spherical 7 14) 11000 1.9 "Data/thesis/test5b")

(test-1 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test1c")
(test-2 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test2c")
(test-3 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test3c")
(test-4 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test4c")
(test-5 '(rectangular 18 9) '(spherical 9 18) 12000 1.9 "Data/thesis/test5c")

(test-1 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test1d")
(test-2 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test2d")
(test-3 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test3d")
(test-4 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test4d")
(test-5 '(rectangular 22 11) '(spherical 11 22) 13000 1.9 "Data/thesis/test5d")

(test-1 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test1e")
(test-2 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test2e")
(test-3 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test3e")
(test-4 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test4e")
(test-5 '(rectangular 26 13) '(spherical 13 26) 14000 1.9 "Data/thesis/test5e")

(test-1 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test1f")
(test-2 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test2f")
(test-3 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test3f")
(test-4 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test4f")
(test-5 '(rectangular 30 15) '(spherical 15 30) 15000 1.9 "Data/thesis/test5f")

(test-1 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test1g")
(test-2 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test2g")
(test-3 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test3g")
(test-4 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test4g")
(test-5 '(rectangular 34 17) '(spherical 17 34) 16000 1.9 "Data/thesis/test5g")

(test-1 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test1h")
(test-2 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test2h")
(test-3 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test3h")
(test-4 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test4h")
(test-5 '(rectangular 38 19) '(spherical 19 38) 17000 1.9 "Data/thesis/test5h")

(test-1 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test1i")
(test-2 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test2i")
(test-3 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test3i")
(test-4 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test4i")
(test-5 '(rectangular 42 21) '(spherical 21 42) 18000 1.9 "Data/thesis/test5i")

(test-1 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test1j")
(test-2 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test2j")
(test-3 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test3j")
(test-4 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test4j")
(test-5 '(rectangular 46 23) '(spherical 23 46) 19000 1.9 "Data/thesis/test5j")

(test-1 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test1k")
(test-2 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test2k")
(test-3 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test3k")
(test-4 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test4k")
(test-5 '(rectangular 50 25) '(spherical 25 50) 20000 1.9 "Data/thesis/test5k")
