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

;;; Even more tests:

(load "pde-test")


;;; Define the test procedures:

(define test-7
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap1))

(define test-8
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap2))


;;; Run the experiments (sorted by size, not test):

(test-7 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test7a")
(test-8 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test8a")

(test-7 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test7b")
(test-8 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test8b")

(test-7 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test7c")
(test-8 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test8c")
