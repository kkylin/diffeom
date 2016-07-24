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

;;; Use this file to collect data for theis work.  Based on pde-test.scm.

(load "pde-test")


;;; Define the test procedures:

(define test-1
  (pde:experiment-too 'pde:make-domain-without-overlaps
		      'combine-equations-without-overlap))

(define test-2
  (pde:experiment-too 'pde:make-domain-with-small-overlaps
		      'combine-equations-without-overlap))

(define test-3
  (pde:experiment-too 'pde:make-domain-with-overlaps
		      'combine-equations-with-overlap))

(define test-4
  (pde:experiment-too 'pde:make-domain-with-larger-overlaps
		      'combine-equations-with-overlap))

(define test-5
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap))

(define test-6
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-using-cmpgrd))


;;; Run the experiments (sorted by size, not test):

(test-1 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test1a")
(test-2 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test2a")
(test-3 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test3a")
(test-4 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test4a")
(test-5 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test5a")
(test-6 '(rectangular 10 5) '(spherical 5 10) 10000 1.7 "Data/disc/test6a")

(test-1 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test1b")
(test-2 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test2b")
(test-3 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test3b")
(test-4 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test4b")
(test-5 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test5b")
(test-6 '(rectangular 20 10) '(spherical 10 20) 10000 1.7 "Data/disc/test6b")

(test-1 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test1c")
(test-2 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test2c")
(test-3 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test3c")
(test-4 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test4c")
(test-5 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test5c")
(test-6 '(rectangular 40 15) '(spherical 20 30) 10000 1.7 "Data/disc/test6c")

(test-1 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test1d")
(test-2 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test2d")
(test-3 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test3d")
(test-4 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test4d")
(test-5 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test5d")
(test-6 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test6d")
