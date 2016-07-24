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

;;; Just the last two tests, which take more memory, from pde-collect.scm.

(load "pde-test")


;;; Define the test procedures:

(define test-7
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap1))

(define test-8
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap2))


;;; Run the experiments (sorted by size, not test):

(test-7 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test7d")
(test-8 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test8d")
