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

;;; Just the last two tests, which take more memory, from pde-collect.scm.

(load "pde-test")


(define test-5
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-with-overlap))

(define test-6
  (pde:experiment-too 'pde:make-simple-domain
		      'combine-equations-using-cmpgrd))

(test-5 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test5d")
(test-6 '(rectangular 60 30) '(spherical 37 61) 10000 1.7 "Data/disc/test6d")
