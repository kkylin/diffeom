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

;;; Figure out if ScmUtils is loaded (by checking if a few key procedures are
;;; defined):

(define *using-scmutils?*
  (let ((procs '(derivative vector:scalar*vector)))
    (and mit-scheme?
	 (not (memq #f (map
			(lambda (proc)
			  (environment-bound? (the-environment) proc))
			procs))))))

(newline)
(display
 (if *using-scmutils?*
     "*** It looks like we're running ScmUtils..."
     "*** ScmUtils doesn't seem to be running.  Going numerical..."))
(newline)

(let ((preload '("misc"
		 "lshared"))

      (core '("charts"
	      "manifold"
	      "vbundle"
	      "product"
	      "boundary"
	      "ranges"
	      "smooth"
	      "spaces"))

      (numdiff '("misc-math"
		 "fem/matlib"
		 "linear"
		 "richardson"))

      (scmutils '("stubs")))

  (for-each load preload)
  (for-each load (if *using-scmutils?* scmutils numdiff))
  (for-each load core))
