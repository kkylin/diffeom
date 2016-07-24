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

;;; Load the PDE stuff separately (so that the parts of the program *not* under
;;; development can still be used without these definitions).

(let ((pde-core
       '("pde-aux"
	 "pde-charts"
	 "pde-cmpgrd"
	 "pde-debug"
	 "pde-elements"
	 "pde-main"
	 "pde-mergers"
	 "pde-nodes"
	 "pde-ops"
	 "pde-tools"
	 "pde-config"

	 "basis-imb"
	 "basis-poly"
	 "basis-real"

	 "pde-examples"))

      (fem-stuff
       '("basis"
	 "2d-poly-basis"
	 "2d-real-basis"
	 "2d-trapezoid"
	 "delaunay"
	 "delaux"
	 "dyntable"
	 "edge"
	 "fem"
	 "relax"
	 "sparse"
	 "util-too"))

      (fem-dir "fem/"))

  (load "load-main")

  (if *using-scmutils?*
      (begin
	(newline)
	(display "*** Warning: PDE code does not work well with ScmUtils!")
	(newline)))

  (for-each (lambda (file)
	      (load (string-append fem-dir file)))
	    fem-stuff)

  (for-each load pde-core))
