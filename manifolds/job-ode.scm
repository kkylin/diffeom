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

(load "load-ode")

(show-time
 (lambda ()

   (let ((segment-size 1.)
	 (count 101)
	 (filename "rigid-man.data")
	 (charts (manifold:get-finite-atlas TSO3)))

     (let loop ((i 0) (x singular-init) (t 0.))
       (if (< i count)
	   (begin

	     (write-line `(step ,i t = ,t))

	     (let ((results
		    (show-time
		     (lambda ()
		       (rigid-body-path x (+ t segment-size) t)))))

	       (let ((port (open-output-file filename #t)))

		 (for-each

		  (lambda (l)
		    (let ((t (car l))
			  (p (cadr l)))

		      (display t port)

		      (let loop ((i 0) (charts charts))
			(if (null? charts)
			    (display "No chart!" port)
			    (let ((chart (car charts)))
			      (if (chart:member? p chart)
				  (begin
				    (display " " port)
				    (display i port)
				    (display " " port)
				    (display (chart:point->coords p chart)
					     port))
				  (loop (+ i 1) (cdr charts))))))

		      (newline port)))

		  (let ((l (sort results (lambda (x y) (< (car x) (car y))))))
		    (if (> t 0)
			(cdr l)
			l)))

		 (close-output-port port))

	       (loop (+ i 1) (cadar results) (caar results)))))))))
