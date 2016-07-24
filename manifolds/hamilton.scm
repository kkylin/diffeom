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

;;; This file implements Hamiltonian mechanics.  Given a sufficiently efficient
;;; implementation of cotangent bundles, this should run faster than using
;;; Lagrangians.

(declare (usual-integrations))


;;; The Hamiltonian should be a smooth map from the cotangent bundle of some
;;; manifold into the real line.

(define (hamiltonian->v.field H)
  (let ((T*M (smooth-map:get-domain H))
	(R (smooth-map:get-range H)))
    (lambda (p)
      (let ((U (manifold:find-best-chart T*M p)))
	(make-tangent U p
		      (hamilton-in-coords
		       (smooth-map:make-transition
			H U (car (manifold:get-finite-atlas R)))
		       (chart:point->coords p U)))))))


;;; Derive Hamilton's equations for f at x (in coordinates):

(define (hamilton-in-coords f x)
  (let* ((2n (vector-length x))
	 (v (make-vector 2n))
	 (n (/ 2n 2)))

    (do ((i n (+ i 1))
	 (j 0 (+ j 1)))
	((>= j n) v)

      (vector-set! v i (- (vector-first (((pdiff j) f) x))))
      (vector-set! v j (vector-first (((pdiff i) f) x))))))
