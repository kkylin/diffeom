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

;;; This file defines some procedures that are useful for working with
;;; imbedding representations of manifolds:

(declare (usual-integrations))


;;; For tangent vectors:

(define (imbedding->tangent M p v)
  (let ((U (manifold:find-best-chart M p)))
    (make-tangent U p (push-forward-in-coords (chart:get-coord-map U) p v))))

(define (tangent->imbedding v)
  (let* ((U (tangent:get-chart v))
	 (p (tangent:get-anchor v)))
    (list p
	  (push-forward-in-coords (chart:get-inverse-map U)
				  (chart:point->coords p U)
				  (tangent:get-coords v)))))

(define tangent->imbedded-velocity
  (compose cadr tangent->imbedding))


;;; For cotangent vectors:

(define (imbedding->cotangent M p v)
  (let ((U (manifold:find-best-chart M p)))
    (make-cotangent U p (pull-back-in-coords (chart:get-inverse-map U)
					     (chart:point->coords p U)
					     v))))

(define (cotangent->imbedding v)
  (let* ((U (cotangent:get-chart v))
	 (p (cotangent:get-anchor v)))
    (list p

	  ;; Need to project the pulled-back functional onto the imbedded
	  ;; surface because it's represented in the standard basis of the
	  ;; ambient space, not the basis of the tangent space to M imbedded
	  ;; inside R^n.

	  (project-onto-basis
	   (make-imbedded-basis U p)
	   (pull-back-in-coords (chart:get-coord-map U)
				p
				(cotangent:get-coords v))))))

(define (make-imbedded-basis chart x)
  (let ((dim (chart:dimension chart)))
    (let loop ((i 0) (vlist '()))
      (if (< i dim)
	  (loop (+ i 1) (cons (tangent->imbedded-velocity
			       (make-tangent chart x (vector:basis dim i 1)))
			      vlist))
	  vlist))))
