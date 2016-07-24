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

;;; This file defines tools for dealing with elements.

(declare (usual-integrations))


;;; A template for element constructors:

(define (pde:element-maker L
			   make-integrator
			   make-basis-function)

  ;; The operator abstraction in the FEM program is a bit artificial (see
  ;; ELEMENT-MAKER).  Note that the interpolation of variable coefficients
  ;; *there* may not be necessary, as the operator is already given a PDE-chart
  ;; along with a list of the nodes.

  (lambda (chart)

    ;; This ugly hack sort of works (for now).  Really ought to just make
    ;; ELEMENT-MAKER to do the right thing (whatever that is).

    (element-maker
     (lambda (nodes)
       (operator:set-context! L chart nodes)
       L)
     make-integrator
     make-basis-function)))


;;; We need to know when a node belongs to an element.  Note that this only
;;; works with simplices.  For more complicated shapes, we would require more
;;; structure (i.e. a list of *faces* of the boundary of the convex element),
;;; which can only be supplied by TESSELATE during domain construction.

(define (element:member? element p)
  (let ((vertices (map node:get-coords (element:get-vertex-nodes element))))
    (in-simplex? p vertices)))


;;; Some procedures that help with ELEMENT:MEMBER?.

(define (in-simplex? point vertices)
  (in-convex-domain?
   point vertices (choose-sublists vertices (vector-length point))))

(define (in-convex-domain? point vertices faces)
  (not (memq #f (map (lambda (face)
		       (same-side? point
				   (find-another-vertex face vertices)
				   (car face)
				   (cdr face)))
		     faces))))

(define (find-another-vertex face vertices)
  (let loop ((vertices vertices))
    (if (null? vertices)
	#f
	(let ((vertex (car vertices)))
	  (if (member vertex face)
	      (loop (cdr vertices))
	      vertex)))))

(define (same-side? p q origin basis)

  ;; See if P and Q lie on the same side of the n-1-dimensional hyperplane
  ;; defined by translating the span of BASIS from 0 to ORIGIN.

  (let ((basis (map (lambda (v) (vector:- v origin)) basis))
	(p (vector:- p origin))
	(q (vector:- q origin)))

    (let ((val (* (det (list->vector (cons p basis)))
		  (det (list->vector (cons q basis))))))
      (or (>= val 0)
	  (almost-zero? val)))))
