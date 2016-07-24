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

;;; This implements the topological and geometric primitives needed for the
;;; Delaunay triangulation algorithm described in:
;;;
;;; "Primitives for the manipulation of general subdivisions and the
;;;  computation of Voronoi diagrams," ACM transactions on graphics, Vol. 4,
;;; No. 2, April 1985, P. 74-123.
;;; Leonidas Guibas and Jorge Stolfi, Xerox PARC and Stanford University.
;;;
;;; This needs the file edge.scm for the edge-algebraic definitions.
;;; The nodal structures are assumed to provide the methods NODE:GET-X and
;;; NODE:GET-Y.

(declare (usual-integrations))


;;; Edge-record methods:

(define (org e-ref)
  (((get-edge-record e-ref) (get-rot-deg e-ref)) 'data))

(define (left e-ref)
  (org (inv-rot e-ref)))

(define (right e-ref)
  (org (rot e-ref)))

(define (dest e-ref)
  (org (sym e-ref)))

(define (set-org! e-ref new)
  ((((get-edge-record e-ref) (get-rot-deg e-ref)) 'set-data!) new))

(define (set-dest! e-ref new)
  (set-org! (sym e-ref) new))


;;; Topological operators:

(define (connect a b)
  ;; Create an edge E that connects A.Dest to B.Org, such that A.Left = E.Left
  ;; = B.Left after the connection is complete.
  (let ((e (make-edge)))
    (set-org! e (dest a))
    (set-dest! e (org b))
    (splice e (lnext a))
    (splice (sym e) b)
    e))

(define (delete-edge e)
  (splice e (oprev e))
  (splice (sym e) (oprev (sym e)))
  (dynamic-table-set! *delaunay-edges* (get-edge-id e) 'deleted))

(define (swap e)
  (let ((a (oprev e))
        (b (oprev (sym e))))
    (splice e a)
    (splice (sym e) b)
    (splice e (lnext a))
    (splice (sym e) (lnext b))
    (set-org! e (dest a))
    (set-dest! e (dest b))))


;;; Geometric primitives:

(define (in-circle a b c d)
  ;; a, b, c, and d should be 2-vectors.
  (let ((m (make-matrix 4 4)))

    (do ((l (map node:get-coords (list a b c d)) (cdr l))
         (i 0 (+ i 1)))
        ((null? l))
      (let* ((p (car l))
             (x (vector-ref p 0))
             (y (vector-ref p 1)))
        (matrix-set! m i 0 x)
        (matrix-set! m i 1 y)
        (matrix-set! m i 2 (+ (square x) (square y)))
        (matrix-set! m i 3 1)))

    (> (det m) 0)))

(define (ccw a b c)
  ;; a, b, and c should be 2-vectors.
  (let ((m (make-matrix 3 3)))

    (do ((l (map node:get-coords (list a b c)) (cdr l))
         (i 0 (+ i 1)))
        ((null? l))
      (let* ((p (car l))
             (x (vector-ref p 0))
             (y (vector-ref p 1)))
        (matrix-set! m i 0 x)
        (matrix-set! m i 1 y)
        (matrix-set! m i 2 1)))

    (> (det m) 0)))

(define (right-of x e)
  (ccw x (dest e) (org e)))

(define (left-of x e)
  (ccw x (org e) (dest e)))


;;; A procedure that lists all of the mesh elements in a Delaunay triangulation
;;; can be very useful, particularly for our FEM applications.

(define (get-edge-mark e-ref)
  (((get-edge-record e-ref) (get-rot-deg e-ref)) 'mark))

(define (set-edge-mark! e-ref val)
  ((((get-edge-record e-ref) (get-rot-deg e-ref)) 'set-mark!) val))

(define (list-faces)
  (let ((edges (list-edges)))
    (if (< (length edges) 3)
	'()
	(let ((faces '()))

	  ;; Reset markings:
	  (for-each
	   (lambda (e)
	     (set-edge-mark! e #f)
	     (set-edge-mark! (sym e) #f))
	   edges)

	  ;; Begin DFS:
	  (let loop ((e (car edges)))
	    (for-each
	     (lambda (a)
	       (if (false? (get-edge-mark a))
		   (let* ((b (lnext a))
			  (c (lnext b)))
		     (set-edge-mark! a #t)
		     (if (node= (dest c) (org a))
			 (begin
			   (set-edge-mark! b #t)
			   (set-edge-mark! c #t)
			   (set! faces (cons (list a b c) faces))))
		     (loop (sym a)))))
	     (get-edge-ring e)))

	  faces))))
