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

;;; This file defines the structures necessary on charts to facilitate PDE
;;; integration.

(declare (usual-integrations))


;;; When these charts are constructed, the space has already been discretized.
;;; It only remains to compute the appropriate coefficients with respect to the
;;; given differential operator.

(define (make-pde-chart chart extra-nodes discretize complex)

  ;; NODES should be the list of nodes used to discretize this chart.
  ;; DISCRETIZE should be a procedure that takes a differential operator and
  ;; whatever extra arguments it requires to produce a list of equations.  Note
  ;; that it doesn't need a list of nodes because that state can already be
  ;; encapsulated in the proceure.

  (chart:install-extra
   chart 'pde-chart
   (vector (complex->vertices complex)
	   extra-nodes
	   discretize
	   complex
	   (concat-node-list extra-nodes)
	   '()))
  chart)

(define (pde-chart? chart)
  (if (chart:get-extra chart 'pde-chart)
      #t
      #f))

(define (chart:get-vertices chart)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-ref result 0)
	#f)))

(define (chart:get-extra-nodes chart)

  ;; This returns a list of lists, where sublists contain extra nodes for
  ;; corresponding faces in the face list of the complex.

  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-ref result 1)
	#f)))

(define (chart:get-nodes chart)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(append (vector-ref result 0) (vector-ref result 4))
	#f)))

(define (chart:get-discretizer chart)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-ref result 2)
	#f)))

(define (chart:set-discretizer! chart discretize)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-set! result 2 discretize)
	#f)))

(define (chart:get-complex chart)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-ref result 3)
	#f)))

(define (chart:get-elements chart)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-ref result 5)
	#f)))

(define (chart:set-elements! chart elements)
  (let ((result (chart:get-extra chart 'pde-chart)))
    (if result
	(vector-set! result 5 elements)
	#f)))

(define (chart:in-an-element? p chart)
  (and (chart:member? p chart)
       (let ((x (chart:point->coords p chart)))
	 (let loop ((elements (chart:get-elements chart)))
	   (if (null? elements)
	       #f
	       (if (element:member? (car elements) x)
		   #t
		   (loop (cdr elements))))))))

(define (chart:coords->elements x chart)
  (if (chart:in-range? x chart)
      (let loop ((elements (chart:get-elements chart)) (result '()))
	(if (null? elements)
	    result
	    (if (element:member? (car elements) x)
		(loop (cdr elements) (cons (car elements) result))
		(loop (cdr elements) result))))
      '()))

(define (chart:point->elements p chart)
  (chart:coords->elements (chart:point->coords p chart) chart))

(define (chart:node->elements node chart)
  (chart:point->elements (node:get-point node) chart))

(define (chart:coords->any-element x chart)
  (if (chart:in-range? x chart)
      (let loop ((elements (chart:get-elements chart)))
	(if (null? elements)
	    #f
	    (if (element:member? (car elements) x)
		(car elements)
		(loop (cdr elements)))))
      #f))

(define (chart:point->any-element p chart)
  (chart:coords->any-element (chart:point->coords p chart) chart))

(define (chart:node->any-element node chart)
  (chart:point->any-element (node:get-point node) chart))


;;; A useful routine:

(define (concat-node-list node-list)
  (let ((all-nodes (apply append node-list)))

    (for-each (lambda (node)
		(node:set-id! node #f))
	      all-nodes)

    (let loop ((l all-nodes) (result '()))
      (if (null? l)
	  result
	  (let ((node (car l)))
	    (if (node:get-id node)
		(loop (cdr l) result)
		(begin
		  (node:set-id! node #t)
		  (loop (cdr l) (cons node result)))))))))


;;; This actually doesn't do anything, except it gives us the flexibility of
;;; using different discretization algorithms on different charts.  For now,
;;; all we have is Galerkin's method.

(define (chart:discretize-pde chart source extra-args)
  (let ((discretize (chart:get-discretizer chart)))
    (if discretize
	(apply discretize `(,chart ,source ,@extra-args))
	(error "Error: Cannot discretize this chart."))))
