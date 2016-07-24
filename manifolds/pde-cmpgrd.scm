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

;;; Following the ideas and algorithms in Chesshire and Henshaw's paper, let's
;;; try to duplicate (as much as possible) what their program, CMPGRD, does.

;;; Note that they use finite differences and this program uses finite
;;; elements.  So their interpolation kluge doesn't seem nearly as natural
;;; here; one is very tempted to interpolate using finite element basis
;;; functions, which doesn't seem to work.

(declare (usual-integrations))


;;; Classify the nodes as interpolation, discretization, or exterior, as done
;;; in the paper.

(define (cmpgrd:combine-equations domain equations)
  (let* ((cv (list->vector (cons 'foo! (manifold:get-finite-atlas domain))))
	 (chart-count (- (vector-length cv) 1))
	 (coords (make-vector (+ chart-count 1) #f)))

    ;; Make the array one longer than it should so that the array indexing
    ;; conforms to the pseudocode in the paper... (Bletch!)

    (vector-set! cv 0 chart-count)
    (vector-set! coords 0 chart-count)

    ;; First, construct and initialize each entry in the array.

    (write-line '(cmpgrd step 0: initializing array...))

    (do ((k 1 (+ k 1)))
	((> k chart-count))
      (vector-set! coords k
		   (make-vector (length (chart:get-nodes (vector-ref cv k)))
				chart-count)))

    ;; Next, follow the steps in the paper.
    ;; Step 1: Assign IDs to *all* nodes.

    (write-line '(cmpgrd step 1: assigning ids to all nodes...))
    (cmpgrd:step-1/assign-ids cv)

    ;; Step 2: Mark exterior nodes.  (Is this really necessary?)

    (write-line '(cmpgrd step 2: marking exterior nodes...))
    (cmpgrd:step-2/mark-exterior-nodes cv coords)

    ;; Step 3: Find interpolating nodes.

    (write-line '(cmpgrd step 3: finding interpolation nodes...))
    (cmpgrd:step-3/find-interpolating-nodes cv coords)

    ;; Step 4: Mark necessary interpolation nodes.

    (write-line '(cmpgrd step 4: marking necessary interpolation nodes...))
    (cmpgrd:step-4/mark-necessary-interpolation-nodes cv coords)

    ;; Step 5: Delete interpolation points.

    (write-line '(cmpgrd step 5: deleting extra interpolation points...))
    (cmpgrd:step-5/delete-interpolation-points cv coords)

    ;; Step 6: Fix the entries in the table.

    (write-line '(cmpgrd step 6: fixing up table...))
    (cmpgrd:step-6/fix-table cv coords)

    ;; Finally, generate the appropriate constraints and produce a matrix:

    (write-line '(cmpgrd final step: generating matrix...))
    (cmpgrd:final-step/generate-matrix cv coords domain equations)))


;;; Step 1: Assign IDs to *all* nodes, sequentially (beginning with 0) within
;;; each chart.

(define (cmpgrd:step-1/assign-ids cv)
  (let ((n (vector-ref cv 0)))
    (do ((k 1 (+ k 1)))
	((> k n))
      (do ((nodes (chart:get-nodes (vector-ref cv k)) (cdr nodes))
	   (i 0 (+ i 1)))
	  ((null? nodes))
	(let ((node (car nodes)))
	  (node:set-id! node i))))))


;;; Step 2: Mark exterior nodes.

(define (cmpgrd:step-2/mark-exterior-nodes cv coords)
  (let ((n (vector-ref cv 0)))
    (do ((k 1 (+ k 1)))
	((> k n))
      (for-each
       (lambda (node)
	 (if (node:boundary? node)
	     (let ((p (node:get-point node)))
	       (do ((k-prime 1 (+ k-prime 1)))
		   ((> k-prime n))
		 (if (not (= k-prime k))
		     (let* ((nodes (chart:get-nodes (vector-ref cv k-prime)))
			    (node (car nodes)))
		       (let loop ((nodes (cdr nodes))
				  (dist (vector:distance
					 p (node:get-point node)))
				  (id (node:get-id node)))
			 (if (null? nodes)
			     (vector-set! (vector-ref coords k-prime) id 0)
			     (let* ((node (car nodes))
				    (new-d (vector:distance
					    p (node:get-point node))))
			       (if (< new-d dist)
				   (loop (cdr nodes) new-d (node:get-id node))
				   (loop (cdr nodes) dist id)))))))))))
       (chart:get-nodes (vector-ref cv k))))))


;;; Step 3: Find interpolating nodes.

(define (cmpgrd:step-3/find-interpolating-nodes cv coords)
  (let ((n (vector-ref cv 0)))
    (let loop ((count 1))
      (let ((change-count 0))

	(do ((k 1 (+ k 1)))
	    ((> k n))
	  (let ((v (vector-ref coords k)))
	    (for-each

	     (lambda (node)
	       (let ((i (node:get-id node))
		     (p (node:get-point node)))

		 ;; In the paper, l = k-prime.

		 (let ((valid-point? #f))
		   (do ((l (vector-ref v i) (- l 1)))
		       ((or (zero? l) valid-point?))
		     (cond ((= k l)
			    (if (node:local-boundary? node)
				(begin
				  (set! change-count (+ change-count 1))
				  (vector-set! v i (- (vector-ref v i) 1)))))

			   ((not (chart:in-an-element? p (vector-ref cv l)))
			    (set! change-count (+ change-count 1))
			    (vector-set! v i (- (vector-ref v i) 1)))

			   (else (set! valid-point? #t)))))))

	     (chart:get-nodes (vector-ref cv k)))))

	(if (> change-count 0)
	    (begin
	      (write-line `(,change-count changes. number of tries = ,count))
	      (loop (+ count 1))))))))


;;; Step 4: Mark necessary interpolation nodes.

(define (cmpgrd:step-4/mark-necessary-interpolation-nodes cv coords)
  (let ((n (vector-ref cv 0)))
    (do ((k 1 (+ k 1)))
	((> k n))
      (let ((v (vector-ref coords k)))
	(for-each
	 (lambda (node)
	   (let* ((i (node:get-id node))
		  (l (vector-ref v i)))
	     (if (and (< l k) (> l 0))
		 (let ((w (vector-ref coords l)))
		   (for-each
		    (lambda (needed)
		      (let ((j (node:get-id needed)))
			(vector-set! w j (- (abs (vector-ref w j))))))
		    (chart:needed-nodes node (vector-ref cv l)))))))
	 (chart:get-nodes (vector-ref cv k)))))))


;;; Step 5: Delete interpolation points.

(define (cmpgrd:step-5/delete-interpolation-points cv coords)
  (let ((n (vector-ref cv 0)))
    (do ((k 1 (+ k 1)))
	((> k n))
      (let* ((v (vector-ref coords k))
	     (i-nodes (cmpgrd:step-5/get-interpolation-nodes cv v k)))

	(for-each
	 (lambda (node)
	   (let ((i (node:get-id node)))
	     (if (> (vector-ref v i) 0)
		 (vector-set! v i 0))))
	 i-nodes)

	(for-each
	 (lambda (node)
	   (if (not (or (node:boundary? node)
			(node:local-boundary? node)))
	       (vector-set! v (node:get-id node) k)))
	 i-nodes)

	(for-each
	 (lambda (node)
	   (let ((l (vector-ref v (node:get-id node))))
	     (if (> l k)
		 (let ((w (vector-ref coords l)))
		   (for-each
		    (lambda (needed)
		      (let ((j (node:get-id needed)))
			(vector-set! w j (- (abs (vector-ref w j))))))
		    (chart:needed-nodes node (vector-ref cv l)))))))
	 i-nodes)))))

(define (cmpgrd:step-5/get-interpolation-nodes cv v k)
  (let loop ((nodes (chart:get-nodes (vector-ref cv k))) (result '()))
    (if (null? nodes)
	result
	(let ((val (vector-ref v (node:get-id (car nodes)))))
	  (if (not (or (= val k) (zero? k)))
	      (loop (cdr nodes) (cons (car nodes) result))
	      (loop (cdr nodes) result))))))


;;; Step 6: Fix the entries in the table.

(define (cmpgrd:step-6/fix-table cv coords)
  (let ((n (vector-ref cv 0)))
    (do ((k 1 (+ k 1)))
	((> k n))
      (let ((v (vector-ref coords k)))
	(for-each
	 (lambda (node)
	   (let* ((i (node:get-id node))
		  (abs-val (abs (vector-ref v i))))
	     (cond ((= abs-val k)
		    (vector-set! v i abs-val))
		   ((> abs-val 0)
		    (vector-set! v i (- abs-val))))))
	 (chart:get-nodes (vector-ref cv k)))))))


;;; Final step: Generate the interpolation equations and produce the final
;;; matrix:

(define (cmpgrd:final-step/generate-matrix cv coords domain equations)

  ;; First, loop through the charts and pick up the nodes.

  (let ((chart-count (vector-ref cv 0))
	(m 0)
	(n 0)
	(mat #f)
	(constraints '()))

    (let next-chart ((k 1) (result '()))
      (if (<= k chart-count)
	  (let ((chart (vector-ref cv k))
		(v (vector-ref coords k)))
	    (let next-node ((nodes (chart:get-nodes chart))
			    (result result))
	      (if (null? nodes)
		  (next-chart (+ k 1) result)
		  (let* ((node (car nodes))
			 (val (vector-ref v (node:get-id node))))
		    (if (< val 0)
			(let* ((other (vector-ref cv (- val)))
			       (eq (chart:pointwise-constraint node other)))
			  (if eq
			      (next-node (cdr nodes) (cons eq result))
			      (next-node (cdr nodes) result)))
			(next-node (cdr nodes) result))))))
	  (set! constraints result)))

    ;; Next, create the matrix.  The IDs need to be reset first:

    (let loop ((nodes (manifold:get-nodes domain)) (count 0))
      (if (null? nodes)
	  (set! m count)
	  (let ((node (car nodes)))
	    (if (node:boundary? node)
		(begin
		  (node:set-id! node 'boundary-node!)
		  (loop (cdr nodes) count))
		(begin
		  (node:set-id! node count)
		  (loop (cdr nodes) (+ count 1)))))))

    (set! n (+ m 1))
    (set! mat (make-sparse-matrix m n))

    ;; Finally, copy the equations into the matrix while replacing equations
    ;; corresponding to interpolation nodes with the corresponding constraint
    ;; equation.

    (let ((ev (make-vector m #f)))

      ;; Need to keep track of equations:

      (for-each
       (lambda (equation)
	 (vector-set! ev
		      (node:get-id (equation:get-node equation))
		      equation))
       equations)

      ;; Constraints can overwrite equations:

      (for-each
       (lambda (constraint)
	 (vector-set! ev
		      (node:get-id (equation:get-node constraint))
		      constraint))
       constraints)

      ;; Now just copy!

      (do ((i 0 (+ i 1)))
	  ((>= i m) mat)
	(let ((eq (vector-ref ev i)))
	  (if eq
	      (begin
		(sparse-matrix-set! mat i m (equation:get-constant eq))
		(let next-term ((terms (equation:get-terms eq)))
		  (if (not (null? terms))
		      (let* ((term (car terms))
			     (j (term:get-id term))
			     (val (term:get-coeff term)))
			(sparse-matrix-set! mat i j val)
			(next-term (cdr terms))))))
	      (write-line `(*** warning: row ,i of matrix is null!))))))))


;;; With the roles of the nodes figured out, here's the real work: Generate the
;;; appropriate constraints.  First, we need to figure out which nodes to
;;; interpolate from:

(define (chart:needed-nodes node chart)

  ;; Instead of using basis functions to interpolate, maybe we should follow
  ;; Chesshire & Henshaw's suggestion and create higher-order *interpolating
  ;; equations* (constraint equations, in our language) by extending the
  ;; constraint to elements neighboring the one containing the given point.

  ;; For now, let's just use the finite element interpolation and test the rest
  ;; of the Chesshire-Henshaw code.

  (element:get-nodes (chart:node->any-element node chart)))
