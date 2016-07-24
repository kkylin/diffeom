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

;;; Manifolds with boundary are probably going to be useful for PDEs:

(declare (usual-integrations))


;;; Boundary charts?  What extra structures are needed?  By convention, a
;;; boundary chart maps the boundary to the half space {x_n >= 0}, so that the
;;; boundary is the space {x_n = 0}.

;;; Of course, any changes made here propagate to tangent and product chart
;;; constructions...    :(

;;; The problem is that the product of two smooth manifolds with boundary will
;;; be a *topological* manifold with boundary, but points (p,q) where p and q
;;; are in the respective boundaries of the component manifolds may not have a
;;; neighborhood that maps to the boundary of a Euclidean half-space.
;;; (Consider the product of the unit interval with itself: There are corners!)

;;; Make these regular domains imbedded inside manifolds and that avoids the
;;; problem -- Can't make those constructions...


(define (add-boundary-to-chart chart i . argl)
  (let ((level 0))

    ;; Locally, the boundary should look like the set of all X such that
    ;; COORD-MAP(X)[i] = LEVEL, where LEVEL is by default 0.

    (if (and (not (null? argl))
	     (real? (car argl)))
	(set! level (car argl)))

    ;; Given the coordinate maps (x_0, ..., x_n), the boundary in the image of
    ;; the chart is the set {x_i = 0}.

    (let ((coord-map (chart:get-coord-map chart))
	  (in-domain? (chart:get-membership-test chart))
	  (in-range? (chart:get-range-test chart)))
      (letrec
	  ((range-boundary?
	    (lambda (x)
	      (and (in-range? x)
		   (almost-equal? level (vector-ref x i)))))

	   (domain-boundary?
	    (lambda (p)
	      (and (in-domain? p)
		   (range-boundary? (coord-map p))))))

	(chart:install-extra chart
			     'boundary-structs
			     (vector i level domain-boundary? range-boundary?))
	chart))))

(define (chart:get-boundary-structs chart)
  (chart:get-extra chart 'boundary-structs))

(define (boundary-chart? chart)
  (if (chart:get-boundary-structs chart)
      #t
      #f))

(define (chart:get-boundary-index chart)
  (let ((result (chart:get-boundary-structs chart)))
    (if result
	(vector-ref result 0)
	#f)))

(define (chart:get-boundary-level chart)
  (let ((result (chart:get-boundary-structs chart)))
    (if result
	(vector-ref result 1)
	#f)))

(define (chart:get-domain-boundary-test chart)
  (let ((result (chart:get-boundary-structs chart)))
    (if result
	(vector-ref result 2)
	#f)))

(define (chart:get-range-boundary-test chart)
  (let ((result (chart:get-boundary-structs chart)))
    (if result
	(vector-ref result 3)
	#f)))

(define (chart:domain-boundary? p chart)
  (let ((boundary? (chart:get-domain-boundary-test chart)))
    (if boundary?
	(boundary? p)
	#f)))

(define (chart:range-boundary? x chart)
  (let ((boundary? (chart:get-range-boundary-test chart)))
    (if boundary?
	(boundary? x)
	#f)))


;;; Make a chart for the boundary out of a chart-with-boundary:

(define (make-boundary-chart chart)
  (let ((boundary-chart (chart:get-extra chart 'boundary-chart)))
    (if boundary-chart
	boundary-chart
	(make-new-boundary-chart chart))))

(define (make-new-boundary-chart chart)
  (let ((in-domain? (chart:get-domain-boundary-test chart))
	(in-range? (chart:get-domain-boundary-test chart))
	(index (chart:get-boundary-index chart))
	(level (chart:get-boundary-level chart))
	(dim (chart:dimension chart)))

    (if (and in-domain? in-range? index level)
	(let ((coord-map (chart:get-coord-map chart))
	      (inverse-map (chart:get-inverse-map chart))

	      (project (lambda (x)
			 (vector:drop-coord x index)))

	      (immerse (lambda (x)
			 (let ((y (vector:add-coord x index)))
			   (vector-set! y index level)
			   y))))

	  (let ((new-coord-map (compose project coord-map))
		(new-inverse-map (compose inverse-map immerse))

		(transition
		 (lambda (Bother)
		   (let ((other (chart:whose-boundary? Bother)))
		     (compose
		      (lambda (x)
			(vector:drop-coord
			 x (chart:get-boundary-index other)))
		      (chart:make-transition-map chart other)
		      immerse)))))

	    (let ((boundary-chart
		   (make-chart (- dim 1) in-domain? in-range?
			       new-coord-map new-inverse-map transition)))

	      (chart:install-extra chart 'boundary-chart boundary-chart)
	      (chart:install-extra
	       boundary-chart 'whose-boundary? (delay chart))
	      boundary-chart)))
	#f)))

(define (chart:whose-boundary? chart)
  (force (chart:get-extra chart 'whose-boundary?)))


;;; Now a manifold with boundary (this may end up being the empty set):

(define (make-boundary-manifold M)
  (let ((charts (manifold:get-finite-atlas M)))
    (if charts

	(let loop ((charts charts) (result '()))
	  (if (null? charts)
	      (if (null? result)
		  #f
		  (charts->manifold result))
	      (let ((boundary-chart (make-boundary-chart (car charts))))
		(if boundary-chart
		    (loop (cdr charts) (cons boundary-chart result))
		    (loop (cdr charts) result)))))

	(let ((find-chart-in-M (manifold:get-general-chart-finder M))
	      (minimize-in-M (manifold:get-general-chart-finder M)))

	  (letrec

	      ((general-find-chart
		(lambda (p . predicates)
		  (call-with-current-continuation
		   (lambda (return)
		     (find-chart-in-M
		      p
		      (lambda (chart)
			(if (chart:domain-boundary? p chart)
			    (let ((new-chart (make-boundary-chart chart)))
			      (let valid? ((predicates predicates))
				(if (null? predicates)
				    (return new-chart)
				    (if ((car predicates) new-chart)
					(valid? (cdr predicates))
					#f))))
			    #f)))))))

	       (find-minimizing-chart
		(lambda (p f <)
		  (cadr (minimize-in-M
			 p
			 (lambda (chart)
			   (if (chart:domain-boundary? p chart)
			       (let ((new-chart (make-boundary-chart chart)))
				 (list new-chart (f new-chart)))
			       #f))
			 (lambda (x y)
			   (or (and x y (< (cadr x) (cadr y)))
			       (and x (false? y))))))))

	       (local-atlas-finder
		(lambda (p)
		  (map (lambda (chart) (make-boundary-chart chart))
		       (manifold:get-local-atlas M p)))))

	    (make-manifold (- (manifold:dimension M) 1)
			   general-find-chart
			   find-minimizing-chart
			   local-atlas-finder))))))
