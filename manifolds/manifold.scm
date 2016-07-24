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

;;; Some obsolete comments:
;;;
;;; How do we represent manifolds as computational objects?  How can we perform
;;; geometric operations in a coordinate-free way?  Can we efficiently encode
;;; charts and mappings?  How can we use the inverse function theorem (or the
;;; implicit function theorem) to automagically construct charts for a
;;; manifold?
;;;
;;; And how much set theory do we need to implement?  This probably depends on
;;; what we want to do with the system.  How can we implement, for example, the
;;; axiom of choice?
;;;
;;; Regarding the question of using the inverse function theorem: We can
;;; probably accomplish this by some kind of numerical differentiation and
;;; using a first-order approximation of the function to define charts.
;;;
;;; Do this quickly and try out ideas.  Don't waste time on completeness or
;;; generality.

;;; December 1996 (from Neal's and Holly's machines):
;;;
;;; Notes on PDEs have been moved to pde.scm, while notes on ODEs are now in
;;; ode.scm.

(declare (usual-integrations))


;;; We need to make some simple Euclidean spaces:

(define (make-new-euclidean-space dim)

  ;; Just need one big happy chart:

  (let* ((test (make-euclidean-test dim))
	 (chart (make-simple-chart dim test test identity identity)))
    (charts->manifold (list chart))))


;;; Will need this a lot in charts, so this will help speed it up a bit:

(define make-euclidean-space (simple-memoize make-new-euclidean-space 26))


;;; Manifolds:

;;; Again, abstract manifolds need only have the right access methods.  This
;;; allows potentially infinite atlases (if the required methods can be
;;; computed efficiently).  It's not even clear the atlases need to be
;;; externally accessible.

;;; To construct a manifold, we need a procedure FIND-CHART that looks up a
;;; chart satisfying some given predicate.  Most other procedures can be
;;; constructed out of this, but it should be arranged so that these procedures
;;; can be replaced, if necessary.

;;; I guess we need to require atlases to be *locally finite*.  Most manifolds
;;; we construct will be compact (or products of compact manifolds with
;;; Euclidean spaces), so it shouldn't matter anyway.

(define (package-manifold-maps dimension
			       general-chart-finder
			       general-chart-minimizer
			       find-chart
			       find-another-chart
			       find-least-distorted
			       get-local-atlas)

  (vector dimension

	  ;; Find a chart containing a given point and satisfying a given list
	  ;; of predicates:
	  general-chart-finder

	  ;; Find a chart containing a given point and minimizing a given
	  ;; function (given an ordering on the function's output):
	  general-chart-minimizer

	  ;; Find a chart containing a given point:
	  find-chart

	  ;; Find a chart containing a given point and not in a given list:
	  find-another-chart

	  ;; Find the least distorted chart containing the given tangent
	  ;; vector:
	  find-least-distorted

	  ;; Find a (finite) set of charts containing a given point.  The
	  ;; procedure is allowed to return (), if p is not in the manifold.
	  ;; Note that everything else is, theoretically, implementable from
	  ;; this.  However, this would be too slow (even for us).
	  get-local-atlas

	  ;; Extra junk:
	  '()))

;;; Here's an easier way to make manifolds: GENERAL-FIND-CHART finds a chart
;;; satisfying a given predicate, and FIND-MINIMIZING-CHART finds a chart that
;;; minimizes a function that always returns either a real *or* #f (#f means
;;; the chart should be thrown out).

;;; As usual, if something cannot be found, #f is returned.

(define (make-manifold dim
		       general-find-chart
		       find-minimizing-chart
		       get-local-atlas)

  (letrec
      ((find-chart
	(lambda (p)
	  (general-find-chart p)))

       (find-another-chart
	(lambda (p charts)
	  (general-find-chart
	   p
	   (lambda (chart)
	     (not (memq chart charts))))))

       (find-least-distorted
	(lambda (tangent)
	  (car (find-minimizing-chart
		(tangent:get-anchor tangent)
		(lambda (chart)
		  (local-distortion chart tangent))
		<)))))

    (package-manifold-maps dim
			   general-find-chart
			   find-minimizing-chart
			   find-chart
			   find-another-chart
			   find-least-distorted
			   get-local-atlas)))


;;; Get the various methods:

(define (manifold:dimension M)
  (vector-ref M 0))

(define (manifold:get-general-chart-finder M)
  (vector-ref M 1))

(define (manifold:get-general-minimizer M)
  (vector-ref M 2))

(define (manifold:get-chart-finder M)
  ;; Return a function that finds a chart containing a given point.
  (vector-ref M 3))

(define (manifold:get-second-opinion M)
  (vector-ref M 4))

(define (manifold:get-least-distorted M)
  (vector-ref M 5))

(define (manifold:get-local-atlas M p)
  ((vector-ref M 6) p))

(define (manifold:install-extra M tag datum)
  (let ((result (assq tag (vector-ref M 7))))
    (if result
	(set-cdr! result datum)
	(vector-set! M 7 (cons (cons tag datum) (vector-ref M 7))))))

(define (manifold:get-extra M tag)
  (let ((result (assq tag (vector-ref M 7))))
    (if result
	(cdr result)
	#f)))

(define (manifold:reset-extra! M)
  (vector-set! M 7 '()))


;;; Some things that are bound to be handy:

(define (manifold:member? M x)
  (if ((manifold:get-chart-finder M) x)
      #t
      #f))

(define (manifold:find-chart M x)

  ;; If the manifold has only one chart, always return it without checking.
  ;; This kludge make smooth functions on tangent bundles of Euclidean spaces
  ;; work with ScmUtils.

  (let ((atlas (manifold:get-finite-atlas M)))
    (if (and atlas (null? (cdr atlas)))
	(car atlas)
	((manifold:get-chart-finder M) x))))

(define (manifold:find-another-chart M x . charts)
  ((manifold:get-second-opinion M) x charts))

(define (manifold:find-least-distorted M tangent)
  ((manifold:get-least-distorted M) tangent))

(define manifold:find-best-chart
  (if *using-scmutils?*
      manifold:find-chart
      (lambda (M x)
	((manifold:get-least-distorted M)
	 (make-tangent (manifold:find-chart M x)
		       x
		       (make-vector (manifold:dimension M) 1))))))


;;; An easy way to construct a large class of manifolds:

(define (charts->manifold charts)

  (if (null? charts)
      (error "No charts given. -- CHARTS->MANIFOLD"))

  (let ((dim (chart:dimension (car charts))))

    (let loop ((charts (cdr charts)))
      (if (not (null? charts))
	  (if (= dim (chart:dimension (car charts)))
	      (loop (cdr charts))
	      (error (string-append "Not all charts have the same dimension!"
				    " -- CHARTS->MANIFOLD")))))

    (letrec
	((general-chart-finder
	  (lambda (p . predicates)
	    (let loop ((charts charts))
	      (if (null? charts)
		  #f
		  (let ((chart (car charts)))
		    (if (chart:member? p chart)
			(let valid? ((predicates predicates))
			  (if (null? predicates)
			      chart
			      (if ((car predicates) chart)
				  (valid? (cdr predicates))
				  (loop (cdr charts)))))
			(loop (cdr charts))))))))

	 (find-minimizing-chart
	  (lambda (p f <)
	    (let loop ((charts charts) (result #f) (min #f))
	      (if (null? charts)
		  (if result
		      (list result min)
		      #f)
		  (let ((chart (car charts)))
		    (if (chart:member? p chart)
			(let ((val (f chart)))
			  (if result
			      (if (< val min)
				  (loop (cdr charts) chart val)
				  (loop (cdr charts) result min))
			      (loop (cdr charts) chart val)))
			(loop (cdr charts) result min)))))))

	 (get-local-atlas
	  (lambda (p)
	    (let loop ((charts charts) (result '()))
	      (if (null? charts)
		  result
		  (let ((chart (car charts)))
		    (if (chart:member? p chart)
			(loop (cdr charts) (cons chart result))
			(loop (cdr charts) result))))))))

      (let ((M (make-manifold dim
			      general-chart-finder
			      find-minimizing-chart
			      get-local-atlas)))
	(manifold:install-extra M 'finite-atlas charts)
	M))))

(define (manifold:get-finite-atlas M)
  (manifold:get-extra M 'finite-atlas))


;;; There are various tools for constructing new manifolds out of old ones,
;;; such as vector bundles and product manifolds.  However, the constructors do
;;; not know about each other: A product of two vector bundles should be a
;;; vector bundle, etc.  I guess if such structures are ever needed, we can
;;; create extra operations.

;;; Note that the product of two manifolds with boundary may have corners on
;;; its boundary, so it may not be a smooth manifold with boundary.  For
;;; example, the product of the unit interval with itself has corners.  (Is
;;; this a problem for manifolds with n > 1?)
