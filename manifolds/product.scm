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

;;; Product manifolds.  (Note that, in this scheme (no pun intended), RxR is
;;; not the same as R^2; they are, of course, diffeomorphic.)

(declare (usual-integrations))


;;; First, we need product charts:

(define (make-product-chart chart-1 chart-2)
  (let* ((dim-1 (chart:dimension chart-1))
	 (dim-2 (chart:dimension chart-2))
	 (dim (+ dim-1 dim-2))

	 (coord-map-1 (chart:get-coord-map chart-1))
	 (coord-map-2 (chart:get-coord-map chart-2))

	 (inverse-map-1 (chart:get-inverse-map chart-1))
	 (inverse-map-2 (chart:get-inverse-map chart-2))

	 (euclidean? (make-euclidean-test dim)))

    (letrec
	((in-domain?
	  (lambda (x)
	    (and (list? x)
		 (not (null? (cdr x)))
		 (chart:member? (car x) chart-1)
		 (chart:member? (cadr x) chart-2))))

	 (in-range?
	  (lambda (x)
	    (and (euclidean? x)
		 (chart:in-range? (vector-head x dim-1) chart-1)
		 (chart:in-range? (vector-end x dim-2) chart-2))))

	 (coord-map
	  (lambda (x)
	    (vector-append (coord-map-1 (product:get-arg-1 x))
			   (coord-map-2 (product:get-arg-2 x)))))

	 (inverse-map
	  (lambda (x)
	    (product:combine (inverse-map-1 (vector-head x dim-1))
			     (inverse-map-2 (vector-end x dim-2)))))

	 (transition
	  (lambda (VxW)
	    (let ((components (chart:get-components VxW))
		  (V (car components))
		  (W (cadr components)))
	      (let ((f (chart:make-transition-map chart-1 V))
		    (g (chart:make-transition-map chart-2 W)))
		(lambda (x)
		  (vector-append (f (vector-head x dim-1))
				 (g (vector-end x dim-2)))))))))

      (let ((new-chart (make-chart
			(+ dim-1 dim-2)
			in-domain? in-range? coord-map inverse-map
			transition)))
	(chart:install-extra new-chart 'product-chart (list chart-1 chart-2))
	new-chart))))

(define (chart:get-components chart)
  (chart:get-extra chart 'product-chart))

(define (chart:first-component chart)
  (car (chart:get-components chart)))

(define (chart:second-component chart)
  (cadr (chart:get-components chart)))


;;; This is slow, but should be sufficient for most examples we construct (such
;;; as the 2-torus and products of Euclidean spaces with other manifolds).

(define (product-manifold M1 M2)
  (let ((atlas-1 (manifold:get-finite-atlas M1))
	(atlas-2 (manifold:get-finite-atlas M2)))
    (let ((M
	   (if (and atlas-1 atlas-2)
	       (charts->manifold (map (lambda (l) (apply make-product-chart l))
				      (all-pairs atlas-1 atlas-2)))
	       (make-general-product-manifold M1 M2))))
      (manifold:install-extra M 'product-manifold (list M1 M2))
      M)))

(define (make-general-product-manifold M1 M2)
  (let ((general-find-1 (manifold:get-general-chart-finder M1))
	(general-find-2 (manifold:get-general-chart-finder M2))
	(find-chart-1 (manifold:get-chart-finder M1))
	(find-chart-2 (manifold:get-chart-finder M2))
	(minimize-1 (manifold:get-general-minimizer M1))
	(minimize-2 (manifold:get-general-minimizer M2)))

    (letrec
	((find-chart
	  (lambda (p . predicates)
	    (if (null? predicates)
		(let ((chart-1 (find-chart-1 (car p)))
		      (chart-2 (find-chart-2 (cadr p))))
		  (if (and chart-1 chart-2)
		      (make-product-chart chart-1 chart-2)
		      #f))
		(call-with-current-continuation
		 (lambda (return)
		   (general-find-1
		    (car p)
		    (lambda (chart-1)
		      (general-find-2
		       (cadr p)
		       (lambda (chart-2)
			 (let ((chart (make-product-chart chart-1 chart-2)))
			   (let valid? ((predicates predicates))
			     (if (null? predicates)
				 (return chart)
				 (if ((car predicates) chart)
				     (valid? (cdr predicates))
				     #f)))))))))))))

	 (minimize-chart
	  (lambda (p f <)
	    (cadr (minimize-1
		   (car p)
		   (lambda (chart-1)
		     (minimize-2
		      (cadr p)
		      (lambda (chart-2)
			(let ((chart (make-product-chart chart-1 chart-2)))
			  (list chart (f chart))))
		      (lambda (x y)
			(< (cadr x) (cadr y)))))
		   (lambda (x y)
		     (< (cadr x) (cadr y)))))))

	 (get-local-atlas
	  (lambda (p)
	    (map make-product-chart
		 (manifold:get-local-atlas M1 (car p))
		 (manifold:get-local-atlas M2 (cadr p))))))

      (make-manifold (apply + (map manifold:dimension (list M1 M2)))
		     find-chart minimize-chart get-local-atlas))))

(define (manifold:get-components M)
  (manifold:get-extra M 'product-manifold))

(define (manifold:first-component M)
  (car (manifold:get-components M)))

(define (manifold:second-component M)
  (cadr (manifold:get-components M)))
