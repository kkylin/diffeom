;;; This file defines tangent bundles as vector bundles (see vbundle.scm).

(declare (usual-integrations))


;;; Make some tangent vectors:

(define (make-tangent chart p v)
  ;; p is the (abstract) point to which v is tangent, and v is the *coordinate
  ;; representation* of the tangent vector in the coordinates provided by the
  ;; given chart.
  (vector 'tangent chart p v))

(define (tangent? x)
  (and (vector? x)
       (> (vector-length x) 0)
       (eq? 'tangent (vector-ref x 0))))

(define (tangent:get-chart v)
  (vector-ref v 1))

(define (tangent:get-anchor v)
  (vector-ref v 2))

(define (tangent:get-coords v)
  (vector-ref v 3))

(define (tangent:dimension v)
  (vector-length (tangent:get-coords v)))

(define (make-binary-tangent-operation op)
  (lambda (v w)
    (let ((p (tangent:get-anchor v))
	  (q (tangent:get-anchor w)))
      (if (equal? p q)
	  (let ((chart (tangent:get-chart v)))
	    (make-tangent chart
			  p
			  (op (tangent:get-coords v)
			      (chart:push-forward w chart))))
	  (error "Cannot add vectors tangent to different points.")))))

(define tangent+ (make-binary-tangent-operation vector:+))
(define tangent- (make-binary-tangent-operation vector:-))

(define (tangent* a v)
  (make-tangent (tangent:get-chart v)
		(tangent:get-anchor v)
		(vector:* a (tangent:get-coords v))))


;;; We can measure the distortion by how close the composition of a coordinate
;;; map with its "inverse" comes to the identity...

(define (local-distortion chart tangent)
  (let ((f (chart:make-transition-map chart chart))
	(x (chart:point->coords (tangent:get-anchor tangent) chart))
	(v (chart:push-forward tangent chart)))
    (vector:distance v (((diff f) x) v))))

(define distorted?
  (let ((close-enuf? (make-comparator 1e-5)))
    (lambda (chart v)
      (not (close-enuf? (local-distortion chart v) 0)))))


;;; Push a tangent vector along a chart:

(define (chart:push-forward tv chart)
  (let ((other (tangent:get-chart tv))
	(v (tangent:get-coords tv)))
    (if (eq? chart other)
	v
	(push-forward-in-coords
	 (chart:make-transition-map other chart)
	 (chart:point->coords (tangent:get-anchor tv) other)
	 v))))

(define (push-forward-in-coords f x v)
  (((diff f) x) v))


;;; Tangent charts:

(define (make-tangent-chart chart)
  (let ((new-chart (chart:get-extra chart 'tangent-chart)))
    (if new-chart
	(force new-chart)
	(make-new-tangent-chart chart))))

(define (make-new-tangent-chart chart)
  (let* ((dim (chart:dimension chart))
	 (2*dim (* 2 dim))

	 (in-M-domain? (chart:get-membership-test chart))
	 (in-M-range? (chart:get-range-test chart))

	 (M-map (chart:get-coord-map chart))
	 (M-inverse (chart:get-inverse-map chart))

	 (dim-vector? (make-euclidean-test dim))
	 (2*dim-vector? (make-euclidean-test 2*dim)))

    (letrec
	((in-domain?
	  (lambda (v)
	    (and (in-M-domain? (tangent:get-anchor v))
		 (dim-vector? (tangent:get-coords v)))))

	 (in-range?
	  (lambda (v)
	    (and (2*dim-vector? v)
		 (in-M-range? (vector-head v dim)))))

	 (coord-map
	  (lambda (v)
	    (vector-append (M-map (tangent:get-anchor v))
			   (chart:push-forward v chart))))

	 (inverse-map
	  (lambda (x)
	    (make-tangent chart
			  (M-inverse (vector-head x dim))
			  (vector-end x dim))))

	 (transition
	  (lambda (Tother)
	    (let* ((other (chart:get-base-chart Tother))
		   (f (chart:make-transition-map chart other)))
	      (lambda (x)
		(let ((anchor (vector-head x dim))
		      (tangent (vector-end x dim)))
		  (vector-append (f anchor)
				 (push-forward-in-coords
				  f anchor tangent))))))))

      (let ((new-chart (make-chart 2*dim in-domain? in-range?
				   coord-map inverse-map transition)))
	(chart:install-extra new-chart 'base-chart (delay chart))
	(chart:install-extra chart 'tangent-chart (delay new-chart))
	new-chart))))

(define (chart:get-base-chart chart)
  (let ((result (chart:get-extra chart 'base-chart)))
    (if result
	(force result)
	#f)))


;;; This is sometimes useful for procedures (such as vector fields) that need
;;; to explicitly manipulate charts on tangent bundles:

(define (make-tangent-chart-finder find-chart-in-M)
  (let ((chart-finder
	 (lambda (x . aux)
	   (let ((chart (apply find-chart-in-M
			       (cons (tangent:get-anchor x) aux))))
	     (if chart
		 (make-tangent-chart chart)
		 #f)))))
    chart-finder))


;;; Here's how we make a tangent bundle:

(define (make-tangent-bundle M)
  (let ((TM (manifold:get-extra M 'tangent-bundle)))
    (if TM
	(force TM)
	(make-new-tangent-bundle M))))

(define (make-new-tangent-bundle M)
  (let ((dim-M (manifold:dimension M)))

    (let ((E
	   (let ((charts (manifold:get-finite-atlas M)))
	     (if charts
		 (charts->manifold (map (lambda (chart)
					  (make-tangent-chart chart))
					charts))

		 (let ((find-chart-in-M (manifold:get-general-chart-finder M))
		       (minimize-in-M (manifold:get-general-minimizer M)))

		   (letrec
		       ((general-find-chart
			 (lambda (p . predicates)
			   (call-with-current-continuation
			    (lambda (return)
			      (find-chart-in-M
			       (tangent:get-anchor p)
			       (lambda (chart)
				 (let ((new-chart (make-tangent-chart chart)))
				   (let valid? ((predicates predicates))
				     (if (null? predicates)
					 (return new-chart)
					 (if ((car predicates) new-chart)
					     (valid? (cdr predicates))
					     #f))))))))))

			(find-minimizing-chart
			 (lambda (p f <)
			   (cadr (minimize-in-M
				  (tangent:get-anchor p)
				  (lambda (chart)
				    (let ((new-chart
					   (make-tangent-chart chart)))
				      (list new-chart (f new-chart))))
				  (lambda (x y)
				    (< (cadr x) (cadr y)))))))

			(local-atlas-finder
			 (lambda (p)
			   (map (lambda (chart) (make-tangent-chart chart))
				(manifold:get-local-atlas
				 M (tangent:get-anchor p))))))

		     (make-manifold (* 2 dim-M)
				    general-find-chart
				    find-minimizing-chart
				    local-atlas-finder))))))

	  (proj tangent:get-anchor)

	  (fiber
	   (lambda (p)
	     (make-fiber tangent+ tangent- tangent*
			 (lambda (v)
			   (equal? p (tangent:get-anchor v)))))))

      (let ((TM (make-vector-bundle M E proj fiber)))
	(manifold:install-extra M 'tangent-bundle (delay TM))
	TM))))
