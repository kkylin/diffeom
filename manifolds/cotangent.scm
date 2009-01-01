;;; This file defines cotangent bundles as vector bundles (see vbundle.scm).

(declare (usual-integrations))


;;; Make some covectors:

(define (make-cotangent chart x v)

  ;; x is a point in the abstract manifold, and chart contains x.  v is an
  ;; element of the dual tangent space at x, represented in the chart as an
  ;; n-vector.

  (vector 'cotangent chart x v))

(define (cotangent? x)
  (and (vector? x)
       (> (vector-length x) 0)
       (eq? 'cotangent (vector-ref x 0))))

(define (cotangent:get-chart v)
  (vector-ref v 1))

(define (cotangent:get-anchor v)
  (vector-ref v 2))

(define (cotangent:get-coords v)
  (vector-ref v 3))

(define (cotangent:dimension v)
  (vector-length (cotangent:get-coords v)))

(define (make-binary-cotangent-operation op)
  (lambda (v w)
    (let ((p (cotangent:get-anchor v))
	  (q (cotangent:get-anchor w)))
      (if (equal? p q)
	  (let ((chart (cotangent:get-chart v)))
	    (make-tangent chart
			  p
			  (op (cotangent:get-coords v)
			      (chart:pull-back w chart))))
	  (error "Cannot add covectors tangent to different points.")))))

(define cotangent+ (make-binary-cotangent-operation vector:+))
(define cotangent- (make-binary-cotangent-operation vector:-))

(define (cotangent* a v)
  (make-tangent (cotangent:get-chart v)
		(cotangent:get-anchor v)
		(vector:* a (cotangent:get-coords v))))

(define (cotangent:act ctv tv)
  (let ((chart (cotangent:get-chart ctv)))
    (vector:dot (cotangent:get-coords ctv) (chart:push-forward tv chart))))


;;; Pull back a covector along a chart:

(define (chart:pull-back ctv chart)
  (let ((orig (cotangent:get-chart ctv))
	(v (cotangent:get-coords ctv)))
    (if (eq? chart orig)
	v
	(pull-back-in-coords
	 (chart:make-transition-map chart orig)
	 (chart:point->coords (cotangent:get-anchor ctv) chart)
	 v))))


;;; Pull back v from T*f(x) to T*x:

(define (pull-back-in-coords f x v)
  (let ((n (vector-length x)))
    (let ((w (make-vector n))
	  (df ((diff f) x)))

      (do ((i 0 (+ i 1)))
	  ((>= i n) w)
	(vector-set! w i (vector:dot v (df (vector:basis n i 1))))))))


;;; Cotangent charts:

(define (make-cotangent-chart chart)
  (let ((new-chart (chart:get-extra chart 'cotangent-chart)))
    (if new-chart
	(force new-chart)
	(make-new-cotangent-chart chart))))

(define (make-new-cotangent-chart chart)
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
	    (and (in-M-domain? (cotangent:get-anchor v))
		 (dim-vector? (cotangent:get-coords v)))))

	 (in-range?
	  (lambda (v)
	    (and (2*dim-vector? v)
		 (in-M-range? (vector-head v dim)))))

	 (coord-map
	  (lambda (v)
	    (vector-append (M-map (cotangent:get-anchor v))
			   (chart:pull-back v chart))))

	 (inverse-map
	  (lambda (x)
	    (make-cotangent chart
			    (M-inverse (vector-head x dim))
			    (vector-end x dim))))

	 (transition
	  (lambda (Tother)
	    (let ((other (chart:get-base-chart Tother)))
	      (let ((f (chart:make-transition-map chart other))
		    (g (chart:make-transition-map other chart)))
		(lambda (x)
		  (let ((anchor (f (vector-head x dim)))
			(cotangent (vector-end x dim)))
		    (vector-append
		     anchor
		     (pull-back-in-coords g anchor cotangent)))))))))

      (let ((new-chart (make-chart 2*dim in-domain? in-range?
				   coord-map inverse-map transition)))
	(chart:install-extra new-chart 'base-chart (delay chart))
	(chart:install-extra chart 'cotangent-chart (delay new-chart))
	new-chart))))

(define (chart:get-base-chart chart)
  (let ((result (chart:get-extra chart 'base-chart)))
    (if result
	(force result)
	#f)))


;;; Here's how we make a cotangent bundle:

(define (make-cotangent-bundle M)
  (let ((T*M (manifold:get-extra M 'cotangent-bundle)))
    (if T*M
	(force T*M)
	(make-new-cotangent-bundle M))))

(define (make-new-cotangent-bundle M)
  (let ((dim-M (manifold:dimension M)))

    (let ((E
	   (let ((charts (manifold:get-finite-atlas M)))
	     (if charts
		 (charts->manifold (map (lambda (chart)
					  (make-cotangent-chart chart))
					charts))

		 (let ((find-chart-in-M (manifold:get-general-chart-finder M))
		       (minimize-in-M (manifold:get-general-minimizer M)))

		   (letrec
		       ((general-find-chart
			 (lambda (p . predicates)
			   (call-with-current-continuation
			    (lambda (return)
			      (find-chart-in-M
			       (cotangent:get-anchor p)
			       (lambda (chart)
				 (let ((new-chart
					(make-cotangent-chart chart)))
				   (let valid? ((predicates predicates))
				     (if (null? predicates)
					 (return new-chart)
					 (if ((car predicates) new-chart)
					     (valid? (cdr predicates))
					     #f))))))))))

			(find-minimizing-chart
			 (lambda (p f <)
			   (cadr (minimize-in-M
				  (cotangent:get-anchor p)
				  (lambda (chart)
				    (let ((new-chart
					   (make-cotangent-chart chart)))
				      (list new-chart (f new-chart))))
				  (lambda (x y)
				    (< (cadr x) (cadr y)))))))

			(local-atlas-finder
			 (lambda (p)
			   (map (lambda (chart) (make-cotangent-chart chart))
				(manifold:get-local-atlas
				 M (cotangent:get-anchor p))))))

		     (make-manifold (* 2 dim-M)
				    general-find-chart
				    find-minimizing-chart
				    local-atlas-finder))))))

	  (proj cotangent:get-anchor)

	  (fiber
	   (lambda (p)
	     (make-fiber cotangent+ cotangent- cotangent*
			 (lambda (v)
			   (equal? p (cotangent:get-anchor v)))))))

      (let ((T*M (make-vector-bundle M E proj fiber)))
	(manifold:install-extra M 'cotangent-bundle (delay T*M))
	T*M))))
