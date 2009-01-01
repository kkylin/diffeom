;;; Charts:

(declare (usual-integrations))


;;; Abstract charts need only contain the right maps.  What they actually do is
;;; up to the particular implementation.

;;; Might have been nice to make charts out of smooth maps, but that might be
;;; more trouble than it's worth.  It's too recursive, and the abstraction has
;;; to bottom out somewhere.  (Why?  Because otherwise it wouldn't be
;;; computable!)  Charts will be made out of structures much like smooth
;;; functions.  We'll try to merge these structures if it appears possible.

(define (make-chart dim in-domain? in-range? coord-map inverse-map transition)

  ;; TRANSITION should be a function that, given another chart, returns a
  ;; transition function to the other chart from this one.  (Within reasonable
  ;; ranges, of course.)

  (vector in-domain? in-range? coord-map inverse-map transition dim '()))

(define (make-simple-chart dim in-domain? in-range? coord-map inverse-map)
  (make-chart dim in-domain? in-range? coord-map inverse-map
	      (lambda (V)
		(compose (chart:get-coord-map V) inverse-map))))


;;; Get the various maps out:

(define (chart:get-membership-test chart)
  ;; Should return a function that tests whether a point is in the chart.
  (vector-ref chart 0))

(define (chart:get-range-test chart)
  (vector-ref chart 1))

(define (chart:get-coord-map chart)
  ;; Should provide the mapping from the manifold to Euclidean coordinates.
  (vector-ref chart 2))

(define (chart:get-inverse-map chart)
  ;; Should provide the inverse of the above.
  (vector-ref chart 3))

(define (chart:get-transition-maker chart)
  (vector-ref chart 4))

(define (chart:dimension chart)
  (vector-ref chart 5))

(define (chart:install-extra chart tag datum)
  (let ((result (assq tag (vector-ref chart 6))))
    (if result
	(set-cdr! result datum)
	(vector-set! chart 6 (cons (cons tag datum) (vector-ref chart 6))))))

(define (chart:get-extra chart tag)
  (let ((result (assq tag (vector-ref chart 6))))
    (if result
	(cdr result)
	#f)))


;;; Some useful wrappers for debugging purposes:

(define (domain-check f chart)
  (lambda (p)
    (if (not (chart:member? p chart))
	(write-line '(warning! stepping out of domain!)))
    (f p)))

(define (range-check g chart)
  (lambda (x)
    (if (not (chart:in-range? x chart))
	(write-line '(warning! stepping out of range!)))
    (g x)))


;;; Some methods that are bound to be handy:

(define (chart:member? x chart)
  ((chart:get-membership-test chart) x))

(define (chart:in-range? x chart)
  ((chart:get-range-test chart) x))

(define (chart:point->coords x chart)
  ((chart:get-coord-map chart) x))

(define (chart:coords->point x chart)
  ((chart:get-inverse-map chart) x))

(define (chart:make-transition-map U V)
  ((chart:get-transition-maker U) V))


;;; Turn the chart maps into smooth maps:

(define (chart:get-range U)
  (make-euclidean-space (chart:dimension U)))

(define (chart:get-domain chart)
  (let ((U (chart:get-extra chart 'chart-as-manifold)))
    (if U
	(force U)
	(let ((U (charts->manifold (list chart))))
	  (chart:install-extra chart 'chart-as-manifold (delay U))
	  U))))

(define (chart:smooth-coord-map chart)
  (make-smooth-map (chart:get-domain chart)
		   (chart:get-range chart)
		   (chart:get-coord-map chart)
		   (lambda (U V) (chart:get-coord-map chart))))

(define (chart:smooth-inverse-map chart)
  (make-smooth-map (chart:get-range chart)
		   (chart:get-domain chart)
		   (chart:get-inverse-map chart)
		   (lambda (U V) (chart:get-inverse-map chart))))


;;; A faster distortion test to compute for the ODE integrator:

(define (chart:stable-coords? x chart)
  (chart:in-range?
   (chart:point->coords (chart:coords->point x chart) chart)
   chart))
