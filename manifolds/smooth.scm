;;; Here is an example to worry about: Let X be a manifold, and let H be a
;;; real-valued smooth function on its cotangent bundle T*X.  Let Y be a
;;; submanifold of X, and let i:Y->X be the inclusion map.  Then we can define
;;; Ti(x,v) = di(x)(v), as usual, and consider H' = Ti*H, the restriction of H
;;; to the submanifold Y.  (Holonomic constraints!)  Now consider this: dH' =
;;; d(Ti*H) = Ti*dH, which basically means (T(Ti))*(dH), right?  And anyway TH'
;;; = T(Ti*H) = T(Ti)*(TH), by covariant functoriality.  Now, the latter
;;; expression is computable under the current system, but does not provide
;;; closure, while the first expression is not even computable, but would
;;; provide closure if it were.  Obviously, this system needs to be
;;; restructured: We should at least be able to express holonomic constraints!

;;; The answer to this is that SMOOTH-MAP:DIFF needs to push the tangent
;;; functor *into* compositions, not pull them out.  So maps are differentiated
;;; in stages rather than as a whole.  This sounds like it stands a chance,
;;; actually.

;;; This file defines the structures for smooth maps between manifolds.  A
;;; smooth map, in addition to a Scheme procedure that computes the point
;;; transformation, should also contain pointers to its domain, range, and
;;; methods for making transition maps.  But then that makes constructing
;;; charts a bit more recursive than necessary, so domains and ranges are not
;;; included (it's not useful information).

(declare (usual-integrations))


;;; How to make one:

(define (make-smooth-map domain range point-function make-transition)

  ;; POINT-FUNCTION should be a scheme procedure that computes the function
  ;; given a point in its domain.  Note that this function is *not* meant to be
  ;; used in manipulations of the function, but only in computing the values.
  ;; In order that the functorial properties are satisfied (whatever that
  ;; means), compositions and exterior differentiation need to be handled
  ;; separately so that the resulting functions are always differentiable.  (In
  ;; particular, this makes ScmUtils work.)

  ;; MAKE-TRANSITION should take as arguments two charts and create a
  ;; transition function between the respective Euclidean spaces.

  (vector domain range point-function make-transition '()))

(define (make-simple-map domain range f)
  (make-smooth-map domain range f (make-simple-transition-maker f)))

(define (make-simple-transition-maker f)
  (lambda (U V)
    (compose (chart:get-coord-map V) f (chart:get-inverse-map U))))


;;; Will probably be useful in our PDE work:

(define R^1 (make-euclidean-space 1))

(define real-line R^1)

(define (make-real-map domain f)
  (make-simple-map domain real-line (compose vector f)))


;;; Accessors:

(define (smooth-map:get-domain f)
  (vector-ref f 0))

(define (smooth-map:get-range f)
  (vector-ref f 1))

(define (smooth-map:get-point-function f)
  (vector-ref f 2))

(define (smooth-map:get-transition-maker f)
  (vector-ref f 3))

(define (smooth-map:make-transition f U V)
  ((smooth-map:get-transition-maker f) U V))

(define (smooth-map:get-extra f tag)
  (let ((result (assq tag (vector-ref f 4))))
    (if result
	(cdr result)
	#f)))

(define (smooth-map:install-extra f tag datum)
  (let ((result (assq tag (vector-ref f 4))))
    (if result
	(set-cdr! result datum)
	(vector-set! f 4 (cons (cons tag datum) (vector-ref f 4))))))


;;; Useful constructs:

(define (apply-smooth-map f p)
  ((smooth-map:get-point-function f) p))

(define (smooth-map:compose f . rest)
  (if (null? rest)
      f
      (let* ((last-guy (car (reverse rest)))
	     (flist (cons f rest))
	     (point-map (apply compose (map smooth-map:get-point-function
					    flist)))
	     (h (make-smooth-map (smooth-map:get-domain last-guy)
				 (smooth-map:get-range f)
				 point-map
				 (make-simple-transition-maker point-map))))
	(smooth-map:install-extra h 'composition flist)
	h)))

(define (smooth-map:decompose f)
  (let ((result (smooth-map:get-extra f 'composition)))
    (if result
	result
	#f)))


;;; It would be nice to make everything else (such as composition and
;;; differentiation) do the right thing with regards to inverses:

(define (make-diffeomorphism f f-inverse)
  (smooth-map:install-extra f 'inverse f-inverse)
  (smooth-map:install-extra f-inverse 'inverse f)
  f)

(define (smooth-map:invert f)
  (let ((result (smooth-map:get-extra f 'inverse)))
    (if result
	result
	#f)))

(define (make-simple-diffeomorphism domain range f g)
  (let ((sf (make-simple-map domain range f))
	(sg (make-simple-map range domain g)))
    (make-diffeomorphism sf sg)))


;;; Some useful (covariant) functors:

;;; This one maps from the category of smooth manifolds into point sets, and
;;; maps smooth maps.

(define forgetful-functor smooth-map:get-point-function)

;;; This uses the differential to map from smooth manifolds into the category
;;; of tangent bundles.

(define (smooth-map:diff f)
  (let ((Tf (smooth-map:get-extra f 'tangent-extension)))
    (if Tf
	Tf
	(let ((flist (smooth-map:decompose f)))
	  (if flist
	      (apply smooth-map:compose (map smooth-map:diff flist))
	      (let ((components (product-map:get-components f)))
		(if components
		    (make-product-map (smooth-map:diff (car component))
				      (smooth-map:diff (cadr component)))
		    (smooth-map:new-diff f))))))))

(define (smooth-map:new-diff smap)
  (let* ((TM (make-tangent-bundle (smooth-map:get-domain smap)))
	 (TN (make-tangent-bundle (smooth-map:get-range smap)))
	 (transit (diff-transition-map smap))
	 (df (diff-point-function smap)))

    (let ((new-map (make-smooth-map TM TN df transit)))
      (smooth-map:install-extra smap 'tangent-extension new-map)
      new-map)))

(define (diff-point-function f)
  (let ((N (smooth-map:get-range f)))
    (lambda (v)
      (let* ((p (tangent:get-anchor v))
	     (q (apply-smooth-map f p))
	     (M-chart (tangent:get-chart v))
	     (N-chart (manifold:find-best-chart N q))
	     (transit (smooth-map:make-transition f M-chart N-chart)))

	(make-tangent N-chart
		      q
		      (push-forward-in-coords transit
					      (chart:point->coords p M-chart)
					      (tangent:get-coords v)))))))


(define (diff-transition-map smap)

  ;; Make a transition map between the tangent charts of two given charts.
  ;; Note that this depends on the fact that SMOOTH-MAP:DIFF decomposes
  ;; compositions of functions into chunks whose transition maps are directly
  ;; differentiable.

  (let ((make-transition-map (smooth-map:get-transition-maker smap)))
    (lambda (TU TV)
      (let* ((U (chart:get-base-chart TU))
	     (V (chart:get-base-chart TV))
	     (f (make-transition-map U V))
	     (dim (chart:dimension U)))

	(lambda (p)
	  (let ((x (vector-head p dim)))
	    (vector-append (f x)
			   (push-forward-in-coords
			    f x (vector-end p dim)))))))))


;;; Another very useul construction:

(define product:combine cons)
(define product:get-arg-1 car)
(define product:get-arg-2 car)

(define (make-product-map f g)
  (let ((fp (smooth-map:get-point-function f))
	(gp (smooth-map:get-point-function g))

	(M-1 (smooth-map:get-domain f))
	(M-2 (smooth-map:get-domain g))

	(N-1 (smooth-map:get-range f))
	(N-2 (smooth-map:get-range g)))

    (let ((point-map
	   (lambda (x)
	     (product:combine (fp (product:get-arg-1 x))
			      (gp (product:get-arg-2 x)))))

	  (make-transition-map
	   (lambda (U V)
	     (let ((dim-1 (chart:dimension U))
		   (dim-2 (chart:dimension V)))
	       (lambda (x)
		 (vector-append (fp (vector-head x dim-1))
				(gp (vector-end x dim-2))))))))

      (let ((f&g (make-smooth-map (product-manifold M1 M2)
				  (product-manifold N1 N2)
				  point-map
				  make-transition-map)))
	(smooth-map:install-extra f&g 'product-map-structs (list f g))
	f&g))))

(define (product-map:get-structs f)
  (smooth-map:get-extra f 'product-map-structs))

(define (product-map:get-components f)
  (let ((result (product-map:get-structs f)))
    (if result
	result
	#f)))


;;; Some useful examples:

(define (make-simple-projection-map n i)

  ;; Make a map from R^n to R^(n-1) by dropping the ith coordinate.

  (make-simple-map (make-euclidean-space n)
		   (make-euclidean-space (- n 1))
		   (lambda (v) (vector:drop-coord v i))))

(define (make-simple-imbedding-map n i)

  ;; Do the opposite:

  (make-simple-map (make-euclidean-space (- n 1))
		   (make-euclidean-space n)
		   (lambda (v) (vector:add-coord v i))))
