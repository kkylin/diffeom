;;; Some manifolds:

;(declare (usual-integrations cos sin acos atan + - * /))


;;; The n-sphere?  What sort of chart should we use?  Stereogrphic projection?
;;; Mercator projection?  Both?  The representation of points isn't so trivial
;;; in this case.  I guess we'll just use the imbedding, since in this context
;;; it's perfectly natural.

(define (make-imbedded-sphere-test dim)
  (let ((euclidean? (make-euclidean-test (+ dim 1))))
    (lambda (v)
      (and (euclidean? v)
	   (almost-equal? (vector:magnitude^2 v) 1)))))


;;; Do the obvious thing: Stereographic projection.

(define (make-stereographic-chart dim pole-dim pole-dir)
  (let* ((ubound 5)
	 (dim+1 (+ dim 1))
	 (pole (vector:basis dim+1 pole-dim pole-dir)))

    (letrec
	((in-domain?
	  (let ((sphere? (make-imbedded-sphere-test dim)))
	    (lambda (v)
	      (and (sphere? v)
		   (not (almost-equal? (vector:distance^2 v pole) 0))
		   (< (- (/ 4 (vector:magnitude^2 (vector:- v pole))) 1)
		      ubound)))))

	 (in-range?
	  (let ((euclidean? (make-euclidean-test dim)))
	    (lambda (v)
	      (and (euclidean? v)
		   (< (vector:magnitude^2 v) ubound)))))

	 (map
	  (lambda (x)
	    (let* ((d (vector:- x pole))
		   (y (vector:* (/ 2 (vector:magnitude^2 d)) d)))
	      (vector:drop-coord (vector:+ y pole) pole-dim))))

	 (inverse
	  (lambda (x)
	    (let* ((d (vector:- (vector:add-coord x pole-dim) pole))
		   (y (vector:* (/ 2 (vector:magnitude^2 d)) d)))
	      (vector:+ y pole)))))

      (let ((chart (make-simple-chart dim in-domain? in-range? map inverse)))
	(make-spherical-range chart (make-vector dim 0) (sqrt ubound))
	chart))))


;;; Of course, in most applications, it's better to have (generalized)
;;; spherical coordinates...

;;; ILIST should be a permutation of 0, 1, ..., dim.  It determines the order
;;; in which the angles are generated.  The singularity is a half-sphere of
;;; dimension (dim - 1), and is orthogonal to the last coordinate in ILIST,
;;; occupying the negative half space with respect to the next-to-last
;;; coordinate.  ROT should be an angle in radians; the final two coordinates
;;; are rotated by this angle before being generated.

(define (make-planar-rotation angle)
  (let* ((cosine (cos angle))
	 (sine (sin angle))
	 (A (list->matrix 2 2 `(,cosine ,sine ,(- sine) ,cosine))))
    (lambda (v)
      (apply-linear-transformation A v))))

(define (make-spherical-chart dim ilist angle)
  (let ((e1 0)
	(e2 (/ pi 9))
	(dim-1 (- dim 1))
	(dim+1 (+ dim 1))
	(rot (make-planar-rotation angle))
	(-rot (make-planar-rotation (- angle))))

    ;; The membership test is rather circular... (No pun intended! :)

    (letrec
	((in-domain?
	  (let ((sphere? (make-imbedded-sphere-test dim)))
	    (lambda (v)
	      (and (sphere? v)
		   (in-range? (map v))))))

	 (in-range?
	  (let ((euclidean? (make-euclidean-test dim)))
	    (lambda (v)
	      (and (euclidean? v)
		   (let valid? ((i 0))
		     (let ((angle (vector-ref v i)))
		       (if (< i dim-1)
			   (and (< e1 angle)
				(< angle (- pi e1))
				(valid? (+ i 1)))
			   (and (< (+ -pi e2) angle)
				(< angle (- pi e2))))))))))

	 (map
	  (lambda (x)
	    (let ((result (make-vector dim 0)))
	      (let loop ((i 0) (ilist ilist) (r 1))
		(if (= i dim-1)
		    (let ((z (rot (vector (vector-ref x (car ilist))
					  (vector-ref x (cadr ilist))))))
		      (vector-set! result i
				   (atan (vector-ref z 1) (vector-ref z 0)))
		      result)
		    (let ((val (/ (vector-ref x (car ilist)) r)))
		      (vector-set! result i (acos val))
		      (loop (+ i 1)
			    (cdr ilist)
			    (* r (sqrt (- 1 (square val)))))))))))

	 (inverse
	  (lambda (x)
	    (let ((p (make-vector dim+1)))
	      (let loop ((i 0) (ilist ilist) (r 1))
		(let ((angle (vector-ref x i)))
		  (if (< i dim-1)
		      (begin
			(vector-set! p (car ilist) (* r (cos angle)))
			(loop (+ i 1) (cdr ilist) (* r (sin angle))))
		      (let ((z (-rot (vector (cos angle) (sin angle)))))
			(vector-set! p (car ilist) (* r (vector-ref z 0)))
			(vector-set! p (cadr ilist) (* r (vector-ref z 1)))
			p))))))))

      (let ((chart (make-simple-chart dim in-domain? in-range? map inverse))
	    (intervals (make-vector dim (make-interval e1 (- pi e1)))))
	(vector-set! intervals dim-1 (make-interval (+ -pi e2) (- pi e2)))
	(make-cell-range chart intervals)
	chart))))


;;; Here's one way to make a sphere; it turns out to be very hard to define
;;; vector fields on its tangent bundle.  (Try the pendulum!)

(define (make-stereographic-sphere dim)
  (charts->manifold (list (make-stereographic-chart dim 0 1.)
			  (make-stereographic-chart dim 0 -1.))))

;;; Another way:

(define (make-spherical-sphere dim)
  (let* ((l1 (list-integers dim))
	 (l2 (reverse l1)))
    (charts->manifold (list (make-spherical-chart dim l1 0)
			    (make-spherical-chart dim l2 pi)))))


;;; Choose a way to make spheres:

(define make-sphere make-spherical-sphere)


;;; The next thing to make is SO(3).  Note that because the inverse function
;;; theorem can't be used to compute coordinate systems directly, it
;;; complicates the creation of charts for Lie subgroups of GL(n).

;;; For now, we won't bother with explicitly representing Lie group structures
;;; computationally.

(define (make-special-orthogonal-group n)

  ;; Don't bother making SO(n) in general:

  (case n

    ;; n = 2 is just the circle group:
    ((2) (make-sphere 1))

    ;; n = 3 is the rotational group in 3-space:
    ((3) (make-rotational-group))

    ;; Otherwise panic:
    (else 
     (error "Sorry! I only know how to make SO(2) and SO(3)! -- MAKE-SO(n)"))))


;;; Make a planar rotation matrix in n-space, in a plane specified by two
;;; canonical coordinate axes:

(define (make-rotation-matrix dim x-axis y-axis theta)
  (let ((rot (make-matrix dim dim)))
    (do ((i 0 (+ i 1)))
	((>= i dim) rot)
      (do ((j 0 (+ j 1)))
	  ((>= j dim))
	(cond ((and (= i x-axis) (= j x-axis))
	       (matrix-set! rot i j (cos theta)))

	      ((and (= i x-axis) (= j y-axis))
	       (matrix-set! rot i j (- (sin theta))))

	      ((and (= i y-axis) (= j y-axis))
	       (matrix-set! rot i j (cos theta)))

	      ((and (= i y-axis) (= j x-axis))
	       (matrix-set! rot i j (sin theta)))

	      ((= i j)
	       (matrix-set! rot i j 1)))))))


;;; Make the rotational group using Euler angles.  This is a nice test because
;;; they contain singularities:

;;; Need support for Lagrangian and Hamiltonian dynamics, too, if you believe
;;; in such things.

(define make-rotational-group
  (let ((result #f))
    (lambda ()
      (if result
	  result
	  (let ((SO3 (charts->manifold
		      (list (make-euler-angles 0 1 0 0)
			    (make-euler-angles 0 1 0 pi)
			    (make-euler-angles 0 2 pi 0)
			    (make-euler-angles 0 2 pi pi)))))

	    ;; This ensures that we get an atlas.  Might be a bit of an
	    ;; overkill, but...

	    (set! result SO3)
	    result)))))

;;; SO(3) = S^2 x S^1?  But it's probably not easier to do, and this (mostly)
;;; works...

(define (make-euler-angles i-axis j-axis r1 r2)

  ;; General strategy: Given a specification of two axes (I and J), we can
  ;; deduce K.  Using I, J, and K, we construct a spherical chart.  We can then
  ;; decompose the rotation R into an S^2 and an S^1 component by its action on
  ;; some fixed vector v: Rv gives the S^2 component, and its action about the
  ;; axis specified by v gives the S^1 component.  We choose v to be the
  ;; standard z-axis.

  (let* ((id (make-identity-matrix 3))
	 (3-vector? (make-euclidean-test 3))
	 (k-axis (- 3 (+ i-axis j-axis)))
	 (S-chart (make-spherical-chart 2 (list k-axis i-axis j-axis) r1))
	 (C-chart (make-spherical-chart 1 (list 0 1) r2))
	 (xv (vector:basis 3 0 1))
	 (zv (vector:basis 3 2 1))
	 (rot (generate-axis-rotation k-axis 2))
	 (inv-rot (transpose rot)))

    (letrec

	((in-domain?
	  (lambda (R)

	    ;; First, check that it's a matrix:
	    (and (matrix? R)

		 ;; Next, check that it's orthogonal:
		 (let ((diff (matrix:- (matrix:* R (transpose R)) id)))
		   (almost-equal? (matrix:max diff) 0))

		 ;; And then check that it works with the S-chart:
		 (chart:member? (apply-linear-transformation R zv) S-chart)

		 ;; Finally, check that it really checks out completely:
		 (in-range? (coord-map R)))))

	 (in-range?
	  (lambda (x)
	    (and (3-vector? x)
		 (chart:in-range? (vector-head x 1) C-chart)
		 (chart:in-range? (vector-tail x 1) S-chart))))

	 (coord-map
	  (let ((i (vector:basis 3 i-axis 1))
		(j (vector:basis 3 j-axis 1))
		(k (vector:basis 3 k-axis 1)))
	    (lambda (R)
	      (let* ((v (apply-linear-transformation R zv))
		     (coords (chart:point->coords v S-chart))
		     (S (generate-rotation coords r1 i-axis j-axis k-axis))
		     (T (matrix:* rot (transpose S) R))
		     (w (apply-linear-transformation T xv))
		     (psi (chart:point->coords (vector-head w 2) C-chart)))
		(vector-append psi coords)))))

	 (inverse-map
	  (lambda (x)
	    (let ((psi (circle->rotation
			(chart:coords->point (vector-head x 1) C-chart)))
		  (S (generate-rotation (vector-tail x 1) r1
					i-axis j-axis k-axis)))
	      (matrix:* S inv-rot psi)))))

      (make-simple-chart 3 in-domain? in-range? coord-map inverse-map))))

(define (circle->rotation v)
  (let ((mat (make-matrix 3 3)))
    (matrix-set! mat 2 2 1)
    (let ((cos (vector-ref v 0))
	  (sin (vector-ref v 1)))
      (matrix-set! mat 0 0 cos)
      (matrix-set! mat 1 1 cos)
      (matrix-set! mat 0 1 (- sin))
      (matrix-set! mat 1 0 sin))
    mat))

(define (generate-rotation coords angle x-axis y-axis z-axis)

  ;; Generate a rotation matrix that takes the z-axis to the coordinates
  ;; specified, minus an extra rotation of ANGLE about the z-axis.

  (matrix:*
   (make-rotation-matrix 3 x-axis y-axis (+ (vector-ref coords 1) angle))
   (make-rotation-matrix 3 z-axis x-axis (vector-ref coords 0))))

(define (generate-axis-rotation axis-1 axis-2)

  ;; Generate an arbitrary (but deterministic) rotation that takes AXIS-1 to
  ;; AXIS-2:

  (if (= axis-1 axis-2)
      (make-identity-matrix 3)
      (let ((i (vector:basis 3 axis-1 1))
	    (j (vector:basis 3 axis-2 1))
	    (R (make-rotation-matrix 3 axis-1 axis-2 (/ pi 2))))
	(if (almost-zero? (vector:distance
			   (apply-linear-transformation R i) j))
	    R
	    (transpose R)))))


;;; We need the unit closed n-ball for testing the PDE solver.

(define (make-ball n . argl)

  ;; It might be easier to do this in terms of the n-sphere.  The unit interval
  ;; (n = 1) case would have to be treated separately, then, since in that case
  ;; the boundary is a 0-dimensional manifold.

  (let ((make-sphere make-sphere))

    ;; Really, no matter how the sphere is made, so long as it has a finite
    ;; atlas, this will work.  But the user might have a preference, e.g. when
    ;; generating the mesh for solving PDEs.

    (if (and (not (null? argl))
	     (procedure? (car argl)))
	(set! make-sphere (car argl)))

    (cond ((= n 1)
	   (error "Sorry!  The unit interval hasn't been implemented yet!"))

	  ((>= n 2)
	   (let ((boundary (make-spherical-sphere (- n 1))))
	     (let* ((n-vector? (make-euclidean-test n))
		    (in? (lambda (p)
			   (and (n-vector? p)
				(< (vector:magnitude p) 2/3))))
		    (center-chart (make-simple-chart
				   n in? in? identity identity)))

	       ;; Define a chart that covers a neighborhood of the center, then
	       ;; define the rest as a deformation retract onto the boundary.

	       (make-spherical-range center-chart (make-vector n 0) 2/3)

	       (let loop ((B-charts (list center-chart))

			  ;; The sphere should have a finite atlas.

			  (S-charts (manifold:get-finite-atlas boundary)))

		 (if (null? S-charts)

		     (charts->manifold B-charts)

		     (let ((S-chart (car S-charts)))
		       (let ((S-coord-map (chart:get-coord-map S-chart))
			     (S-inverse (chart:get-inverse-map S-chart))
			     (in-S-domain? (chart:get-membership-test S-chart))
			     (in-S-range? (chart:get-range-test S-chart)))

			 ;; Construct a chart on the ball out of a chart on the
			 ;; sphere:

			 (let ((coord-map
				(lambda (p)
				  (let* ((len (vector:magnitude p))
					 (x (vector:add-coord
					     (S-coord-map (vector:* (/ len) p))
					     0)))
				    (vector-set! x 0 len)
				    x)))

			       (inverse-map
				(lambda (x)
				  (vector:* (vector-ref x 0)
					    (S-inverse (vector-tail x 1)))))

			       (in-domain?
				(lambda (p)
				  (and (n-vector? p)
				       (let ((len (vector:magnitude p)))
					 (and (< 1/3 len)
					      (<= len 1)
					      (in-S-domain?
					       (vector:* (/ len) p)))))))

			       (in-range?
				(lambda (x)
				  (and (n-vector? x)
				       (let ((len (vector-ref x 0)))
					 (and (< 1/3 len) (<= len 1)))
				       (in-S-range? (vector-tail x 1))))))

			   (let ((B-chart (make-simple-chart
					   n in-domain? in-range?
					   coord-map inverse-map)))

			     ;; This thing has a boundary:

			     (add-boundary-to-chart B-chart 0 1)

			     ;; Check the range type; this is useful for
			     ;; meshing.  Note that we only handle n-cells
			     ;; for now.

			     (if (chart:cell-range? S-chart)
				 (let ((int (cell-range:get-interval-list
					     S-chart)))
				   (make-cell-range
				    B-chart (cons (make-interval 1/3 1) int))))

			     ;; Keep going:

			     (loop (cons B-chart B-charts)
				   (cdr S-charts)))))))))))

	  (else (error "Error: Invalid argument. -- MAKE-BALL")))))
